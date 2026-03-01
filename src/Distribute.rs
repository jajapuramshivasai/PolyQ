use crate::Z8::engine::PhasePolynomialZ8;
use num_complex::Complex;
use rayon::prelude::*;
use std::f64::consts::PI;

/// Compute statevector for a Z8 phase polynomial with some fixed input bits.
///
/// explanation (step by step):
/// - Let `t = poly.num_vars` and `n = poly.num_qubits`.
/// - We give some input bits in `input_bitstring`. These bits are fixed.
/// - The remaining `free_bits = t - n` bits we will sum over (all possible assignments).
/// - For each assignment of the free bits we make the full boolean vector `x`.
/// - We read the output logical bits from `poly.wire_array` (use the last element of each wire).
/// - We pack the output logical bits into an integer `chosenbits` using MSB-first order.
///   That means the first element of `poly.wire_array` becomes the most-significant bit.
/// - For this `x` we compute the Z8 phase residue `val_out` by evaluating the terms in `poly.terms`.
/// - We count how many times each pair `(chosenbits, val_out)` happens across all free-bit assignments.
/// - After counting, amplitude for index `chosenbits` = sum over `r=0..7` of
///   `count( chosenbits, r ) * exp(i * r * PI/4)`, then multiply by normalization
///   `2^{-0.5 * free_bits}`.
///
/// Memory note:
/// - To save memory, each worker only keeps counts for its own output range `start..end`.
///   This keeps memory per worker low when we run many workers.
///
/// Indexing note:
/// - Indexing is MSB-first as described above.
pub fn distributed_statevector_z8(
    poly: &PhasePolynomialZ8,
    input_bitstring: &[bool],
    num_workers: usize,
) -> Vec<Complex<f64>> {
    let n = poly.num_qubits;
    let dim = 1usize << n;
    let t = poly.num_vars;
    assert!(t >= n, "num_vars must be >= num_qubits");
    let free_bits = t - n;
    let x_range = 1usize << free_bits;

    // Precompute wire outputs (one per logical output bit) and sanity-check
    let ovs: Vec<usize> = poly.wire_array.iter().map(|w| *w.last().unwrap()).collect();
    assert_eq!(ovs.len(), n, "wire_array must contain one output per qubit");

    // Map each free variable to the terms that depend on it for incremental updates
    let terms_by_free_bit: Vec<Vec<usize>> = if free_bits > 0 {
        let mut map = vec![Vec::new(); free_bits];
        for (term_idx, (_weight, idxs)) in poly.terms.iter().enumerate() {
            for &var in idxs {
                if let Some(free_idx) = var.checked_sub(n) {
                    if free_idx < free_bits {
                        map[free_idx].push(term_idx);
                    }
                }
            }
        }
        map
    } else {
        Vec::new()
    };

    // Precompute constant factors
    let norm = (2f64).powf(-0.5 * (free_bits as f64));
    let phase_weights: [Complex<f64>; 8] = [
        Complex::from_polar(1.0, 0.0),
        Complex::from_polar(1.0, PI / 4.0),
        Complex::from_polar(1.0, PI / 2.0),
        Complex::from_polar(1.0, 3.0 * PI / 4.0),
        Complex::from_polar(1.0, PI),
        Complex::from_polar(1.0, 5.0 * PI / 4.0),
        Complex::from_polar(1.0, 3.0 * PI / 2.0),
        Complex::from_polar(1.0, 7.0 * PI / 4.0),
    ];

    // split the free-variable space instead of output indices
    let y_chunk = (x_range + num_workers - 1) / num_workers;

    // Each worker returns a partial amplitude vector for the full output dimension
    let partials: Vec<Vec<Complex<f64>>> = (0..num_workers)
        .into_par_iter()
        .map(|worker_id| {
            let y_start = worker_id * y_chunk;
            let y_end = ((worker_id + 1) * y_chunk).min(x_range);

            // create per-worker partial state vector
            let mut partial = vec![Complex::new(0.0, 0.0); dim];

            // prepare assignment vector and fill input bits
            let mut x = vec![false; t];
            for (j, &b) in input_bitstring.iter().enumerate() {
                x[j] = b;
            }

            // initialize free bits to y_start's Gray code, and compute initial term activity/phase
            let mut term_active = vec![false; poly.terms.len()];
            let mut current_phase = 0i32;
            if y_start < y_end {
                // set free bits according to the Gray code of y_start
                let mut gray = y_start ^ (y_start >> 1);
                for ind in 0..free_bits {
                    let bit = (gray >> (free_bits - 1 - ind)) & 1;
                    x[n + ind] = bit == 1;
                }
                // compute activity & phase
                for (term_idx, (weight, idxs)) in poly.terms.iter().enumerate() {
                    let active = idxs.iter().all(|&j| x[j]);
                    term_active[term_idx] = active;
                    if active {
                        current_phase = (current_phase + *weight).rem_euclid(8);
                    }
                }

                for y in y_start..y_end {
                    // compute output index
                    let mut chosenbits = 0usize;
                    for (j, &ov) in ovs.iter().enumerate() {
                        if x[ov] {
                            chosenbits |= 1 << (n - 1 - j);
                        }
                    }
                    let phase_idx = current_phase as usize;
                    partial[chosenbits] += phase_weights[phase_idx];

                    if y + 1 < y_end {
                        let next_gray = (y + 1) ^ ((y + 1) >> 1);
                        let changed = gray ^ next_gray;
                        let flipped_bit = changed.trailing_zeros() as usize;
                        debug_assert!(flipped_bit < free_bits);
                        let free_idx = free_bits - 1 - flipped_bit;
                        let new_bit = ((next_gray >> flipped_bit) & 1) == 1;
                        x[n + free_idx] = new_bit;
                        for &term_idx in &terms_by_free_bit[free_idx] {
                            let term = &poly.terms[term_idx];
                            let new_active = term.1.iter().all(|&j| x[j]);
                            if new_active != term_active[term_idx] {
                                term_active[term_idx] = new_active;
                                let delta = if new_active { term.0 } else { -term.0 };
                                current_phase = (current_phase + delta).rem_euclid(8);
                            }
                        }
                        gray = next_gray;
                    }
                }
            }

            partial
        })
        .collect();

    // reduce partials into a single statevector
    let mut state = vec![Complex::new(0.0, 0.0); dim];
    for p in partials {
        for (i, &c) in p.iter().enumerate() {
            state[i] += c;
        }
    }
    // apply normalization
    for amp in &mut state {
        *amp *= norm;
    }
    state
}

/// Evaluate a function `f` over all bitstrings of length `num_vars`, distributing the work across `num_workers`.
/// Each worker evaluates a chunk of the statevector and the results are combined at the end.
pub fn distributed_evaluate<F>(num_vars: usize, num_workers: usize, f: F) -> Vec<Complex<f64>>
where
    F: Fn(usize) -> Complex<f64> + Sync,
{
    let dim = 1 << num_vars;
    let chunk_size = (dim + num_workers - 1) / num_workers;
    // Parallelize over chunks
    (0..num_workers)
        .into_par_iter()
        .flat_map(|worker_id| {
            let start = worker_id * chunk_size;
            let end = ((worker_id + 1) * chunk_size).min(dim);
            (start..end).map(|i| f(i)).collect::<Vec<_>>()
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Z8::engine::phase_polynomial_z8;
    use crate::quantum_circuit::{Gate, QuantumCircuit};
    use num_complex::Complex;
    use std::time::Instant;
    #[test]
    fn test_distributed_evaluate() {
        let num_vars = 4;
        let num_workers = 4;
        // Example: f(i) = Complex::new(i as f64, 0.0)
        let result = distributed_evaluate(num_vars, num_workers, |i| Complex::new(i as f64, 0.0));
        assert_eq!(result.len(), 1 << num_vars);
        assert_eq!(result[0], Complex::new(0.0, 0.0));
        assert_eq!(result[15], Complex::new(15.0, 0.0));
    }

    #[test]
    fn test_distributed_statevector_bell() {
        // Bell state: H(0); CZ(0,1)
        let circuit = QuantumCircuit {
            num_qubits: 2,
            gates: vec![
                Gate {
                    name: "h".to_string(),
                    qubits: vec![0],
                    params: vec![],
                },
                Gate {
                    name: "cz".to_string(),
                    qubits: vec![0, 1],
                    params: vec![],
                },
            ],
        };
        let poly = phase_polynomial_z8(&circuit);
        let input = vec![false, false];
        // 1 processor
        let start1 = Instant::now();
        let state1 = distributed_statevector_z8(&poly, &input, 1);
        let dur1 = start1.elapsed();
        println!("Bell state (1 proc): {:?} in {:?}", state1, dur1);
        // 2 processors
        let start2 = Instant::now();
        let state2 = distributed_statevector_z8(&poly, &input, 2);
        let dur2 = start2.elapsed();
        println!("Bell state (2 proc): {:?} in {:?}", state2, dur2);
        // Both should match expected: H on q0 then CZ(0,1) gives (|00> + |10>)/sqrt(2)
        let norm = (0.5f64).sqrt();
        let expected = vec![
            Complex::new(norm, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(norm, 0.0),
            Complex::new(0.0, 0.0),
        ];
        for (a, b) in state1.iter().zip(expected.iter()) {
            assert!((a - b).norm() < 1e-10);
        }
        for (a, b) in state2.iter().zip(expected.iter()) {
            assert!((a - b).norm() < 1e-10);
        }
    }

    #[test]
    fn test_distributed_statevector_large_qasm() {
        // parse the 21-qubit QASM from benchmarks and compare with baseline
        use std::path::PathBuf;
        let qasm_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("Benchmark/qasm2_exports/clifford_t_21q_11.qasm2");
        let circuit = QuantumCircuit::from_qasm_file(&qasm_path)
            .expect("failed to parse qasm export");
        let class = circuit.classify();
        let transpiled = crate::Transpile_circuit::transpile_to_gateset(&circuit, class);
        let poly = phase_polynomial_z8(&transpiled);
        let input = vec![false; transpiled.num_qubits];
        // compute distributed with several workers
        // compute baseline once
        let base = crate::Z8::engine::simulate_phase_polynomial_z8(&poly, &input);
        // run with varying worker counts and record timings
        for &workers in &[1usize, 2, 4] {
            let start = Instant::now();
            let dist = distributed_statevector_z8(&poly, &input, workers);
            let dur = start.elapsed();
            println!("workers {} -> duration {:?}", workers, dur);
            assert_eq!(dist.len(), base.len());
            for (a, b) in dist.iter().zip(base.iter()) {
                assert!((a - b).norm() < 1e-10);
            }
            let norm_sq: f64 = dist.iter().map(|c| c.norm_sqr()).sum();
            assert!((norm_sq - 1.0).abs() < 1e-10);
        }
    }
}

/*
Example: 4 qubits, 16 basis states, 4 workers
Spawn all workers with f(input,x_i),g(x_i) ,output_bitstring range specific to worker



Worker 1        || {0000 to 0011} \.  ||
Worker 2        || {0100 to 0111} /\. || final state vector
Worker 2        || {1000 to 1011} \/. ||
Worker 2        || {1100 to 1111} /.  ||

### Algorithm: Distributed Phase Polynomial Simulation with Gray Code
**Input**: PhasePolynomialZ8 (terms, wire_array), input_bitstring, num_workers
**Output**: Statevector (Vec<Complex>)

1.  **Preprocessing & Setup**
    * Set `n` = number of qubits; `t` = total variables.
    * Set `free_bits` = `t - n`.
    * Set `x_range` = $2^{free\_bits}$ (Use `u128` to prevent overflow).
    * **Bit-to-Term Mapping**: 
        * Create a map where each free bit index `i` points to a list of terms containing it.
        * Store each term as `(weight, bitmask_of_other_bits)` for fast dependency checking.
    * **Partitioning**: Divide the output statevector range ($2^n$) into `num_workers` chunks.

2.  **Parallel Execution (For each Worker)**
    * Define `start` and `end` indices for the worker's assigned chunk.
    * Allocate `counts[local_dim][8]` initialized to zero (stores frequency of phase residues 0-7).

3.  **Sum-over-Paths with Gray Code**
    * Initialize `gray_code` = 0.
    * Initialize `current_phase` = evaluate_polynomial(input_bitstring, all_free_bits_zero).
    
    * **Loop** `y` from 0 to `x_range - 1`:
        * **Step A: Determine Output Index**
            * Calculate `chosenbits` based on `gray_code` and `wire_array` (Logical mapping).
        
        * **Step B: Update Counts**
            * If `chosenbits` is within `[start, end)`:
                * Increment `counts[chosenbits - start][current_phase]`.

        * **Step C: Incremental Transition** (Skip on last iteration)
            * Identify the bit `i` that flips in the next Gray Code: `next_y = (y+1) ^ ((y+1) >> 1)`.
            * `flipped_bit_idx = trailing_zeros(y ^ next_y)`.
            * For each `term` in `bit_to_terms[flipped_bit_idx]`:
                * If all *other* bits in the term mask are `1` in the current `gray_code`:
                    * If bit `i` is flipping `0 -> 1`: `current_phase = (current_phase + weight) mod 8`.
                    * Else (flipping `1 -> 0`): `current_phase = (current_phase - weight) mod 8`.
            * `gray_code = next_y`.

4.  **Amplitude Synthesis**
    * For each `index` in worker's range:
        * `amplitude = 0`.
        * For `r` from 0 to 7:
            * If `counts[index][r] > 0`:
                * `amplitude += counts[index][r] * exp(i * r * PI / 4)`.
        * `statevector[index] = amplitude * 2^(-0.5 * free_bits)`.

5.  **Collect and Return** combined statevector from all workers.

Complexity Reduction: Standard evaluation is O(Terms⋅2 ^free_bits). 
Gray Code updates reduce this to O(AvgTermsPerBit⋅2 ^free_bits).

Memory Efficiency: By using 8-slot integer counters per basis state instead of complex numbers, we significantly reduce cache misses and memory pressure during the hot loop.

Branch Pruning: The bitmask check (gray_code & mask) == mask allows the CPU to skip entire groups of terms that are currently "inactive" (evaluating to zero)
*/

/*
TODO:
1. when we distribute across a lot of threads there are redundent variables to loop over for each thred
like few MSB so for each worker first pre evaluate that and do {optimise_polynomial}
for each worker to get a smaller polynomial with fewer variables to loop over
gray code to cache non changing variables across iterations of the inner loop
2. Integrate with other simulation techniques e.g weight cycle, sympletic etc
3. Implement rsmpi support for distributed execution across machines
4. NTT based acceleration
*/
