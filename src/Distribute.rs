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
    let chunk_size = (dim + num_workers - 1) / num_workers;

    // Precompute wire outputs (one per logical output bit) and sanity-check
    let ovs: Vec<usize> = poly.wire_array.iter().map(|w| *w.last().unwrap()).collect();
    assert_eq!(ovs.len(), n, "wire_array must contain one output per qubit");

    // Precompute constant factors
    let norm = (2f64).powf(-0.5 * (free_bits as f64));
    let phase_angles: [f64; 8] = [
        0.0,
        PI / 4.0,
        PI / 2.0,
        3.0 * PI / 4.0,
        PI,
        5.0 * PI / 4.0,
        3.0 * PI / 2.0,
        7.0 * PI / 4.0,
    ];

    (0..num_workers)
        .into_par_iter()
        .flat_map_iter(|worker_id| {
            let start = worker_id * chunk_size;
            let end = ((worker_id + 1) * chunk_size).min(dim);

            // Only allocate storage for the worker's chunk to save memory.
            // This uses about (end-start) * 8 * 4 bytes per worker instead of dim * 8 * 4.
            let local_dim = end - start;
            let mut counts = vec![[0u32; 8]; local_dim];

            // Reusable buffer for full variable assignment
            let mut x = vec![false; t];
            for (j, &b) in input_bitstring.iter().enumerate() {
                x[j] = b;
            }

            // Iterate over all assignments to free variables once
            for y in 0..x_range {
                // Fill free variables. Preserve previous MSB-first mapping used in original code:
                for ind in 0..free_bits {
                    let bit = (y >> (free_bits - 1 - ind)) & 1;
                    x[n + ind] = bit == 1;
                }

                // Compute output index (chosenbits) using MSB-first mapping to match prior behavior
                let mut chosenbits = 0usize;
                for (j, &ov) in ovs.iter().enumerate() {
                    if x[ov] {
                        chosenbits |= 1 << (n - 1 - j);
                    }
                }

                // Compute phase residue (0..7)
                let mut val_out: u8 = 0;
                for (weight, idxs) in &poly.terms {
                    let mut v = true;
                    for &j in idxs {
                        v &= x[j];
                    }
                    if v {
                        val_out = (val_out + (*weight as u8)) % 8;
                    }
                }

                // Only record counts for indices that belong to this worker's chunk.
                if chosenbits >= start && chosenbits < end {
                    counts[chosenbits - start][val_out as usize] += 1;
                }
            }

            // Now convert counts for this worker's range into amplitudes
            (start..end).map(move |i| {
                let arr = &counts[i - start];
                // If all zeros, amplitude is zero
                if arr.iter().all(|&c| c == 0) {
                    return Complex::new(0.0, 0.0);
                }
                let mut amp = Complex::new(0.0, 0.0);
                for j in 0..8 {
                    // magnitude = arr[j], angle = phase_angles[j]
                    amp += Complex::from_polar(arr[j] as f64, phase_angles[j]);
                }
                amp * norm
            })
        })
        .collect()
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
*/
