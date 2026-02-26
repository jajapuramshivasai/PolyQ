#[cfg(test)]
mod phase_poly_tests {
    use super::*;
    use approx::assert_relative_eq;
    use num_complex::Complex;

    // Helper to verify state vectors with high precision
    fn assert_state_eq(actual: &[Complex<f64>], expected: &[Complex<f64>]) {
        for (a, b) in actual.iter().zip(expected.iter()) {
            assert_relative_eq!(a.re, b.re, epsilon = 1e-10);
            assert_relative_eq!(a.im, b.im, epsilon = 1e-10);
        }
    }

    #[test]
    fn test_clifford_plus_s_circuit() {
        // Test Clifford+S circuit: H(0), S(0), CNOT(0,1)
        let circuit = QuantumCircuit {
            num_qubits: 2,
            gates: vec![
                Gate {
                    name: "h".to_string(),
                    qubits: vec![0],
                    params: vec![],
                },
                Gate {
                    name: "s".to_string(),
                    qubits: vec![0],
                    params: vec![],
                },
                Gate {
                    name: "h".to_string(),
                    qubits: vec![1],
                    params: vec![],
                },
                Gate {
                    name: "cz".to_string(),
                    qubits: vec![0, 1],
                    params: vec![],
                },
                Gate {
                    name: "h".to_string(),
                    qubits: vec![1],
                    params: vec![],
                },
            ],
        };

        let poly = phase_polynomial_z8(&circuit);
        let input = vec![false, false]; // |00>
        let state = simulate_phase_polynomial_z8(&poly, &input);

        // Expected: 1/sqrt(2) (|00> + i|11>)
        let norm = (0.5f64).sqrt();
        let expected = vec![
            Complex::new(norm, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, norm), // This is the 'i' component
        ];
        assert_state_eq(&state, &expected);
    }

    #[test]
    fn test_clifford_plus_t_relative_phase() {
        // Test Clifford+T circuit: H(0), T(0), CNOT(0,1)
        let circuit = QuantumCircuit {
            num_qubits: 2,
            gates: vec![
                Gate {
                    name: "h".to_string(),
                    qubits: vec![0],
                    params: vec![],
                },
                Gate {
                    name: "t".to_string(),
                    qubits: vec![0],
                    params: vec![],
                },
                Gate {
                    name: "h".to_string(),
                    qubits: vec![1],
                    params: vec![],
                },
                Gate {
                    name: "cz".to_string(),
                    qubits: vec![0, 1],
                    params: vec![],
                },
                Gate {
                    name: "h".to_string(),
                    qubits: vec![1],
                    params: vec![],
                },
            ],
        };

        let poly = phase_polynomial_z8(&circuit);
        let state = simulate_phase_polynomial_z8(&poly, &vec![false, false]);

        // Expected: 1/sqrt(2) |00> + (1+i)/2 |11>
        // Note: (1+i)/2 has magnitude 1/sqrt(2)
        let expected = vec![
            Complex::new((0.5f64).sqrt(), 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.5, 0.5), // (1+i)/2
        ];
        assert_state_eq(&state, &expected);
    }

    #[test]
    fn test_ry_approximation_sequence() {
        // Test Ry(pi/4) approximation using H-T-H-S sequence
        let circuit = QuantumCircuit {
            num_qubits: 1,
            gates: vec![
                Gate {
                    name: "h".to_string(),
                    qubits: vec![0],
                    params: vec![],
                },
                Gate {
                    name: "t".to_string(),
                    qubits: vec![0],
                    params: vec![],
                },
                Gate {
                    name: "h".to_string(),
                    qubits: vec![0],
                    params: vec![],
                },
                Gate {
                    name: "s".to_string(),
                    qubits: vec![0],
                    params: vec![],
                },
            ],
        };

        let poly = phase_polynomial_z8(&circuit);
        let state = simulate_phase_polynomial_z8(&poly, &vec![false]);

        // Manual calculation: H-T-H-S |0> = [(1 + e^{iπ/4})/2, i(1 - e^{iπ/4})/2]
        let e_ipi4 = Complex::new((2.0f64).sqrt() / 2.0, (2.0f64).sqrt() / 2.0);
        let expected = vec![
            (Complex::new(1.0, 0.0) + e_ipi4) * 0.5,
            (Complex::new(1.0, 0.0) - e_ipi4) * Complex::new(0.0, 0.5),
        ];

        assert_state_eq(&state, &expected);
    }
}

use crate::quantum_circuit::{Gate, QuantumCircuit};
use num_complex::Complex;
use std::collections::HashMap;

// PolyQ-style phase polynomial for Z8 (Clifford+T+S) circuits
#[derive(Debug, Clone)]
pub struct PhasePolynomialZ8 {
    pub terms: Vec<(i32, Vec<usize>)>, // (weight, [var indices])
    pub wire_array: Vec<Vec<usize>>,   // wire_array[q] = [var indices on q]
    pub num_qubits: usize,
    pub num_vars: usize, // total number of variables (n + #Hadamards)
}

// Build PolyQ-style phase polynomial for Clifford+T+S circuits (Z8)
pub fn phase_polynomial_z8(circuit: &QuantumCircuit) -> PhasePolynomialZ8 {
    let n = circuit.num_qubits;
    let mut wire_array: Vec<Vec<usize>> = (0..n).map(|i| vec![i]).collect();
    let mut max_new_var = n;
    let mut terms = Vec::new();
    for gate in &circuit.gates {
        let name = gate.name.to_lowercase();
        let qubits = &gate.qubits;
        match name.as_str() {
            "h" => {
                let q = qubits[0];
                wire_array[q].push(max_new_var);
                // Hadamard: add quadratic term between last two variables on this wire, weight=4
                let last = wire_array[q].len();
                terms.push((4, vec![wire_array[q][last - 2], wire_array[q][last - 1]]));
                max_new_var += 1;
            }
            "cz" => {
                let q0 = qubits[0];
                let q1 = qubits[1];
                terms.push((
                    4,
                    vec![
                        *wire_array[q0].last().unwrap(),
                        *wire_array[q1].last().unwrap(),
                    ],
                ));
            }
            "s" => {
                let q = qubits[0];
                terms.push((2, vec![*wire_array[q].last().unwrap()]));
            }
            "t" => {
                let q = qubits[0];
                terms.push((1, vec![*wire_array[q].last().unwrap()]));
            }
            "sdg" => {
                let q = qubits[0];
                terms.push((6, vec![*wire_array[q].last().unwrap()]));
            }
            "tdg" => {
                let q = qubits[0];
                terms.push((7, vec![*wire_array[q].last().unwrap()]));
            }
            _ => {}
        }
    }
    PhasePolynomialZ8 {
        terms,
        wire_array,
        num_qubits: n,
        num_vars: max_new_var,
    }
}

#[cfg(test)]
mod bell_extra {
    use super::*;
    use approx::assert_relative_eq;
    use num_complex::Complex;
    #[test]
    fn test_basic_bell_state() {
        // Bell state |Φ+> = (|00> + |11>)/sqrt(2)
        // Circuit: H(0); CZ(0,1)
        let circuit = QuantumCircuit {
            num_qubits: 2,
            gates: vec![
                Gate {
                    name: "h".to_string(),
                    qubits: vec![0],
                    params: vec![],
                },
                Gate {
                    name: "h".to_string(),
                    qubits: vec![1],
                    params: vec![],
                },
                Gate {
                    name: "cz".to_string(),
                    qubits: vec![0, 1],
                    params: vec![],
                },
                Gate {
                    name: "h".to_string(),
                    qubits: vec![1],
                    params: vec![],
                },
            ],
        };
        let poly = phase_polynomial_z8(&circuit);
        let input = vec![false, false];
        let state = simulate_phase_polynomial_z8(&poly, &input);
        let norm = (0.5f64).sqrt();
        let expected = vec![
            Complex::new(norm, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(norm, 0.0),
        ];
        for (a, b) in state.iter().zip(expected.iter()) {
            assert_relative_eq!(a.re, b.re, epsilon = 1e-10);
            assert_relative_eq!(a.im, b.im, epsilon = 1e-10);
        }
    }
    // End of bell_extra mod
}

// Simulate PolyQ-style phase polynomial for Z8 (Clifford+T+S) circuits
pub fn simulate_phase_polynomial_z8(
    poly: &PhasePolynomialZ8,
    input_bitstring: &[bool],
) -> Vec<Complex<f64>> {
    let n = poly.num_qubits;
    let t = poly.num_vars;
    let ovs: Vec<usize> = poly.wire_array.iter().map(|w| *w.last().unwrap()).collect();
    let dim = 1 << n;
    let mut s_ldic: HashMap<usize, [u32; 8]> = HashMap::new();
    let x_range = 1 << (t - n);
    let mut x = vec![false; t];
    // Fill input variables
    for (i, &b) in input_bitstring.iter().enumerate() {
        x[i] = b;
    }
    // For each assignment to intermediate variables
    for i in 0..x_range {
        // Set intermediate variables
        for (ind, bit) in (0..(t - n))
            .map(|j| ((i >> (t - n - 1 - j)) & 1))
            .enumerate()
        {
            x[n + ind] = bit == 1;
        }
        // Evaluate phase polynomial mod 8
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
        // Output bitstring (output variables, big-endian order, fix variable order)
        let mut chosenbits = 0usize;
        for (j, &ov) in ovs.iter().rev().enumerate() {
            if x[ov] {
                chosenbits |= 1 << j;
            }
        }
        let arr = s_ldic.entry(chosenbits).or_insert([0u32; 8]);
        arr[val_out as usize] += 1;
    }
    // Now, for each output, compute amplitude using FFT hardcoded for Z8
    let mut state = vec![Complex::new(0.0, 0.0); dim];
    let norm = (2f64).powf(-0.5 * ((t - n) as f64));
    println!("[DEBUG] n = {}, t = {}, norm = {}", n, t, norm);
    for (k, arr) in s_ldic.iter() {
        // Z8 FFT: sum_{j=0}^7 arr[j] * exp(i*pi/4*j)
        let mut amp = Complex::new(0.0, 0.0);
        for j in 0..8 {
            let theta = (std::f64::consts::PI / 4.0) * (j as f64);
            let c = Complex::from_polar(arr[j] as f64, theta);
            amp += c;
        }
        state[*k] = amp * norm;
    }
    println!("[DEBUG] state = {:?}", state);
    state
}
