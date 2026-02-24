#[cfg(test)]
mod additional_z4_tests {
    use super::*;
    use approx::assert_relative_eq;
    use num_complex::Complex;

    // Helper to compare state vectors
    fn assert_state_eq(actual: &[Complex<f64>], expected: &[Complex<f64>]) {
        for (a, b) in actual.iter().zip(expected.iter()) {
            assert_relative_eq!(a.re, b.re, epsilon = 1e-10);
            assert_relative_eq!(a.im, b.im, epsilon = 1e-10);
        }
    }

    #[test]
    fn test_sdag_gate() {
        // Test Sdg gate: should produce |0> - i|1> / sqrt(2)
        let circuit = QuantumCircuit {
            num_qubits: 1,
            gates: vec![
                Gate { name: "h".to_string(), qubits: vec![0], params: vec![] },
                Gate { name: "sdg".to_string(), qubits: vec![0], params: vec![] },
            ],
        };
        let poly = phase_polynomial_z4(&circuit);
        let state = simulate_phase_polynomial_z4(&poly, &[false]);
        
        let norm = (0.5f64).sqrt();
        let expected = vec![
            Complex::new(norm, 0.0),
            Complex::new(0.0, -norm), // -i component
        ];
        assert_state_eq(&state, &expected);
    }

    #[test]
    fn test_z_gate() {
        // Test Z gate: should produce |0> - |1> / sqrt(2)
        let circuit = QuantumCircuit {
            num_qubits: 1,
            gates: vec![
                Gate { name: "h".to_string(), qubits: vec![0], params: vec![] },
                Gate { name: "z".to_string(), qubits: vec![0], params: vec![] },
            ],
        };
        let poly = phase_polynomial_z4(&circuit);
        let state = simulate_phase_polynomial_z4(&poly, &[false]);
        
        let norm = (0.5f64).sqrt();
        let expected = vec![
            Complex::new(norm, 0.0),
            Complex::new(-norm, 0.0),
        ];
        assert_state_eq(&state, &expected);
    }

    #[test]
    fn test_double_hadamard_identity() {
        // Double Hadamard should return to initial state (up to global phase)
        let circuit = QuantumCircuit {
            num_qubits: 1,
            gates: vec![
                Gate { name: "h".to_string(), qubits: vec![0], params: vec![] },
                Gate { name: "h".to_string(), qubits: vec![0], params: vec![] },
            ],
        };
        let poly = phase_polynomial_z4(&circuit);
        let state = simulate_phase_polynomial_z4(&poly, &[false]); // Input |0>
        
        // Expected: |0> (magnitudes 1, 0)
        assert_relative_eq!(state[0].norm(), 1.0, epsilon = 1e-10);
        assert_relative_eq!(state[1].norm(), 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_cz_entanglement() {
        // Test CZ entanglement: produces 1/2(|00> + |01> + |10> - |11>)
        let circuit = QuantumCircuit {
            num_qubits: 2,
            gates: vec![
                Gate { name: "h".to_string(), qubits: vec![0], params: vec![] },
                Gate { name: "h".to_string(), qubits: vec![1], params: vec![] },
                Gate { name: "cz".to_string(), qubits: vec![0, 1], params: vec![] },
            ],
        };
        let poly = phase_polynomial_z4(&circuit);
        let state = simulate_phase_polynomial_z4(&poly, &[false, false]);
        
        let expected = vec![
            Complex::new(0.5, 0.0),
            Complex::new(0.5, 0.0),
            Complex::new(0.5, 0.0),
            Complex::new(-0.5, 0.0),
        ];
        assert_state_eq(&state, &expected);
    }
}


use crate::quantum_circuit::{QuantumCircuit, Gate};
use num_complex::Complex;
use std::collections::HashMap;

// PolyQ-style phase polynomial for Z4 (Clifford+S) circuits
#[derive(Debug, Clone)]
pub struct PhasePolynomialZ4 {
    pub terms: Vec<(i32, Vec<usize>)>, // (weight, [var indices])
    pub wire_array: Vec<Vec<usize>>,   // wire_array[q] = [var indices on q]
    pub num_qubits: usize,
    pub num_vars: usize,               // total number of variables (n + #Hadamards)
}

// Build PolyQ-style phase polynomial for Clifford+S circuits (Z4)
pub fn phase_polynomial_z4(circuit: &QuantumCircuit) -> PhasePolynomialZ4 {
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
				// Hadamard: add quadratic term between last two variables on this wire, weight=2 (Z4)
				let last = wire_array[q].len();
				terms.push((2, vec![wire_array[q][last-2], wire_array[q][last-1]]));
				max_new_var += 1;
			},
			"z" => {
				let q = qubits[0];
				terms.push((2, vec![*wire_array[q].last().unwrap()]));
			},
			"cz" => {
				let q0 = qubits[0];
				let q1 = qubits[1];
				terms.push((2, vec![*wire_array[q0].last().unwrap(), *wire_array[q1].last().unwrap()]));
			},
			"s" => {
				let q = qubits[0];
				terms.push((1, vec![*wire_array[q].last().unwrap()]));
			},
			"sdg" | "sdag" => {
				let q = qubits[0];
				terms.push((3, vec![*wire_array[q].last().unwrap()]));
			},
			_ => {}
		}
	}
	PhasePolynomialZ4 {
		terms,
		wire_array,
		num_qubits: n,
		num_vars: max_new_var,
	}
}

// Simulate PolyQ-style phase polynomial for Z4 (Clifford+S) circuits
pub fn simulate_phase_polynomial_z4(poly: &PhasePolynomialZ4, input_bitstring: &[bool]) -> Vec<Complex<f64>> {
	let n = poly.num_qubits;
	let t = poly.num_vars;
	let ovs: Vec<usize> = poly.wire_array.iter().map(|w| *w.last().unwrap()).collect();
	let dim = 1 << n;
	let mut s_ldic: HashMap<usize, [u32; 4]> = HashMap::new();
	let x_range = 1 << (t - n);
	let mut x = vec![false; t];
	// Fill input variables
	for (i, &b) in input_bitstring.iter().enumerate() { x[i] = b; }
	// For each assignment to intermediate variables
	for i in 0..x_range {
		// Set intermediate variables
		for (ind, bit) in (0..(t-n)).map(|j| ((i >> (t-n-1-j)) & 1)) .enumerate() {
			x[n+ind] = bit == 1;
		}
		// Evaluate phase polynomial mod 4
		let mut val_out: u8 = 0;
		for (weight, idxs) in &poly.terms {
			let mut v = true;
			for &j in idxs {
				v &= x[j];
			}
			if v {
				val_out = (val_out + (*weight as u8)) % 4;
			}
		}
		// Output bitstring (output variables, big-endian order)
		let mut chosenbits = 0usize;
		for (j, &ov) in ovs.iter().rev().enumerate() {
			if x[ov] { chosenbits |= 1 << j; }
		}
		let arr = s_ldic.entry(chosenbits).or_insert([0u32; 4]);
		arr[val_out as usize] += 1;
	}
	// Now, for each output, compute amplitude using FFT hardcoded for Z4
	let mut state = vec![Complex::new(0.0, 0.0); dim];
	let norm = (2f64).powf(-0.5 * ((t-n) as f64));
	for (k, arr) in s_ldic.iter() {
		// Z4 FFT: sum_{j=0}^3 arr[j] * exp(i*pi/2*j)
		let mut amp = Complex::new(0.0, 0.0);
		for j in 0..4 {
			let theta = (std::f64::consts::PI / 2.0) * (j as f64);
			let c = Complex::from_polar(arr[j] as f64, theta);
			amp += c;
		}
		state[*k] = amp * norm;
	}
	state
}
