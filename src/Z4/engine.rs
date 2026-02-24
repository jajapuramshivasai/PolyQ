use crate::quantum_circuit::{QuantumCircuit, Gate};
use num_complex::Complex;

/// Represents a quadratic phase polynomial for a Clifford+S circuit (Z4 phases)
/// Only linear (Z, S) and quadratic (CZ) terms are allowed. Phases are multiples of pi/2.
#[derive(Debug, Clone)]
pub struct PhasePolynomialZ4 {
	/// Each term is a tuple: (qubit indices, phase)
	/// - Linear: vec![q] for Z(q) or S(q)
	/// - Quadratic: vec![q0, q1] for CZ(q0, q1)
	/// phase is in units of pi/2 (i.e., 1 means pi/2, 2 means pi, 3 means 3pi/2, 4 means 2pi, etc.)
	pub terms: Vec<(Vec<usize>, i32)>,
	pub num_qubits: usize,
}

/// Construct the quadratic phase polynomial for a Clifford+S circuit
/// Only linear (Z, S) and quadratic (CZ) terms are included.
pub fn phase_polynomial_z4(circuit: &QuantumCircuit) -> PhasePolynomialZ4 {
	let mut terms = Vec::new();
	let n = circuit.num_qubits;
	for gate in &circuit.gates {
		match gate.name.to_lowercase().as_str() {
			"z" => {
				// Z(q): phase pi if q=1 (2 * pi/2)
				let q = gate.qubits[0];
				terms.push((vec![q], 2));
			},
			"s" => {
				// S(q): phase pi/2 if q=1
				let q = gate.qubits[0];
				terms.push((vec![q], 1));
			},
			"cz" => {
				// CZ(q0, q1): phase pi if q0=1 and q1=1 (2 * pi/2)
				let q0 = gate.qubits[0];
				let q1 = gate.qubits[1];
				terms.push((vec![q0, q1], 2));
			},
			"sdg" | "sdag" => {
				// Sdg(q): phase -pi/2 if q=1, which is 3 * pi/2 mod 2pi
				let q = gate.qubits[0];
				terms.push((vec![q], 3));
			},
			_ => {
				// Ignore all other gates for phase polynomial (no cubic or higher terms)
			},
		}
	}
	PhasePolynomialZ4 { terms, num_qubits: n }
}

/// Compute the final statevector from a Z4 phase polynomial and input bitstring
pub fn simulate_phase_polynomial_z4(poly: &PhasePolynomialZ4, input_bitstring: &[bool]) -> Vec<Complex<f64>> {
	let n = poly.num_qubits;
	let dim = 1 << n;
	let mut state = vec![Complex::new(0.0, 0.0); dim];
	// For quadratic phase polynomials, the statevector is:
	// |psi> = 1/sqrt(2^n) sum_x exp(i (pi/2) P(x)) |x>
	// where P(x) = sum_j a_j x_j + sum_{j<k} b_{jk} x_j x_k
	let norm = (1.0 / (dim as f64).sqrt()) as f64;
	for i in 0..dim {
		let x: Vec<bool> = (0..n).map(|j| ((i >> (n - 1 - j)) & 1) == 1).collect();
		let mut phase = 0i32;
		for (qs, p) in &poly.terms {
			if qs.len() == 1 {
				if x[qs[0]] { phase += p; }
			} else if qs.len() == 2 {
				if x[qs[0]] && x[qs[1]] { phase += p; }
			}
		}
		let theta = (std::f64::consts::PI / 2.0) * (phase as f64);
		state[i] = Complex::from_polar(norm, theta);
	}
	state
}
