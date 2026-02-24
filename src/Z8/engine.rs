use crate::quantum_circuit::{QuantumCircuit, Gate};
use num_complex::Complex;

/// Represents a quadratic phase polynomial for a Clifford+T+S circuit (Z8 phases)
/// Only linear (Z, S, Sdg, T, Tdg) and quadratic (CZ) terms are allowed. Phases are multiples of pi/4.
#[derive(Debug, Clone)]
pub struct PhasePolynomialZ8 {
	/// Each term is a tuple: (qubit indices, phase)
	/// - Linear: vec![q] for Z(q), S(q), Sdg(q), T(q), Tdg(q)
	/// - Quadratic: vec![q0, q1] for CZ(q0, q1)
	/// phase is in units of pi/4 (i.e., 1 means pi/4, 2 means pi/2, ..., 8 means 2pi, etc.)
	pub terms: Vec<(Vec<usize>, i32)>,
	pub num_qubits: usize,
}

/// Construct the quadratic phase polynomial for a Clifford+T+S circuit
/// Only linear (Z, S, Sdg, T, Tdg) and quadratic (CZ) terms are included.
pub fn phase_polynomial_z8(circuit: &QuantumCircuit) -> PhasePolynomialZ8 {
	let mut terms = Vec::new();
	let n = circuit.num_qubits;
	for gate in &circuit.gates {
		match gate.name.to_lowercase().as_str() {
			"z" => {
				// Z(q): phase pi if q=1 (4 * pi/4)
				let q = gate.qubits[0];
				terms.push((vec![q], 4));
			},
			"s" => {
				// S(q): phase pi/2 if q=1 (2 * pi/4)
				let q = gate.qubits[0];
				terms.push((vec![q], 2));
			},
			"sdg" | "sdag" => {
				// Sdg(q): phase -pi/2 if q=1 (6 * pi/4 mod 2pi)
				let q = gate.qubits[0];
				terms.push((vec![q], 6));
			},
			"t" => {
				// T(q): phase pi/4 if q=1
				let q = gate.qubits[0];
				terms.push((vec![q], 1));
			},
			"tdg" | "tdag" => {
				// Tdg(q): phase -pi/4 if q=1 (7 * pi/4 mod 2pi)
				let q = gate.qubits[0];
				terms.push((vec![q], 7));
			},
			"cz" => {
				// CZ(q0, q1): phase pi if q0=1 and q1=1 (4 * pi/4)
				let q0 = gate.qubits[0];
				let q1 = gate.qubits[1];
				terms.push((vec![q0, q1], 4));
			},
			_ => {
				// Ignore all other gates for phase polynomial (no cubic or higher terms)
			},
		}
	}
	PhasePolynomialZ8 { terms, num_qubits: n }
}

/// Compute the final statevector from a Z8 phase polynomial and input bitstring
pub fn simulate_phase_polynomial_z8(poly: &PhasePolynomialZ8, input_bitstring: &[bool]) -> Vec<Complex<f64>> {
	let n = poly.num_qubits;
	let dim = 1 << n;
	let mut state = vec![Complex::new(0.0, 0.0); dim];
	// For quadratic phase polynomials, the statevector is:
	// |psi> = 1/sqrt(2^n) sum_x exp(i (pi/4) P(x)) |x>
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
		let theta = (std::f64::consts::PI / 4.0) * (phase as f64);
		state[i] = Complex::from_polar(norm, theta);
	}
	state
}
