#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_bell_states_phase_polynomial() {
        // Bell state circuit: H(0); CZ(0,1)
        let circuit = QuantumCircuit {
            num_qubits: 2,
            gates: vec![
                Gate { name: "h".to_string(), qubits: vec![0], params: vec![] },
                Gate { name: "cz".to_string(), qubits: vec![0, 1], params: vec![] },
            ],
        };
        let poly = phase_polynomial(&circuit);
        let norm = (0.5f64).sqrt();

        // |Φ+> = (|00> + |11>)/sqrt(2)
        let input = vec![false, false];
        let state = simulate_phase_polynomial(&poly, &input);
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

        // |Φ-> = (|00> - |11>)/sqrt(2)
        let input = vec![false, true];
        let state = simulate_phase_polynomial(&poly, &input);
        let expected = vec![
            Complex::new(norm, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(-norm, 0.0),
        ];
        for (a, b) in state.iter().zip(expected.iter()) {
            assert_relative_eq!(a.re, b.re, epsilon = 1e-10);
            assert_relative_eq!(a.im, b.im, epsilon = 1e-10);
        }

        // |Ψ+> = (|01> + |10>)/sqrt(2)
        let input = vec![true, false];
        let state = simulate_phase_polynomial(&poly, &input);
        let expected = vec![
            Complex::new(0.0, 0.0),
            Complex::new(norm, 0.0),
            Complex::new(norm, 0.0),
            Complex::new(0.0, 0.0),
        ];
        for (a, b) in state.iter().zip(expected.iter()) {
            assert_relative_eq!(a.re, b.re, epsilon = 1e-10);
            assert_relative_eq!(a.im, b.im, epsilon = 1e-10);
        }

        // |Ψ-> = (|01> - |10>)/sqrt(2)
        let input = vec![true, true];
        let state = simulate_phase_polynomial(&poly, &input);
        let expected = vec![
            Complex::new(0.0, 0.0),
            Complex::new(norm, 0.0),
            Complex::new(-norm, 0.0),
            Complex::new(0.0, 0.0),
        ];
        for (a, b) in state.iter().zip(expected.iter()) {
            assert_relative_eq!(a.re, b.re, epsilon = 1e-10);
            assert_relative_eq!(a.im, b.im, epsilon = 1e-10);
        }
    }
}

use crate::quantum_circuit::{QuantumCircuit, Gate};

/// Represents a quadratic phase polynomial for a Clifford circuit
/// Only linear (Z) and quadratic (CZ) terms are allowed.
#[derive(Debug, Clone)]
pub struct PhasePolynomial {
    /// Each term is a tuple: (qubit indices, phase)
    /// - Linear: vec![q] for Z(q)
    /// - Quadratic: vec![q0, q1] for CZ(q0, q1)
    /// phase is in units of pi (i.e., 1 means pi, 2 means 2pi, etc.)
    pub terms: Vec<(Vec<usize>, i32)>,
    pub num_qubits: usize,
}

/// Construct the quadratic phase polynomial for a Clifford circuit
/// Only linear (Z) and quadratic (CZ) terms are included.
pub fn phase_polynomial(circuit: &QuantumCircuit) -> PhasePolynomial {
    let mut terms = Vec::new();
    let n = circuit.num_qubits;
    for gate in &circuit.gates {
        match gate.name.to_lowercase().as_str() {
            "z" => {
                // Linear term: Z(q)
                let q = gate.qubits[0];
                terms.push((vec![q], 1));
            },
            "cz" => {
                // Quadratic term: CZ(q0, q1)
                let q0 = gate.qubits[0];
                let q1 = gate.qubits[1];
                terms.push((vec![q0, q1], 1));
            },
            _ => {
                // Ignore all other gates for phase polynomial (no cubic or higher terms)
            },
        }
    }
    PhasePolynomial { terms, num_qubits: n }
}

use num_complex::Complex;

/// Compute the final statevector from a phase polynomial and input bitstring
pub fn simulate_phase_polynomial(poly: &PhasePolynomial, input_bitstring: &[bool]) -> Vec<Complex<f64>> {
    let n = poly.num_qubits;
    let dim = 1 << n;
    let mut state = vec![Complex::new(0.0, 0.0); dim];
    // For quadratic phase polynomials, the statevector is:
    // |psi> = 1/sqrt(2^n) sum_x exp(i pi P(x)) |x>
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
        let theta = std::f64::consts::PI * (phase as f64);
        state[i] = Complex::from_polar(norm, theta);
    }
    state
}

