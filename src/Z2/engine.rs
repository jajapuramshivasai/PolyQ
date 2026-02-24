use crate::quantum_circuit::QuantumCircuit;

#[derive(Debug, Clone)]
pub struct PhasePolynomial {
    pub terms: Vec<(i32, Vec<usize>)>,
    pub num_qubits: usize,
    pub num_vars: usize, // total variables (input + intermediate)
    pub output_vars: Vec<usize>, // variable index for each output qubit
}

pub fn phase_polynomial(circuit: &QuantumCircuit) -> PhasePolynomial {
    let n = circuit.num_qubits;
    let mut wire_array: Vec<Vec<usize>> = (0..n).map(|i| vec![i]).collect();
    let mut max_new_var = n;
    let mut terms = Vec::new();
    for gate in &circuit.gates {
        match gate.name.to_lowercase().as_str() {
            "h" => {
                let q = gate.qubits[0];
                wire_array[q].push(max_new_var);
                terms.push((1, vec![wire_array[q][wire_array[q].len()-2], wire_array[q][wire_array[q].len()-1]]));
                max_new_var += 1;
            },
            "z" => {
                let q = gate.qubits[0];
                terms.push((1, vec![*wire_array[q].last().unwrap()]));
            },
            "cz" => {
                let q0 = gate.qubits[0];
                let q1 = gate.qubits[1];
                terms.push((1, vec![*wire_array[q0].last().unwrap(), *wire_array[q1].last().unwrap()]));
            },
            _ => {},
        }
    }
    // Output variable for each qubit is the last variable on its wire
    let output_vars = wire_array.iter().map(|w| *w.last().unwrap()).collect();
    PhasePolynomial {
        terms,
        num_qubits: n,
        num_vars: max_new_var,
        output_vars,
    }
}

use num_complex::Complex;

pub fn simulate_phase_polynomial(poly: &PhasePolynomial, input_bitstring: &[bool]) -> Vec<Complex<f64>> {
    use std::collections::HashSet;
    let n = poly.num_qubits;
    let t = poly.num_vars;
    let dim = 1 << n;
    let mut state = vec![Complex::new(0.0, 0.0); dim];

    // PolyQ: input variables are 0..n, output variables are poly.output_vars
    let input_vars: HashSet<usize> = (0..n).collect();
    let output_vars: HashSet<usize> = poly.output_vars.iter().copied().collect();
    let all_vars: HashSet<usize> = (0..t).collect();
    let intermediate_vars: Vec<usize> = all_vars.difference(&input_vars).filter(|v| !output_vars.contains(v)).copied().collect();
    let num_intermediate = intermediate_vars.len();

    // For each output bitstring (basis state)
    for i in 0..dim {
        let mut amp = Complex::new(0.0, 0.0);
        let output_bits: Vec<bool> = (0..n).map(|j| ((i >> (n - 1 - j)) & 1) == 1).collect();
        for y in 0..(1 << num_intermediate) {
            let mut x = vec![false; t];
            for (j, &val) in input_bitstring.iter().enumerate() {
                x[j] = val;
            }
            for (q, &var_idx) in poly.output_vars.iter().enumerate() {
                x[var_idx] = output_bits[q];
            }
            for (k, &var_idx) in intermediate_vars.iter().enumerate() {
                x[var_idx] = ((y >> (num_intermediate - 1 - k)) & 1) == 1;
            }
            // PolyQ: phase = sum_j weight_j * prod(x_{i in term_j}) mod 2
            let mut phase = 0i32;
            for (weight, indices) in &poly.terms {
                let mut v = true;
                for &j in indices {
                    v &= x[j];
                }
                if v {
                    phase += *weight;
                }
            }
            let sign = if ((phase % 2 + 2) % 2) == 0 { 1.0 } else { -1.0 };
            amp += Complex::new(sign, 0.0);
        }
        state[i] = amp;
    }
    // Normalize the statevector to unit norm (PolyQ does this after summing)
    let norm: f64 = state.iter().map(|c| c.norm_sqr()).sum::<f64>().sqrt();
    if norm > 0.0 {
        for v in &mut state {
            *v /= norm;
        }
    }
    state
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::quantum_circuit::Gate;
    use approx::assert_relative_eq;

    #[test]
    fn test_bell_states_phase_polynomial() {
        // Bell state circuit: H(0); H(1); CZ(0,1); H(1)
        let circuit = QuantumCircuit {
            num_qubits: 2,
            gates: vec![
                Gate { name: "h".to_string(), qubits: vec![0], params: vec![] },
                Gate { name: "h".to_string(), qubits: vec![1], params: vec![] },
                Gate { name: "cz".to_string(), qubits: vec![0, 1], params: vec![] },
                Gate { name: "h".to_string(), qubits: vec![1], params: vec![] }
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
    }
}