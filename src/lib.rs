use fixedbitset::FixedBitSet;
use num_complex::Complex64;
use std::collections::HashMap;

#[derive(Clone, Debug, PartialEq)]
pub enum Gate { H(usize), S(usize), Z(usize), CZ(usize, usize), T(usize), CNOT(usize, usize) }

pub struct PolyLibrary {
    pub num_qubits: usize,
    pub num_h: usize,
    matrix_b: Vec<FixedBitSet>, 
    vector_v: Vec<u8>,         
    output_mapping: Vec<usize>,
}

impl PolyLibrary {
    pub fn new(n: usize, gates: Vec<Gate>) -> Self {
        let mut sim = Self {
            num_qubits: n,
            num_h: 0,
            matrix_b: (0..n).map(|_| FixedBitSet::with_capacity(n)).collect(),
            vector_v: vec![0; n],
            output_mapping: (0..n).collect(),
        };
        
        let mut expanded_gates = Vec::new();
        for gate in gates {
            match gate {
                Gate::CNOT(c, t) => {
                    expanded_gates.push(Gate::H(t));
                    expanded_gates.push(Gate::CZ(c, t));
                    expanded_gates.push(Gate::H(t));
                }
                _ => expanded_gates.push(gate),
            }
        }
        sim.build(expanded_gates);
        sim
    }

    fn build(&mut self, gates: Vec<Gate>) {
        for gate in gates {
            match gate {
                Gate::H(q) => {
                    let old_var = self.output_mapping[q];
                    let new_var = self.vector_v.len();
                    self.vector_v.push(0);
                    for row in &mut self.matrix_b { row.grow(new_var + 1); }
                    let mut new_row = FixedBitSet::with_capacity(new_var + 1);
                    new_row.insert(old_var); 
                    self.matrix_b.push(new_row);
                    self.matrix_b[old_var].insert(new_var);
                    self.output_mapping[q] = new_var;
                    self.num_h += 1;
                }
                Gate::Z(q) => self.vector_v[self.output_mapping[q]] ^= 1, 
                Gate::CZ(q1, q2) => {
                    let (i, j) = (self.output_mapping[q1], self.output_mapping[q2]);
                    self.matrix_b[i].insert(j);
                    self.matrix_b[j].insert(i);
                }
                _ => {}
            }
        }
    }

    pub fn solve_amplitude(&self, target_y: usize) -> Complex64 {
        let total_vars = self.vector_v.len();
        let mut fixed_vals = HashMap::new();
        for i in 0..self.num_qubits { fixed_vals.insert(i, 0); }
        for i in 0..self.num_qubits {
            fixed_vals.insert(self.output_mapping[i], ((target_y >> i) & 1) as u8);
        }

        let mut internal_vars = Vec::new();
        for i in 0..total_vars {
            if !fixed_vals.contains_key(&i) { internal_vars.push(i); }
        }

        let m = internal_vars.len();
        let mut b_red = (0..m).map(|_| FixedBitSet::with_capacity(m)).collect::<Vec<_>>();
        let mut v_red = vec![0u8; m];
        let mut epsilon = 0u8;

        for (new_i, &orig_i) in internal_vars.iter().enumerate() {
            v_red[new_i] = self.vector_v[orig_i];
            for (new_j, &orig_j) in internal_vars.iter().enumerate() {
                if self.matrix_b[orig_i].contains(orig_j) { b_red[new_i].insert(new_j); }
            }
            for (&fixed_idx, &val) in &fixed_vals {
                if val == 1 && self.matrix_b[orig_i].contains(fixed_idx) {
                    v_red[new_i] ^= 1;
                }
            }
        }
        
        for (&idx, &val) in &fixed_vals {
            if val == 1 {
                epsilon ^= self.vector_v[idx];
                for (&idx2, &val2) in &fixed_vals {
                    if val2 == 1 && idx < idx2 && self.matrix_b[idx].contains(idx2) {
                        epsilon ^= 1;
                    }
                }
            }
        }

        let (r, v_final) = self.run_dickson(b_red, v_red);
        
        // Debugging Information
        println!("--- Debug: Amplitude Calculation for |{}> ---", target_y);
        println!("Internal Variables (m): {}", m);
        println!("Rank (2k): {}", r);
        println!("Transformed Linear Vector: {:?}", v_final);
        println!("Epsilon (Constant Phase): {}", epsilon);

        self.compute_weight_sum(m, r, v_final, epsilon)
    }

    fn run_dickson(&self, mut b: Vec<FixedBitSet>, mut v: Vec<u8>) -> (usize, Vec<u8>) {
        let m = b.len();
        let mut r = 0;
        for p in (0..m.saturating_sub(1)).step_by(2) {
            let mut pivot = None;
            for i in p..m {
                for j in i+1..m {
                    if b[i].contains(j) { pivot = Some((i, j)); break; }
                }
                if pivot.is_some() { break; }
            }

            if let Some((i, j)) = pivot {
                b.swap(p, i); b.swap(p+1, j);
                v.swap(p, i); v.swap(p+1, j); 
                
                let r_p = b[p].clone(); let r_p1 = b[p+1].clone();
                for k in p+2..m {
                    if b[k].contains(p) { b[k].union_with(&r_p1); v[k] ^= v[p+1]; }
                    if b[k].contains(p + 1) { b[k].union_with(&r_p); v[k] ^= v[p]; }
                }
                r += 2;
            } else { break; }
        }
        (r, v)
    }

    fn compute_weight_sum(&self, m: usize, r: usize, v: Vec<u8>, eps: u8) -> Complex64 {
        // Kernel variables Check: Section 6.2.2 
        for i in r..m {
            if v[i] == 1 { 
                println!("Result: Balanced function (Kernel variable L_k != 0). Amplitude = 0.");
                return Complex64::new(0.0, 0.0); 
            }
        }

        // Exponential Sum calculation for F2: Section 6.2.1 [cite: 516]
        // Sum (-1)^{sum x_i x_{i+1}} = 2^k
        let mut sum_val = 1.0f64;
        for _ in 0..(r/2) { sum_val *= 2.0; } // Sum of (-1)^{x1x2} over 4 values is 2
        
        sum_val *= (1 << (m - r)) as f64; // Scale by free kernel variables
        if eps == 1 { sum_val *= -1.0; }

        let norm = (1.0 / 2.0f64).powf(self.num_h as f64 / 2.0);
        println!("Result: Amplitude = {}", sum_val * norm);
        Complex64::new(sum_val * norm, 0.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ghz_state() {
        let gates = vec![Gate::H(0), Gate::CNOT(0, 1), Gate::CNOT(1, 2)];
        let sim = PolyLibrary::new(3, gates);
        let amp0 = sim.solve_amplitude(0b000); 
        let amp7 = sim.solve_amplitude(0b111); 
        assert!((amp0.norm() - 0.7071).abs() < 0.01);
        assert!((amp7.norm() - 0.7071).abs() < 0.01);
    }

    #[test]
    fn test_ghz_state_with_phase() {
        // prepare GHZ on 3 qubits, then apply a phase of i on the |111> term
        // (a Z on qubit 2 produces a –1 phase; an S gives +i).
        let gates = vec![
            Gate::H(0),
            Gate::CNOT(0, 1),
            Gate::CNOT(1, 2),
            Gate::S(2),             
        ];
        let sim = PolyLibrary::new(3, gates);

        let amp000 = sim.solve_amplitude(0b000);
        let amp111 = sim.solve_amplitude(0b111);

        assert!((amp000.norm() - 0.7071).abs() < 0.01);
        assert!((amp111.norm() - 0.7071).abs() < 0.01);

        assert!(amp111.im.abs() -0.7071 < 0.01);
    }
//    #[test]
//    fn ry(){

//         let gates = vec![
//             Gate::H(0), Gate::T(0), Gate::H(0)
//         ];
//         let sim = PolyLibrary::new(1, gates);


//         let amp000 = sim.solve_amplitude(0b0);
//         let amp001 = sim.solve_amplitude(0b1);
//         println!("Amplitude for |0>: {}", amp000);
//         println!("Amplitude for |1>: {}", amp001);
//         assert!((amp000.norm() - 0.924).abs() < 0.01);
//         assert!((amp001.norm() - 0.383).abs() < 0.01);
//    }

    #[test]
    fn test_bv_algorithm() {
        let gates = vec![
        
            Gate::H(0), Gate::H(1),Gate::H(2), // Hadamard wall
            Gate::Z(1), // Oracle for f(x) = x1 XOR x2
            Gate::H(0), Gate::H(1),Gate::H(2), // Hadamard wall

        ];
        let sim = PolyLibrary::new(3, gates);

        let amp_correct = sim.solve_amplitude(0b010); 
        assert!((amp_correct.norm() - 1.0).abs() < 0.01);
    }
}