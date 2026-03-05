pub mod z2 {
    use crate::quantum_circuit::QuantumCircuit;

    #[derive(Debug, Clone)]
    pub struct PhasePolynomial {
        pub terms: Vec<(i32, Vec<usize>)>,
        pub num_qubits: usize,
        pub num_vars: usize,
        pub output_vars: Vec<usize>,
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
                    terms.push((
                        1,
                        vec![
                            wire_array[q][wire_array[q].len() - 2],
                            wire_array[q][wire_array[q].len() - 1],
                        ],
                    ));
                    max_new_var += 1;
                }
                "z" => {
                    let q = gate.qubits[0];
                    terms.push((1, vec![*wire_array[q].last().unwrap()]));
                }
                "cz" => {
                    let q0 = gate.qubits[0];
                    let q1 = gate.qubits[1];
                    terms.push((
                        1,
                        vec![
                            *wire_array[q0].last().unwrap(),
                            *wire_array[q1].last().unwrap(),
                        ],
                    ));
                }
                _ => {}
            }
        }
        let output_vars = wire_array.iter().map(|w| *w.last().unwrap()).collect();
        PhasePolynomial {
            terms,
            num_qubits: n,
            num_vars: max_new_var,
            output_vars,
        }
    }

    use num_complex::Complex;

    pub fn simulate_phase_polynomial(
        poly: &PhasePolynomial,
        input_bitstring: &[bool],
    ) -> Vec<Complex<f64>> {
        use std::collections::HashSet;
        let n = poly.num_qubits;
        let t = poly.num_vars;
        let dim = 1 << n;
        let mut state = vec![Complex::new(0.0, 0.0); dim];

        let input_vars: HashSet<usize> = (0..n).collect();
        let output_vars: HashSet<usize> = poly.output_vars.iter().copied().collect();
        let all_vars: HashSet<usize> = (0..t).collect();
        let intermediate_vars: Vec<usize> = all_vars
            .difference(&input_vars)
            .filter(|v| !output_vars.contains(v))
            .copied()
            .collect();
        let num_intermediate = intermediate_vars.len();

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
        let norm: f64 = state.iter().map(|c| c.norm_sqr()).sum::<f64>().sqrt();
        if norm > 0.0 {
            for v in &mut state {
                *v /= norm;
            }
        }
        state
    }
}

pub mod z4 {
    use crate::quantum_circuit::QuantumCircuit;
    use num_complex::Complex;
    use std::collections::HashMap;

    #[derive(Debug, Clone)]
    pub struct PhasePolynomialZ4 {
        pub terms: Vec<(i32, Vec<usize>)>, 
        pub wire_array: Vec<Vec<usize>>,   
        pub num_qubits: usize,
        pub num_vars: usize, 
    }

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
                    let last = wire_array[q].len();
                    terms.push((2, vec![wire_array[q][last - 2], wire_array[q][last - 1]]));
                    max_new_var += 1;
                }
                "z" => {
                    let q = qubits[0];
                    terms.push((2, vec![*wire_array[q].last().unwrap()]));
                }
                "cz" => {
                    let q0 = qubits[0];
                    let q1 = qubits[1];
                    terms.push((
                        2,
                        vec![
                            *wire_array[q0].last().unwrap(),
                            *wire_array[q1].last().unwrap(),
                        ],
                    ));
                }
                "s" => {
                    let q = qubits[0];
                    terms.push((1, vec![*wire_array[q].last().unwrap()]));
                }
                "sdg" | "sdag" => {
                    let q = qubits[0];
                    terms.push((3, vec![*wire_array[q].last().unwrap()]));
                }
                _ => {}
            }
        }
        PhasePolynomialZ4 { terms, wire_array, num_qubits: n, num_vars: max_new_var }
    }

    pub fn simulate_phase_polynomial_z4(
        poly: &PhasePolynomialZ4,
        input_bitstring: &[bool],
    ) -> Vec<Complex<f64>> {
        let n = poly.num_qubits;
        let t = poly.num_vars;
        let ovs: Vec<usize> = poly.wire_array.iter().map(|w| *w.last().unwrap()).collect();
        let dim = 1 << n;
        let mut s_ldic: HashMap<usize, [u32; 4]> = HashMap::new();
        let x_range = 1 << (t - n);
        let mut x = vec![false; t];
        
        for (i, &b) in input_bitstring.iter().enumerate() {
            x[i] = b;
        }
        
        for i in 0..x_range {
            for (ind, bit) in (0..(t - n))
                .map(|j| ((i >> (t - n - 1 - j)) & 1))
                .enumerate()
            {
                x[n + ind] = bit == 1;
            }
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
            let mut chosenbits = 0usize;
            for (j, &ov) in ovs.iter().rev().enumerate() {
                if x[ov] {
                    chosenbits |= 1 << j;
                }
            }
            let arr = s_ldic.entry(chosenbits).or_insert([0u32; 4]);
            arr[val_out as usize] += 1;
        }
        
        let mut state = vec![Complex::new(0.0, 0.0); dim];
        let norm = (2f64).powf(-0.5 * ((t - n) as f64));
        for (k, arr) in s_ldic.iter() {
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
}

pub mod symplectic {
    use crate::quantum_circuit::QuantumCircuit;
    use num_complex::Complex;
    use std::collections::{HashMap, HashSet};

    pub type Term = (i32, Vec<usize>);

    // Bit-Packed GF(2) Row Vector
    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    pub struct VarMask(pub [u64; 4]);

    impl VarMask {
        #[inline(always)]
        pub fn new() -> Self { VarMask([0; 4]) }

        #[inline(always)]
        pub fn basis(index: usize) -> Self {
            let mut m = [0; 4];
            m[index / 64] = 1 << (index % 64);
            VarMask(m)
        }

        #[inline(always)]
        pub fn get(&self, index: usize) -> u8 {
            ((self.0[index / 64] >> (index % 64)) & 1) as u8
        }

        #[inline(always)]
        pub fn set(&mut self, index: usize) {
            self.0[index / 64] |= 1 << (index % 64);
        }

        #[inline(always)]
        pub fn dot(&self, other: &VarMask) -> u8 {
            let mut p = 0;
            for i in 0..4 { p ^= (self.0[i] & other.0[i]).count_ones(); }
            (p & 1) as u8
        }
    }

    impl std::ops::BitXorAssign for VarMask {
        #[inline(always)]
        fn bitxor_assign(&mut self, rhs: Self) {
            for i in 0..4 { self.0[i] ^= rhs.0[i]; }
        }
    }

    impl std::ops::BitXor for VarMask {
        type Output = Self;
        #[inline(always)]
        fn bitxor(self, rhs: Self) -> Self::Output {
            VarMask([
                self.0[0] ^ rhs.0[0],
                self.0[1] ^ rhs.0[1],
                self.0[2] ^ rhs.0[2],
                self.0[3] ^ rhs.0[3],
            ])
        }
    }

    pub fn circuit_to_z8_polynomial(circuit: &QuantumCircuit) -> (Vec<Term>, Vec<usize>, usize, usize) {
        let n = circuit.num_qubits;
        let mut wirearray: Vec<Vec<usize>> = (0..n).map(|i| vec![i]).collect();
        let mut maxnewvar = n;
        let mut terms: Vec<Term> = Vec::new();
        let mut hcount: usize = 0;
        for gate in &circuit.gates {
            let name = gate.name.to_lowercase();
            match name.as_str() {
                "h" => {
                    hcount += 1;
                    let q = gate.qubits[0];
                    let oldv = *wirearray[q].last().unwrap();
                    wirearray[q].push(maxnewvar);
                    terms.push((4, vec![oldv, maxnewvar]));
                    maxnewvar += 1;
                }
                "z" => {
                    let v = *wirearray[gate.qubits[0]].last().unwrap();
                    terms.push((4, vec![v]));
                }
                "s" => {
                    let v = *wirearray[gate.qubits[0]].last().unwrap();
                    terms.push((2, vec![v]));
                }
                "t" => {
                    let v = *wirearray[gate.qubits[0]].last().unwrap();
                    terms.push((1, vec![v]));
                }
                "cz" => {
                    let a = *wirearray[gate.qubits[0]].last().unwrap();
                    let b = *wirearray[gate.qubits[1]].last().unwrap();
                    terms.push((4, vec![a, b]));
                }
                "ccz" => {
                    let a = *wirearray[gate.qubits[0]].last().unwrap();
                    let b = *wirearray[gate.qubits[1]].last().unwrap();
                    let c = *wirearray[gate.qubits[2]].last().unwrap();
                    terms.push((4, vec![a, b, c]));
                }
                _ => {}
            }
        }
        let y_idx: Vec<usize> = wirearray.iter().map(|w| *w.last().unwrap()).collect();
        (terms, y_idx, maxnewvar, hcount)
    }

    fn mod8(x: i32) -> i32 {
        ((x % 8) + 8) % 8
    }

    fn extract_components_z8(terms: &[Term], m: usize) -> (Vec<Vec<i32>>, Vec<i32>, i32, Vec<Term>) {
        let mut q = vec![vec![0i32; m]; m];
        let mut l = vec![0i32; m];
        let mut eps8 = 0i32;
        let mut cubic: Vec<Term> = Vec::new();
        for (c0, idx) in terms.iter() {
            let c = mod8(*c0);
            if c == 0 { continue; }
            match idx.len() {
                0 => eps8 = (eps8 + c) & 7,
                1 => l[idx[0]] = (l[idx[0]] + c) & 7,
                2 => {
                    let mut a = idx[0];
                    let mut b = idx[1];
                    if a == b {
                        l[a] = (l[a] + c) & 7;
                    } else {
                        if a > b { std::mem::swap(&mut a, &mut b); }
                        q[a][b] = (q[a][b] + c) & 7;
                    }
                }
                3 => {
                    let mut v = idx.clone();
                    v.sort();
                    cubic.push((c, v));
                }
                _ => {}
            }
        }
        (q, l, eps8, cubic)
    }

    pub fn restrict_and_split_to_z4(
        terms: &[Term],
        y_idx: &[usize],
        n_qubits: usize,
        x_bits: usize,
        y_bits: usize,
    ) -> Option<(Vec<VarMask>, Vec<i8>, i32, i32, Vec<Term>, bool, Vec<usize>)> {
        let mut fixed: HashMap<usize, i32> = HashMap::new();
        
        for i in 0..n_qubits {
            fixed.insert(i, ((x_bits >> (n_qubits - 1 - i)) & 1) as i32);
        }
        for q in 0..n_qubits {
            let outv = y_idx[q];
            let want = ((y_bits >> (n_qubits - 1 - q)) & 1) as i32;
            if let Some(&v) = fixed.get(&outv) {
                if v != want { return None; }
            } else {
                fixed.insert(outv, want);
            }
        }
        
        let mut subst: Vec<Term> = Vec::new();
        for (c0, idx) in terms.iter() {
            let mut prod_fixed = 1i32;
            let mut rem: Vec<usize> = Vec::new();
            for &v in idx.iter() {
                if let Some(&val) = fixed.get(&v) {
                    prod_fixed &= val;
                    if prod_fixed == 0 { break; }
                } else {
                    rem.push(v);
                }
            }
            if prod_fixed == 1 {
                subst.push((mod8(*c0), rem));
            }
        }
        
        let mut allu: HashSet<usize> = HashSet::new();
        for (_, ids) in subst.iter() { for &g in ids { allu.insert(g); } }
        let mut u_vars: Vec<usize> = allu.into_iter().collect();
        u_vars.sort();
        let mu = u_vars.len();
        
        let mut vmap: HashMap<usize, usize> = HashMap::new();
        for (i, &g) in u_vars.iter().enumerate() { vmap.insert(g, i); }
        
        let mut renum: Vec<Term> = Vec::new();
        let mut eps8_extra = 0i32;
        for (c, ids) in subst.iter() {
            if ids.is_empty() {
                eps8_extra = (eps8_extra + c) & 7;
            } else {
                let mut r = Vec::new();
                for &v in ids.iter() { r.push(*vmap.get(&v).unwrap()); }
                renum.push((*c, r));
            }
        }
        
        let (q, lz8, eps8, cubic) = extract_components_z8(&renum, mu);
        let eps8 = (eps8 + eps8_extra) & 7;
        let eps_odd = eps8 & 1;
        let eps4 = (eps8 >> 1) & 3;
        
        if mu == 0 {
            return Some((vec![], vec![], eps4, eps_odd, cubic, true, vec![]));
        }
        
        let mut b_rows = vec![VarMask::new(); mu];
        for i in 0..mu {
            for j in (i+1)..mu {
                if q[i][j] % 8 == 4 { 
                    b_rows[i].set(j);
                    b_rows[j].set(i);
                }
            }
        }
        
        let mut l4 = vec![0i8; mu];
        for i in 0..mu { l4[i] = ((lz8[i] & 7) >> 1) as i8; }
        
        let mut odd_linear = Vec::new();
        for i in 0..mu { if (lz8[i] & 1) == 1 { odd_linear.push(i); } }
        Some((b_rows, l4, eps4, eps_odd, cubic, false, odd_linear))
    }

    pub fn symplectic_gram_schmidt(b_rows: &[VarMask]) -> (Vec<VarMask>, usize) {
        let m = b_rows.len();
        if m == 0 { return (vec![], 0); }
        let mut w = vec![VarMask::new(); m];
        for i in 0..m { w[i] = VarMask::basis(i); }
        
        let mut active: Vec<usize> = (0..m).collect();
        let mut t_cols = Vec::with_capacity(m);
        let mut k = 0usize;
        
        loop {
            let mut pair: Option<(usize, usize)> = None;
            for ia in 0..active.len() {
                let a = active[ia];
                let u = &w[a];
                
                for ib in (ia+1)..active.len() {
                    let b = active[ib];
                    let v = &w[b];
                    
                    let mut b_v = VarMask::new();
                    for j in 0..m {
                        if v.get(j) == 1 {
                            b_v ^= b_rows[j];
                        }
                    }
                    if u.dot(&b_v) == 1 { 
                        pair = Some((ia, ib)); 
                        break; 
                    }
                }
                if pair.is_some() { break; }
            }
            if pair.is_none() { break; }
            
            let (ia, ib) = pair.unwrap();
            let a = active[ia];
            let b = active[ib];
            let u = w[a];
            let v = w[b];
            t_cols.push(u);
            t_cols.push(v);
            k += 1;
            
            let mut b_u = VarMask::new();
            let mut b_v = VarMask::new();
            for j in 0..m {
                if u.get(j) == 1 { b_u ^= b_rows[j]; }
                if v.get(j) == 1 { b_v ^= b_rows[j]; }
            }
            
            let mut new_active = Vec::with_capacity(active.len() - 2);
            for &t in active.iter() {
                if t == a || t == b { continue; }
                let mut w_t = w[t];
                
                let cu = w_t.dot(&b_u);
                let cv = w_t.dot(&b_v);
                
                if cu == 1 { w_t ^= v; }
                if cv == 1 { w_t ^= u; }
                
                w[t] = w_t;
                new_active.push(t);
            }
            active = new_active;
        }
        for &t in active.iter() { t_cols.push(w[t]); }
        (t_cols, k)
    }

    fn i_pow(e: i32) -> Complex<f64> {
        let e = ((e % 8) + 8) % 8;
        match e {
            0 => Complex::new(1.0, 0.0),
            1 => Complex::new(2f64.sqrt()/2.0, 2f64.sqrt()/2.0),
            2 => Complex::new(0.0, 1.0),
            3 => Complex::new(-2f64.sqrt()/2.0, 2f64.sqrt()/2.0),
            4 => Complex::new(-1.0, 0.0),
            5 => Complex::new(-2f64.sqrt()/2.0, -2f64.sqrt()/2.0),
            6 => Complex::new(0.0, -1.0),
            7 => Complex::new(2f64.sqrt()/2.0, -2f64.sqrt()/2.0),
            _ => Complex::new(1.0, 0.0),
        }
    }

    fn pair_sum_clifford(alpha4: i32, beta4: i32) -> Complex<f64> {
        let a = i_pow(2 * beta4);
        let b = i_pow(2 * alpha4);
        let c = i_pow(2 * ((alpha4 + beta4) & 3));
        Complex::new(1.0, 0.0) + a + b - c
    }

    fn kernel_factor_clifford(t4: i32) -> Complex<f64> {
        Complex::new(1.0, 0.0) + i_pow(2 * (t4 & 3))
    }

    pub fn amplitude_clifford_t(
        terms: &[Term],
        y_idx: &[usize],
        n_qubits: usize,
        h_count: usize,
        x_bits: usize,
        y_bits: usize,
    ) -> Complex<f64> {
        let out = restrict_and_split_to_z4(terms, y_idx, n_qubits, x_bits, y_bits);
        if out.is_none() { return Complex::new(0.0, 0.0); }
        let (B, L4, eps4, eps_odd, cubic_terms, trivial, odd_linear_vars) = out.unwrap();
        if trivial || B.len() == 0 { 
            return i_pow(2 * eps4 + eps_odd) * Complex::new(2.0f64.powf(-0.5 * h_count as f64), 0.0); 
        }
        
        let mu = B.len();
        let mut cubic_vars: Vec<usize> = Vec::new();
        for (_, idxs) in cubic_terms.iter() { for &v in idxs.iter() { cubic_vars.push(v); } }
        for &v in odd_linear_vars.iter() { cubic_vars.push(v); }
        cubic_vars.sort(); cubic_vars.dedup();
        let t = cubic_vars.len();
        
        if t == 0 {
            let (t_cols, k) = symplectic_gram_schmidt(&B);
            let mut lp = vec![0i32; mu];
            for j in 0..mu {
                let mut sum = 0i32;
                for i in 0..mu {
                    if t_cols[j].get(i) == 1 { sum += L4[i] as i32; }
                }
                lp[j] = (sum % 4 + 4) % 4;
            }
            let mut s = i_pow(2 * eps4 + eps_odd);
            for p in 0..k { s *= pair_sum_clifford(lp[2*p], lp[2*p+1]); }
            for j in (2*k)..mu { s *= kernel_factor_clifford(lp[j]); }
            return s * Complex::new(2.0f64.powf(-0.5 * h_count as f64), 0.0);
        }
        
        let mut total = Complex::new(0.0, 0.0);
        let mut cubic_vars_sorted = cubic_vars.clone();
        cubic_vars_sorted.sort();
        
        let slices = 1usize << t;
        for assign in 0..slices {
            let mut assign_map: HashSet<usize> = HashSet::new();
            for i in 0..t { 
                if ((assign >> i) & 1) == 1 { assign_map.insert(cubic_vars_sorted[i]); } 
            }
            
            let mut subst: Vec<Term> = Vec::new();
            for (c0, idx) in terms.iter() {
                let mut prod_fixed = 1i32;
                let mut rem: Vec<usize> = Vec::new();
                for &v in idx.iter() {
                    if assign_map.contains(&v) { 
                        prod_fixed &= 1;
                    } else if cubic_vars_sorted.binary_search(&v).is_ok() {
                        prod_fixed = 0; 
                        break;
                    } else { 
                        rem.push(v); 
                    }
                }
                if prod_fixed == 1 { subst.push((mod8(*c0), rem)); }
            }
            
            let new_out = restrict_and_split_to_z4(&subst, y_idx, n_qubits, x_bits, y_bits);
            if new_out.is_none() { continue; }
            let (newB, newL4, neweps4, neweps_odd, newcubic, newtriv, _) = new_out.unwrap();
            if !newcubic.is_empty() { continue; }
            
            let slice_val = if newtriv || newB.len() == 0 {
                i_pow(2 * neweps4 + neweps_odd)
            } else {
                let (t_cols, k) = symplectic_gram_schmidt(&newB);
                let mut lp = vec![0i32; newB.len()];
                for j in 0..newB.len() {
                    let mut sum = 0i32;
                    for i in 0..newB.len() { if t_cols[j].get(i) == 1 { sum += newL4[i] as i32; } }
                    lp[j] = (sum % 4 + 4) % 4;
                }
                let mut s = i_pow(2 * neweps4 + neweps_odd);
                for p in 0..k { s *= pair_sum_clifford(lp[2*p], lp[2*p+1]); }
                for j in (2*k)..newB.len() { s *= kernel_factor_clifford(lp[j]); }
                s
            };
            total += slice_val;
        }
        total * Complex::new(2.0f64.powf(-0.5 * h_count as f64), 0.0)
    }

    pub fn simulate_statevector_z8(circuit: &QuantumCircuit, x_bits: usize) -> Vec<Complex<f64>> {
        let (terms, y_idx, _, hcount) = circuit_to_z8_polynomial(circuit);
        let n = circuit.num_qubits;
        let dim = 1usize << n;
        let mut psi = vec![Complex::new(0.0, 0.0); dim];
        for y in 0..dim {
            psi[y] = amplitude_clifford_t(&terms, &y_idx, n, hcount, x_bits, y);
        }
        psi
    }
}

// #[cfg(test)]
// mod bench_tests {
//     use super::{symplectic, z4};
//     use crate::quantum_circuit::Gate;
//     use std::time::Instant;
//     use num_complex::Complex;

//     fn make_test_circuit(n: usize, depth: usize) -> crate::quantum_circuit::QuantumCircuit {
//         let mut gates = Vec::new();
//         // Add Hadamards to create intermediate variables (superposition)
//         for q in 0..n {
//             gates.push(Gate { name: "h".to_string(), qubits: vec![q], params: vec![] });
//         }
//         // Add entangling gates
//         for _ in 0..depth {
//             for q in 0..(n.saturating_sub(1)) {
//                 gates.push(Gate { name: "cz".to_string(), qubits: vec![q, q + 1], params: vec![] });
//             }
//         }
//         crate::quantum_circuit::QuantumCircuit { num_qubits: n, gates }
//     }

//     /// A lean Z4 solver strictly for a single amplitude (Standard Sum-Over-Paths)
//     /// Time Complexity: O(2^m) where m is the number of intermediate variables
//     fn z4_single_amplitude(
//         poly: &z4::PhasePolynomialZ4,
//         target_y: usize,
//     ) -> Complex<f64> {
//         let n = poly.num_qubits;
//         let t = poly.num_vars;
//         let m = t - n; // Number of intermediate variables
//         let ovs: Vec<usize> = poly.wire_array.iter().map(|w| *w.last().unwrap()).collect();
        
//         let mut sum_amp = Complex::new(0.0, 0.0);
//         let x_range = 1 << m;
//         let mut x = vec![false; t];
        
//         // Input variables are fixed to 0 (since input is |00...0>)
//         for i in 0..n { x[i] = false; }

//         for i in 0..x_range {
//             // 1. Assign intermediate variables
//             for j in 0..m {
//                 x[n + j] = ((i >> (m - 1 - j)) & 1) == 1;
//             }

//             // 2. Check if this path naturally leads to our target output state
//             let mut chosenbits = 0usize;
//             for (j, &ov) in ovs.iter().rev().enumerate() {
//                 if x[ov] { chosenbits |= 1 << j; }
//             }
            
//             // If it doesn't map to our target state, discard this path
//             if chosenbits != target_y { continue; }

//             // 3. Evaluate the Z4 phase polynomial for this path
//             let mut val_out: u8 = 0;
//             for (weight, idxs) in &poly.terms {
//                 let mut v = true;
//                 for &idx in idxs { v &= x[idx]; }
//                 if v { val_out = (val_out + (*weight as u8)) % 4; }
//             }

//             let theta = (std::f64::consts::PI / 2.0) * (val_out as f64);
//             sum_amp += Complex::from_polar(1.0, theta);
//         }
        
//         sum_amp * (2.0f64).powf(-0.5 * m as f64)
//     }

//     #[test]
//     fn bench_single_amplitude_comparison() {
//         println!("\n--- Single Amplitude Comparison (|00...0>) ---");
//         println!("{:<5} | {:<15} | {:<15}", "n", "Z4 (Sum-over-Paths)", "Symplectic (GF2)");
//         println!("------------------------------------------------");

//         // We run up to n=20 to clearly show the crossover point where O(m^3) beats O(2^m)
//         for n in 10..=20 {
//             let circuit = make_test_circuit(n, 1);
            
//             // --- Z4 Evaluation ---
//             // We cap Z4 at n=20 because the O(2^m) scaling makes it too slow to benchmark efficiently above that
//             let dur_z4 = if n <= 20 {
//                 let poly = z4::phase_polynomial_z4(&circuit);
//                 let start = Instant::now();
//                 let _amp = z4_single_amplitude(&poly, 0); // target = |00...0>
//                 format!("{:?}", start.elapsed())
//             } else {
//                 "DNT (Too Slow)".to_string()
//             };

//             // --- Symplectic Evaluation ---
//             // Extract the polynomial and run the GF(2) Gram-Schmidt solver
//             let (terms, y_idx, _, hcount) = symplectic::circuit_to_z8_polynomial(&circuit);
//             let start = Instant::now();
//             let _amp = symplectic::amplitude_clifford_t(&terms, &y_idx, n, hcount, 0, 0);
//             let dur_sym = format!("{:?}", start.elapsed());

//             println!("{:<5} | {:<15} | {:<15}", n, dur_z4, dur_sym);
//         }
//         println!("------------------------------------------------\n");
//     }
// }

#[cfg(test)]
mod bench_tests {
    use super::{symplectic, z4};
    use crate::quantum_circuit::Gate;
    use std::time::Instant;
    use num_complex::Complex;

    fn make_test_circuit(n: usize, depth: usize) -> crate::quantum_circuit::QuantumCircuit {
        let mut gates = Vec::new();
        // Add Hadamards to create intermediate variables (superposition)
        for q in 0..n {
            gates.push(Gate { name: "h".to_string(), qubits: vec![q], params: vec![] });
        }
        // Add entangling gates
        for _ in 0..depth {
            for q in 0..(n.saturating_sub(1)) {
                gates.push(Gate { name: "cz".to_string(), qubits: vec![q, q + 1], params: vec![] });
            }
        }
        crate::quantum_circuit::QuantumCircuit { num_qubits: n, gates }
    }

    /// A lean Z4 solver strictly for a single amplitude (Standard Sum-Over-Paths)
    fn z4_single_amplitude(
        poly: &z4::PhasePolynomialZ4,
        target_y: usize,
    ) -> Complex<f64> {
        let n = poly.num_qubits;
        let t = poly.num_vars;
        let m = t - n; 
        let ovs: Vec<usize> = poly.wire_array.iter().map(|w| *w.last().unwrap()).collect();
        
        let mut sum_amp = Complex::new(0.0, 0.0);
        let x_range = 1 << m;
        let mut x = vec![false; t];
        
        for i in 0..n { x[i] = false; }

        for i in 0..x_range {
            for j in 0..m { x[n + j] = ((i >> (m - 1 - j)) & 1) == 1; }

            let mut chosenbits = 0usize;
            for (j, &ov) in ovs.iter().rev().enumerate() {
                if x[ov] { chosenbits |= 1 << j; }
            }
            
            if chosenbits != target_y { continue; }

            let mut val_out: u8 = 0;
            for (weight, idxs) in &poly.terms {
                let mut v = true;
                for &idx in idxs { v &= x[idx]; }
                if v { val_out = (val_out + (*weight as u8)) % 4; }
            }

            let theta = (std::f64::consts::PI / 2.0) * (val_out as f64);
            sum_amp += Complex::from_polar(1.0, theta);
        }
        
        sum_amp * (2.0f64).powf(-0.5 * m as f64)
    }

    #[test]
    fn test_single_amplitude_correctness() {
        println!("\n--- Single Amplitude Correctness Check (|00...0>) ---");
        println!("{:<4} | {:<22} | {:<12} | {:<12} | Match?", "n", "Calculated Amp", "Z4 Time", "Sym Time");
        println!("------------------------------------------------------------------------");

        // Test up to n=12 to keep Z4 execution time reasonable during tests
        for n in 2..=12 {
            let circuit = make_test_circuit(n, 1);
            
            // --- Z4 Evaluation ---
            let poly = z4::phase_polynomial_z4(&circuit);
            let start_z4 = Instant::now();
            let amp_z4 = z4_single_amplitude(&poly, 0); 
            let dur_z4 = format!("{:?}", start_z4.elapsed());

            // --- Symplectic Evaluation ---
            let (terms, y_idx, _, hcount) = symplectic::circuit_to_z8_polynomial(&circuit);
            let start_sym = Instant::now();
            let amp_sym = symplectic::amplitude_clifford_t(&terms, &y_idx, n, hcount, 0, 0);
            let dur_sym = format!("{:?}", start_sym.elapsed());

            // --- Correctness Assertion ---
            // Floating point math requires a small epsilon for equality checks
            let is_match = (amp_z4 - amp_sym).norm() < 1e-8;
            
            assert!(
                is_match, 
                "Amplitude mismatch at n={}! Z4: {:?}, Sym: {:?}", 
                n, amp_z4, amp_sym
            );

            let amp_str = format!("{:.4} + {:.4}i", amp_sym.re, amp_sym.im);
            println!("{:<4} | {:<22} | {:<12} | {:<12} | {}", n, amp_str, dur_z4, dur_sym, if is_match { "Yes" } else { "NO" });
        }
        println!("------------------------------------------------------------------------\n");
    }
}