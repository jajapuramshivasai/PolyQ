// Root cause (why norm_sq became 2):
//
// In `amplitude_clifford_t_accel` we "fixed" BOTH:
//   - input vars 0..n-1
//   - output vars for each qubit (poly.output_vars[i])
//
// But if a qubit has *no Hadamard*, then output_vars[i] == i.
// That means we accidentally forced the *same* variable to be both input and output,
// which is correct physically (output must equal input), BUT we enforced it by simply
// overwriting in the HashMap:
//
//   fixed.insert(i, input[i]);
//   fixed.insert(i, output_bit);
//
// If those disagree, we should return amplitude 0.
// If they agree, fine.
//
// In the failing 3-qubit test, qubit 1 likely had output_var==1 (no H on that wire),
// so outputs varied over y but we never enforced the constraint "y_bit == input_bit"
// as a delta; we overwrote fixed values in a way that effectively doubled total
// probability mass across outputs, producing norm_sq=2.
//
// Fix: when inserting into `fixed`, check for conflicts.
// If conflict: amplitude = 0 immediately.
// Additionally: if output_var is same as input var, it should not increase the number
// of free variables, and the algorithm remains correct once conflicts are handled.
//
// Below is a minimal patch: add `insert_fixed_checked` and use it when populating fixed.

use fixedbitset::FixedBitSet;
use num_complex::Complex64;
use std::collections::{HashMap, HashSet};

#[derive(Clone, Debug, PartialEq)]
pub enum Gate {
    H(usize),
    Z(usize),
    S(usize),
    T(usize),
    CZ(usize, usize),
}

#[derive(Clone, Debug)]
pub struct Z8Term {
    pub weight: u8,
    pub vars: Vec<usize>,
}

#[derive(Clone, Debug)]
pub struct CompiledPhasePoly {
    pub num_qubits: usize,
    pub num_vars: usize,
    pub num_h: usize,
    pub output_vars: Vec<usize>,
    pub b4: Vec<FixedBitSet>,
    pub v4: Vec<u8>,
    pub eps4: u8,
    pub rem: Vec<Z8Term>,
}

pub fn compile_clifford_t(num_qubits: usize, gates: &[Gate]) -> CompiledPhasePoly {
    let n = num_qubits;
    let mut wire: Vec<Vec<usize>> = (0..n).map(|i| vec![i]).collect();
    let mut next_var = n;
    let mut num_h = 0usize;

    let mut b4: Vec<FixedBitSet> = (0..n).map(|_| FixedBitSet::with_capacity(n)).collect();
    let mut v4: Vec<u8> = vec![0; n];
    let eps4 = 0u8;

    let mut rem: Vec<Z8Term> = Vec::new();

    let mut grow_to =
        |new_t: usize, b4: &mut Vec<FixedBitSet>, v4: &mut Vec<u8>| {
            while v4.len() < new_t {
                v4.push(0);
            }
            while b4.len() < new_t {
                b4.push(FixedBitSet::with_capacity(new_t));
            }
            for row in b4.iter_mut() {
                row.grow(new_t);
            }
        };

    for g in gates {
        match *g {
            Gate::H(q) => {
                let prev = *wire[q].last().unwrap();
                let cur = next_var;
                next_var += 1;
                num_h += 1;

                wire[q].push(cur);
                grow_to(next_var, &mut b4, &mut v4);

                b4[prev].insert(cur);
                b4[cur].insert(prev);
            }
            Gate::CZ(a, b) => {
                let va = *wire[a].last().unwrap();
                let vb = *wire[b].last().unwrap();
                grow_to(next_var, &mut b4, &mut v4);
                b4[va].insert(vb);
                b4[vb].insert(va);
            }
            Gate::Z(q) => {
                let v = *wire[q].last().unwrap();
                grow_to(next_var, &mut b4, &mut v4);
                v4[v] ^= 1;
            }
            Gate::S(q) => {
                let v = *wire[q].last().unwrap();
                rem.push(Z8Term {
                    weight: 2,
                    vars: vec![v],
                });
            }
            Gate::T(q) => {
                let v = *wire[q].last().unwrap();
                rem.push(Z8Term {
                    weight: 1,
                    vars: vec![v],
                });
            }
        }
    }

    let output_vars = wire.iter().map(|w| *w.last().unwrap()).collect::<Vec<_>>();

    CompiledPhasePoly {
        num_qubits: n,
        num_vars: next_var,
        num_h,
        output_vars,
        b4,
        v4,
        eps4,
        rem,
    }
}

fn eval_rem_mod8(rem: &[Z8Term], x: &[u8]) -> u8 {
    let mut acc = 0u8;
    for t in rem {
        let mut v = 1u8;
        for &idx in &t.vars {
            v &= x[idx];
            if v == 0 {
                break;
            }
        }
        if v == 1 {
            acc = (acc + (t.weight % 8)) % 8;
        }
    }
    acc
}

fn dickson_reduce(mut b: Vec<FixedBitSet>, mut v: Vec<u8>) -> (usize, Vec<u8>) {
    let m = b.len();
    let mut r = 0usize;

    let mut p = 0usize;
    while p + 1 < m {
        let mut pivot: Option<(usize, usize)> = None;
        'outer: for i in p..m {
            for j in (i + 1)..m {
                if b[i].contains(j) {
                    pivot = Some((i, j));
                    break 'outer;
                }
            }
        }

        let Some((i, j)) = pivot else {
            break;
        };

        b.swap(p, i);
        v.swap(p, i);

        let mut j2 = j;
        if j2 == p {
            j2 = i;
        }
        b.swap(p + 1, j2);
        v.swap(p + 1, j2);

        let rp = b[p].clone();
        let rp1 = b[p + 1].clone();

        for k in (p + 2)..m {
            if b[k].contains(p) {
                b[k].symmetric_difference_with(&rp1);
                v[k] ^= v[p + 1];
            }
            if b[k].contains(p + 1) {
                b[k].symmetric_difference_with(&rp);
                v[k] ^= v[p];
            }
        }

        r += 2;
        p += 2;
    }

    (r, v)
}

fn z2_quadratic_exponential_sum(b: Vec<FixedBitSet>, v: Vec<u8>, eps: u8) -> f64 {
    let m = v.len();
    let (r, v2) = dickson_reduce(b, v);

    for i in r..m {
        if v2[i] == 1 {
            return 0.0;
        }
    }

    let mag = (1u128 << (m - r / 2)) as f64;
    if eps == 1 { -mag } else { mag }
}

fn insert_fixed_checked(fixed: &mut HashMap<usize, u8>, idx: usize, val: u8) -> bool {
    match fixed.get(&idx) {
        None => {
            fixed.insert(idx, val & 1);
            true
        }
        Some(&old) => old == (val & 1),
    }
}

pub fn amplitude_clifford_t_accel(
    poly: &CompiledPhasePoly,
    input: &[u8],
    target_y: usize,
) -> Complex64 {
    let n = poly.num_qubits;
    let t = poly.num_vars;
    assert_eq!(input.len(), n);

    let mut fixed: HashMap<usize, u8> = HashMap::new();

    // Insert input constraints
    for i in 0..n {
        if !insert_fixed_checked(&mut fixed, i, input[i]) {
            return Complex64::new(0.0, 0.0);
        }
    }
    // Insert output constraints (may collide with input vars if output_vars[i]==i)
    for i in 0..n {
        let ov = poly.output_vars[i];
        let bit = ((target_y >> i) & 1) as u8;
        if !insert_fixed_checked(&mut fixed, ov, bit) {
            return Complex64::new(0.0, 0.0);
        }
    }

    let mut vset: HashSet<usize> = HashSet::new();
    for term in &poly.rem {
        for &idx in &term.vars {
            if !fixed.contains_key(&idx) {
                vset.insert(idx);
            }
        }
    }

    let internal: Vec<usize> = (0..t).filter(|i| !fixed.contains_key(i)).collect();
    let vvars: Vec<usize> = internal.iter().cloned().filter(|i| vset.contains(i)).collect();
    let uvars: Vec<usize> = internal.iter().cloned().filter(|i| !vset.contains(i)).collect();

    let nv = vvars.len();
    let nu = uvars.len();

    let mut x_full = vec![0u8; t];
    for (&k, &val) in fixed.iter() {
        x_full[k] = val;
    }

    let mut eps_base = poly.eps4 & 1;

    for (&idx, &val) in fixed.iter() {
        if val == 1 && (poly.v4[idx] & 1) == 1 {
            eps_base ^= 1;
        }
    }
    let fixed_keys: Vec<usize> = fixed.keys().cloned().collect();
    for a in 0..fixed_keys.len() {
        for b in (a + 1)..fixed_keys.len() {
            let i = fixed_keys[a];
            let j = fixed_keys[b];
            if fixed[&i] == 1 && fixed[&j] == 1 && poly.b4[i].contains(j) {
                eps_base ^= 1;
            }
        }
    }

    let mut bu: Vec<FixedBitSet> = (0..nu).map(|_| FixedBitSet::with_capacity(nu)).collect();
    let mut vu_base: Vec<u8> = vec![0u8; nu];
    for (ui, &orig_u) in uvars.iter().enumerate() {
        vu_base[ui] = poly.v4[orig_u] & 1;
    }

    for (ui, &orig_u) in uvars.iter().enumerate() {
        for (uj, &orig_uj) in uvars.iter().enumerate() {
            if ui < uj && poly.b4[orig_u].contains(orig_uj) {
                bu[ui].insert(uj);
                bu[uj].insert(ui);
            }
        }
        for (&fidx, &fval) in fixed.iter() {
            if fval == 1 && poly.b4[orig_u].contains(fidx) {
                vu_base[ui] ^= 1;
            }
        }
    }

    let mut cross: Vec<Vec<usize>> = vec![Vec::new(); nu];
    for (ui, &orig_u) in uvars.iter().enumerate() {
        for (vj_pos, &orig_v) in vvars.iter().enumerate() {
            if poly.b4[orig_u].contains(orig_v) {
                cross[ui].push(vj_pos);
            }
        }
    }

    let mut amp = Complex64::new(0.0, 0.0);
    let total_v = 1usize << nv;

    for mask in 0..total_v {
        for (j, &orig_v) in vvars.iter().enumerate() {
            x_full[orig_v] = ((mask >> j) & 1) as u8;
        }

        let r = eval_rem_mod8(&poly.rem, &x_full);
        let theta = std::f64::consts::PI * (r as f64) / 4.0;
        let phase = Complex64::from_polar(1.0, theta);

        let mut vu = vu_base.clone();
        for ui in 0..nu {
            let mut togg = 0u8;
            for &vj_pos in &cross[ui] {
                togg ^= ((mask >> vj_pos) & 1) as u8;
            }
            vu[ui] ^= togg;
        }

        let mut eps = eps_base;

        for (j, &orig_v) in vvars.iter().enumerate() {
            if ((mask >> j) & 1) == 1 && (poly.v4[orig_v] & 1) == 1 {
                eps ^= 1;
            }
        }
        for (&fidx, &fval) in fixed.iter() {
            if fval == 1 {
                for (j, &orig_v) in vvars.iter().enumerate() {
                    if ((mask >> j) & 1) == 1 && poly.b4[fidx].contains(orig_v) {
                        eps ^= 1;
                    }
                }
            }
        }
        for a in 0..nv {
            if ((mask >> a) & 1) == 0 {
                continue;
            }
            for b in (a + 1)..nv {
                if ((mask >> b) & 1) == 0 {
                    continue;
                }
                if poly.b4[vvars[a]].contains(vvars[b]) {
                    eps ^= 1;
                }
            }
        }

        let inner = z2_quadratic_exponential_sum(bu.clone(), vu, eps);
        amp += phase * inner;
    }

    let norm = (2f64).powf(-(poly.num_h as f64) / 2.0);
    amp * norm
}

#[cfg(test)]
mod tests {
    use rayon::vec;

    use super::*;

    fn c(re: f64, im: f64) -> Complex64 {
        Complex64::new(re, im)
    }

    fn approx_eq(a: Complex64, b: Complex64, eps: f64) {
        assert!(
            (a - b).norm() <= eps,
            "a={:?} b={:?} |a-b|={}",
            a,
            b,
            (a - b).norm()
        );
    }

    #[test]
    fn test_bell_state_clifford_amplitudes() {
        let gates = vec![Gate::H(0), Gate::H(1), Gate::CZ(0, 1), Gate::H(1)];
        let poly = compile_clifford_t(2, &gates);

        let input = vec![0u8, 0u8];
        let a00 = amplitude_clifford_t_accel(&poly, &input, 0);
        let a11 = amplitude_clifford_t_accel(&poly, &input, 3);
        let a01 = amplitude_clifford_t_accel(&poly, &input, 1);
        let a10 = amplitude_clifford_t_accel(&poly, &input, 2);

        let norm = (0.5f64).sqrt();
        approx_eq(a00, c(norm, 0.0), 1e-10);
        approx_eq(a11, c(norm, 0.0), 1e-10);
        approx_eq(a01, c(0.0, 0.0), 1e-10);
        approx_eq(a10, c(0.0, 0.0), 1e-10);
    }

    #[test]
    fn test_clifford_plus_s_phase() {
        let gates = vec![
            Gate::H(0),
            Gate::S(0),
            Gate::H(1),
            Gate::CZ(0, 1),
            Gate::H(1),
        ];
        let poly = compile_clifford_t(2, &gates);

        let input = vec![0u8, 0u8];
        let a00 = amplitude_clifford_t_accel(&poly, &input, 0);
        let a11 = amplitude_clifford_t_accel(&poly, &input, 3);

        let norm = (0.5f64).sqrt();
        approx_eq(a00, c(norm, 0.0), 1e-10);
        approx_eq(a11, c(0.0, norm), 1e-10);
    }

    #[test]
    fn test_clifford_plus_t_phase() {
        let gates = vec![
            Gate::H(0),
            Gate::T(0),
            Gate::H(1),
            Gate::CZ(0, 1),
            Gate::H(1),
        ];
        let poly = compile_clifford_t(2, &gates);

        let input = vec![0u8, 0u8];
        let a00 = amplitude_clifford_t_accel(&poly, &input, 0);
        let a11 = amplitude_clifford_t_accel(&poly, &input, 3);

        approx_eq(a00, c((0.5f64).sqrt(), 0.0), 1e-10);
        approx_eq(a11, c(0.5, 0.5), 1e-10);
    }

    #[test]
    fn test_three_qubit_mixed_norm_sanity() {
        let gates = vec![
            Gate::H(0),
            Gate::CZ(0, 1),
            Gate::T(1),
            Gate::H(2),
            Gate::CZ(1, 2),
            Gate::S(0),
            Gate::Z(2),
        ];
        let poly = compile_clifford_t(3, &gates);
        let input = vec![0u8, 0u8, 0u8];

        let mut norm_sq = 0.0f64;
        for y in 0..(1usize << 3) {
            let a = amplitude_clifford_t_accel(&poly, &input, y);
            norm_sq += a.norm_sqr();
        }
        assert!(
            (norm_sq - 1.0).abs() < 1e-8,
            "norm_sq={}",
            norm_sq
        );
    }

    #[test]
    fn ry(){
        let gates = vec![Gate::H(0), Gate::T(0), Gate::H(0)];
        let poly = compile_clifford_t(1, &gates);
        let input = vec![0u8];
        let a0 = amplitude_clifford_t_accel(&poly, &input, 0);
        let a1 = amplitude_clifford_t_accel(&poly, &input, 1);

        //composere result [ 0.854+0.354j, 0.354+0.146j ]
        println!("a0={:?} a1={:?}", a0, a1);
        approx_eq(a0, c(0.8535533905932737, 0.35355339059327373), 1e-3);
        approx_eq(a1, c(0.1464466094067262, -0.35355339059327373), 1e-3);

        
    }

    #[test]
    fn benchmark_like_clifford_vs_nonclifford_amplitude() {
        use std::time::Instant;

        let n = 12;

        let mut g1 = Vec::new();
        for q in 0..n {
            g1.push(Gate::H(q));
        }
        for q in 0..(n - 1) {
            g1.push(Gate::CZ(q, q + 1));
        }
        for q in (0..n).rev() {
            g1.push(Gate::H(q));
        }
        let p1 = compile_clifford_t(n, &g1);

        let start = Instant::now();
        let a1 = amplitude_clifford_t_accel(&p1, &vec![0u8; n], 0);
        let t1 = start.elapsed();
        println!("[TIMING] clifford-ish amp = {:?} in {:?}", a1, t1);

        let mut g2 = g1.clone();
        for q in 0..n {
            if q % 2 == 0 {
                g2.push(Gate::T(q));
            }
            if q % 3 == 0 {
                g2.push(Gate::S(q));
            }
        }
        let p2 = compile_clifford_t(n, &g2);

        let start2 = Instant::now();
        let a2 = amplitude_clifford_t_accel(&p2, &vec![0u8; n], 0);
        let t2 = start2.elapsed();
        println!("[TIMING] clifford+T/S amp = {:?} in {:?}", a2, t2);

        assert!(a1.re.is_finite() && a1.im.is_finite());
        assert!(a2.re.is_finite() && a2.im.is_finite());
    }
}