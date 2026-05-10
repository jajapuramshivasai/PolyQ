use num_complex::Complex64;
use std::collections::{HashMap, HashSet};
use std::f64::consts::{FRAC_1_SQRT_2, PI};

// ========== Circuit Gate Definitions =============

#[derive(Debug, Clone)]
enum Gate {
    H(usize),           // Hadamard
    CZ(usize, usize),   // Controlled-Z
    CX(usize, usize),   // CNOT
    SWAP(usize, usize), // Swap
    Z(usize),
    S(usize),
    SDG(usize),
    RZ(usize, f64), // RZ rotation (angle in radians)
    T(usize),
    ID(usize), // Identity (used as placeholder in branching)
}

// ============ Branch Tree for Non-Clifford Gates ============

#[derive(Clone)]
struct BranchNode {
    weight: Complex64,
    gates: Vec<Gate>,
    children: Vec<BranchNode>,
    label: String,
}

impl BranchNode {
    fn leaf(weight: Complex64, gates: Vec<Gate>, label: &str) -> Self {
        BranchNode {
            weight,
            gates,
            children: vec![],
            label: label.to_string(),
        }
    }
}

// =========== Step 1: Circuit Partition ===========

struct CircuitPartition {
    e_left: Vec<Gate>,
    central: Vec<Gate>,
    e_right: Vec<Gate>,
}

fn partition_circuit(gates: &[Gate]) -> CircuitPartition {
    let first_h = gates
        .iter()
        .position(|g| matches!(g, Gate::H(_)))
        .unwrap_or(0);
    let last_h = gates
        .iter()
        .rposition(|g| matches!(g, Gate::H(_)))
        .unwrap_or(gates.len() - 1);
    CircuitPartition {
        e_left: gates[..first_h].to_vec(),
        central: gates[first_h..=last_h].to_vec(),
        e_right: if last_h + 1 < gates.len() {
            gates[last_h + 1..].to_vec()
        } else {
            vec![]
        },
    }
}

// =========== Step 2: Affine block as matrix (CX/SWAP/ID) ===========

type BitMatrix = Vec<Vec<u8>>;
type BitVec = Vec<u8>;

fn identity_matrix(n: usize) -> BitMatrix {
    (0..n)
        .map(|i| {
            let mut v = vec![0; n];
            v[i] = 1;
            v
        })
        .collect()
}

/// Compose two F2 bit matrices: t2 * t1
fn matmul_mod2(a: &BitMatrix, b: &BitMatrix) -> BitMatrix {
    let n = a.len();
    let mut res = vec![vec![0; n]; n];
    for i in 0..n {
        for j in 0..n {
            let mut sum = 0;
            for k in 0..n {
                sum ^= a[i][k] & b[k][j];
            }
            res[i][j] = sum;
        }
    }
    res
}

/// Affine transform struct: y' = T*y, x' = T*x
#[derive(Clone)]
struct Affine {
    t_matrix: BitMatrix,
    t_vec: BitVec, // not used for CNOT/SWAP only
}
impl Affine {
    fn identity(n: usize) -> Self {
        Affine {
            t_matrix: identity_matrix(n),
            t_vec: vec![0; n],
        }
    }

    fn apply(&self, bs: &BitVec) -> BitVec {
        let n = self.t_matrix.len();
        let mut res = vec![0; n];
        for i in 0..n {
            res[i] = self.t_vec[i];
            for j in 0..n {
                res[i] ^= self.t_matrix[i][j] & bs[j];
            }
            res[i] %= 2;
        }
        res
    }

    // Compose: result = self after other (i.e. self * other)
    fn compose(&self, other: &Affine) -> Affine {
        let t_matrix = matmul_mod2(&self.t_matrix, &other.t_matrix);
        let n = t_matrix.len();
        let mut t_vec = vec![0; n];
        for i in 0..n {
            for j in 0..n {
                t_vec[i] ^= self.t_matrix[i][j] & other.t_vec[j];
            }
            t_vec[i] ^= self.t_vec[i];
            t_vec[i] %= 2;
        }
        Affine { t_matrix, t_vec }
    }
}

/// Build F2 affine for input gates (just CNOT/SWAP)
fn build_affine(n: usize, gates: &[Gate]) -> Affine {
    let mut aff = Affine::identity(n);
    for g in gates {
        match *g {
            Gate::CX(c, t) => {
                if c != t {
                    let (a, b) = aff.t_matrix.split_at_mut(std::cmp::max(c, t));
                    let (row_c, row_t) = if c < t {
                        (&a[c], &mut b[0])
                    } else {
                        (&b[0], &mut a[t])
                    };
                    for i in 0..row_t.len() {
                        row_t[i] ^= row_c[i];
                    }
                }
            }
            Gate::SWAP(a, b) => {
                aff.t_matrix.swap(a, b);
                for row in &mut aff.t_matrix {
                    row.swap(a, b);
                }
            }
            _ => {},
        }
    }
    aff
}

// ========== Step 3: Branch Tree Construction and Z4 Polynomial ============

/// For any non-Clifford gate (RZ/T), branch
/// Weight for RZ(w): cos(w/2)*I + -i sin(w/2)*Z
fn build_branch_tree(gates: &[Gate], n_qubits: usize) -> BranchNode {
    fn _grow(
        idx: usize,
        gates: &[Gate],
        weight: Complex64,
        hist: Vec<Gate>,
        label: String,
    ) -> BranchNode {
        if idx == gates.len() {
            return BranchNode::leaf(weight, hist, &label);
        }
        match gates[idx] {
            Gate::RZ(q, theta) => {
                let mut hist_id = hist.clone();
                hist_id.push(Gate::ID(q));
                let mut hist_z = hist.clone();
                hist_z.push(Gate::Z(q));
                let c = weight * Complex64::new((theta / 2.0).cos(), 0.0);
                let s = weight * Complex64::new(0.0, -(theta / 2.0).sin());
                let c_node = _grow(idx + 1, gates, c, hist_id, format!("{label}->I"));
                let z_node = _grow(idx + 1, gates, s, hist_z, format!("{label}->Z"));
                let mut node = BranchNode::leaf(weight, vec![], &label);
                node.children = vec![c_node, z_node];
                node
            }
            Gate::T(q) => {
                let theta = PI / 4.0;
                let mut hist_id = hist.clone();
                hist_id.push(Gate::ID(q));
                let mut hist_z = hist.clone();
                hist_z.push(Gate::Z(q));
                let c = weight * Complex64::new((theta / 2.0).cos(), 0.0);
                let s = weight * Complex64::new(0.0, -(theta / 2.0).sin());
                let c_node = _grow(idx + 1, gates, c, hist_id, format!("{label}->I"));
                let z_node = _grow(idx + 1, gates, s, hist_z, format!("{label}->Z"));
                let mut node = BranchNode::leaf(weight, vec![], &label);
                node.children = vec![c_node, z_node];
                node
            }
            _ => {
                let mut new_hist = hist.clone();
                new_hist.push(gates[idx].clone());
                _grow(idx + 1, gates, weight, new_hist, label)
            }
        }
    }
    _grow(
        0,
        gates,
        Complex64::new(1.0, 0.0),
        vec![],
        "ROOT".to_string(),
    )
}

// ========== Step 4: Block Diagonalization (F2, symmetric) ============

/// Block-diagonalize quadratic form over F2.
/// Input: symmetric F2 square matrix. Output: ops to block-diagonalize.
/// Returns (ops, reduced_form, rank)
fn block_diagonalize_f2(mat: &BitMatrix) -> (Vec<(usize, usize)>, BitMatrix, usize) {
    let n = mat.len();
    let mut mat = mat.clone();
    let mut ops = vec![];
    let mut p = 0;
    let mut rank = 0;
    while p + 1 < n {
        // Find pivot (i, j)
        let mut pivot = None;
        for i in p..n {
            for j in i + 1..n {
                if mat[i][j] == 1 {
                    pivot = Some((i, j));
                    break;
                }
            }
            if pivot.is_some() {
                break;
            }
        }
        if pivot.is_none() {
            break;
        }
        let (i, j) = pivot.unwrap();
        if i != p {
            mat.swap(p, i);
            for row in &mut mat {
                row.swap(p, i);
            }
            ops.push((p, i)); // SWAP rows/cols p <-> i
        }
        let j_act = if j == p { i } else { j };
        if j_act != p + 1 {
            mat.swap(p + 1, j_act);
            for row in &mut mat {
                row.swap(p + 1, j_act);
            }
            ops.push((p + 1, j_act));
        }
        // Eliminate
        let (rp, rp1) = (mat[p].clone(), mat[p + 1].clone());
        for k in (p + 2)..n {
            if mat[k][p] == 1 {
                for m in 0..n {
                    mat[k][m] ^= rp1[m];
                    mat[m][k] ^= rp1[m];
                }
            }
            if mat[k][p + 1] == 1 {
                for m in 0..n {
                    mat[k][m] ^= rp[m];
                    mat[m][k] ^= rp[m];
                }
            }
        }
        rank += 2;
        p += 2;
    }
    (ops, mat, rank)
}

// ========== Step 5: Clifford Amplitude/Statevector Computation ===========

// Get vector dot product mod 2 (for F2)
fn vec_dot2(a: &BitVec, b: &BitVec) -> u8 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum::<u8>() % 2
}

/// Amplitude for a (Clifford) quadratic form, optionally with (branch) phase vector
fn clifford_amplitude(
    x: &BitVec,
    y: &BitVec,
    qmat: &BitMatrix,       // quadratic form (u^T Q u)
    lvec: &BitVec,          // linear form (l^T u)
    v0: u8,                 // global phase offset
    ops: &[(usize, usize)], // block-diag ops (for lvec)
    rank: usize,
    hadamard_count: usize,
) -> Complex64 {
    // Form: exp(i*v0*π/2) * ... * sum_z2(u) (-1)^(Q(u)+l(u)) / 2^{h/2}
    // Here, we return norm if lvec only zero, else 0 if kernel collision
    for k in rank..lvec.len() {
        if lvec[k] != 0 {
            return Complex64::new(0.0, 0.0);
        }
    }
    // For demonstration, just Clifford sum: norm
    let norm = (2.0f64).powf(-(hadamard_count as f64) / 2.0);
    // Apply any simple phase
    let phases = [1.0, 1.0, -1.0, -1.0]; // "Z4 phases": only 0/2 for Clifford
    Complex64::new(norm * phases[v0 as usize], 0.0)
}

// ============= Util: Generate all n-bit BitVecs ==============

fn all_bitstrings(n: usize) -> Vec<BitVec> {
    let mut res = Vec::with_capacity(1 << n);
    for bits in 0..(1 << n) {
        let mut v = vec![0; n];
        for i in 0..n {
            v[i] = ((bits >> i) & 1) as u8;
        }
        res.push(v);
    }
    res
}

// =============== DEMO / MAIN ===============

fn main() {
    // Toy example 2-qubit Clifford+T circuit: H(0) -- CZ(0,1) -- T(1) -- H(1)
    let q = 2;
    let circuit = vec![Gate::H(0), Gate::CZ(0, 1), Gate::T(1), Gate::H(1)];

    // Step 1: Partition
    let CircuitPartition {
        e_left,
        central,
        e_right,
    } = partition_circuit(&circuit);
    println!(
        "E_left: {:?}\nCentral: {:?}\nE_right: {:?}",
        e_left, central, e_right
    );

    // Step 2: Affine (not relevant here: only CX/SWAP in E left/right)
    let aff_left = build_affine(q, &e_left);
    let aff_right = build_affine(q, &e_right);

    // Step 3: Branch out non-Clifford (T) in central, get all leaves
    let branch_root = build_branch_tree(&central, q);

    // For each branch, reconstruct Clifford part and compute amplitude matrix
    // For demo: only works fully for Clifford circuits (no T/RZ).
    let leaves = {
        let mut out = vec![];
        fn collect_leaves(
            node: &BranchNode,
            list: &mut Vec<(&Vec<Gate>, Complex64, String)>,
            label: String,
        ) {
            if node.children.is_empty() {
                list.push((&node.gates, node.weight, label));
            } else {
                for (i, c) in node.children.iter().enumerate() {
                    collect_leaves(c, list, format!("{label}->{}", c.label));
                }
            }
        }
        collect_leaves(&branch_root, &mut out, "ROOT".to_string());
        out
    };

    let all_x = all_bitstrings(q);
    let all_y = all_bitstrings(q);

    // For each branch: build quadratic form Q, L for Clifford gates
    for (gates, weight, branch_label) in leaves {
        // Only support H and CZ for Clifford (for demo)
        let mut hadamard_count = 0;
        // Quadratic form: Q(u) = u^T B u (B: F2 symmetric matrix)
        let mut bmat = vec![vec![0; q]; q];
        for g in gates.iter() {
            match *g {
                Gate::H(idx) => hadamard_count += 1,
                Gate::CZ(i, j) => {
                    bmat[i][j] ^= 1;
                    bmat[j][i] ^= 1;
                }
                _ => {}
            }
        }
        // Block-diagonalize
        let (ops, bmat_red, rank) = block_diagonalize_f2(&bmat);

        // Print analytic formula for this branch
        println!("BRANCH: {}", branch_label);
        println!("  Weight: {}", weight);
        println!("  Quadratic matrix:\n    {:?}", bmat);
        println!("  Reduced form: {:?}, rank: {}", bmat_red, rank);

        // For demo: Linear (L) zero, global phase 0 for Clifford
        let lvec = vec![0; q];
        let v0 = 0;

        // Demo: Print amplitudes (for input |00⟩..|11⟩)
        println!("  Statevector:");
        for (yy, y) in all_y.iter().enumerate() {
            let aff_y = aff_right.apply(y);
            let mut total = Complex64::new(0.0, 0.0);
            for (xx, x) in all_x.iter().enumerate() {
                let aff_x = aff_left.apply(x);
                let amp = clifford_amplitude(
                    &aff_x,
                    &aff_y,
                    &bmat,
                    &lvec,
                    v0,
                    &ops,
                    rank,
                    hadamard_count,
                );
                total += amp;
                println!(
                    "    T({}, {}|{}) = {:.3} + {:.3}i",
                    xx, yy, branch_label, amp.re, amp.im
                );
            }
            println!(
                "    Norm[{}|{}] = {:.3}",
                y.iter().map(|bi| bi.to_string()).collect::<String>(),
                branch_label,
                total.norm()
            );
        }
    }
}
