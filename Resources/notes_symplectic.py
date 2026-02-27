"""
Quantum circuit simulator for Clifford+T circuits using Z8 phase polynomial representation.
This module implements a simulator for circuits over H, Z, S, CZ, T by converting
them to a Z8 phase polynomial, extracting Z4 quadratic components and Z8 cubic terms, and
evaluating amplitudes via symplectic Gram–Schmidt on the Clifford part plus enumeration
over cubic slices for T/CCZ contributions.
Dependencies
------------
- numpy: Matrix and vector operations.
- qiskit: Circuit representation and reference Statevector for benchmarking.
- time: Benchmark timing utilities.
Limitations
-----------
- Supports only H, Z, S, CZ, T, other gates raise NotImplementedError.
- T/CCZ introduce cubic/odd-linear slicing, which increases complexity.
- Memory and runtime scale with qubits and number of T/CCZ gates.
"""
from __future__ import annotations
import time
from typing import List, Tuple, Optional
import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector

# =========================
# Circuit → Z8 polynomial
# =========================
def circuit_to_z8_polynomial(
    qc: QuantumCircuit,
) -> Tuple[List[Tuple[int, List[int]]], List[int], int, int]:
    """
    Convert a Qiskit circuit with H, Z, S, CZ, T, CCZ into a Z8 phase polynomial.
    Each gate contributes terms over Z8: H introduces an auxiliary variable with a 4·old·new
    bilinear term; Z/S/T contribute linear phases with coefficients 4/2/1; CZ adds a 4·u·v
    bilinear; CCZ adds a 4·u·v·w cubic term.
    Parameters
    ----------
    qc : QuantumCircuit
        Circuit composed of H, Z, S, CZ, T, CCZ; non-unitary ops like measure/reset are skipped.
    Returns
    -------
    terms : list[tuple[int, list[int]]]
        List of (coefficient mod 8, variable index list) for phase-polynomial terms.
    y_idx : list[int]
        Final variable index on each output wire after H-induced wire splits.
    maxnewvar : int
        Total variable count including auxiliaries introduced by H gates.
    h_count : int
        Total number of Hadamard gates for overall normalization factor \(2^{-h/2}\).
    Notes
    -----
    Variable indices start from input bits 0..n-1 and increase as H gates introduce new symbols.
    """
    n = qc.num_qubits
    wirearray: List[List[int]] = [[i] for i in range(n)]  # initial variables x_0..x_{n-1}
    maxnewvar = n
    terms: List[Tuple[int, List[int]]] = []  # (coeff mod 8, [indices])
    h_count = 0
    for instruction in qc.data:
        name = instruction.operation.name.lower()
        qids = [qc.find_bit(q).index for q in instruction.qubits]
        if name == "h":
            h_count += 1
            q = qids[0]
            oldv = wirearray[q][-1]
            wirearray[q].append(maxnewvar)
            # H contributes 4 * old * new (introduces a new summation variable)
            terms.append((4, [oldv, maxnewvar]))
            maxnewvar += 1
        elif name == "z":
            v = wirearray[qids[0]][-1]
            terms.append((4, [v]))
        elif name == "s":
            v = wirearray[qids[0]][-1]
            terms.append((2, [v]))
        elif name == "t":
            v = wirearray[qids[0]][-1]
            terms.append((1, [v]))
        elif name == "cz":
            c = wirearray[qids[0]][-1]
            t = wirearray[qids[1]][-1]
            terms.append((4, [c, t]))
        elif name == "ccz":
            a = wirearray[qids[0]][-1]
            b = wirearray[qids[1]][-1]
            c = wirearray[qids[2]][-1]
            # CCZ contributes a cubic phase (-1 on |111>), i.e., coefficient 4 in Z8 on xyz
            terms.append((4, [a, b, c]))
        elif name in ("id", "barrier", "delay", "measure", "reset"):
            # Non-unitary/visualization ops are ignored for phase-polynomial simulation
            continue
        else:
            raise NotImplementedError(
                f"Only H, Z, S, CZ, T, CCZ are supported, found: {name}"
            )
    y_idx = [wirearray[q][-1] for q in range(n)]
    return terms, y_idx, maxnewvar, h_count

# =========================
# Extract Z4 quadratic and Z8 cubic components
# =========================
def extract_components_z8(
    terms: List[Tuple[int, List[int]]],
    m: int,
) -> Tuple[np.ndarray, np.ndarray, int, List[Tuple[int, List[int]]]]:
    """
    Split Z8 phase-polynomial terms into quadratic matrix, linear vector, constant, and cubic list.
    Terms with degree 0/1/2 are accumulated into a constant eps8, a Z8-linear vector L, and an
    upper-triangular Z8 matrix Q on indices i<j; degree-3 terms are collected for later slicing.
    Parameters
    ----------
    terms : list[tuple[int, list[int]]]
        Z8 terms as (coefficient, indices) with degrees 0..3 expected.
    m : int
        Number of variables in scope for Q and L allocation.
    Returns
    -------
    Q : np.ndarray
        m x m integer matrix (stored mod 8 logically) containing quadratic coefficients for i < j.
    L : np.ndarray
        Length-m integer vector (mod 8) of linear coefficients.
    eps8 : int
        Z8 constant term (mod 8).
    cubic_terms : list[tuple[int, list[int]]]
        Degree-3 terms (coefficient mod 8, sorted indices).
    Notes
    -----
    Degree > 3 is not expected for this gate set; such input raises a ValueError.
    """
    Q = np.zeros((m, m), dtype=int)  # only i<j used
    L = np.zeros(m, dtype=int)
    eps8 = 0
    cubic_terms: List[Tuple[int, List[int]]] = []
    for c, idx in terms:
        c &= 7
        if c == 0:
            continue
        d = len(idx)
        if d == 0:
            eps8 = (eps8 + c) & 7
        elif d == 1:
            L[idx[0]] = (L[idx[0]] + c) & 7
        elif d == 2:
            i, j = sorted(idx)
            if i == j:
                L[i] = (L[i] + c) & 7
            else:
                Q[i, j] = (Q[i, j] + c) & 7
        elif d == 3:
            cubic_terms.append((c, sorted(idx)))
        else:
            raise ValueError("Degree > 3 not expected")
    return Q, L, eps8, cubic_terms

def restrict_and_split_to_z4(
    terms: List[Tuple[int, List[int]]],
    y_idx: List[int],
    n_qubits: int,
    x_bits: int,
    y_bits: int,
) -> Optional[
    Tuple[
        np.ndarray,
        np.ndarray,
        int,
        int,
        List[Tuple[int, List[int]]],
        bool,
        List[int],
    ]
]:
    """
    Apply input/output bit constraints, eliminate fixed vars, and split Z8 to Z4 components.
    Fixed variables are substituted using the input x_bits and output y_bits consistency on
    the final wire variables y_idx; remaining internal vars are renumbered to 0..mu-1.
    The Z8 components reduce to Z4 via L4 = floor(L/2) mod 4 and eps4 = floor(eps8/2) mod 4,
    and B is the GF(2) adjacency extracted from Q where entries are 4 mod 8.
    Parameters
    ----------
    terms : list[tuple[int, list[int]]]
        Z8 polynomial terms prior to restriction.
    y_idx : list[int]
        Final variable index per output wire for enforcing y_bits.
    n_qubits : int
        Number of circuit qubits used to interpret x_bits and y_bits.
    x_bits : int
        Input computational basis state as integer bitstring (LSB = qubit 0).
    y_bits : int
        Output computational basis state as integer bitstring (LSB = qubit 0).
    Returns
    -------
    B : np.ndarray | None
        mu x mu binary symmetric matrix over GF(2) with 1 where Q[i,j] == 4 mod 8; None if inconsistent I/O.
    L4 : np.ndarray
        Length-mu Z4 linear coefficients (mod 4).
    eps4 : int
        Z4 constant term (mod 4).
    eps_odd : int
        The odd bit of eps8 (eps8 mod 2), to carry remaining i-phase.
    cubic_terms : list[tuple[int, list[int]]]
        Any remaining cubic terms after restriction (should be empty for a slice).
    trivial : bool
        True if mu == 0 after restriction (no internal summation variables).
    odd_linear_vars : list[int]
        Indices with odd Z8 linear coefficients to include in slicing set.
    Notes
    -----
    Returns None if x_bits/y_bits assignments are inconsistent under y_idx mapping (amplitude = 0).
    """
    fixed = {i: ((x_bits >> i) & 1) for i in range(n_qubits)}
    for q in range(n_qubits):
        outv = y_idx[q]
        want = (y_bits >> q) & 1
        if outv in fixed:
            if fixed[outv] != want:
                return None
        else:
            fixed[outv] = want
    # Substitute fixed variables
    subst: List[Tuple[int, List[int]]] = []
    for c, idx in terms:
        prod_fixed = 1
        rem: List[int] = []
        for v in idx:
            if v in fixed:
                prod_fixed &= fixed[v]
                if prod_fixed == 0:
                    break
            else:
                rem.append(v)
        if prod_fixed == 1:
            subst.append((c & 7, rem))
    # Collect remaining internal vars and renumber to 0..mu-1
    allu = set()
    for _, ids in subst:
        allu.update(ids)
    u_vars = sorted(allu)
    mu = len(u_vars)
    vmap = {g: i for i, g in enumerate(u_vars)}
    renum: List[Tuple[int, List[int]]] = []
    eps8_extra = 0
    for c, ids in subst:
        if len(ids) == 0:
            eps8_extra = (eps8_extra + c) & 7
        else:
            renum.append((c, [vmap[v] for v in ids]))
    # Build Q, L (Z8) and split to Z4
    Q, Lz8, eps8, cubic_terms = extract_components_z8(renum, mu)
    eps8 = (eps8 + eps8_extra) & 7
    eps_odd = eps8 & 1
    eps4 = (eps8 >> 1) & 3
    if mu == 0:
        B = np.zeros((0, 0), dtype=np.uint8)
        L4 = np.zeros((0,), dtype=np.int8)
        return B, L4, int(eps4), int(eps_odd), cubic_terms, True, []
    # Bilinear B from Q: edges where Q[i,j] == 4 mod 8
    B = np.zeros((mu, mu), dtype=np.uint8)
    mask = (Q % 8) == 4
    idxs = np.argwhere(mask)
    for i, j in idxs:
        if i < j:
            B[i, j] = 1
            B[j, i] = 1
    # Z4 linear vector
    L4 = np.zeros(mu, dtype=np.int8)
    for i in range(mu):
        c = int(Lz8[i] & 7)
        L4[i] = (c >> 1) & 3
    odd_linear_vars = [i for i in range(mu) if (Lz8[i] & 1) == 1]
    return B, L4, int(eps4), int(eps_odd), cubic_terms, False, odd_linear_vars

# =========================
# Dickson reduction (B)
# =========================
def symplectic_gram_schmidt(B: np.ndarray) -> Tuple[np.ndarray, int]:
    """
    Perform symplectic Gram–Schmidt over GF(2) to block-diagonalize a symmetric binary matrix.
    Produces a change-of-basis matrix T such that T^T B T is block-diagonal with k hyperbolic
    2x2 blocks; remaining columns span the kernel part (no pairs).
    Parameters
    ----------
    B : np.ndarray
        m x m symmetric binary matrix over GF(2) encoding the quadratic form’s bilinear part.
    Returns
    -------
    T : np.ndarray
        m x m' binary matrix (m' ≥ rank) giving the new basis vectors as columns.
    k : int
        Number of symplectic pairs (2x2 hyperbolic blocks) extracted from B.
    Notes
    -----
    The routine greedily finds pairs with nonzero bilinear coupling and orthogonalizes others.
    """
    m = B.shape[0]
    W = np.eye(m, dtype=np.uint8)
    active = list(range(m))
    T_cols: List[np.ndarray] = []
    k = 0
    def bilinear(u: np.ndarray, v: np.ndarray) -> int:
        return int(((u @ B) @ v) & 1)
    while True:
        pair = None
        for ia in range(len(active)):
            a = active[ia]
            u = W[:, a]
            uTB = (u @ B) & 1
            for ib in range(ia + 1, len(active)):
                b = active[ib]
                v = W[:, b]
                if int(uTB @ v) & 1:
                    pair = (ia, ib)
                    break
            if pair:
                break
        if not pair:
            break
        ia, ib = pair
        a = active[ia]
        b = active[ib]
        u = W[:, a].copy()
        v = W[:, b].copy()
        T_cols.append(u)
        T_cols.append(v)
        k += 1
        new_active = []
        for t in active:
            if t == a or t == b:
                continue
            w = W[:, t]
            cu = bilinear(u, w)
            cv = bilinear(v, w)
            if cu:
                w = w ^ v
            if cv:
                w = w ^ u
            W[:, t] = w
            new_active.append(t)
        active = new_active
    for t in active:
        T_cols.append(W[:, t].copy())
    if len(T_cols) == 0:
        return np.zeros((0, 0), dtype=np.uint8), 0
    T = np.stack(T_cols, axis=1).astype(np.uint8)
    return T, k

# =========================
# Closed-form sums
# =========================
def i_pow(e: int) -> complex:
    """
    Compute the 8th-root-of-unity power i^(e mod 8) used for phase factors.
    Parameters
    ----------
    e : int
        Exponent taken modulo 8 before evaluation.
    Returns
    -------
    complex
        Value of i^e for e in {0,..,7} on the unit circle.
    Notes
    -----
    Implemented as a small lookup table for numerical stability and speed.
    """
    e &= 7
    return (
        1 + 0j,
        np.exp(1j * np.pi / 4),
        1j,
        np.exp(3j * np.pi / 4),
        -1 + 0j,
        np.exp(5j * np.pi / 4),
        -1j,
        np.exp(7j * np.pi / 4),
    )[e]

def pair_sum_clifford(alpha4: int, beta4: int) -> complex:
    """
    Closed-form amplitude factor for a symplectic pair in the Clifford sum.
    Parameters
    ----------
    alpha4 : int
        Z4 linear coefficient of the first variable in the pair.
    beta4 : int
        Z4 linear coefficient of the second variable in the pair.
    Returns
    -------
    complex
        Complex factor multiplying the global amplitude from summing a 2-variable block.
    """
    return 1.0 + i_pow(2 * beta4) + i_pow(2 * alpha4) - i_pow(2 * ((alpha4 + beta4) & 3))

def kernel_factor_clifford(t4: int) -> complex:
    """
    Closed-form amplitude factor for a kernel (unpaired) variable in the Clifford sum.
    Parameters
    ----------
    t4 : int
        Z4 linear coefficient (mod 4) for the kernel variable.
    Returns
    -------
    complex
        Complex factor from summing a single unpaired variable.
    """
    return 1.0 + i_pow(2 * (t4 & 3))

# =========================
# Amplitude computation
# =========================
def amplitude_clifford_t(
    terms: List[Tuple[int, List[int]]],
    y_idx: List[int],
    n_qubits: int,
    h_count: int,
    x_bits: int,
    y_bits: int,
) -> complex:
    """
    Compute ⟨y_bits|U|x_bits⟩ for a Clifford+T circuit encoded as Z8 terms.
    The routine restricts by x/y, splits to Z4, and if cubic/odd-linear variables exist,
    enumerates their 0/1 assignments; each slice is then a pure Clifford sum evaluated by
    symplectic Gram–Schmidt, and all slices are accumulated with the global i-phase and
    \(2^{-h/2}\) normalization.
    Parameters
    ----------
    terms : list[tuple[int, list[int]]]
        Z8 phase-polynomial terms for the circuit.
    y_idx : list[int]
        Final variable indices of output wires.
    n_qubits : int
        Number of qubits in the circuit.
    h_count : int
        Total H gate count for normalization \(2^{-h/2}\).
    x_bits : int
        Input computational basis index (integer).
    y_bits : int
        Output computational basis index (integer).
    Returns
    -------
    complex
        Transition amplitude from |x_bits⟩ to |y_bits⟩.
    Notes
    -----
    Slicing set is the union of cubic-term variables and odd-Z8-linear variables to remove
    non-Clifford contributions before the Clifford evaluation step.
    """
    out = restrict_and_split_to_z4(terms, y_idx, n_qubits, x_bits, y_bits)
    if out is None:
        return 0.0 + 0.0j
    B, L4, eps4, eps_odd, cubic_terms, trivial, odd_linear_vars = out
    if trivial or B.shape[0] == 0:
        return i_pow(2 * eps4 + eps_odd) * (2.0 ** (-0.5 * h_count))
    mu = B.shape[0]
    # Collect slicing vars: cubic and odd-linear Z8 variables
    cubic_vars = set()
    for _, idx in cubic_terms:
        cubic_vars.update(idx)
    cubic_vars.update(odd_linear_vars)
    t = len(cubic_vars)
    # Pure Clifford case: no slicing needed
    if t == 0:
        T, k = symplectic_gram_schmidt(B)
        Lp = np.zeros_like(L4)
        for j in range(mu):
            mask_bool = (T[:, j] != 0)
            if mask_bool.any():
                Lp[j] = int(np.sum(L4[mask_bool]) % 4)
            else:
                Lp[j] = 0
        S = i_pow(2 * eps4 + eps_odd)
        for p in range(k):
            a4 = int(Lp[2 * p] & 3)
            b4 = int(Lp[2 * p + 1] & 3)
            S *= pair_sum_clifford(a4, b4)
        for j in range(2 * k, mu):
            S *= kernel_factor_clifford(int(Lp[j] & 3))
        return S * (2.0 ** (-0.5 * h_count))
    # Clifford+T/CCZ: enumerate slicing assignments
    total = 0.0 + 0.0j
    cubic_vars = sorted(cubic_vars)
    vmap = {cubic_vars[i]: i for i in range(t)}
    for assign_int in range(1 << t):
        assign = [(assign_int >> i) & 1 for i in range(t)]
        var_map = {cubic_vars[i]: assign[i] for i in range(t)}
        # Substitute slicing vars
        subst: List[Tuple[int, List[int]]] = []
        for c, idx in terms:
            prod_fixed = 1
            rem: List[int] = []
            for v in idx:
                if v in var_map:
                    prod_fixed &= var_map[v]
                    if prod_fixed == 0:
                        break
                else:
                    rem.append(v)
            if prod_fixed == 1:
                subst.append((c & 7, rem))
        # Restrict again with the same x/y
        new_out = restrict_and_split_to_z4(subst, y_idx, n_qubits, x_bits, y_bits)
        if new_out is None:
            continue
        new_B, new_L4, new_eps4, new_eps_odd, new_cubic_terms, new_trivial, _ = new_out
        # After slicing, no cubic terms should remain in a valid slice
        if new_cubic_terms:
            continue
        if new_trivial or new_B.shape[0] == 0:
            slice_val = i_pow(2 * new_eps4 + new_eps_odd)
        else:
            T, k = symplectic_gram_schmidt(new_B)
            Lp = np.zeros_like(new_L4)
            for j in range(new_B.shape[0]):
                mask_bool = (T[:, j] != 0)
                if mask_bool.any():
                    Lp[j] = int(np.sum(new_L4[mask_bool]) % 4)
                else:
                    Lp[j] = 0
            S = i_pow(2 * new_eps4 + new_eps_odd)
            for p in range(k):
                a4 = int(Lp[2 * p] & 3)
                b4 = int(Lp[2 * p + 1] & 3)
                S *= pair_sum_clifford(a4, b4)
            for j in range(2 * k, new_B.shape[0]):
                S *= kernel_factor_clifford(int(Lp[j] & 3))
            slice_val = S
        total += slice_val
    return total * (2.0 ** (-0.5 * h_count))

# =========================
# Statevector simulation
# =========================
def simulate_statevector_z8(qc: QuantumCircuit, x_bits: int = 0) -> np.ndarray:
    """
    Compute the full output statevector for input |x_bits⟩ using the Z8 simulator.
    Parameters
    ----------
    qc : QuantumCircuit
        Circuit to simulate (H, Z, S, CZ, T, CCZ).
    x_bits : int, optional
        Input basis state index as integer; default 0 for |0…0⟩.
    Returns
    -------
    np.ndarray
        Complex statevector of length 2^n ordered by integer y in little-endian convention.
    Notes
    -----
    The amplitude ⟨y|U|x⟩ is computed for all y via amplitude_clifford_t and collected into a vector.
    """
    terms, y_idx, _, h_count = circuit_to_z8_polynomial(qc)
    n = qc.num_qubits
    dim = 1 << n
    psi = np.empty(dim, dtype=np.complex128)
    for y in range(dim):
        psi[y] = amplitude_clifford_t(terms, y_idx, n, h_count, x_bits, y)
    return psi