import numpy as np
from typing import List, Dict, Set, Tuple, Optional, Any
from qiskit import QuantumCircuit

class DicksonOp:
    def __init__(self, op_type: str, a: int, b: Optional[int] = None):
        self.type, self.a, self.b = op_type, a, b

class DicksonEngine:
    """
    Optimized F2/Z4 Quadratic Form Engine.
    Handles structural caching and high-speed statevector evaluation via Gray code.
    """
    def __init__(self, circuit: QuantumCircuit):
        self.circuit = circuit
        self.num_qubits = circuit.num_qubits
        self.gates: List[Tuple[str, List[int]]] = []

        # Structure (Cached)
        self.b4: Optional[np.ndarray] = None
        self.output_vars: List[int] = []
        self.num_h = 0
        self.n_vars = 0
        self.rank = 0
        self.uvars_skeleton: List[int] = []
        self.ops: List[DicksonOp] = []
        self.b_reduced: Optional[np.ndarray] = None
        self.v4 = np.zeros(0)
        self._passthrough_bits: int = 0

        self._translate_circuit()
        self.compile_structure()

    def _translate_circuit(self):
        for instr in self.circuit.data:
            name = instr.operation.name.lower()
            idxs = [self.circuit.find_bit(q).index for q in instr.qubits]
            if name in ['h', 'z', 's', 'sdg', 'cz', 'id']:
                self.gates.append((name.upper(), idxs))

    def compile_structure(self):
        wires = [[i] for i in range(self.num_qubits)]
        next_v, self.num_h = self.num_qubits, 0
        adj: Dict[int, Set[int]] = {i: set() for i in range(self.num_qubits)}
        for g_type, idxs in self.gates:
            if g_type == 'H':
                q = idxs[0]
                prev, cur = wires[q][-1], next_v
                next_v += 1
                self.num_h += 1
                wires[q].append(cur)
                adj.setdefault(cur, set()).add(prev)
                adj[prev].add(cur)
            elif g_type == 'CZ':
                v1, v2 = wires[idxs[0]][-1], wires[idxs[1]][-1]
                adj[v1].add(v2)
                adj[v2].add(v1)
        self.n_vars = next_v
        self.b4 = np.zeros((self.n_vars, self.n_vars), dtype=np.int8)
        for i in range(self.n_vars):
            for j in adj.get(i, []):
                self.b4[i, j] = 1
        self.output_vars = [wires[q][-1] for q in range(self.num_qubits)]
        self.uvars_skeleton = [
            i for i in range(self.n_vars)
            if i >= self.num_qubits and i not in self.output_vars
        ]
        nu = len(self.uvars_skeleton)
        b_u = np.zeros((nu, nu), dtype=np.int8)
        for ui, u_o in enumerate(self.uvars_skeleton):
            for uj, uj_o in enumerate(self.uvars_skeleton):
                if ui < uj and self.b4[u_o, uj_o]:
                    b_u[ui, uj] = b_u[uj, ui] = 1
        self.ops, self.b_reduced, self.rank = self._plan_dickson(b_u, nu)
        self.v4 = np.zeros(self.n_vars, dtype=np.int8)

        self._passthrough_bits = 0
        for i in range(self.num_qubits):
            if self.output_vars[i] < self.num_qubits:
                self._passthrough_bits |= (1 << i)

    def _plan_dickson(self, b: np.ndarray, n: int):
        b_w, ops, r, p = np.copy(b), [], 0, 0
        while p + 1 < n:
            pivot = next(((i, j) for i in range(p, n) for j in range(i + 1, n) if b_w[i, j] == 1), None)
            if not pivot: break
            i, j = pivot
            if i != p:
                b_w[[p, i]] = b_w[[i, p]]
                b_w[:, [p, i]] = b_w[:, [i, p]]
                ops.append(DicksonOp('SWAP', p, i))
                j_act = p if j == i else (i if j == p else j)
            else: j_act = j
            if j_act != p + 1:
                b_w[[p + 1, j_act]] = b_w[[j_act, p + 1]]
                b_w[:, [p + 1, j_act]] = b_w[:, [j_act, p + 1]]
                ops.append(DicksonOp('SWAP', p + 1, j_act))
            rp, rp1 = np.copy(b_w[p, :]), np.copy(b_w[p + 1, :])
            for k in range(p + 2, n):
                if b_w[k, p]:
                    b_w[k, :] ^= rp1
                    b_w[:, k] ^= rp1
                    ops.append(DicksonOp('ADD', p + 1, k))
                if b_w[k, p + 1]:
                    b_w[k, :] ^= rp
                    b_w[:, k] ^= rp
                    ops.append(DicksonOp('ADD', p, k))
            r += 2
            p += 2
        return ops, b_w, r

    def print_analytic_formula(self, transition_mode=True, weight=1.0, branch_label=""):
        nu = len(self.uvars_skeleton)
        n_cols = 2 * self.num_qubits + 1 if transition_mode else self.num_qubits + 1
        coeff_matrix = np.zeros((nu, n_cols), dtype=np.int8)
        for ui, orig_u in enumerate(self.uvars_skeleton):
            coeff_matrix[ui, -1] = self.v4[orig_u] % 4
            if transition_mode:
                for xi in range(self.num_qubits):
                    if self.b4[orig_u, xi]: coeff_matrix[ui, xi] = 2 
            for yi in range(self.num_qubits):
                if self.b4[orig_u, self.output_vars[yi]]:
                    col_idx = self.num_qubits + yi if transition_mode else yi
                    coeff_matrix[ui, col_idx] = (coeff_matrix[ui, col_idx] + 2) % 4
        for op in self.ops:
            if op.type == 'SWAP': coeff_matrix[[op.a, op.b]] = coeff_matrix[[op.b, op.a]]
            elif op.type == 'ADD': coeff_matrix[op.b] = (coeff_matrix[op.b] + coeff_matrix[op.a]) % 4

        def build_xor_expr(j):
            terms = ["1"] if coeff_matrix[j, -1] in (2, 3) else []
            if transition_mode:
                terms += [f"x_{xi}" for xi in range(self.num_qubits) if coeff_matrix[j, xi] == 2]
                terms += [f"y_{yi}" for yi in range(self.num_qubits) if coeff_matrix[j, self.num_qubits + yi] == 2]
            else: terms += [f"y_{yi}" for yi in range(self.num_qubits) if coeff_matrix[j, yi] == 2]
            return " XOR ".join(terms) if terms else "0"

        header = f"BRANCH: {branch_label} (Weight: {weight})" if branch_label else "CLIFFORD ANALYTIC FORM"
        print(f"{'-'*40}\n{header}\n{'-'*40}")
        n1_terms = [f"({build_xor_expr(j)})" for j in range(self.rank)]
        print(f"1. Exponential Sum Variables (N1):")
        if n1_terms:
            for idx, term in enumerate(n1_terms): print(f"   v_{idx}: {term}")
        else: print("   None")
        zero_conditions = [build_xor_expr(k) for k in range(self.rank, nu)]
        print("\n2. Zero Amplitude (Kernel) Conditions:")
        if zero_conditions:
            for idx, cond in enumerate(zero_conditions): print(f"   K_{idx}: ({cond}) == 0")
        else: print("   No kernel constraints.\n")

    def set_phases(self, circuit: QuantumCircuit):
        self.v4.fill(0)
        wires = [[i] for i in range(self.num_qubits)]
        nv = self.num_qubits
        for instr in circuit.data:
            name = instr.operation.name.lower()
            idxs = [circuit.find_bit(q).index for q in instr.qubits]
            q = idxs[0]
            if name == 'h':
                wires[q].append(nv)
                nv += 1
            elif name == 'z': self.v4[wires[q][-1]] = (self.v4[wires[q][-1]] + 2) % 4
            elif name == 's': self.v4[wires[q][-1]] = (self.v4[wires[q][-1]] + 1) % 4
            elif name == 'sdg': self.v4[wires[q][-1]] = (self.v4[wires[q][-1]] + 3) % 4

    def get_amplitude(self, y: int, x: int = 0) -> complex:
        fixed = [None] * self.n_vars
        for i in range(self.num_qubits): fixed[i] = (x >> i) & 1
        for i in range(self.num_qubits):
            bit, ov = (y >> i) & 1, self.output_vars[i]
            if fixed[ov] is not None and fixed[ov] != bit: return 0j
            fixed[ov] = bit
        eps = self._calc_eps_from_fixed(fixed)
        nu = len(self.uvars_skeleton)
        vu = np.zeros(nu, dtype=np.int8)
        for ui, u_o in enumerate(self.uvars_skeleton):
            vu[ui] = self.v4[u_o] % 4
            for v_idx, val in enumerate(fixed):
                if val == 1 and self.b4[u_o, v_idx]: vu[ui] = (vu[ui] + 2) % 4
        for op in self.ops:
            if op.type == 'SWAP': vu[op.a], vu[op.b] = vu[op.b], vu[op.a]
            elif op.type == 'ADD': vu[op.b] = (vu[op.b] + vu[op.a]) % 4
        return ([1, 1j, -1, -1j][eps % 4] * self._eval_canonical_sum(vu, nu) * (2 ** (-self.num_h / 2)))

    def _calc_eps_from_fixed(self, fixed):
        eps, f_list = 0, [v for v, val in enumerate(fixed) if val == 1]
        for i, f in enumerate(f_list):
            eps = (eps + self.v4[f]) % 4
            for f2 in f_list[i + 1:]:
                if self.b4[f, f2]: eps = (eps + 2) % 4
        return eps

    def _eval_canonical_sum(self, vu: np.ndarray, nu: int) -> complex:
        s, p, phases = 1.0, 0, [1, 1j, -1, -1j]
        while p < self.rank:
            p_s = 0j
            for x1, x2 in [(0, 0), (0, 1), (1, 0), (1, 1)]:
                ph = ((2 if self.b_reduced[p, p + 1] and x1 and x2 else 0) + (vu[p] if x1 else 0) + (vu[p + 1] if x2 else 0)) % 4
                p_s += phases[ph]
            s *= p_s
            p += 2
        for k in range(self.rank, nu):
            v = vu[k] % 4
            if v == 0: s *= 2.0
            elif v == 1: s *= (1 + 1j)
            elif v == 2: return 0j
            elif v == 3: s *= (1 - 1j)
        return s

    def get_statevector_gray(self, weight, total_sv: np.ndarray, x: int = 0) -> np.ndarray:
        nu, dim = len(self.uvars_skeleton), 2 ** self.num_qubits
        phases, passthrough = [1, 1j, -1, -1j], self._passthrough_bits
        m_mat = np.zeros((nu, self.num_qubits), dtype=np.int8)
        for ui, u_o in enumerate(self.uvars_skeleton):
            for yi in range(self.num_qubits):
                if self.b4[u_o, self.output_vars[yi]]: m_mat[ui, yi] = 2
        for op in self.ops:
            if op.type == 'SWAP': m_mat[[op.a, op.b]] = m_mat[[op.b, op.a]]
            elif op.type == 'ADD': m_mat[op.b] = (m_mat[op.b] + m_mat[op.a]) % 4
        vu = np.zeros(nu, dtype=np.int8)
        for ui, u_o in enumerate(self.uvars_skeleton):
            vu[ui] = self.v4[u_o] % 4
            for xi in range(self.num_qubits):
                if (x >> xi) & 1 and self.b4[u_o, xi]: vu[ui] = (vu[ui] + 2) % 4
        for op in self.ops:
            if op.type == 'SWAP': vu[op.a], vu[op.b] = vu[op.b], vu[op.a]
            elif op.type == 'ADD': vu[op.b] = (vu[op.b] + vu[op.a]) % 4
        norm, gray = 2 ** (-self.num_h / 2), 0
        total_sv[0] += weight * self.get_amplitude(0, x)
        for i in range(1, dim):
            ng = i ^ (i >> 1)
            bit = (gray ^ ng).bit_length() - 1
            if (ng >> bit) & 1: vu = (vu + m_mat[:, bit]) % 4
            else: vu = (vu - m_mat[:, bit]) % 4
            gray = ng
            if gray & passthrough: continue
            total_sv[gray] += (weight * phases[self._calc_eps(gray) % 4] * self._eval_canonical_sum(vu, nu) * norm)
        return total_sv

    def _calc_eps(self, y):
        f_idxs, eps = [self.output_vars[j] for j in range(self.num_qubits) if (y >> j) & 1], 0
        for i, f in enumerate(f_idxs):
            eps = (eps + self.v4[f]) % 4
            for f2 in f_idxs[i + 1:]:
                if self.b4[f, f2]: eps = (eps + 2) % 4
        return eps

    def get_transition_matrix(self) -> np.ndarray:
        dim = 2 ** self.num_qubits
        tm = np.zeros((dim, dim), dtype=np.complex128)
        for x in range(dim):
            for y in range(dim): tm[y, x] = self.get_amplitude(y, x)
        return tm

class BranchNode:
    def __init__(self, weight, history=None, label="ROOT"):
        self.weight, self.history, self.label, self.children = weight, history if history is not None else [], label, []

class UniversalQC:
    def __init__(self, circuit: QuantumCircuit):
        self.circuit, self.num_qubits = circuit, circuit.num_qubits
        self.input_phases, self.output_phases, self.global_phase = np.zeros(self.num_qubits), np.zeros(self.num_qubits), 0.0
        self.root, self.engine, self._skeleton = None, None, None

    def build_tree(self):
        data = self.circuit.data
        f_h = [next((i for i, ins in enumerate(data) if ins.operation.name.lower() == 'h' and self.circuit.find_bit(ins.qubits[0]).index == q), None) for q in range(self.num_qubits)]
        l_h = [next((i for i in range(len(data)-1, -1, -1) if data[i].operation.name.lower() == 'h' and self.circuit.find_bit(data[i].qubits[0]).index == q), None) for q in range(self.num_qubits)]
        branching_data, self.global_phase = [], 0.0
        self.input_phases.fill(0.0); self.output_phases.fill(0.0)

        for i, ins in enumerate(data):
            name = ins.operation.name.lower()
            idxs = [self.circuit.find_bit(q).index for q in ins.qubits]
            if name in ['z', 's', 'sdg', 't', 'rz']:
                t = {'z': np.pi, 's': np.pi/2, 'sdg': -np.pi/2, 't': np.pi/4}.get(name, ins.operation.params[0] if name == 'rz' else 0)
                q0 = idxs[0]
                g_phase_offset = -t/2 if name == 'rz' else (np.pi/8 if name == 't' else 0)
                if f_h[q0] is None or i < f_h[q0]:
                    self.input_phases[q0] += t
                    self.global_phase += g_phase_offset
                elif i > l_h[q0]:
                    self.output_phases[q0] += t
                    self.global_phase += g_phase_offset
                else:
                    if name == 't': self.global_phase += np.pi / 8
                    branching_data.append(ins)
            else: branching_data.append(ins)

        def _grow(idx: int, node: BranchNode):
            if idx == len(branching_data): return
            ins = branching_data[idx]
            name, idxs = ins.operation.name.lower(), [self.circuit.find_bit(q).index for q in ins.qubits]
            if name in ['h', 'cz', 'z', 's', 'sdg', 'id']:
                node.history.append((name, idxs))
                _grow(idx + 1, node)
            elif name in ('rz', 't'):
                theta = np.pi/4 if name == 't' else ins.operation.params[0]
                c1 = BranchNode(node.weight * np.cos(theta/2), node.history + [('id', idxs)], "I")
                c2 = BranchNode(node.weight * (-1j * np.sin(theta/2)), node.history + [('z', idxs)], "Z")
                node.children += [c1, c2]; _grow(idx + 1, c1); _grow(idx + 1, c2)

        self.root = BranchNode(1.0 + 0j)
        _grow(0, self.root)
        self._skeleton = QuantumCircuit(self.num_qubits)
        for ins in self.circuit.data:
            if ins.operation.name.lower() in ['h', 'cz']: self._skeleton.append(ins)
        self.engine = DicksonEngine(self._skeleton)

    def get_statevector(self, x: int = 0) -> np.ndarray:
        if self.root is None: raise RuntimeError("Call build_tree() before get_statevector().")
        in_ph = sum(self.input_phases[q] for q in range(self.num_qubits) if (x >> q) & 1)
        sv = np.zeros(2 ** self.num_qubits, dtype=np.complex128)
        stack = [self.root]
        while stack:
            curr = stack.pop()
            if curr.children:
                stack.extend(curr.children); continue
            
            # Use fresh circuit, node.history already contains all Clifford gates
            b_qc = QuantumCircuit(self.num_qubits)
            for g, idxs in curr.history:
                if g != 'id': getattr(b_qc, g)(*idxs)
            self.engine.set_phases(b_qc)
            self.engine.get_statevector_gray(curr.weight, sv, x=x)
            
        final_sv = sv * np.exp(1j * (self.global_phase + in_ph))
        for q in range(self.num_qubits):
            if self.output_phases[q]:
                final_sv[(np.arange(2 ** self.num_qubits) >> q) & 1 == 1] *= np.exp(1j * self.output_phases[q])
        return final_sv

    def print_full_analytic_decomposition(self, transition_mode=False):
        if self.root is None: print("Tree not built."); return
        print(f"{'='*60}\nUNIVERSAL QC ANALYTIC DECOMPOSITION\n{'='*60}")
        p_terms = ([f"{self.global_phase:.6f}"] if abs(self.global_phase) > 1e-9 else []) + \
                  [f"({ph:.6f} * x_{q})" for q, ph in enumerate(self.input_phases) if abs(ph) > 1e-9] + \
                  [f"({ph:.6f} * y_{q})" for q, ph in enumerate(self.output_phases) if abs(ph) > 1e-9]
        print(f"Phase Dependence: phase = exp(i * [{' + '.join(p_terms) if p_terms else '0'}])\n")
        leaves, stack = [], [(self.root, "ROOT")]
        while stack:
            node, path = stack.pop()
            if not node.children: leaves.append((node, path))
            else:
                for i, child in enumerate(node.children): stack.append((child, f"{path} -> {child.label}"))
        
        for node, path in leaves:
            # Use fresh circuit, node.history already contains all Clifford gates
            b_qc = QuantumCircuit(self.num_qubits)
            for g, idxs in node.history:
                if g != 'id': getattr(b_qc, g)(*idxs)
            self.engine.set_phases(b_qc)
            self.engine.print_analytic_formula(transition_mode=transition_mode, weight=node.weight, branch_label=path)