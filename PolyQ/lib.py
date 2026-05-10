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
        eps, f_list = [0, [v for v, val in enumerate(fixed) if val == 1]]
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
    def __init__(self, weight, label="ROOT"):
        self.weight = weight
        self.label = label
        self.children = []
        # Pre-cached Z4 state for O(1) state generation
        self.v4_state: Optional[np.ndarray] = None 

class UniversalQC:
    def __init__(self, circuit: QuantumCircuit):
        self.circuit = circuit
        self.num_qubits = circuit.num_qubits
        self.global_phase = 0.0
        self.EL_gates = []
        self.core_gates = []
        self.ER_gates = []
        self.root = None
        self.engine = None
        self._skeleton = None

    def _split_circuit(self):
        """Splits circuit into E_L (left fringes), Core U, and E_R (right fringes)"""
        data = self.circuit.data
        n = self.num_qubits
        
        # 1. Forward pass: Collect E_L gates 
        in_EL = [True] * n
        EL_gates, core_gates = [], []
        
        for ins in data:
            name = ins.operation.name.lower()
            idxs = [self.circuit.find_bit(q).index for q in ins.qubits]
            
            if name in ['h', 'rx', 'ry']:
                for q in idxs: in_EL[q] = False
                core_gates.append(ins)
            elif name in ['cx', 'swap']:
                if in_EL[idxs[0]] and in_EL[idxs[1]]: EL_gates.append(ins)
                else:
                    for q in idxs: in_EL[q] = False
                    core_gates.append(ins)
            else:
                if all(in_EL[q] for q in idxs): EL_gates.append(ins)
                else: core_gates.append(ins)

        # 2. Reverse pass: Collect E_R gates from the core remainder
        in_ER = [True] * n
        ER_gates, final_core = [], []
        
        for ins in reversed(core_gates):
            name = ins.operation.name.lower()
            idxs = [self.circuit.find_bit(q).index for q in ins.qubits]
            
            if name in ['h', 'rx', 'ry']:
                for q in idxs: in_ER[q] = False
                final_core.insert(0, ins)
            elif name in ['cx', 'swap']:
                if in_ER[idxs[0]] and in_ER[idxs[1]]: ER_gates.insert(0, ins)
                else:
                    for q in idxs: in_ER[q] = False
                    final_core.insert(0, ins)
            else:
                if all(in_ER[q] for q in idxs): ER_gates.insert(0, ins)
                else: final_core.insert(0, ins)

        self.EL_gates = EL_gates
        self.core_gates = final_core
        self.ER_gates = ER_gates

    def _evaluate_EL(self, x: int) -> Tuple[int, complex]:
        """Evaluates Left fringe affine transformation to provide analytical start state"""
        bits = [(x >> i) & 1 for i in range(self.num_qubits)]
        phase = 0j
        for ins in self.EL_gates:
            name = ins.operation.name.lower()
            idxs = [self.circuit.find_bit(q).index for q in ins.qubits]
            
            if name == 'cx': bits[idxs[1]] ^= bits[idxs[0]]
            elif name == 'swap': bits[idxs[0]], bits[idxs[1]] = bits[idxs[1]], bits[idxs[0]]
            elif name == 'z' and bits[idxs[0]]: phase += np.pi
            elif name == 's' and bits[idxs[0]]: phase += np.pi / 2
            elif name == 'sdg' and bits[idxs[0]]: phase -= np.pi / 2
            elif name == 't' and bits[idxs[0]]: phase += np.pi / 4
            elif name == 'rz':
                theta = ins.operation.params[0]
                if bits[idxs[0]]: phase += theta / 2
                else: phase -= theta / 2
                
        new_x = sum(b << i for i, b in enumerate(bits))
        return new_x, np.exp(1j * phase)

    def _apply_ER(self, sv: np.ndarray):
        """Applies Right fringe transformations in-place analytically to avoid core inflation"""
        idx = np.arange(len(sv))
        for ins in self.ER_gates:
            name = ins.operation.name.lower()
            idxs = [self.circuit.find_bit(q).index for q in ins.qubits]
            
            if name == 'cx':
                c, t = idxs[0], idxs[1]
                idx_c1 = idx[(idx & (1 << c)) != 0]
                idx_c1_t0 = idx_c1[(idx_c1 & (1 << t)) == 0]
                idx_c1_t1 = idx_c1_t0 | (1 << t)
                sv[idx_c1_t0], sv[idx_c1_t1] = sv[idx_c1_t1], sv[idx_c1_t0]
            elif name == 'swap':
                b0, b1 = (idx >> idxs[0]) & 1, (idx >> idxs[1]) & 1
                diff = b0 != b1
                swapped_idx = idx.copy()
                swapped_idx[diff] ^= (1 << idxs[0]) | (1 << idxs[1])
                sv[:] = sv[swapped_idx]
            elif name in ['z', 's', 'sdg', 't', 'rz']:
                mask = 1 << idxs[0]
                idx_1 = (idx & mask) != 0
                idx_0 = ~idx_1
                if name == 'z': sv[idx_1] *= -1
                elif name == 's': sv[idx_1] *= 1j
                elif name == 'sdg': sv[idx_1] *= -1j
                elif name == 't': sv[idx_1] *= np.exp(1j * np.pi / 4)
                elif name == 'rz':
                    theta = ins.operation.params[0]
                    sv[idx_1] *= np.exp(1j * theta / 2)
                    sv[idx_0] *= np.exp(-1j * theta / 2)

    def build_tree(self):
        self._split_circuit()
        
        # Calculate Phase Accumulation
        core_t_count = sum(1 for ins in self.core_gates if ins.operation.name.lower() == 't')
        self.global_phase = core_t_count * (np.pi / 8)

        # Transpile Core skeleton: Support inner CX and SWAP natively in b_reduced form
        self._skeleton = QuantumCircuit(self.num_qubits)
        for ins in self.core_gates:
            name = ins.operation.name.lower()
            idxs = [self.circuit.find_bit(q).index for q in ins.qubits]
            if name in ['h', 'cz']: self._skeleton.append(ins)
            elif name == 'cx':
                self._skeleton.h(idxs[1]); self._skeleton.cz(idxs[0], idxs[1]); self._skeleton.h(idxs[1])
            elif name == 'swap':
                a, b = idxs[0], idxs[1]
                self._skeleton.h(b); self._skeleton.cz(a, b); self._skeleton.h(b)
                self._skeleton.h(a); self._skeleton.cz(b, a); self._skeleton.h(a)
                self._skeleton.h(b); self._skeleton.cz(a, b); self._skeleton.h(b)

        self.engine = DicksonEngine(self._skeleton)

        # Pre-calculate Variable Tracking for O(1) Linear lookups
        gate_to_vars = []
        wires = [[i] for i in range(self.num_qubits)]
        nv = self.num_qubits
        
        for ins in self.core_gates:
            name = ins.operation.name.lower()
            idxs = [self.circuit.find_bit(q).index for q in ins.qubits]
            gate_to_vars.append([wires[q][-1] for q in idxs])
            
            if name == 'h':
                wires[idxs[0]].append(nv); nv += 1
            elif name == 'cx':
                wires[idxs[1]].append(nv); nv += 1
                wires[idxs[1]].append(nv); nv += 1
            elif name == 'swap':
                a, b = idxs[0], idxs[1]
                wires[b].append(nv); nv += 1; wires[b].append(nv); nv += 1
                wires[a].append(nv); nv += 1; wires[a].append(nv); nv += 1
                wires[b].append(nv); nv += 1; wires[b].append(nv); nv += 1

        self.root = BranchNode(1.0 + 0j)
        self.root.v4_state = np.zeros(self.engine.n_vars, dtype=np.int8)

        def _grow(idx: int, node: BranchNode):
            if idx == len(self.core_gates): return
            ins = self.core_gates[idx]
            name = ins.operation.name.lower()
            q_vars = gate_to_vars[idx]

            if name in ['h', 'cz', 'cx', 'swap', 'id']:
                _grow(idx + 1, node)
            elif name == 'z':
                node.v4_state[q_vars[0]] = (node.v4_state[q_vars[0]] + 2) % 4
                _grow(idx + 1, node)
            elif name == 's':
                node.v4_state[q_vars[0]] = (node.v4_state[q_vars[0]] + 1) % 4
                _grow(idx + 1, node)
            elif name == 'sdg':
                node.v4_state[q_vars[0]] = (node.v4_state[q_vars[0]] + 3) % 4
                _grow(idx + 1, node)
            elif name in ('rz', 't'):
                theta = np.pi/4 if name == 't' else ins.operation.params[0]
                c1 = BranchNode(node.weight * np.cos(theta/2), label=node.label+"_I")
                c1.v4_state = node.v4_state.copy()
                
                c2 = BranchNode(node.weight * (-1j * np.sin(theta/2)), label=node.label+"_Z")
                c2.v4_state = node.v4_state.copy()
                c2.v4_state[q_vars[0]] = (c2.v4_state[q_vars[0]] + 2) % 4

                node.children.extend([c1, c2])
                _grow(idx + 1, c1)
                _grow(idx + 1, c2)

        _grow(0, self.root)

    def get_statevector(self, x: int = 0) -> np.ndarray:
        if self.root is None: raise RuntimeError("Call build_tree() before get_statevector().")
        
        # Step 1: Pre-process E_L Transform
        x_core, input_phase = self._evaluate_EL(x)
        sv = np.zeros(2 ** self.num_qubits, dtype=np.complex128)
        
        stack = [self.root]
        while stack:
            curr = stack.pop()
            if curr.children:
                stack.extend(curr.children); continue
            
            # Fast O(1) linear term transfer (No QuantumCircuit reconstruction)
            self.engine.v4 = curr.v4_state.copy()
            self.engine.get_statevector_gray(curr.weight, sv, x=x_core)
            
        final_sv = sv * np.exp(1j * self.global_phase) * input_phase
        
        # Step 2: Post-process E_R Permutations
        self._apply_ER(final_sv)
        
        return final_sv

    def print_full_analytic_decomposition(self, transition_mode=False):
        if self.root is None: print("Tree not built."); return
        print(f"{'='*60}\nUNIVERSAL QC ANALYTIC DECOMPOSITION\n{'='*60}")
        print(f"E_L Gates: {len(self.EL_gates)} (Input Affine Transform & Phases)")
        print(f"E_R Gates: {len(self.ER_gates)} (Output Affine Permutations)")
        print(f"Core Global Phase from T-gates: {self.global_phase:.6f}\n")
        
        leaves, stack = [], [(self.root, "ROOT")]
        while stack:
            node, path = stack.pop()
            if not node.children: leaves.append((node, path))
            else:
                for i, child in enumerate(node.children): stack.append((child, f"{path} -> {child.label}"))
        
        for node, path in leaves:
            self.engine.v4 = node.v4_state.copy()
            self.engine.print_analytic_formula(transition_mode=transition_mode, weight=node.weight, branch_label=path)