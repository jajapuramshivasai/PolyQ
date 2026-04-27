import numpy as np
from typing import List, Dict, Set, Tuple, Optional, Any
from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector

class DicksonOp:
    """Represents an F2 linear transformation step during Dickson decomposition."""
    def __init__(self, op_type: str, a: int, b: Optional[int] = None):
        self.type = op_type
        self.a = a  # Pivot for ADD, first index for SWAP
        self.b = b  # Target for ADD, second index for SWAP

class QC:
    """
    A Quantum Circuit simulator using Dickson decomposition for F2 quadratic forms.
    Supports H, Z, S, and CZ gates via Qiskit QuantumCircuit.
    """
    def __init__(self, circuit: QuantumCircuit):
        self.circuit = circuit
        self.num_qubits = circuit.num_qubits
        self.gates: List[Tuple[Any, ...]] = []
        self.compiled = False
        
        # Internal state for the graph-based representation
        self.b4: Optional[np.ndarray] = None
        self.v4: Optional[np.ndarray] = None
        self.output_vars: List[int] = []
        self.num_h = 0
        self.n_vars = 0
        self.rank = 0
        self.uvars_skeleton: List[int] = []
        self.ops: List[DicksonOp] = []
        self.b_reduced: Optional[np.ndarray] = None
        
        self._translate_qiskit_circuit()

    def _translate_qiskit_circuit(self) -> None:
        """Extracts supported gates from the Qiskit circuit."""
        for instruction in self.circuit.data:
            gate = instruction.operation
            qargs = instruction.qubits
            # Map qubit objects to their integer indices in the circuit
            indices = [self.circuit.find_bit(q).index for q in qargs]
            
            name = gate.name.lower()
            if name == 'h':
                self.gates.append(('H', indices[0]))
            elif name == 'z':
                self.gates.append(('Z', indices[0]))
            elif name == 's':
                self.gates.append(('S', indices[0]))
            elif name == 'sdg':
                self.gates.append(('SDG', indices[0]))
            elif name == 'cz':
                self.gates.append(('CZ', indices[0], indices[1]))
            elif name == 'barrier':
                continue
            else:
                raise ValueError(f"Gate '{name}' is not supported by this simulator. "
                                 "Supported gates: H, Z, S, Sdg, CZ.")

    def compile(self) -> None:
        """Translates the gate sequence into an adjacency matrix over F2."""
        wires = [[i] for i in range(self.num_qubits)]
        next_var = self.num_qubits
        self.num_h = 0
        b4_adj: Dict[int, Set[int]] = {}
        v4_dict: Dict[int, int] = {}
        
        def ensure_var(v: int):
            b4_adj.setdefault(v, set())
            v4_dict.setdefault(v, 0)

        for i in range(self.num_qubits): ensure_var(i)

        for gate in self.gates:
            g_type = gate[0]
            if g_type == 'H':
                q = gate[1]
                prev, cur = wires[q][-1], next_var
                next_var += 1
                self.num_h += 1
                wires[q].append(cur)
                ensure_var(cur)
                b4_adj[prev].add(cur); b4_adj[cur].add(prev)
            elif g_type == 'CZ':
                v1, v2 = wires[gate[1]][-1], wires[gate[2]][-1]
                b4_adj[v1].add(v2); b4_adj[v2].add(v1)
            elif g_type == 'Z':
                v = wires[gate[1]][-1]
                v4_dict[v] = (v4_dict[v] + 2) % 4
            elif g_type == 'S':
                v = wires[gate[1]][-1]
                v4_dict[v] = (v4_dict[v] + 1) % 4
            elif g_type == 'SDG':
                v = wires[gate[1]][-1]
                v4_dict[v] = (v4_dict[v] + 3) % 4

        self.n_vars = next_var
        self.b4 = np.zeros((self.n_vars, self.n_vars), dtype=np.int8)
        self.v4 = np.zeros(self.n_vars, dtype=np.int8)
        
        for i in range(self.n_vars):
            self.v4[i] = v4_dict.get(i, 0)
            for j in b4_adj.get(i, []): self.b4[i, j] = 1

        self.output_vars = [wires[q][-1] for q in range(self.num_qubits)]
        self.uvars_skeleton = [i for i in range(self.n_vars) if i >= self.num_qubits and i not in self.output_vars]
        nu = len(self.uvars_skeleton)
        
        b_u_skel = np.zeros((nu, nu), dtype=np.int8)
        for ui, orig_u in enumerate(self.uvars_skeleton):
            for uj, orig_uj in enumerate(self.uvars_skeleton):
                 if ui < uj and self.b4[orig_u, orig_uj]:
                     b_u_skel[ui, uj] = b_u_skel[uj, ui] = 1
                     
        self.ops, self.b_reduced, self.rank = self._plan_dickson_f2(b_u_skel, nu)
        self.compiled = True

    def _plan_dickson_f2(self, b_matrix: np.ndarray, n_vars: int) -> Tuple[List[DicksonOp], np.ndarray, int]:
        """Performs Dickson's decomposition to simplify the quadratic form."""
        b, ops, r, p = np.copy(b_matrix), [], 0, 0
        while p + 1 < n_vars:
            pivot = next(((i, j) for i in range(p, n_vars) for j in range(i + 1, n_vars) if b[i, j] == 1), None)
            if not pivot: break
            i, j = pivot
            
            if i != p:
                b[[p, i]] = b[[i, p]]
                b[:, [p, i]] = b[:, [i, p]]
                ops.append(DicksonOp('SWAP', p, i))
                
            j_act = i if j == p else j
            if j_act != p + 1:
                b[[p + 1, j_act]] = b[[j_act, p + 1]]
                b[:, [p + 1, j_act]] = b[:, [j_act, p + 1]]
                ops.append(DicksonOp('SWAP', p + 1, j_act))
                
            rp, rp1 = np.copy(b[p, :]), np.copy(b[p + 1, :])
            
            for k in range(p + 2, n_vars):
                if b[k, p] == 1:
                    b[k, :] ^= rp1
                    b[:, k] ^= rp1
                    ops.append(DicksonOp('ADD', p + 1, k))
                if b[k, p + 1] == 1:
                    b[k, :] ^= rp
                    b[:, k] ^= rp
                    ops.append(DicksonOp('ADD', p, k))
            r += 2; p += 2
        return ops, b, r

    def get_amplitude(self, y_val: int, x_val: int = 0) -> complex:
        """Evaluates the transition amplitude <y|U|x>."""
        if not self.compiled: raise RuntimeError("Circuit must be compiled first.")
            
        fixed: List[Optional[int]] = [None] * self.n_vars
        for i in range(self.num_qubits): fixed[i] = (x_val >> i) & 1
            
        for i in range(self.num_qubits):
            bit, ov = (y_val >> i) & 1, self.output_vars[i]
            if fixed[ov] is None: fixed[ov] = bit
            elif fixed[ov] != bit: return 0j
                
        eps_base = 0
        f_list = [v for v, val in enumerate(fixed) if val == 1]
        
        for f in f_list:
            eps_base = (eps_base + self.v4[f]) % 4
            for f2 in f_list:
                if f < f2 and self.b4[f, f2]: eps_base = (eps_base + 2) % 4

        nu = len(self.uvars_skeleton)
        vu_base = np.zeros(nu, dtype=np.int8)
        
        for ui, orig_u in enumerate(self.uvars_skeleton):
            vu_base[ui] = self.v4[orig_u] % 4
            for f in f_list:
                if self.b4[orig_u, f]: vu_base[ui] = (vu_base[ui] + 2) % 4
                    
        for op in self.ops:
            if op.type == 'SWAP': 
                vu_base[op.a], vu_base[op.b] = vu_base[op.b], vu_base[op.a]
            elif op.type == 'ADD': 
                vu_base[op.b] = (vu_base[op.b] + vu_base[op.a]) % 4

        phase_map = {0: 1.0+0j, 1: 1j, 2: -1.0+0j, 3: -1j}
        return phase_map[eps_base] * self._eval_canonical_sum(vu_base, nu) * (2.0 ** (-self.num_h / 2.0))

    def _eval_canonical_sum(self, vu_base: np.ndarray, nu: int) -> complex:
        sum_val, p = 1.0 + 0j, 0
        while p < self.rank:
            has_edge = self.b_reduced[p, p + 1]
            pair_sum = 0j
            for x1 in (0, 1):
                for x2 in (0, 1):
                    ph = ((2 if has_edge and x1 and x2 else 0) + (vu_base[p] if x1 else 0) + (vu_base[p + 1] if x2 else 0)) % 4
                    if ph == 0: pair_sum += 1
                    elif ph == 1: pair_sum += 1j
                    elif ph == 2: pair_sum += -1
                    elif ph == 3: pair_sum += -1j
            sum_val *= pair_sum
            p += 2
            
        for k in range(self.rank, nu):
            val = vu_base[k] % 4
            if val == 0: sum_val *= 2.0
            elif val == 1: sum_val *= (1.0 + 1j)
            elif val == 2: return 0j 
            elif val == 3: sum_val *= (1.0 - 1j)
        return sum_val

    def get_transition_matrix(self) -> np.ndarray:
        dim = 2 ** self.num_qubits
        return np.array([[self.get_amplitude(y, x) for x in range(dim)] for y in range(dim)], dtype=np.complex128)
    
    def print_circuit_parameters(self):
        if not self.compiled: raise RuntimeError("Circuit must be compiled first.")
        nu = len(self.uvars_skeleton)
        m = self.rank // 2
        print(f"\n{'='*60}\n   CIRCUIT METRICS & SIMULATION PARAMETERS\n{'='*60}")
        print(f"Total Qubits (n)           : {self.num_qubits}")
        print(f"Total F2 Variables         : {self.n_vars}")
        print(f"Hadamard Count (h)         : {self.num_h}")
        print(f"Normalization Factor       : 2^(-{self.num_h}/2)")
        print(f"Unbound Internal Vars (nu) : {nu}")
        print(f"Dickson Rank (2m)          : {self.rank} (m = {m})")
        print(f"Kernel Size (Zero Conds.)  : {nu - self.rank}")
        print(f"{'='*60}\n")
        
    def print_analytic_formula(self, transition_mode=True):
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
            else:
                terms += [f"y_{yi}" for yi in range(self.num_qubits) if coeff_matrix[j, yi] == 2]
            return " XOR ".join(terms) if terms else "0"

        title = "<y|U|x>" if transition_mode else "Statevector <y|U|0>"
        print(f"\n{'='*60}\n   ANALYTICAL FUNCTIONS FOR {title}\n{'='*60}")
        n1_terms = [f"({build_xor_expr(j)})" for j in range(self.rank)]
        print(f"1. EXPONENTIAL SUM VARIABLES\n   N1 = {' + '.join(n1_terms) if n1_terms else '0'}\n   N0 = {self.rank} - N1\n")
        zero_conditions = [build_xor_expr(k) for k in range(self.rank, nu)]
        print("2. ZERO AMPLITUDE CONDITIONS (Balanced Function)")
        if zero_conditions:
            for idx, cond in enumerate(zero_conditions): print(f"   Kernel_{idx}: ({cond}) == 1")
        else: print("   No kernel variables exist.")
        print("="*60 + "\n")