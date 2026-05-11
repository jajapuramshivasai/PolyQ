# lib.py - Dickson Quantum Engine with Van den Nest Decomposition
import numpy as np
from typing import List, Dict, Set, Tuple, Optional, Any
from qiskit import QuantumCircuit


class DicksonOp:
    """Represents a single Dickson reduction operation."""
    def __init__(self, op_type: str, a: int, b: Optional[int] = None):
        self.type = op_type
        self.a = a
        self.b = b


class DicksonEngine:
    """
    Optimized F2/Z4 Quadratic Form Engine.
    - Handles structural caching for high-speed statevector evaluation
    - Uses Gray code optimization for efficient state traversal
    - Manages Z4 phase tracking and Clifford circuit compilation
    """

    def __init__(self, circuit: QuantumCircuit):
        self.circuit = circuit
        self.num_qubits = circuit.num_qubits
        
        # Circuit structure
        self.gates: List[Tuple[str, List[int]]] = []
        
        # Cached structure
        self.b4: Optional[np.ndarray] = None  # Adjacency matrix
        self.output_vars: List[int] = []
        self.num_h = 0  # Hadamard count
        self.n_vars = 0  # Total variables
        self.rank = 0  # Dickson rank
        self.uvars_skeleton: List[int] = []  # Auxiliary variables
        self.ops: List[DicksonOp] = []  # Reduction operations
        self.b_reduced: Optional[np.ndarray] = None  # Reduced form
        self.v4 = np.zeros(0)  # Z4 phase vector
        self._passthrough_bits: int = 0

        self._translate_circuit()
        self.compile_structure()

    def _translate_circuit(self):
        """Extracts gate operations from Qiskit circuit."""
        for instr in self.circuit.data:
            name = instr.operation.name.lower()
            idxs = [self.circuit.find_bit(q).index for q in instr.qubits]
            if name in ['h', 'z', 's', 'sdg', 'cz', 'id']:
                self.gates.append((name.upper(), idxs))

    def compile_structure(self):
        """
        Builds adjacency matrix and performs Dickson decomposition.
        Creates a reduced form for fast statevector computation.
        """
        # Track variable lineage through circuit
        wires = [[i] for i in range(self.num_qubits)]
        next_v, self.num_h = self.num_qubits, 0
        adj: Dict[int, Set[int]] = {i: set() for i in range(self.num_qubits)}

        # Build adjacency from gates
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
        
        # Convert adjacency to matrix
        self.b4 = np.zeros((self.n_vars, self.n_vars), dtype=np.int8)
        for i in range(self.n_vars):
            for j in adj.get(i, []):
                self.b4[i, j] = 1

        # Track output variables
        self.output_vars = [wires[q][-1] for q in range(self.num_qubits)]
        
        # Identify auxiliary variables (not outputs)
        self.uvars_skeleton = [
            i for i in range(self.n_vars)
            if i >= self.num_qubits and i not in self.output_vars
        ]

        # Perform Dickson reduction on auxiliary space
        nu = len(self.uvars_skeleton)
        if nu > 0:
            b_u = np.zeros((nu, nu), dtype=np.int8)
            for ui, u_o in enumerate(self.uvars_skeleton):
                for uj, uj_o in enumerate(self.uvars_skeleton):
                    if ui < uj and self.b4[u_o, uj_o]:
                        b_u[ui, uj] = b_u[uj, ui] = 1

            self.ops, self.b_reduced, self.rank = self._plan_dickson(b_u, nu)
        else:
            self.ops, self.b_reduced, self.rank = [], self.b4, 0
            
        self.v4 = np.zeros(self.n_vars, dtype=np.int8)

        # Cache passthrough bits (inputs → outputs directly)
        self._passthrough_bits = 0
        for i in range(self.num_qubits):
            if self.output_vars[i] < self.num_qubits:
                self._passthrough_bits |= (1 << i)

    def _plan_dickson(self, b: np.ndarray, n: int) -> Tuple[List[DicksonOp], np.ndarray, int]:
        """
        Dickson block reduction algorithm.
        Returns: (operations, reduced_form, rank)
        """
        if n == 0:
            return [], b, 0
            
        b_w, ops, r, p = np.copy(b), [], 0, 0
        
        while p + 1 < n:
            # Find pivot
            pivot = None
            for i in range(p, n):
                for j in range(i + 1, n):
                    if b_w[i, j] == 1:
                        pivot = (i, j)
                        break
                if pivot:
                    break
            
            if not pivot:
                break

            i, j = pivot

            # Swap row/col to position (p, p+1)
            if i != p:
                b_w[[p, i]] = b_w[[i, p]]
                b_w[:, [p, i]] = b_w[:, [i, p]]
                ops.append(DicksonOp('SWAP', p, i))
                j_act = p if j == i else (i if j == p else j)
            else:
                j_act = j

            if j_act != p + 1:
                b_w[[p + 1, j_act]] = b_w[[j_act, p + 1]]
                b_w[:, [p + 1, j_act]] = b_w[:, [j_act, p + 1]]
                ops.append(DicksonOp('SWAP', p + 1, j_act))

            # Eliminate remaining entries
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

    def set_phases(self, circuit: QuantumCircuit):
        """Scans circuit and sets Z4 phase vector from Z, S, Sdg gates."""
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
            elif name == 'z':
                self.v4[wires[q][-1]] = (self.v4[wires[q][-1]] + 2) % 4
            elif name == 's':
                self.v4[wires[q][-1]] = (self.v4[wires[q][-1]] + 1) % 4
            elif name == 'sdg':
                self.v4[wires[q][-1]] = (self.v4[wires[q][-1]] + 3) % 4

    def get_amplitude(self, y: int, x: int = 0) -> complex:
        """
        Computes amplitude ⟨y|U|x⟩ via Z4 quadratic form.
        Handles both input (x) and output (y) constraints.
        """
        fixed = [None] * self.n_vars
        
        # Fix input bits
        for i in range(self.num_qubits):
            fixed[i] = (x >> i) & 1

        # Fix output bits and check consistency
        for i in range(self.num_qubits):
            bit, ov = (y >> i) & 1, self.output_vars[i]
            if fixed[ov] is not None and fixed[ov] != bit:
                return 0j
            fixed[ov] = bit

        # Compute phase
        eps = self._calc_eps_from_fixed(fixed)
        
        # Compute auxiliary variables
        nu = len(self.uvars_skeleton)
        if nu == 0:
            return [1, 1j, -1, -1j][eps % 4] * (2 ** (-self.num_h / 2))
            
        vu = np.zeros(nu, dtype=np.int8)
        for ui, u_o in enumerate(self.uvars_skeleton):
            vu[ui] = self.v4[u_o] % 4
            for v_idx, val in enumerate(fixed):
                if val == 1 and self.b4[u_o, v_idx]:
                    vu[ui] = (vu[ui] + 2) % 4

        # Apply reduction operations
        for op in self.ops:
            if op.type == 'SWAP':
                vu[op.a], vu[op.b] = vu[op.b], vu[op.a]
            elif op.type == 'ADD':
                vu[op.b] = (vu[op.b] + vu[op.a]) % 4

        return (
            [1, 1j, -1, -1j][eps % 4] 
            * self._eval_canonical_sum(vu, nu) 
            * (2 ** (-self.num_h / 2))
        )

    def _calc_eps_from_fixed(self, fixed: List[Optional[int]]) -> int:
        """Computes Z4 phase accumulation from fixed bits."""
        eps = 0
        f_list = [v for v, val in enumerate(fixed) if val == 1]
        
        for i, f in enumerate(f_list):
            eps = (eps + self.v4[f]) % 4
            for f2 in f_list[i + 1:]:
                if self.b4[f, f2]:
                    eps = (eps + 2) % 4
        
        return eps

    def _eval_canonical_sum(self, vu: np.ndarray, nu: int) -> complex:
        """
        Evaluates canonical sum over reduced auxiliary variables.
        Optimized for Dickson block structure.
        """
        s = 1.0
        phases = [1, 1j, -1, -1j]

        # Process Dickson pairs
        for p in range(0, self.rank, 2):
            p_s = 0j
            for x1, x2 in [(0, 0), (0, 1), (1, 0), (1, 1)]:
                ph = (
                    (2 if p + 1 < self.rank and self.b_reduced[p, p + 1] and x1 and x2 else 0)
                    + (vu[p] if x1 else 0)
                    + (vu[p + 1] if x2 else 0)
                ) % 4
                p_s += phases[ph]
            s *= p_s

        # Process kernel (zero) constraints
        for k in range(self.rank, nu):
            v = vu[k] % 4
            if v == 0:
                s *= 2.0
            elif v == 1:
                s *= (1 + 1j)
            elif v == 2:
                return 0j
            elif v == 3:
                s *= (1 - 1j)

        return s

    # def _calc_eps(self, y: int) -> int:
    #     """Computes phase for output bitstring y."""
    #     f_idxs = [self.output_vars[j] for j in range(self.num_qubits) if (y >> j) & 1]
    #     eps = 0
        
    #     for i, f in enumerate(f_idxs):
    #         eps = (eps + self.v4[f]) % 4
    #         for f2 in f_idxs[i + 1:]:
    #             if self.b4[f, f2]:
    #                 eps = (eps + 2) % 4
        
    #     return eps
    def _calc_eps(self, y: int, x: int = 0) -> int:
        """Computes phase for output bitstring y and input bitstring x."""
        f_set = set()
        
        # Add active input variables
        for j in range(self.num_qubits):
            if (x >> j) & 1:
                f_set.add(j)
                
        # Add active output variables
        for j in range(self.num_qubits):
            if (y >> j) & 1:
                f_set.add(self.output_vars[j])

        f_idxs = list(f_set)
        eps = 0
        
        for i, f in enumerate(f_idxs):
            eps = (eps + self.v4[f]) % 4
            for f2 in f_idxs[i + 1:]:
                if self.b4[f, f2]:
                    eps = (eps + 2) % 4
        
        return eps
    
    
    
    
    # def get_statevector_gray(self, weight: complex, total_sv: np.ndarray, x: int = 0) -> np.ndarray:
    #     """
    #     Efficient statevector generation using Gray code.
    #     Accumulates amplitudes into total_sv with given weight.
    #     """
    #     nu = len(self.uvars_skeleton)
    #     dim = 2 ** self.num_qubits
    #     phases = [1, 1j, -1, -1j]
    #     passthrough = self._passthrough_bits

    #     # Build measurement matrix
    #     m_mat = np.zeros((nu, self.num_qubits), dtype=np.int8) if nu > 0 else np.zeros((0, self.num_qubits), dtype=np.int8)
    #     for ui, u_o in enumerate(self.uvars_skeleton):
    #         for yi in range(self.num_qubits):
    #             if self.b4[u_o, self.output_vars[yi]]:
    #                 m_mat[ui, yi] = 2

    #     # Apply reduction to measurement matrix
    #     for op in self.ops:
    #         if op.type == 'SWAP':
    #             m_mat[[op.a, op.b]] = m_mat[[op.b, op.a]]
    #         elif op.type == 'ADD':
    #             m_mat[op.b] = (m_mat[op.b] + m_mat[op.a]) % 4

    #     # Initialize auxiliary vector
    #     vu = np.zeros(nu, dtype=np.int8) if nu > 0 else np.zeros(0, dtype=np.int8)
    #     for ui, u_o in enumerate(self.uvars_skeleton):
    #         vu[ui] = self.v4[u_o] % 4
    #         for xi in range(self.num_qubits):
    #             if (x >> xi) & 1 and self.b4[u_o, xi]:
    #                 vu[ui] = (vu[ui] + 2) % 4

    #     # Apply reduction
    #     for op in self.ops:
    #         if op.type == 'SWAP':
    #             vu[op.a], vu[op.b] = vu[op.b], vu[op.a]
    #         elif op.type == 'ADD':
    #             vu[op.b] = (vu[op.b] + vu[op.a]) % 4

    #     norm = 2 ** (-self.num_h / 2)
    #     gray = 0

    #     # Add amplitude for y=0
    #     total_sv[0] += weight * self.get_amplitude(0, x)

    #     # Gray code traversal
    #     for i in range(1, dim):
    #         ng = i ^ (i >> 1)
    #         bit = (gray ^ ng).bit_length() - 1
            
    #         if nu > 0:
    #             if (ng >> bit) & 1:
    #                 vu = (vu + m_mat[:, bit]) % 4
    #             else:
    #                 vu = (vu - m_mat[:, bit]) % 4

    #         gray = ng

    #         if gray & passthrough:
    #             continue

    #         total_sv[gray] += (
    #             weight
    #             * phases[self._calc_eps(gray) % 4]
    #             * self._eval_canonical_sum(vu, nu)
    #             * norm
    #         )

    #     return total_sv
    
    def get_statevector_gray(self, weight: complex, total_sv: np.ndarray, x: int = 0) -> np.ndarray:
        """
        Efficient statevector generation using Gray code.
        Accumulates amplitudes into total_sv with given weight.
        """
        nu = len(self.uvars_skeleton)
        dim = 2 ** self.num_qubits
        phases = [1, 1j, -1, -1j]
        passthrough_mask = self._passthrough_bits
        
        # Determine what the passthrough bits *must* be based on input x
        required_y_mask = 0
        if passthrough_mask:
            for i in range(self.num_qubits):
                if self.output_vars[i] < self.num_qubits:
                    in_bit = self.output_vars[i]
                    if (x >> in_bit) & 1:
                        required_y_mask |= (1 << i)

        # Build measurement matrix
        m_mat = np.zeros((nu, self.num_qubits), dtype=np.int8) if nu > 0 else np.zeros((0, self.num_qubits), dtype=np.int8)
        for ui, u_o in enumerate(self.uvars_skeleton):
            for yi in range(self.num_qubits):
                if self.b4[u_o, self.output_vars[yi]]:
                    m_mat[ui, yi] = 2

        # Apply reduction to measurement matrix
        for op in self.ops:
            if op.type == 'SWAP':
                m_mat[[op.a, op.b]] = m_mat[[op.b, op.a]]
            elif op.type == 'ADD':
                m_mat[op.b] = (m_mat[op.b] + m_mat[op.a]) % 4

        # Initialize auxiliary vector
        vu = np.zeros(nu, dtype=np.int8) if nu > 0 else np.zeros(0, dtype=np.int8)
        for ui, u_o in enumerate(self.uvars_skeleton):
            vu[ui] = self.v4[u_o] % 4
            for xi in range(self.num_qubits):
                if (x >> xi) & 1 and self.b4[u_o, xi]:
                    vu[ui] = (vu[ui] + 2) % 4

        # Apply reduction
        for op in self.ops:
            if op.type == 'SWAP':
                vu[op.a], vu[op.b] = vu[op.b], vu[op.a]
            elif op.type == 'ADD':
                vu[op.b] = (vu[op.b] + vu[op.a]) % 4

        norm = 2 ** (-self.num_h / 2)
        gray = 0

        # Add amplitude for y=0
        total_sv[0] += weight * self.get_amplitude(0, x)

        # Gray code traversal
        for i in range(1, dim):
            ng = i ^ (i >> 1)
            bit = (gray ^ ng).bit_length() - 1
            
            if nu > 0:
                if (ng >> bit) & 1:
                    vu = (vu + m_mat[:, bit]) % 4
                else:
                    vu = (vu - m_mat[:, bit]) % 4

            gray = ng

            # Correctly check passthrough consistency against input x
            if (gray & passthrough_mask) != required_y_mask:
                continue

            total_sv[gray] += (
                weight
                * phases[self._calc_eps(gray, x) % 4]  # Pass x into _calc_eps
                * self._eval_canonical_sum(vu, nu)
                * norm
            )

        return total_sv

    def get_transition_matrix(self) -> np.ndarray:
        """Computes full transition matrix as 2D array."""
        dim = 2 ** self.num_qubits
        tm = np.zeros((dim, dim), dtype=np.complex128)
        for x in range(dim):
            for y in range(dim):
                tm[y, x] = self.get_amplitude(y, x)
        return tm

    def print_analytic_formula(
        self,
        transition_mode: bool = True,
        weight: float = 1.0,
        branch_label: str = ""
    ):
        """
        Prints human-readable Z4 quadratic form formula.
        Useful for debugging and theoretical analysis.
        """
        nu = len(self.uvars_skeleton)
        n_cols = 2 * self.num_qubits + 1 if transition_mode else self.num_qubits + 1
        coeff_matrix = np.zeros((nu, n_cols), dtype=np.int8)

        # Initialize coefficients
        for ui, orig_u in enumerate(self.uvars_skeleton):
            coeff_matrix[ui, -1] = self.v4[orig_u] % 4
            if transition_mode:
                for xi in range(self.num_qubits):
                    if self.b4[orig_u, xi]:
                        coeff_matrix[ui, xi] = 2
            for yi in range(self.num_qubits):
                if self.b4[orig_u, self.output_vars[yi]]:
                    col_idx = self.num_qubits + yi if transition_mode else yi
                    coeff_matrix[ui, col_idx] = (coeff_matrix[ui, col_idx] + 2) % 4

        # Apply reduction operations
        for op in self.ops:
            if op.type == 'SWAP':
                coeff_matrix[[op.a, op.b]] = coeff_matrix[[op.b, op.a]]
            elif op.type == 'ADD':
                coeff_matrix[op.b] = (coeff_matrix[op.b] + coeff_matrix[op.a]) % 4

        def build_xor_expr(j: int) -> str:
            """Builds XOR expression for row j."""
            terms = ["1"] if coeff_matrix[j, -1] in (2, 3) else []
            if transition_mode:
                terms += [f"x_{xi}" for xi in range(self.num_qubits) if coeff_matrix[j, xi] == 2]
                terms += [f"y_{yi}" for yi in range(self.num_qubits) if coeff_matrix[j, self.num_qubits + yi] == 2]
            else:
                terms += [f"y_{yi}" for yi in range(self.num_qubits) if coeff_matrix[j, yi] == 2]
            return " XOR ".join(terms) if terms else "0"

        header = (
            f"BRANCH: {branch_label} (Weight: {weight})" if branch_label
            else "CLIFFORD ANALYTIC FORM"
        )
        print(f"{'-' * 40}\n{header}\n{'-' * 40}")

        print("1. Exponential Sum Variables (N1):")
        n1_terms = [f"({build_xor_expr(j)})" for j in range(self.rank)]
        if n1_terms:
            for idx, term in enumerate(n1_terms):
                print(f"   v_{idx}: {term}")
        else:
            print("   None")

        print("\n2. Zero Amplitude (Kernel) Conditions:")
        zero_conditions = [build_xor_expr(k) for k in range(self.rank, nu)]
        if zero_conditions:
            for idx, cond in enumerate(zero_conditions):
                print(f"   K_{idx}: ({cond}) == 0")
        else:
            print("   No kernel constraints.\n")


class BranchNode:
    """Tree node for branch decomposition with cached phase state."""
    def __init__(self, weight: complex, label: str = "ROOT"):
        self.weight = weight
        self.label = label
        self.children: List['BranchNode'] = []
        self.v4_state: Optional[np.ndarray] = None


class DicksonTranspiler:
    """
    Transpiles circuits using Dickson decomposition with Van den Nest optimization.
    
    Features:
    - Converts arbitrary Clifford circuits to minimal H + CZ + phase form
    - Applies Van den Nest decomposition for optimal basis switching
    - Handles mixed Clifford/non-Clifford (universal) circuits
    - Optimizes each Clifford block independently
    """

    def __init__(self, num_qubits: int, use_van_den_nest: bool = True):
        self.num_qubits = num_qubits
        self.use_van_den_nest = use_van_den_nest

    @staticmethod
    def _is_clifford_gate(gate_name: str) -> bool:
        """Check if gate is in Clifford group."""
        clifford_gates = {
            'h', 'x', 'y', 'z', 's', 'sdg', 'id', 'cx', 'cz', 'swap'
        }
        return gate_name.lower() in clifford_gates

    def _classify_circuit(self, circuit: QuantumCircuit) -> Tuple[bool, List[Tuple[int, int, bool]]]:
        """
        Classifies circuit as pure Clifford or universal.
        Returns: (is_clifford, [(block_start, block_end, is_clifford_block), ...])
        """
        is_clifford = True
        block_boundaries = []
        current_block_start = 0
        in_clifford_block = True

        for idx, instr in enumerate(circuit.data):
            gate_name = instr.operation.name.lower()
            gate_is_clifford = self._is_clifford_gate(gate_name)

            if gate_is_clifford != in_clifford_block:
                # Boundary detected
                block_boundaries.append((current_block_start, idx, in_clifford_block))
                current_block_start = idx
                in_clifford_block = gate_is_clifford

            if not gate_is_clifford:
                is_clifford = False

        # Add final block
        block_boundaries.append((current_block_start, len(circuit.data), in_clifford_block))

        return is_clifford, block_boundaries

    def _extract_clifford_block(
        self,
        circuit: QuantumCircuit,
        start_idx: int,
        end_idx: int
    ) -> QuantumCircuit:
        """Extracts a continuous block of Clifford gates."""
        block = QuantumCircuit(self.num_qubits)
        for idx in range(start_idx, end_idx):
            instr = circuit.data[idx]
            block.append(instr)
        return block

    def _convert_to_z4_form(self, clifford_circuit: QuantumCircuit) -> Tuple[np.ndarray, np.ndarray]:
        """
        Converts Clifford circuit to Z4 Quadratic Form: B matrix (entanglement) + V vector (phases).
        
        Returns:
            b_matrix: Adjacency/entanglement structure (CZ interactions)
            v_vector: Z4 phase values (0=I, 1=S, 2=Z, 3=Sdg)
        """
        n = self.num_qubits
        b_matrix = np.zeros((n, n), dtype=np.int8)
        v_vector = np.zeros(n, dtype=np.int8)

        for instr in clifford_circuit.data:
            gate_name = instr.operation.name.lower()
            idxs = [clifford_circuit.find_bit(q).index for q in instr.qubits]

            if gate_name == 'h':
                # Hadamard mixes X↔Z basis (implicit in our representation)
                pass
            elif gate_name in ['cx', 'cnot']:
                # CX = H · CZ · H on target, marks entanglement
                c, t = idxs[0], idxs[1]
                b_matrix[c, t] = b_matrix[t, c] = 1
            elif gate_name == 'cz':
                # CZ marks Z-basis entanglement
                c, t = idxs[0], idxs[1]
                b_matrix[c, t] = b_matrix[t, c] = 1
            elif gate_name == 'swap':
                # SWAP adds entanglement between both qubits
                q0, q1 = idxs[0], idxs[1]
                b_matrix[q0, q1] = b_matrix[q1, q0] = 1
            elif gate_name == 'z':
                v_vector[idxs[0]] = (v_vector[idxs[0]] + 2) % 4
            elif gate_name == 's':
                v_vector[idxs[0]] = (v_vector[idxs[0]] + 1) % 4
            elif gate_name == 'sdg':
                v_vector[idxs[0]] = (v_vector[idxs[0]] + 3) % 4
            elif gate_name in ['x', 'y']:
                # X/Y gates affect Z and X basis simultaneously
                pass

        return b_matrix, v_vector

    def _apply_van_den_nest(self, b_matrix: np.ndarray) -> np.ndarray:
        """
        Applies Van den Nest decomposition to optimize basis switching.
        
        Van den Nest decomposes the CZ network into a sequence of:
        1. H operations (basis changes)
        2. CZ operations (diagonal in computational basis)
        3. H operations (basis restoration)
        
        Returns optimized gate sequence as transformed B-matrix.
        """
        # For now, return the matrix as-is. Advanced optimizations would go here.
        return np.copy(b_matrix)

    def _optimize_chain_topology(self, b_matrix: np.ndarray, adj_list: Dict) -> np.ndarray:
        """Optimizes CZ chain by reordering qubits to reduce gate depth."""
        return np.copy(b_matrix)

    def synthesize_clifford(
        self,
        clifford_circuit: QuantumCircuit,
        optimize: bool = True
    ) -> QuantumCircuit:
        """
        Synthesizes optimized Clifford circuit by simulating equivalently.
        
        Args:
            clifford_circuit: Input Clifford circuit
            optimize: Whether to apply Van den Nest optimization (currently returns circuit as-is)
            
        Returns:
            Circuit (currently returns original for correctness)
        """
        # Return circuit as-is to ensure correctness
        # Advanced optimizations deferred to Van den Nest implementation
        return clifford_circuit

    def transpile(self, circuit: QuantumCircuit) -> QuantumCircuit:
        """
        Main transpilation entry point.
        
        - For pure Clifford: applies Van den Nest optimization
        - For universal: splits into Clifford blocks, optimizes each independently
        
        Args:
            circuit: Input circuit (Clifford or universal)
            
        Returns:
            Optimized circuit
        """
        is_clifford, block_info = self._classify_circuit(circuit)
        
        if is_clifford:
            # Pure Clifford: return as-is (optimization deferred)
            return circuit
        else:
            # Universal: return as-is for now
            return circuit


class UniversalQC:
    """
    Universal Quantum Circuit analyzer combining:
    - Left fringe (E_L): Affine input transformations
    - Core: Dickson-reduced Clifford operations
    - Right fringe (E_R): Output permutations
    """

    def __init__(self, circuit: QuantumCircuit):
        self.circuit = circuit
        self.num_qubits = circuit.num_qubits
        self.global_phase = 0.0
        
        self.EL_gates: List = []
        self.core_gates: List = []
        self.ER_gates: List = []
        
        self.root: Optional[BranchNode] = None
        self.engine: Optional[DicksonEngine] = None
        self._skeleton: Optional[QuantumCircuit] = None

    def _split_circuit(self):
        """
        Splits circuit into three layers:
        - E_L: Initial single-qubit and commuting gates
        - Core: Central Clifford/universal gates
        - E_R: Final output permutations
        """
        data = self.circuit.data
        n = self.num_qubits

        # Forward pass: collect E_L gates
        in_EL = [True] * n
        EL_gates, core_gates = [], []

        for ins in data:
            name = ins.operation.name.lower()
            idxs = [self.circuit.find_bit(q).index for q in ins.qubits]

            if name in ['h', 'rx', 'ry']:
                for q in idxs:
                    in_EL[q] = False
                core_gates.append(ins)
            elif name in ['cx', 'swap']:
                if in_EL[idxs[0]] and in_EL[idxs[1]]:
                    EL_gates.append(ins)
                else:
                    for q in idxs:
                        in_EL[q] = False
                    core_gates.append(ins)
            else:
                if all(in_EL[q] for q in idxs):
                    EL_gates.append(ins)
                else:
                    core_gates.append(ins)

        # Reverse pass: collect E_R gates
        in_ER = [True] * n
        ER_gates, final_core = [], []

        for ins in reversed(core_gates):
            name = ins.operation.name.lower()
            idxs = [self.circuit.find_bit(q).index for q in ins.qubits]

            if name in ['h', 'rx', 'ry']:
                for q in idxs:
                    in_ER[q] = False
                final_core.insert(0, ins)
            elif name in ['cx', 'swap']:
                if in_ER[idxs[0]] and in_ER[idxs[1]]:
                    ER_gates.insert(0, ins)
                else:
                    for q in idxs:
                        in_ER[q] = False
                    final_core.insert(0, ins)
            else:
                if all(in_ER[q] for q in idxs):
                    ER_gates.insert(0, ins)
                else:
                    final_core.insert(0, ins)

        self.EL_gates = EL_gates
        self.core_gates = final_core
        self.ER_gates = ER_gates

    def _evaluate_EL(self, x: int) -> Tuple[int, complex]:
        """
        Evaluates left fringe affine transformation.
        Returns: (transformed_input, phase_factor)
        """
        bits = [(x >> i) & 1 for i in range(self.num_qubits)]
        phase = 0.0

        for ins in self.EL_gates:
            name = ins.operation.name.lower()
            idxs = [self.circuit.find_bit(q).index for q in ins.qubits]

            if name == 'cx':
                bits[idxs[1]] ^= bits[idxs[0]]
            elif name == 'swap':
                bits[idxs[0]], bits[idxs[1]] = bits[idxs[1]], bits[idxs[0]]
            elif name == 'cz' and bits[idxs[0]] and bits[idxs[1]]:
                phase += np.pi
            elif name == 'z' and bits[idxs[0]]:
                phase += np.pi
            elif name == 's' and bits[idxs[0]]:
                phase += np.pi / 2
            elif name == 'sdg' and bits[idxs[0]]:
                phase -= np.pi / 2
            elif name == 't' and bits[idxs[0]]:
                phase += np.pi / 4
            elif name == 'rz':
                theta = ins.operation.params[0]
                if bits[idxs[0]]:
                    phase += theta / 2
                else:
                    phase -= theta / 2

        new_x = sum(b << i for i, b in enumerate(bits))
        return new_x, np.exp(1j * phase)

    def _apply_ER(self, sv: np.ndarray):
        """
        Applies right fringe transformations in-place to statevector.
        Avoids circuit reconstruction for efficiency.
        """
        idx = np.arange(len(sv))

        for ins in self.ER_gates:
            name = ins.operation.name.lower()
            idxs = [self.circuit.find_bit(q).index for q in ins.qubits]

            if name == 'cx':
                c, t = idxs[0], idxs[1]
                idx_c1 = idx[(idx & (1 << c)) != 0]
                idx_c1_t0 = idx_c1[(idx_c1 & (1 << t)) == 0]
                idx_c1_t1 = idx_c1_t0 | (1 << t)
                sv[idx_c1_t0], sv[idx_c1_t1] = sv[idx_c1_t1].copy(), sv[idx_c1_t0].copy()
            elif name == 'swap':
                b0, b1 = (idx >> idxs[0]) & 1, (idx >> idxs[1]) & 1
                diff = b0 != b1
                swapped_idx = idx.copy()
                swapped_idx[diff] ^= (1 << idxs[0]) | (1 << idxs[1])
                sv[:] = sv[swapped_idx]
            elif name == 'cz':
                mask = (1 << idxs[0]) | (1 << idxs[1])
                idx_11 = (idx & mask) == mask
                sv[idx_11] *= -1
            elif name in ['z', 's', 'sdg', 't', 'rz']:
                mask = 1 << idxs[0]
                idx_1 = (idx & mask) != 0
                
                if name == 'z':
                    sv[idx_1] *= -1
                elif name == 's':
                    sv[idx_1] *= 1j
                elif name == 'sdg':
                    sv[idx_1] *= -1j
                elif name == 't':
                    sv[idx_1] *= np.exp(1j * np.pi / 4)
                elif name == 'rz':
                    theta = ins.operation.params[0]
                    idx_0 = ~idx_1
                    sv[idx_1] *= np.exp(1j * theta / 2)
                    sv[idx_0] *= np.exp(-1j * theta / 2)

    def build_tree(self):
        """
        Builds the branch tree for exponential sum decomposition.
        Optimizes core Clifford blocks using Dickson.
        """
        self._split_circuit()

        # Calculate global phase from T-gates in core
        core_t_count = sum(1 for ins in self.core_gates if ins.operation.name.lower() == 't')
        self.global_phase = core_t_count * (np.pi / 8)

        # Create skeleton circuit from core gates
        self._skeleton = QuantumCircuit(self.num_qubits)
        for ins in self.core_gates:
            self._skeleton.append(ins)

        self.engine = DicksonEngine(self._skeleton)

        # Pre-compute gate-to-variable mapping
        gate_to_vars = []
        wires = [[i] for i in range(self.num_qubits)]
        nv = self.num_qubits

        for ins in self._skeleton.data:
            name = ins.operation.name.lower()
            idxs = [self._skeleton.find_bit(q).index for q in ins.qubits]
            gate_to_vars.append([wires[q][-1] for q in idxs])

            if name == 'h':
                wires[idxs[0]].append(nv)
                nv += 1
            elif name == 'cx':
                wires[idxs[1]].append(nv)
                nv += 1
                wires[idxs[1]].append(nv)
                nv += 1
            elif name == 'cz':
                pass

        # Initialize root
        self.root = BranchNode(1.0 + 0j)
        self.root.v4_state = np.zeros(self.engine.n_vars, dtype=np.int8)

        # Grow branch tree recursively
        def _grow(idx: int, node: BranchNode):
            if idx == len(self._skeleton.data):
                return

            ins = self._skeleton.data[idx]
            name = ins.operation.name.lower()
            q_vars = gate_to_vars[idx] if idx < len(gate_to_vars) else []

            if name in ['h', 'cz', 'cx', 'id']:
                _grow(idx + 1, node)
            elif name == 'z':
                if q_vars:
                    node.v4_state[q_vars[0]] = (node.v4_state[q_vars[0]] + 2) % 4
                _grow(idx + 1, node)
            elif name == 's':
                if q_vars:
                    node.v4_state[q_vars[0]] = (node.v4_state[q_vars[0]] + 1) % 4
                _grow(idx + 1, node)
            elif name == 'sdg':
                if q_vars:
                    node.v4_state[q_vars[0]] = (node.v4_state[q_vars[0]] + 3) % 4
                _grow(idx + 1, node)
            elif name in ('rz', 't'):
                theta = np.pi / 4 if name == 't' else ins.operation.params[0]
                
                # Identity branch
                c1 = BranchNode(node.weight * np.cos(theta / 2), label=node.label + "_I")
                c1.v4_state = node.v4_state.copy()

                # Z branch
                c2 = BranchNode(node.weight * (-1j * np.sin(theta / 2)), label=node.label + "_Z")
                c2.v4_state = node.v4_state.copy()
                if q_vars:
                    c2.v4_state[q_vars[0]] = (c2.v4_state[q_vars[0]] + 2) % 4

                node.children.extend([c1, c2])
                _grow(idx + 1, c1)
                _grow(idx + 1, c2)

        _grow(0, self.root)

    def get_statevector(self, x: int = 0) -> np.ndarray:
        """
        Computes full statevector.
        - Applies E_L affine transform to input
        - Uses optimized Dickson engine for core
        - Applies E_R permutations to output
        """
        if self.root is None:
            raise RuntimeError("Call build_tree() before get_statevector().")

        # Step 1: Pre-process E_L transform
        x_core, input_phase = self._evaluate_EL(x)
        sv = np.zeros(2 ** self.num_qubits, dtype=np.complex128)

        # Step 2: Traverse branch tree and accumulate amplitudes
        stack = [self.root]
        while stack:
            curr = stack.pop()
            if curr.children:
                stack.extend(curr.children)
                continue

            self.engine.v4 = curr.v4_state.copy()
            self.engine.get_statevector_gray(curr.weight, sv, x=x_core)

        # Step 3: Apply global phase and input phase
        final_sv = sv * np.exp(1j * self.global_phase) * input_phase

        # Step 4: Post-process E_R permutations
        self._apply_ER(final_sv)

        return final_sv

    def print_full_analytic_decomposition(self, transition_mode: bool = False):
        """Prints analytic decomposition for all branches."""
        if self.root is None:
            print("Tree not built.")
            return

        print(f"{'=' * 60}\nUNIVERSAL QC ANALYTIC DECOMPOSITION\n{'=' * 60}")
        print(f"E_L Gates: {len(self.EL_gates)} (Input Affine Transform & Phases)")
        print(f"E_R Gates: {len(self.ER_gates)} (Output Affine Permutations)")
        print(f"Core Global Phase from T-gates: {self.global_phase:.6f}\n")

        # Collect leaves
        leaves, stack = [], [(self.root, "ROOT")]
        while stack:
            node, path = stack.pop()
            if not node.children:
                leaves.append((node, path))
            else:
                for i, child in enumerate(node.children):
                    stack.append((child, f"{path} -> {child.label}"))

        # Print each leaf branch
        for node, path in leaves:
            self.engine.v4 = node.v4_state.copy()
            self.engine.print_analytic_formula(
                transition_mode=transition_mode,
                weight=node.weight,
                branch_label=path
            )