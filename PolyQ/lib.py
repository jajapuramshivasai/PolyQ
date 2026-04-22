import numpy as np
from dataclasses import dataclass
from typing import List, Tuple, Optional, Dict
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
import os

# Optional verification dependency (only used if you call verify_against_qiskit)
try:
    from qiskit import QuantumCircuit  # type: ignore
    from qiskit.quantum_info import Statevector  # type: ignore
    _HAS_QISKIT = True
except Exception:
    _HAS_QISKIT = False


# ============================================================
# Dickson reduction ops (for the Clifford-only amplitude engine)
# ============================================================
class DicksonOp:
    def __init__(self, op_type, a, b=None, target=None, pivot=None):
        self.type = op_type
        # SWAP uses a,b; ADD uses pivot,target (we store as a,b)
        if op_type == "SWAP":
            self.a, self.b = a, b
        elif op_type == "ADD":
            self.a, self.b = pivot, target  # pivot -> a, target -> b
        else:
            self.a, self.b = a, b


# ============================================================
# QC: Clifford amplitude simulator + branching wrapper for Rz
# ============================================================
class QC:
    """
    Simulator for circuits built from H, S, Z, CZ plus branched Rz(omega).

    - Clifford-only circuits: analytic Dickson/F2 method (fast).
    - With Rz: sum-over-Cliffords branching:
        Rz(ω) = cos(ω) I  +  (-i*sin(ω)) Z
      (as requested)

    Optimizations:
      - Cache compiled Clifford branches (each branch compiled once).
      - Cache branch list for a given threshold.
      - Parallelized statevector generation (optional; see get_statevector_0()).
    """

    def __init__(self, num_qubits: int):
        self.num_qubits = num_qubits
        self.gates: List[Tuple] = []

        # Clifford compilation artifacts (for *this* circuit only, Clifford-only case)
        self.compiled = False
        self.b4: Optional[np.ndarray] = None
        self.v4: Optional[np.ndarray] = None
        self.output_vars: Optional[List[int]] = None
        self.num_h, self.n_vars, self.rank = 0, 0, 0
        self.uvars_skeleton: Optional[List[int]] = None
        self.ops: List[DicksonOp] = []
        self.b_reduced: Optional[np.ndarray] = None

        # Branching caches
        # threshold -> list of _CompiledBranch
        self._branch_cache: Dict[float, List["QC._CompiledBranch"]] = {}

        # Clifford-only amplitude cache for x=0 (optional; most helpful when repeatedly querying few y's)
        self._amp0_cache: Dict[int, complex] = {}

    # -----------------
    # Gate constructors
    # -----------------
    def h(self, q: int): self.gates.append(("H", q))
    def z(self, q: int): self.gates.append(("Z", q))
    def s(self, q: int): self.gates.append(("S", q))
    def cz(self, q1: int, q2: int): self.gates.append(("CZ", q1, q2))

    def rz(self, q: int, omega: float):
        """
        Non-Clifford (in general) rotation about Z.
        Branching rule used:
            Rz(ω) = cos(ω) I - i sin(ω) Z
        (as requested)
        """
        self.gates.append(("RZ", q, float(omega)))

    # ----------------------------------------
    # Internal: compile Clifford-only structure
    # ----------------------------------------
    def compile(self):
        """
        Compile the circuit under the assumption all gates are Clifford.
        If RZ gates exist, compilation will fail; use branching entry points.
        """
        if any(g[0] == "RZ" for g in self.gates):
            raise RuntimeError(
                "Circuit contains RZ gates; use branching evaluation (get_amplitude/get_transition_matrix), "
                "or remove RZ before calling compile()."
            )

        wires = [[i] for i in range(self.num_qubits)]
        next_var = self.num_qubits
        self.num_h = 0
        b4_adj, v4_dict = {}, {}

        def ensure_var(v):
            b4_adj.setdefault(v, set())
            v4_dict.setdefault(v, 0)

        for i in range(self.num_qubits):
            ensure_var(i)

        for gate in self.gates:
            g_type = gate[0]
            if g_type == "H":
                q = gate[1]
                prev, cur = wires[q][-1], next_var
                next_var += 1
                self.num_h += 1
                wires[q].append(cur)
                ensure_var(cur)
                b4_adj[prev].add(cur)
                b4_adj[cur].add(prev)
            elif g_type == "CZ":
                v1, v2 = wires[gate[1]][-1], wires[gate[2]][-1]
                b4_adj[v1].add(v2)
                b4_adj[v2].add(v1)
            elif g_type == "Z":
                v = wires[gate[1]][-1]
                v4_dict[v] = (v4_dict[v] + 2) % 4
            elif g_type == "S":
                v = wires[gate[1]][-1]
                v4_dict[v] = (v4_dict[v] + 1) % 4
            else:
                raise ValueError(f"Unsupported gate in Clifford compile: {gate}")

        self.n_vars = next_var
        self.b4 = np.zeros((self.n_vars, self.n_vars), dtype=np.int8)
        self.v4 = np.zeros(self.n_vars, dtype=np.int8)

        for i in range(self.n_vars):
            self.v4[i] = v4_dict.get(i, 0)
            for j in b4_adj.get(i, []):
                self.b4[i, j] = 1

        self.output_vars = [wires[q][-1] for q in range(self.num_qubits)]
        self.uvars_skeleton = [
            i for i in range(self.n_vars)
            if i >= self.num_qubits and i not in self.output_vars
        ]
        nu = len(self.uvars_skeleton)

        b_u_skel = np.zeros((nu, nu), dtype=np.int8)
        for ui, orig_u in enumerate(self.uvars_skeleton):
            for uj, orig_uj in enumerate(self.uvars_skeleton):
                if ui < uj and self.b4[orig_u, orig_uj]:
                    b_u_skel[ui, uj] = b_u_skel[uj, ui] = 1

        self.ops, self.b_reduced, self.rank = self._plan_dickson_f2(b_u_skel, nu)
        self.compiled = True

        # Clifford-only: clear amplitude cache since structure changed
        self._amp0_cache.clear()

    def _plan_dickson_f2(self, b_matrix, n_vars):
        b, ops, r, p = np.copy(b_matrix), [], 0, 0
        while p + 1 < n_vars:
            pivot = next(((i, j)
                          for i in range(p, n_vars)
                          for j in range(i + 1, n_vars)
                          if b[i, j] == 1), None)
            if not pivot:
                break
            i, j = pivot

            if i != p:
                b[[p, i]] = b[[i, p]]
                b[:, [p, i]] = b[:, [i, p]]
                ops.append(DicksonOp("SWAP", p, i))

            j_act = i if j == p else j
            if j_act != p + 1:
                b[[p + 1, j_act]] = b[[j_act, p + 1]]
                b[:, [p + 1, j_act]] = b[:, [j_act, p + 1]]
                ops.append(DicksonOp("SWAP", p + 1, j_act))

            rp, rp1 = np.copy(b[p, :]), np.copy(b[p + 1, :])

            for k in range(p + 2, n_vars):
                if b[k, p] == 1:
                    b[k, :] ^= rp1
                    b[:, k] ^= rp1
                    ops.append(DicksonOp("ADD", a=None, b=None, pivot=p + 1, target=k))
                if b[k, p + 1] == 1:
                    b[k, :] ^= rp
                    b[:, k] ^= rp
                    ops.append(DicksonOp("ADD", a=None, b=None, pivot=p, target=k))
            r += 2
            p += 2
        return ops, b, r

    def _eval_canonical_sum(self, vu_base, nu):
        sum_val, p = 1.0 + 0j, 0
        while p < self.rank:
            has_edge = self.b_reduced[p, p + 1]  # type: ignore[index]
            pair_sum = 0j
            for x1 in (0, 1):
                for x2 in (0, 1):
                    ph = (
                        (2 if has_edge and x1 and x2 else 0)
                        + (vu_base[p] if x1 else 0)
                        + (vu_base[p + 1] if x2 else 0)
                    ) % 4
                    if ph == 0:
                        pair_sum += 1
                    elif ph == 1:
                        pair_sum += 1j
                    elif ph == 2:
                        pair_sum += -1
                    elif ph == 3:
                        pair_sum += -1j
            sum_val *= pair_sum
            p += 2

        for k in range(self.rank, nu):
            val = int(vu_base[k]) % 4
            if val == 0:
                sum_val *= 2.0
            elif val == 1:
                sum_val *= (1.0 + 1j)
            elif val == 2:
                return 0j
            elif val == 3:
                sum_val *= (1.0 - 1j)
        return sum_val

    # ============================================================
    # Branching (sum-over-Cliffords) for RZ
    # ============================================================
    @dataclass(frozen=True)
    class _Branch:
        weight: complex
        gates: Tuple[Tuple, ...]  # Clifford-only gate list (immutable)

    @dataclass
    class _CompiledBranch:
        weight: complex
        qc: "QC"  # Clifford-only compiled QC

    def _decompose_gate_for_branching(self, gate: Tuple) -> List[Tuple[complex, List[Tuple]]]:
        g = gate[0]
        if g in ("H", "S", "Z", "CZ"):
            return [(1.0 + 0j, [gate])]
        if g == "RZ":
            _, q, omega = gate
            w0 = complex(np.cos(omega), 0.0)
            w1 = (-1j) * complex(np.sin(omega), 0.0)
            return [
                (w0, []),            # I
                (w1, [("Z", q)]),    # Z
            ]
        raise ValueError(f"Unsupported gate: {gate}")

    def _build_clifford_branches(self, threshold: float) -> List["_Branch"]:
        """
        Expand full circuit into Clifford-only branches, pruning by |weight| < threshold.
        """
        branches: List[QC._Branch] = [QC._Branch(weight=1.0 + 0j, gates=tuple())]

        for gate in self.gates:
            decomposed = self._decompose_gate_for_branching(gate)

            if len(decomposed) == 1:
                w, gg = decomposed[0]
                if gg:
                    gg_t = tuple(gg)
                    branches = [QC._Branch(weight=b.weight * w, gates=b.gates + gg_t) for b in branches
                                if abs(b.weight * w) >= threshold]
                else:
                    branches = [QC._Branch(weight=b.weight * w, gates=b.gates) for b in branches
                                if abs(b.weight * w) >= threshold]
                if not branches:
                    break
                continue

            new_branches: List[QC._Branch] = []
            for b in branches:
                for w, gg in decomposed:
                    new_w = b.weight * w
                    if abs(new_w) < threshold:
                        continue
                    new_branches.append(QC._Branch(weight=new_w, gates=b.gates + tuple(gg)))
            branches = new_branches
            if not branches:
                break

        return branches

    def _get_compiled_branches(self, branch_threshold: float) -> List["_CompiledBranch"]:
        """
        Build (or fetch from cache) compiled Clifford branches for the given threshold.
        Each branch is compiled once, then reused for many amplitudes/statevectors.
        """
        cached = self._branch_cache.get(branch_threshold)
        if cached is not None:
            return cached

        branches = self._build_clifford_branches(threshold=branch_threshold)
        compiled: List[QC._CompiledBranch] = []

        for br in branches:
            tmp = QC(self.num_qubits)
            tmp.gates = list(br.gates)
            tmp.compile()
            compiled.append(QC._CompiledBranch(weight=br.weight, qc=tmp))

        self._branch_cache[branch_threshold] = compiled
        return compiled

    # ============================================================
    # Clifford amplitude core (compiled)
    # ============================================================
    def _get_amplitude_clifford_compiled(self, y_val: int, x_val: int = 0) -> complex:
        if not self.compiled:
            raise RuntimeError("Circuit must be compiled first.")

        fixed = [None] * self.n_vars  # type: ignore[arg-type]
        for i in range(self.num_qubits):
            fixed[i] = (x_val >> i) & 1

        for i in range(self.num_qubits):
            bit = (y_val >> i) & 1
            ov = self.output_vars[i]  # type: ignore[index]
            if fixed[ov] is None:
                fixed[ov] = bit
            elif fixed[ov] != bit:
                return 0j

        eps_base = 0
        f_list = [v for v, val in enumerate(fixed) if val == 1]

        # phase from fixed vars
        v4 = self.v4
        b4 = self.b4
        for f in f_list:
            eps_base = (eps_base + int(v4[f])) % 4  # type: ignore[index]
            for f2 in f_list:
                if f < f2 and b4[f, f2]:  # type: ignore[index]
                    eps_base = (eps_base + 2) % 4

        nu = len(self.uvars_skeleton)  # type: ignore[arg-type]
        vu_base = np.zeros(nu, dtype=np.int8)

        uvars = self.uvars_skeleton  # type: ignore[assignment]
        for ui, orig_u in enumerate(uvars):
            vu = int(v4[orig_u]) % 4  # type: ignore[index]
            for f in f_list:
                if b4[orig_u, f]:  # type: ignore[index]
                    vu = (vu + 2) % 4
            vu_base[ui] = vu

        for op in self.ops:
            if op.type == "SWAP":
                vu_base[op.a], vu_base[op.b] = vu_base[op.b], vu_base[op.a]
            elif op.type == "ADD":
                pivot = op.a
                target = op.b
                vu_base[target] = (vu_base[target] + vu_base[pivot]) % 4

        phase_map = {0: 1.0 + 0j, 1: 1j, 2: -1.0 + 0j, 3: -1j}
        return phase_map[eps_base] * self._eval_canonical_sum(vu_base, nu) * (2.0 ** (-self.num_h / 2.0))

    # ============================================================
    # Public API (general amplitude / matrix)
    # ============================================================
    def get_amplitude(self, y_val: int, x_val: int = 0, *, branch_threshold: float = 1e-12) -> complex:
        if any(g[0] == "RZ" for g in self.gates):
            total = 0j
            compiled = self._get_compiled_branches(branch_threshold)
            for br in compiled:
                total += br.weight * br.qc._get_amplitude_clifford_compiled(y_val=y_val, x_val=x_val)
            return total

        if not self.compiled:
            self.compile()
        return self._get_amplitude_clifford_compiled(y_val=y_val, x_val=x_val)

    def get_transition_matrix(self, *, branch_threshold: float = 1e-12) -> np.ndarray:
        dim = 2 ** self.num_qubits
        return np.array(
            [[self.get_amplitude(y, x, branch_threshold=branch_threshold) for x in range(dim)] for y in range(dim)],
            dtype=np.complex128,
        )

    # ============================================================
    # NEW (requested): optimized amplitude + statevector for input |0...0>
    # ============================================================
    def get_output_amplitude_0(self, y_val: int, *, branch_threshold: float = 1e-12, cache: bool = True) -> complex:
        """
        Optimized evaluation of amplitude <y|U|0...0>.

        - Uses x_val = 0.
        - For Clifford-only: compiles once and can cache amplitudes per y.
        - For RZ branching: caches compiled branches and evaluates each compiled branch.
        """
        if any(g[0] == "RZ" for g in self.gates):
            # Branching mode: caching is handled at branch compilation level.
            compiled = self._get_compiled_branches(branch_threshold)
            # (Optional) could add per-y caching here too; often less useful unless repeatedly querying same y.
            total = 0j
            for br in compiled:
                total += br.weight * br.qc._get_amplitude_clifford_compiled(y_val=y_val, x_val=0)
            return total

        # Clifford-only fast path
        if not self.compiled:
            self.compile()

        if cache:
            hit = self._amp0_cache.get(y_val)
            if hit is not None:
                return hit
            amp = self._get_amplitude_clifford_compiled(y_val=y_val, x_val=0)
            self._amp0_cache[y_val] = amp
            return amp

        return self._get_amplitude_clifford_compiled(y_val=y_val, x_val=0)

    def get_statevector_0(
        self,
        *,
        branch_threshold: float = 1e-12,
        parallel: bool = True,
        backend: str = "thread",
        workers: Optional[int] = None,
        chunk_size: int = 256,
        cache_amplitudes: bool = False,
    ) -> np.ndarray:
        """
        Optimized full statevector for input |0...0>:
            sv[y] = <y|U|0...0>

        Performance notes:
          - In branching mode, the expensive part is compiling each branch; we cache that.
          - After that, computing sv is embarrassingly parallel over y.
          - Default parallel backend is threads (works well because numpy + python overhead mixes;
            process backend can be faster for heavier workloads but has pickling overhead).

        Params:
          - parallel: enable parallel evaluation over y
          - backend: "thread" or "process"
          - workers: number of worker threads/processes (default: os.cpu_count()).
          - chunk_size: number of y's per submitted task (reduces overhead).
          - cache_amplitudes: if True and Clifford-only, fills internal cache for subsequent single-amplitude queries.
        """
        dim = 2 ** self.num_qubits
        sv = np.empty(dim, dtype=np.complex128)

        # Pre-compile / pre-fetch compiled branches once
        has_rz = any(g[0] == "RZ" for g in self.gates)
        if has_rz:
            compiled_branches = self._get_compiled_branches(branch_threshold)

            # Define a local worker that only closes over compiled_branches (threads ok; processes will need pickling)
            def eval_chunk(y_start: int, y_end: int):
                out = np.empty(y_end - y_start, dtype=np.complex128)
                for idx, y in enumerate(range(y_start, y_end)):
                    total = 0j
                    for br in compiled_branches:
                        total += br.weight * br.qc._get_amplitude_clifford_compiled(y_val=y, x_val=0)
                    out[idx] = total
                return y_start, out
        else:
            if not self.compiled:
                self.compile()

            def eval_chunk(y_start: int, y_end: int):
                out = np.empty(y_end - y_start, dtype=np.complex128)
                for idx, y in enumerate(range(y_start, y_end)):
                    out[idx] = self._get_amplitude_clifford_compiled(y_val=y, x_val=0)
                return y_start, out

        if not parallel or dim <= chunk_size:
            # Sequential
            y0, out = eval_chunk(0, dim)
            sv[y0:y0 + len(out)] = out
        else:
            if workers is None:
                workers = os.cpu_count() or 1

            Executor = ThreadPoolExecutor if backend.lower().startswith("thread") else ProcessPoolExecutor

            # For process backend, compiled branches/QC objects are not reliably picklable.
            # So we strongly recommend threads unless you redesign serialization.
            if Executor is ProcessPoolExecutor and has_rz:
                raise RuntimeError(
                    "backend='process' is not supported in branching mode because compiled QC branches aren't picklable. "
                    "Use backend='thread' or parallel=False."
                )

            futures = []
            with Executor(max_workers=workers) as ex:
                for y_start in range(0, dim, chunk_size):
                    y_end = min(dim, y_start + chunk_size)
                    futures.append(ex.submit(eval_chunk, y_start, y_end))

                for fut in as_completed(futures):
                    y_start, out = fut.result()
                    sv[y_start:y_start + len(out)] = out

        if (not has_rz) and cache_amplitudes:
            # Fill cache for later single amplitude queries
            self._amp0_cache = {y: sv[y] for y in range(dim)}

        return sv

    # ============================================================
    # Printing / diagnostics
    # ============================================================
    def print_circuit_parameters(self):
        if any(g[0] == "RZ" for g in self.gates):
            branches = self._build_clifford_branches(threshold=0.0)
            print(f"\n{'='*60}\n   CIRCUIT METRICS (BRANCHING MODE)\n{'='*60}")
            print(f"Total Qubits (n)           : {self.num_qubits}")
            print(f"Total Gates                : {len(self.gates)}")
            print(f"RZ Gates                   : {sum(1 for g in self.gates if g[0]=='RZ')}")
            print(f"Total Clifford Branches    : {len(branches)} (before pruning)")
            print(f"{'='*60}\n")
            return

        if not self.compiled:
            raise RuntimeError("Circuit must be compiled first.")

        nu = len(self.uvars_skeleton)  # type: ignore[arg-type]
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

    def verify_against_qiskit(self, *, branch_threshold: float = 1e-12):
        if not _HAS_QISKIT:
            raise RuntimeError("qiskit is not installed; cannot verify against qiskit.")

        dim = 2 ** self.num_qubits
        custom_sv = self.get_statevector_0(branch_threshold=branch_threshold, parallel=False)

        qiskit_qc = QuantumCircuit(self.num_qubits)
        for g in self.gates:
            gt = g[0]
            if gt == "H":
                qiskit_qc.h(g[1])
            elif gt == "S":
                qiskit_qc.s(g[1])
            elif gt == "Z":
                qiskit_qc.z(g[1])
            elif gt == "CZ":
                qiskit_qc.cz(g[1], g[2])
            elif gt == "RZ":
                # Requested Rz(ω)=exp(-i ω Z) whereas qiskit rz(theta)=exp(-i theta/2 Z).
                qiskit_qc.rz(2.0 * g[2], g[1])
            else:
                raise ValueError(f"Unsupported gate for qiskit verify: {g}")

        max_diff = np.max(np.abs(Statevector(qiskit_qc).data - custom_sv))
        print(f"{'✅ PASS!' if max_diff < 1e-10 else '❌ FAIL!'} Max deviation: {max_diff:.2e}")


# ============================================================
# Example usage
# ============================================================
if __name__ == "__main__":
    qc = QC(num_qubits=3)

    qc.h(0)
    qc.cz(0, 1)
    qc.rz(0, np.pi / 8)
    qc.s(2)
    qc.h(0)

    qc.print_circuit_parameters()

    # Single amplitude for |0...0> input
    amp = qc.get_output_amplitude_0(y_val=5, branch_threshold=1e-12)
    print("Amplitude <y=5|U|0> =", amp)

    # Full statevector for |0...0> input (parallel)
    sv = qc.get_statevector_0(branch_threshold=1e-12, parallel=True, backend="thread", chunk_size=256)
    print("Statevector norm:", np.linalg.norm(sv))

    # Optional: verify vs qiskit if installed
    qc.verify_against_qiskit()