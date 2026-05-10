import numpy as np

from scipy.linalg import expm

from qiskit import QuantumCircuit
from qiskit.circuit import Parameter
from qiskit.quantum_info import Operator, SparsePauliOp
from qiskit.circuit.library import PauliEvolutionGate

# ============================================================
# Helper: Cancel adjacent H gates
# ============================================================

def cancel_h_pairs(qc):
    """
    Scans through the QuantumCircuit and cancels adjacent H gates.
    H . H = I, so any pair of H gates on the same qubit with no
    intervening operations on that qubit will be removed.
    """
    new_qc = QuantumCircuit(*qc.qregs, *qc.cregs, name=qc.name)
    pending_h = set()
    
    for inst in qc.data:
        # Compatibility across different Qiskit versions
        op = getattr(inst, 'operation', inst[0])
        qubits = getattr(inst, 'qubits', inst[1])
        
        if op.name == 'h':
            q = qubits[0]
            if q in pending_h:
                pending_h.remove(q) # H . H cancels out
            else:
                pending_h.add(q)    # Store pending H
        else:
            # Flush pending H gates for qubits involved in this operation
            for q in qubits:
                if q in pending_h:
                    new_qc.h(q)
                    pending_h.remove(q)
            
            # Append the current operation
            if hasattr(inst, 'operation'):
                new_qc.append(inst)
            else:
                new_qc.append(inst[0], inst[1], inst[2])
                
    # Flush any remaining H gates at the end of the circuit
    for q in qc.qubits:
        if q in pending_h:
            new_qc.h(q)
            
    return new_qc

# ============================================================
# Helper: CZ-based CNOT
# ============================================================

def cx_via_cz(qc, control, target):
    """
    Implements CX(control,target)
    using only H + CZ.
    """
    qc.h(target)
    qc.cz(control, target)
    qc.h(target)

# ============================================================
# Main transpiler
# ============================================================

def transpile_ppr(
    pauli_string,
    theta,
    rz_location=None,
    label=None
):
    """
    Transpile:

        exp(-i theta P / 2)

    into:

        {H, S, Z, CZ, RZ}
    """

    n = len(pauli_string)

    if rz_location is None:
        rz_location = n - 1

    target = rz_location

    qc = QuantumCircuit(n, name=label)

    # ========================================================
    # 1. PAULI BASIS TRANSFORMS
    # ========================================================

    for q, p in enumerate(pauli_string):

        if p == 'X':
            qc.h(q)

        elif p == 'Y':
            qc.sdg(q)
            qc.h(q)

        elif p == 'Z':
            pass

        elif p == 'I':
            pass
            
        else:
            raise ValueError(f'Unsupported Pauli: {p}')

    # ========================================================
    # 2. PARITY ACCUMULATION
    # ========================================================

    # Build CNOT ladder using CZ+H.

    for q in range(n):

        if q == target or pauli_string[q] == 'I':
            continue

        cx_via_cz(qc, q, target)

    # ========================================================
    # 3. CENTRAL PHASE
    # ========================================================

    qc.rz(theta, target)

    # ========================================================
    # 4. UNCOMPUTE
    # ========================================================

    for q in reversed(range(n)):

        if q == target or pauli_string[q] == 'I':
            continue

        cx_via_cz(qc, q, target)

    # ========================================================
    # 5. UNDO PAULI BASIS TRANSFORMS
    # ========================================================

    for q, p in reversed(list(enumerate(pauli_string))):

        if p == 'X':
            qc.h(q)

        elif p == 'Y':
            qc.h(q)
            qc.s(q)

    # Optimize the circuit by canceling identical self-inverse H sequences
    return cancel_h_pairs(qc)

# ============================================================
# Exact reference matrices
# ============================================================

PAULI = {
    'I': np.eye(2, dtype=complex),
    'X': np.array([[0, 1], [1, 0]], dtype=complex),
    'Y': np.array([[0, -1j], [1j, 0]], dtype=complex),
    'Z': np.array([[1, 0], [0, -1]], dtype=complex),
}

def pauli_matrix(pauli_string):
    M = PAULI[pauli_string[0]]
    for p in pauli_string[1:]:
        M = np.kron(M, PAULI[p])
    return M

def exact_ppr(pauli_string, theta):
    P = pauli_matrix(pauli_string)
    return expm(-1j * theta * P / 2)

# ============================================================
# Verification against exact matrix
# ============================================================

def verify_exact(
    pauli_string,
    theta_value=np.pi / 5,
    rz_location=None
):
    theta = Parameter('θ')
    qc = transpile_ppr(
        pauli_string,
        theta,
        rz_location=rz_location
    )
    qc_bound = qc.assign_parameters({theta: theta_value})
    U_circuit = Operator(qc_bound).data
    U_exact = exact_ppr(pauli_string, theta_value)

    ok = np.allclose(U_circuit, U_exact)

    print('\n=====================================')
    print('Exact Matrix Verification')
    print('=====================================')
    print('Pauli String :', pauli_string)
    print('RZ location  :', rz_location)
    print('PASS' if ok else 'FAIL')
    return ok

# ============================================================
# Verification against Qiskit PauliEvolutionGate
# ============================================================

def verify_qiskit(
    pauli_string,
    theta_value=np.pi / 5,
    rz_location=None
):
    theta = Parameter('θ')
    qc_custom = transpile_ppr(
        pauli_string,
        theta,
        rz_location=rz_location
    )
    qc_custom = qc_custom.assign_parameters(
        {theta: theta_value}
    )

    # --------------------------------------------------------
    # Reference implementation
    # --------------------------------------------------------
    op = SparsePauliOp(pauli_string)
    evo = PauliEvolutionGate(
        op,
        time=theta_value / 2
    )
    qc_ref = QuantumCircuit(len(pauli_string))
    qc_ref.append(evo, range(len(pauli_string)))

    U_custom = Operator(qc_custom).data
    U_ref = Operator(qc_ref).data
    ok = np.allclose(U_custom, U_ref)

    print('\n=====================================')
    print('Qiskit Verification')
    print('=====================================')
    print('Pauli String :', pauli_string)
    print('RZ location  :', rz_location)
    print('PASS' if ok else 'FAIL')
    return ok

# ============================================================
# Gate counting
# ============================================================

def gate_summary(qc):
    ops = qc.count_ops()
    print('\nGate Counts:')
    for k, v in ops.items():
        print(f'{k:>5s} : {v}')
    print('\nH-count:', ops.get('h', 0))

# ============================================================
# Example circuits
# ============================================================

if __name__ == "__main__":
    theta = Parameter('θ')
    
    qc1 = transpile_ppr('XZZX', theta, rz_location=0)
    print('\n=== XZZX ===\n')
    print(qc1.draw('text'))
    gate_summary(qc1)
    
    qc2 = transpile_ppr('YZZY', theta, rz_location=0)
    print('\n=== YZZY ===\n')
    print(qc2.draw('text'))
    gate_summary(qc2)
    
    # ============================================================
    # Verification tests
    # ============================================================
    verify_exact('XZZX', rz_location=0)
    verify_exact('XZZX', rz_location=2)
    verify_exact('YZZZZY', rz_location=3)
    
    verify_qiskit('XZZX', rz_location=0)
    verify_qiskit('XZZX', rz_location=1)
    verify_qiskit('YZZZZY', rz_location=4)