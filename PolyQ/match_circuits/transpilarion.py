import numpy as np

from qiskit import QuantumCircuit
from qiskit.circuit import Parameter
from qiskit.quantum_info import Operator


# ============================================================
# Minimal-H CZ-only Pauli Product Rotation transpiler
#
# Target gateset:
#   {H, Z, S, CZ, RZ}
#
# Supported Pauli string form:
#   A0 Z Z Z ... Z A{n-1}
#
# where endpoints A_i ∈ {X,Y,Z}
#
# Implements:
#
#   exp(-i theta P / 2)
#
# using:
#   - minimal Hadamards
#   - CZ-only entanglers
#   - one RZ
#
# ============================================================


def transpile_ppr(pauli_string, theta, label=None):
    """
    Transpile Pauli Product Rotation into
    {H, Z, S, CZ, RZ} gateset.

    Parameters
    ----------
    pauli_string : str
        Example:
            "XZZX"
            "YZZY"
            "XZZZZY"

    theta : float or Parameter
        Rotation angle.

    Returns
    -------
    QuantumCircuit
    """

    n = len(pauli_string)

    qc = QuantumCircuit(n, name=label)

    target = n - 1

    # ========================================================
    # 1. BASIS TRANSFORMS
    # ========================================================

    #
    # X -> H Z H
    #
    # Y -> Sdg H Z H S
    #
    # Minimal H count:
    #   exactly one H per X/Y qubit
    #

    for q, p in enumerate(pauli_string):

        if p == "X":
            qc.h(q)

        elif p == "Y":
            qc.sdg(q)
            qc.h(q)

        elif p == "Z":
            pass

        else:
            raise ValueError(f"Unsupported Pauli: {p}")

    # ========================================================
    # 2. CZ PARITY ENCODING
    # ========================================================

    for q in range(n - 1):
        qc.cz(q, target)

    # ========================================================
    # 3. CENTRAL PHASE ROTATION
    # ========================================================

    qc.rz(theta, target)

    # ========================================================
    # 4. UNCOMPUTE PARITY
    # ========================================================

    for q in reversed(range(n - 1)):
        qc.cz(q, target)

    # ========================================================
    # 5. UNDO BASIS TRANSFORMS
    # ========================================================

    for q, p in reversed(list(enumerate(pauli_string))):

        if p == "X":
            qc.h(q)

        elif p == "Y":
            qc.h(q)
            qc.s(q)

    return qc


# ============================================================
# Utility:
# Count H gates
# ============================================================

def count_h_gates(qc):
    return qc.count_ops().get("h", 0)


# ============================================================
# Example Usage
# ============================================================

theta = Parameter("θ")

qc1 = transpile_ppr("XZZX", theta)

print("\n=== XZZX PPR ===\n")
print(qc1.draw("text"))

print("\nGate counts:")
print(qc1.count_ops())

print("\nH-count:", count_h_gates(qc1))


# ============================================================
# Another Example
# ============================================================

qc2 = transpile_ppr("YZZZZY", theta)

print("\n=== YZZZZY PPR ===\n")
print(qc2.draw("text"))

print("\nGate counts:")
print(qc2.count_ops())

print("\nH-count:", count_h_gates(qc2))


# ============================================================
# Verify unitary correctness
# ============================================================

#
# Compare against exact matrix:
#
#   exp(-i θ P / 2)
#

from scipy.linalg import expm


PAULI = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.array([[1, 0], [0, -1]], dtype=complex),
}


def pauli_matrix(pauli_string):

    M = PAULI[pauli_string[0]]

    for p in pauli_string[1:]:
        M = np.kron(M, PAULI[p])

    return M


def exact_ppr(pauli_string, theta_value):

    P = pauli_matrix(pauli_string)

    return expm(-1j * theta_value * P / 2)


def verify(pauli_string, theta_value=np.pi / 5):

    theta = Parameter("θ")

    qc = transpile_ppr(pauli_string, theta)

    qc_bound = qc.assign_parameters({theta: theta_value})

    U_circuit = Operator(qc_bound).data

    U_exact = exact_ppr(pauli_string, theta_value)

    ok = np.allclose(U_circuit, U_exact)

    print(f"\nVerification for {pauli_string}: {'PASS' if ok else 'FAIL'}")

    return ok


verify("XZZX")
verify("YZZY")
verify("XZZZZY")


# ============================================================
# OPTIONAL:
# Transpile a list of commuting PPRs
# ============================================================

def transpile_ppr_layer(pprs):
    """
    Parameters
    ----------
    pprs : list of tuples

        Example:
        [
            ("XZZX", theta1),
            ("YZZY", theta2),
        ]
    """

    n = len(pprs[0][0])

    qc = QuantumCircuit(n)

    for pauli_string, theta in pprs:

        qc.compose(
            transpile_ppr(pauli_string, theta),
            inplace=True
        )

    return qc


# ============================================================
# Example layered circuit
# ============================================================

theta1 = Parameter("θ1")
theta2 = Parameter("θ2")

layer = transpile_ppr_layer([
    ("XZZX", theta1),
    ("YZZY", theta2),
])

print("\n=== Layered PPR Circuit ===\n")
print(layer.draw("text"))