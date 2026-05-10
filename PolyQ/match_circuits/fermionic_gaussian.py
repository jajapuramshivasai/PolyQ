import numpy as np
import re

def is_matchgate_matrix(B, tol=1e-8):
    """
    Helper: Tests if a single 4x4 matrix B is a 2-qubit Matchgate.
    """
    if B.shape != (4, 4):
        return False

    # Matchgate Identities (M1 to M5)
    M1 = (B[0,0]*B[3,3] - B[2,2]*B[1,1] - B[0,3]*B[3,0] + B[2,1]*B[1,2])
    M2 = (B[2,0]*B[3,3] - B[2,2]*B[3,1] - B[3,0]*B[2,3] + B[2,1]*B[3,2])
    M3 = (B[1,0]*B[3,3] + B[1,1]*B[3,2] - B[3,0]*B[1,3] - B[1,2]*B[3,1])
    M4 = (B[0,1]*B[3,3] + B[1,1]*B[2,3] - B[0,3]*B[3,1] - B[2,1]*B[1,3])
    M5 = (B[0,2]*B[3,3] - B[2,2]*B[1,3] - B[0,3]*B[3,2] + B[1,2]*B[2,3])

    identities_zero = all(np.isclose(m, 0, atol=tol) for m in [M1, M2, M3, M4, M5])
    
    b44_nonzero = not np.isclose(B[3,3], 0, atol=tol)
    is_diagonal = np.allclose(B, np.diag(np.diagonal(B)), atol=tol)
    
    return identities_zero and (b44_nonzero or is_diagonal)


def is_flo_pauli_generator(pauli_string):
    """
    Helper: Tests if a single Pauli string belongs to the L_2 Lie Algebra.
    """
    pauli_string = pauli_string.upper()
    if not set(pauli_string).issubset({'I', 'X', 'Y', 'Z'}):
        return False

    # L_2 valid patterns: Single Paulis or U(Z...Z)V where U,V in {X,Y}
    pattern = r'^I*([XYZ]|[XY]Z*[XY])?I*$'
    return bool(re.match(pattern, pauli_string))


# ==========================================
# Circuit-Level Testers
# ==========================================

def is_match_circuit(circuit_matrices, tol=1e-8):
    """
    Tests if an entire circuit (list of 4x4 matrices) is a valid Matchcircuit.
    Returns True if every matrix in the circuit is a valid matchgate.
    """
    invalid_indices = []
    for idx, matrix in enumerate(circuit_matrices):
        if not is_matchgate_matrix(matrix, tol):
            invalid_indices.append(idx)
            
    if invalid_indices:
        print(f"Circuit is NOT a Matchcircuit. Invalid gates at indices: {invalid_indices}")
        return False
    
    print("Circuit IS a valid Matchcircuit.")
    return True


def is_flo_circuit(circuit_ppr_strings):
    """
    Tests if an entire circuit (list of Pauli strings) is a valid FLO circuit.
    Returns True if every PPR gate in the circuit belongs to the L_2 Lie algebra.
    """
    invalid_gates = []
    for idx, p_string in enumerate(circuit_ppr_strings):
        if not is_flo_pauli_generator(p_string):
            invalid_gates.append((idx, p_string))
            
    if invalid_gates:
        print("Circuit is NOT a valid FLO circuit. Invalid gates found:")
        for idx, p_string in invalid_gates:
            print(f"  - Gate {idx}: '{p_string}'")
        return False
        
    print("Circuit IS a valid FLO circuit.")
    return True


# ==========================================
# Example Usage
# ==========================================
if __name__ == "__main__":
    
    print("--- Testing PPR (Pauli Product Rotation) Circuits ---")
    
    # Circuit 1: All valid FLO operators
    valid_ppr_circuit = [
        "IXXI",    # X_2 X_3
        "IXZZYI",  # X_2 Z_3 Z_4 Y_5
        "IZII",    # Z_2
        "YI",      # Y_1
    ]
    is_flo_circuit(valid_ppr_circuit)
    # Output: Circuit IS a valid FLO circuit.

    print("\n")
    
    # Circuit 2: Contains some invalid FLO operators
    invalid_ppr_circuit = [
        "IXXI",    # Valid
        "IXYXI",   # Invalid: Y inside the boundary
        "IXXZI",   # Invalid: Does not end in X or Y
        # "II"       # Valid
        "IIZZI"
    ]
    is_flo_circuit(invalid_ppr_circuit)
    # Output: Circuit is NOT a valid FLO circuit. Invalid gates found: 1 and 2.

    print("\n--- Testing Matrix Circuits ---")
    
    # A list of matrices (e.g., the Identity and a random matrix)
    valid_matrix = np.eye(4)
    invalid_matrix = np.random.rand(4, 4)
    
    circuit_matrices = [valid_matrix, valid_matrix, invalid_matrix]
    
    is_match_circuit(circuit_matrices)
    # Output: Circuit is NOT a Matchcircuit. Invalid gates at indices: [2]