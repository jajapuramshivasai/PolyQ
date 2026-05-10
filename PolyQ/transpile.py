import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import Clifford
from qiskit.synthesis.clifford import synth_clifford_layers

def van_den_nest_transpile(circuit: QuantumCircuit):
    """
    Decomposes a Clifford circuit into Van den Nest normal form:
    C1 (H-free) -> H_A (Hadamard layer) -> C2 (H-free)
    """
    # 1. Convert input circuit to a Clifford operator object
    # This automatically simplifies the representation to its functional core
    cliff = Clifford.from_circuit(circuit)
    
    # 2. Synthesize using the layered approach
    # This produces a circuit consisting of: 
    # [Pauli, S, CZ, CX, H, S, CZ, CX] layers
    decomposed_qc = synth_clifford_layers(cliff)
    
    return decomposed_qc

# --- Example Usage ---

# Create a random/complex Clifford circuit
n_qubits = 3
qc = QuantumCircuit(n_qubits)
qc.h(0)
qc.cx(0, 1)
qc.s(1)
qc.h(2)
qc.cx(1, 2)
qc.h(1)
qc.cz(0, 2)

print("Original Circuit Depth:", qc.depth())
print(qc.draw(output='text'))
# Transpile to Normal Form
normal_form_qc = van_den_nest_transpile(qc)

print("\nNormal Form Circuit (Van den Nest Structure):")
print(normal_form_qc.decompose().draw(output='text'))

# Verify Equivalence
original_cliff = Clifford(qc)
new_cliff = Clifford(normal_form_qc)
print("\nIs the decomposition equivalent?", original_cliff == new_cliff)