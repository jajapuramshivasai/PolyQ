import pytest
import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector
from PolyQ.engine import QC

def get_qiskit_sv(circuit):
    """Helper to get ground truth from Qiskit."""
    return Statevector(circuit).data

@pytest.mark.parametrize("num_qubits", [1, 2, 3])
def test_hadamard_gate(num_qubits):
    """Test single and multi-qubit Hadamard gates."""
    qiskit_qc = QuantumCircuit(num_qubits)
    for i in range(num_qubits):
        qiskit_qc.h(i)
    
    sim = QC(qiskit_qc)
    sim.compile()
    
    custom_sv = np.array([sim.get_amplitude(y) for y in range(2**num_qubits)])
    qiskit_sv = get_qiskit_sv(qiskit_qc)
    
    assert np.allclose(custom_sv, qiskit_sv, atol=1e-10)

def test_cz_gate():
    """Test Entanglement via CZ gate."""
    qiskit_qc = QuantumCircuit(2)
    qiskit_qc.h(0)
    qiskit_qc.cz(0, 1)
    
    sim = QC(qiskit_qc)
    sim.compile()
    
    custom_sv = np.array([sim.get_amplitude(y) for y in range(4)])
    qiskit_sv = get_qiskit_sv(qiskit_qc)
    
    assert np.allclose(custom_sv, qiskit_sv, atol=1e-10)

def test_phase_gates():
    """Test Z, S, and Sdg phase gates."""
    qiskit_qc = QuantumCircuit(1)
    qiskit_qc.h(0)
    qiskit_qc.s(0)
    qiskit_qc.sdg(0)
    qiskit_qc.z(0)
    
    sim = QC(qiskit_qc)
    sim.compile()
    
    custom_sv = np.array([sim.get_amplitude(y) for y in range(2)])
    qiskit_sv = get_qiskit_sv(qiskit_qc)
    
    assert np.allclose(custom_sv, qiskit_sv, atol=1e-10)

def test_unitarity():
    """Verify that the generated transition matrix is unitary (U* U = I)."""
    qiskit_qc = QuantumCircuit(2)
    qiskit_qc.h(0)
    qiskit_qc.cz(0, 1)
    qiskit_qc.h(1)
    qiskit_qc.s(1)
    
    sim = QC(qiskit_qc)
    sim.compile()
    
    U = sim.get_transition_matrix()
    identity_check = U.conj().T @ U
    assert np.allclose(identity_check, np.eye(4), atol=1e-10)

def test_zero_amplitude_condition():
    """Test circuits that result in destructive interference (0 amplitude)."""
    qiskit_qc = QuantumCircuit(1)
    qiskit_qc.h(0)
    qiskit_qc.z(0)
    qiskit_qc.h(0)
    # H Z H is equivalent to X gate. <0|X|0> should be 0.
    
    sim = QC(qiskit_qc)
    sim.compile()
    
    amplitude = sim.get_amplitude(0, 0)
    assert np.abs(amplitude) < 1e-10

def test_ghz():
    """Test generation of GHZ state."""
    num_qubits = 3
    qiskit_qc = QuantumCircuit(num_qubits)
    qiskit_qc.h(0)
    for i in range(1, num_qubits):
        qiskit_qc.cz(0, i)
        
    sim = QC(qiskit_qc)
    sim.compile()
    
    custom_sv = np.array([sim.get_amplitude(y) for y in range(2**num_qubits)])
    qiskit_sv = get_qiskit_sv(qiskit_qc)
    
    assert np.allclose(custom_sv, qiskit_sv, atol=1e-10)

def test_bell_state_11():
    """Test generation of specific Bell state combinations."""
    qiskit_qc = QuantumCircuit(2)
    qiskit_qc.h(0)
    qiskit_qc.h(1)
    qiskit_qc.cz(0, 1)
    qiskit_qc.z(1)
    qiskit_qc.h(1)
    qiskit_qc.h(0)
    qiskit_qc.z(0)
    qiskit_qc.h(0)
    
    sim = QC(qiskit_qc)
    sim.compile()
    
    custom_sv = np.array([sim.get_amplitude(y) for y in range(4)])
    qiskit_sv = get_qiskit_sv(qiskit_qc)
    
    assert np.allclose(custom_sv, qiskit_sv, atol=1e-10)