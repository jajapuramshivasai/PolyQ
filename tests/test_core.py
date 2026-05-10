
import pytest
import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector
from PolyQ.lib import UniversalQC

def get_qiskit_sv(circuit):
    """Helper to get ground truth from Qiskit."""
    return Statevector(circuit).data

# --- TEST PART 4: FRINGE SPLITTING & AFFINE TRANSFORMS ---

def test_universal_qc_splitting():
    """
    Tests if the UniversalQC correctly partitions gates into 
    E_L (Left Fringe), Core, and E_R (Right Fringe).
    """
    qc = QuantumCircuit(3)
    # E_L Gates (No superpositions yet)
    qc.cx(0, 1)
    qc.swap(1, 2)
    qc.t(0)
    
    # Core Gates (Superposition / Branching)
    qc.h(0)
    qc.cz(1, 2) # note this is pulled into E_l rightfully so
    qc.h(1)
    
    # E_R Gates (Fringe after core)
    qc.cx(1, 2)
    qc.s(0)

    sim = UniversalQC(qc)
    sim._split_circuit()

    # Verify counts based on the splitting logic
    assert len(sim.EL_gates) == 4
    assert len(sim.core_gates) == 2
    assert len(sim.ER_gates) == 2

def test_el_evaluation():
    """Tests the exact calculation of the E_L affine state and phase shift."""
    qc = QuantumCircuit(2)
    qc.cx(0, 1)        # x_1 ^= x_0
    qc.z(1)            # Phase flip if x_1 is 1
    qc.rz(np.pi, 0)    # Phase shift based on x_0
    
    sim = UniversalQC(qc)
    sim._split_circuit()

    # Input x = 1 (|01> in integer representation, meaning q0=1, q1=0)
    # After CX(0, 1): q0=1, q1=1 -> x_core = 3
    # Phase: Z(1) adds pi. Rz(pi, 0) adds pi/2 (since q0=1). Total phase = 1.5pi -> -1j
    
    x_core, phase = sim._evaluate_EL(x=1)
    
    assert x_core == 3
    # Check phase (exp(1.5j * pi) is approx -1j)
    assert np.allclose(phase, -1j, atol=1e-10)

def test_er_application():
    """Tests the right fringe output permutation and phase adjustments."""
    qc = QuantumCircuit(2)
    # Give it an empty core to only test E_R
    qc.cx(0, 1)
    qc.swap(0, 1)
    qc.s(0)
    
    sim = UniversalQC(qc)
    # Forcibly set ER gates to bypass splitting logic for unit testing
    sim.ER_gates = qc.data 
    
    # Start with a dummy statevector: [1, 2, 3, 4]
    sv = np.array([1.0, 2.0, 3.0, 4.0], dtype=np.complex128)
    sim._apply_ER(sv)
    
    # Ground truth simulation
    expected_sv = np.array([1.0, 2.0, 3.0, 4.0], dtype=np.complex128)
    # Apply CX(0, 1): |01> (1) <-> |11> (3) -> [1, 4, 3, 2]
    # Apply SWAP(0, 1): |01> (1) <-> |10> (2) -> [1, 3, 4, 2]
    # Apply S(0): multiply indices where q0=1 (idx 1 and 3) by 1j -> [1, 3j, 4, 2j]
    expected_sv = np.array([1.0, 3.0j, 4.0, 2.0j], dtype=np.complex128)
    
    assert np.allclose(sv, expected_sv, atol=1e-10)

def test_full_pipeline_with_fringes():
    """
    Integration test based on df.py to ensure E_L, Core, and E_R 
    work together smoothly in the full statevector generation.
    """
    qc = QuantumCircuit(4)
    # E_L
    qc.cx(0, 1)
    qc.swap(1, 2)
    qc.t(0)
    qc.rz(np.pi / 3, 2)
    # Core U
    qc.h(0); qc.h(1); qc.h(2)
    qc.cz(0, 1)
    qc.cx(1, 2) 
    qc.rz(np.pi / 4, 1) 
    qc.h(0); qc.h(1); qc.h(2)
    # E_R
    qc.cx(2, 3)
    qc.swap(0, 1)
    qc.s(3)

    sim = UniversalQC(qc)
    sim.build_tree()
    sv_custom = sim.get_statevector(x=0)
    sv_qiskit = get_qiskit_sv(qc)
    
    assert np.allclose(sv_custom, sv_qiskit, atol=1e-8)

def test_global_phase_accumulation():
    """Verifies that internal core T-gates accurately accumulate global phase."""
    qc = QuantumCircuit(1)
    qc.h(0)
    qc.t(0) # T gate in the core creates global phase tracking
    qc.h(0)
    
    sim = UniversalQC(qc)
    sim.build_tree()
    
    # T gate should add a global phase of pi/8
    assert np.isclose(sim.global_phase, np.pi / 8, atol=1e-10)
    
    sv_custom = sim.get_statevector()
    sv_qiskit = get_qiskit_sv(qc)
    
    assert np.allclose(sv_custom, sv_qiskit, atol=1e-10)

# --- TEST PART 5: ERROR HANDLING ---

def test_get_statevector_without_build_tree():
    """Ensures get_statevector raises an error if the tree isn't built."""
    qc = QuantumCircuit(1)
    qc.h(0)
    
    sim = UniversalQC(qc)
    
    with pytest.raises(RuntimeError, match="Call build_tree.*"):
        sim.get_statevector()