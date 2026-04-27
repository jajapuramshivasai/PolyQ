import pytest
import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector
# Ensure your merged classes are accessible (e.g., from optimized_sim.py)
from lib import DicksonEngine, UniversalQC 

def get_qiskit_sv(circuit):
    """Helper to get ground truth from Qiskit."""
    return Statevector(circuit).data

# --- TEST PART 1: CORE ENGINE & GRAY CODE OPTIMIZATION ---

@pytest.mark.parametrize("num_qubits", [2, 3, 4])
def test_engine_consistency(num_qubits):
    """
    Verifies that the Gray code statevector evaluation matches 
    the point-wise get_amplitude results and Qiskit.
    """
    qc = QuantumCircuit(num_qubits)
    for i in range(num_qubits): qc.h(i)
    qc.cz(0, 1)
    if num_qubits > 2: qc.s(2)
    
    engine = DicksonEngine(qc)
    engine.set_phases(qc)
    
    # 1. Gray code statevector (Optimized)
    sv_gray = engine.get_statevector_gray(np.zeros(2**num_qubits, dtype=complex))
    
    # 2. Manual amplitude collection
    sv_manual = np.array([engine.get_amplitude(y) for y in range(2**num_qubits)])
    
    # 3. Ground Truth
    sv_qiskit = get_qiskit_sv(qc)
    
    assert np.allclose(sv_gray, sv_manual, atol=1e-10)
    assert np.allclose(sv_gray, sv_qiskit, atol=1e-10)

def test_engine_structural_cache():
    """Verifies that the engine can be reused with different phase gates (v4)."""
    skeleton = QuantumCircuit(2)
    skeleton.h(0)
    skeleton.cz(0, 1)
    
    engine = DicksonEngine(skeleton)
    
    # Test Case 1: No phase gates
    qc1 = skeleton.copy()
    engine.set_phases(qc1)
    sv1 = engine.get_statevector_gray(np.zeros(4, dtype=complex))
    assert np.allclose(sv1, get_qiskit_sv(qc1))
    
    # Test Case 2: Update phases (S gate) without recompiling structure
    qc2 = skeleton.copy(); qc2.s(0)
    engine.set_phases(qc2) # Only updates v4
    sv2 = engine.get_statevector_gray(np.zeros(4, dtype=complex))
    assert np.allclose(sv2, get_qiskit_sv(qc2))

# --- TEST PART 2: UNIVERSAL SIM & PHASE CONDENSATION ---

def test_phase_condensation():
    """
    Tests if input/output zone phases are correctly condensed 
    and applied outside the branching sum.
    """
    qc = QuantumCircuit(2)
    qc.rz(np.pi/4, 0) # Input Zone (Condensed)
    qc.h(0)
    qc.cz(0, 1)
    qc.h(1)
    qc.s(1)           # Output Zone (Condensed)
    
    sim = UniversalQC(qc)
    sim.build_tree()
    
    # Verify tree depth: Since no internal non-Clifford gates exist, 
    # the tree should have only 1 branch (root).
    assert len(sim.root.children) == 0
    
    sv_custom = sim.get_statevector()
    sv_qiskit = get_qiskit_sv(qc)
    assert np.allclose(sv_custom, sv_qiskit, atol=1e-10)

def test_magic_state_branching():
    """Test full branching with internal T gates and Rz rotations."""
    qc = QuantumCircuit(2)
    qc.h(0)
    qc.rz(np.pi/3, 0) # Internal Branching (between Hadamards)
    qc.h(0)
    
    sim = UniversalQC(qc)
    sim.build_tree()
    
    # Should have 2 branches (I and Z) for the Rz gate
    assert len(sim.root.children) == 2
    
    sv_custom = sim.get_statevector()
    sv_qiskit = get_qiskit_sv(qc)
    assert np.allclose(sv_custom, sv_qiskit, atol=1e-10)

# --- TEST PART 3: DIAGNOSTICS ---

def test_diagnostic_prints():
    """Ensures analytical and parameter printing works for branched circuits."""
    qc = QuantumCircuit(2)
    qc.h(0)
    qc.t(0) # Non-Clifford
    qc.h(0)
    
    sim = UniversalQC(qc)
    sim.build_tree()
    
    # These should execute without error
    sim.print_analytic_formula()
    
    # Test transition matrix with condensation
    tm = sim.engine.get_transition_matrix()
    assert tm.shape == (4, 4)
    
    
def test_bell_11():
    """Test Bell state |11> creation with internal T gate."""
    qc = QuantumCircuit(2)
    qc.h(0)
    
    qc.h(1)
    qc.cz(0, 1)
    qc.z(0)
    qc.h(1)
    
    qc.h(0)
    qc.z(0)
    qc.h(0)
    
    
    sim = UniversalQC(qc)
    sim.build_tree()
    
    sv_custom = sim.get_statevector()
    sv_qiskit = get_qiskit_sv(qc)
    # sim.print_analytic_formula()
    # sim.transition_matrix_print()

    assert np.allclose(sv_custom, sv_qiskit, atol=1e-10)

if __name__ == "__main__":
    # If running manually without pytest
    print("Running Manual Verification...")
    test_phase_condensation()
    test_magic_state_branching()
    print("Verification Successful!")
    
    test_bell_11()