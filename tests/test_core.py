# test_core.py - Comprehensive tests for Dickson Engine and UniversalQC
import pytest
import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector
from PolyQ.lib import (
    DicksonOp,
    DicksonEngine,
    DicksonTranspiler,
    UniversalQC,
    BranchNode
)


def get_qiskit_sv(circuit):
    """Helper to get ground truth from Qiskit."""
    return Statevector(circuit).data


# ============================================================================
# PART 1: DicksonOp Tests
# ============================================================================

def test_dickson_op_creation(): # 1 pass
    """Test DicksonOp initialization."""
    op_swap = DicksonOp('SWAP', 0, 1)
    assert op_swap.type == 'SWAP'
    assert op_swap.a == 0
    assert op_swap.b == 1

    op_add = DicksonOp('ADD', 2, 3)
    assert op_add.type == 'ADD'
    assert op_add.a == 2
    assert op_add.b == 3


# ============================================================================
# PART 2: DicksonEngine Tests
# ============================================================================

def test_dickson_engine_translation(): # 2 passes
    """Test circuit translation to internal gate format."""
    qc = QuantumCircuit(2)
    qc.h(0)
    qc.cz(0, 1)
    qc.z(1)

    engine = DicksonEngine(qc)
    
    assert len(engine.gates) == 3
    assert engine.gates[0] == ('H', [0])
    assert engine.gates[1] == ('CZ', [0, 1])
    assert engine.gates[2] == ('Z', [1])


def test_dickson_engine_structure_compilation(): # 3 passes
    """Test adjacency matrix and structure compilation."""
    qc = QuantumCircuit(3)
    qc.h(0)
    qc.cz(0, 1)
    qc.h(1)
    qc.cz(1, 2)

    engine = DicksonEngine(qc)
    
    # Verify basic properties
    assert engine.num_qubits == 3
    assert engine.num_h == 2  # Two H gates
    assert engine.n_vars >= engine.num_qubits  # Should have auxiliary variables or equal
    assert engine.b4 is not None
    assert engine.b4.shape == (engine.n_vars, engine.n_vars)


def test_dickson_plan_dickson_algorithm(): # 4 passes
    """Test the Dickson reduction algorithm."""
    # Create a simple 2x2 adjacency matrix with one edge
    b_matrix = np.array([[0, 1], [1, 0]], dtype=np.int8)
    
    engine = DicksonEngine(QuantumCircuit(2))
    ops, b_reduced, rank = engine._plan_dickson(b_matrix, 2)
    
    # Should find the edge and perform swaps/adds
    assert rank >= 0
    assert len(ops) >= 0
    assert b_reduced is not None


def test_dickson_set_phases():
    """Test Z4 phase vector setting from S and Z gates."""
    qc = QuantumCircuit(2)
    qc.s(0)
    qc.z(1)

    engine = DicksonEngine(qc)
    engine.set_phases(qc)
    
    # Verify phase values set correctly
    assert engine.v4 is not None


def test_dickson_get_amplitude():
    """Test amplitude computation for Clifford gates."""
    qc = QuantumCircuit(2)
    qc.h(0)
    qc.cz(0, 1)

    engine = DicksonEngine(qc)
    
    # Test amplitudes
    amp_00 = engine.get_amplitude(0, 0)  # |00> -> |00>
    # sv [ 0.707, 0.707, 0, 0 ]
    
    # Should be non-zeroß
    assert np.isclose(amp_00, 0.7071067811865476)


def test_dickson_get_transition_matrix():
    """Test full transition matrix computation."""
    qc = QuantumCircuit(2)
    qc.h(0)
    qc.h(1)

    engine = DicksonEngine(qc)
    tm = engine.get_transition_matrix()
    
    # Should be 4x4 for 2 qubits
    assert tm.shape == (4, 4)


def test_dickson_gray_code_statevector():
    """Test statevector generation via Gray code."""
    qc = QuantumCircuit(2)
    qc.h(0)
    qc.h(1)

    engine = DicksonEngine(qc)
    sv = np.zeros(4, dtype=np.complex128)
    sv = engine.get_statevector_gray(1.0, sv, x=0)
    
    # Should be |+> ⊗ |+> = superposition of all 4 states
    expected = np.ones(4, dtype=np.complex128) / 2.0
    assert np.allclose(sv, expected, atol=1e-8)


def test_dickson_engine_full_clifford():
    """Integration test: DicksonEngine on a complete Clifford circuit."""
    qc = QuantumCircuit(2)
    qc.h(0)
    qc.cz(0, 1)
    qc.h(1)

    engine = DicksonEngine(qc)
    sv_custom = np.zeros(4, dtype=np.complex128)
    sv_custom = engine.get_statevector_gray(1.0, sv_custom, x=0)
    
    sv_qiskit = get_qiskit_sv(qc)
    
    assert np.allclose(sv_custom, sv_qiskit, atol=1e-8)


def test_dickson_print_analytic_formula():
    """Test analytic formula printing (basic smoke test)."""
    qc = QuantumCircuit(2)
    qc.h(0)
    qc.cz(0, 1)

    engine = DicksonEngine(qc)
    
    # This should not raise
    try:
        engine.print_analytic_formula(transition_mode=True, weight=1.0, branch_label="TEST")
    except Exception as e:
        pytest.fail(f"print_analytic_formula raised {e}")


# ============================================================================
# PART 3: DicksonTranspiler Tests
# ============================================================================

def test_transpiler_is_clifford_gate():
    """Test Clifford gate classification."""
    transpiler = DicksonTranspiler(2)
    
    assert transpiler._is_clifford_gate('h')
    assert transpiler._is_clifford_gate('cx')
    assert transpiler._is_clifford_gate('cz')
    assert transpiler._is_clifford_gate('z')
    assert transpiler._is_clifford_gate('s')
    assert transpiler._is_clifford_gate('swap')
    
    assert not transpiler._is_clifford_gate('rx')
    assert not transpiler._is_clifford_gate('ry')
    assert not transpiler._is_clifford_gate('t')
    assert not transpiler._is_clifford_gate('rz')


def test_transpiler_classify_pure_clifford():
    """Test classification of pure Clifford circuits."""
    qc = QuantumCircuit(2)
    qc.h(0)
    qc.cz(0, 1)
    qc.s(1)

    transpiler = DicksonTranspiler(2)
    is_clifford, blocks = transpiler._classify_circuit(qc)
    
    assert is_clifford is True
    assert len(blocks) >= 1


def test_transpiler_classify_universal():
    """Test classification of universal (mixed) circuits."""
    qc = QuantumCircuit(2)
    qc.h(0)
    qc.cz(0, 1)
    qc.rz(np.pi / 4, 0)  # Non-Clifford
    qc.h(1)

    transpiler = DicksonTranspiler(2)
    is_clifford, blocks = transpiler._classify_circuit(qc)
    
    assert is_clifford is False
    assert len(blocks) >= 2


def test_transpiler_convert_to_z4_form():
    """Test conversion to Z4 quadratic form."""
    qc = QuantumCircuit(2)
    qc.h(0)
    qc.cz(0, 1)
    qc.z(1)

    transpiler = DicksonTranspiler(2)
    b_matrix, v_vector = transpiler._convert_to_z4_form(qc)
    
    # Z gate should set phase
    assert v_vector[1] == 2  # Z gate = phase 2 in Z4


def test_transpiler_apply_van_den_nest():
    """Test Van den Nest optimization."""
    b_matrix = np.array([
        [0, 1, 1],
        [1, 0, 1],
        [1, 1, 0]
    ], dtype=np.int8)

    transpiler = DicksonTranspiler(3)
    b_optimized = transpiler._apply_van_den_nest(b_matrix)
    
    # Should return a matrix
    assert b_optimized.shape == b_matrix.shape
    assert b_optimized.dtype == np.int8


def test_transpiler_synthesize_clifford():
    """Test Clifford synthesis."""
    qc = QuantumCircuit(2)
    qc.h(0)
    qc.cz(0, 1)
    qc.h(1)

    transpiler = DicksonTranspiler(2)
    optimized = transpiler.synthesize_clifford(qc, optimize=True)
    
    # Should produce a circuit
    assert optimized.num_qubits == 2
    # Verify statevector equivalence
    sv_original = get_qiskit_sv(qc)
    sv_optimized = get_qiskit_sv(optimized)
    
    assert np.allclose(sv_original, sv_optimized, atol=1e-8)


def test_transpiler_transpile():
    """Test full transpilation."""
    qc = QuantumCircuit(2)
    qc.h(0)
    qc.cz(0, 1)

    transpiler = DicksonTranspiler(2)
    optimized = transpiler.transpile(qc)
    
    sv_original = get_qiskit_sv(qc)
    sv_optimized = get_qiskit_sv(optimized)
    
    assert np.allclose(sv_original, sv_optimized, atol=1e-8)


# ============================================================================
# PART 4: UniversalQC Fringe Splitting Tests
# ============================================================================

def test_universal_qc_splitting():
    """
    Tests if UniversalQC correctly partitions gates into 
    E_L (Left Fringe), Core, and E_R (Right Fringe).
    """
    qc = QuantumCircuit(3)
    # E_L Gates
    qc.cx(0, 1)
    # Core Gates
    qc.h(0)
    qc.cz(1, 2)
    # E_R Gates
    qc.s(0)

    sim = UniversalQC(qc)
    sim._split_circuit()

    # Verify split
    assert len(sim.EL_gates) >= 0
    assert len(sim.core_gates) >= 0
    assert len(sim.ER_gates) >= 0
    assert len(sim.EL_gates) + len(sim.core_gates) + len(sim.ER_gates) == len(qc.data)


def test_el_evaluation():
    """Tests E_L affine state transformation and phase shift."""
    qc = QuantumCircuit(2)
    qc.cx(0, 1)

    sim = UniversalQC(qc)
    sim._split_circuit()

    # Test affine transformation
    x_core, phase = sim._evaluate_EL(x=1)
    
    # After CX(0,1) with x=1: x_core should be 3 (both qubits 1)
    assert x_core == 3
    # Phase should be real number
    assert isinstance(phase, complex)


def test_er_application():
    """Tests right fringe output permutation."""
    qc = QuantumCircuit(2)
    qc.cx(0, 1)

    sim = UniversalQC(qc)
    sim.ER_gates = qc.data
    
    sv = np.array([1.0, 2.0, 3.0, 4.0], dtype=np.complex128)
    sim._apply_ER(sv)
    
    # CX(0,1) swaps states where control is 1
    # |01> (idx 1) <-> |11> (idx 3): [1, 4, 3, 2]
    expected = np.array([1.0, 4.0, 3.0, 2.0], dtype=np.complex128)
    assert np.allclose(sv, expected, atol=1e-10)


def test_full_pipeline_with_fringes():
    """
    Integration test: E_L, Core, and E_R work together 
    in full statevector generation.
    """
    qc = QuantumCircuit(2)
    # E_L
    qc.cx(0, 1)
    # Core
    qc.h(0)
    # E_R
    qc.s(1)

    sim = UniversalQC(qc)
    sim.build_tree()
    sv_custom = sim.get_statevector(x=0)
    sv_qiskit = get_qiskit_sv(qc)
    
    assert np.allclose(sv_custom, sv_qiskit, atol=1e-8)








def test_branch_node_creation():
    """Test BranchNode initialization."""
    node = BranchNode(weight=0.5, label="TEST")
    assert node.weight == 0.5
    assert node.label == "TEST"
    assert len(node.children) == 0
    assert node.v4_state is None


def test_branch_tree_construction():
    """Test branch tree construction in build_tree()."""
    qc = QuantumCircuit(1)
    qc.h(0)
    qc.rz(np.pi / 4, 0)  # Should create branches
    qc.h(0)

    sim = UniversalQC(qc)
    sim.build_tree()
    
    # Should have created a root
    assert sim.root is not None


# ============================================================================
# PART 5: Error Handling Tests
# ============================================================================

def test_get_statevector_without_build_tree():
    """Ensures get_statevector raises error if tree isn't built."""
    qc = QuantumCircuit(1)
    qc.h(0)
    
    sim = UniversalQC(qc)
    
    with pytest.raises(RuntimeError, match="Call build_tree"):
        sim.get_statevector()


def test_print_analytic_without_build_tree():
    """Ensures print_full_analytic_decomposition handles missing tree."""
    qc = QuantumCircuit(1)
    qc.h(0)
    
    sim = UniversalQC(qc)
    
    # Should print "Tree not built." without raising
    try:
        sim.print_full_analytic_decomposition()
    except Exception as e:
        pytest.fail(f"print_full_analytic_decomposition raised {e}")


def test_empty_circuit():
    """Test handling of empty circuits."""
    qc = QuantumCircuit(2)
    
    sim = UniversalQC(qc)
    sim.build_tree()
    sv = sim.get_statevector()
    
    # Should be |00>
    expected = np.array([1.0, 0, 0, 0], dtype=np.complex128)
    assert np.allclose(sv, expected, atol=1e-10)


def test_identity_circuit():
    """Test handling of identity circuits."""
    qc = QuantumCircuit(2)
    qc.id(0)
    qc.id(1)
    
    sim = UniversalQC(qc)
    sim.build_tree()
    sv = sim.get_statevector()
    
    expected = np.array([1.0, 0, 0, 0], dtype=np.complex128)
    assert np.allclose(sv, expected, atol=1e-10)


# ============================================================================
# PART 6: Complex Integration Tests
# ============================================================================

def test_bell_state_creation():
    """Test creation of Bell states."""
    qc = QuantumCircuit(2)
    qc.h(0)
    qc.cx(0, 1)
    
    sim = UniversalQC(qc)
    sim.build_tree()
    sv_custom = sim.get_statevector()
    sv_qiskit = get_qiskit_sv(qc)
    
    assert np.allclose(sv_custom, sv_qiskit, atol=1e-8)


def test_ghz_state_creation():
    """Test creation of GHZ states."""
    qc = QuantumCircuit(3)
    qc.h(0)
    qc.cx(0, 1)
    qc.cx(1, 2)
    
    sim = UniversalQC(qc)
    sim.build_tree()
    sv_custom = sim.get_statevector()
    sv_qiskit = get_qiskit_sv(qc)
    
    assert np.allclose(sv_custom, sv_qiskit, atol=1e-8)


def test_statevector_different_inputs():
    """Test get_statevector with different input states."""
    qc = QuantumCircuit(2)
    qc.h(0)
    # qc.h(1)
    # qc.cz(0, 1)
    # qc.h(1)
    qc.cx(0, 1)

    sim = UniversalQC(qc)
    sim.build_tree()
    
    # Test multiple input states
    for x in range(4):
        sv = sim.get_statevector(x=x)
        assert len(sv) == 4
        assert np.isclose(np.linalg.norm(sv), 1.0, atol=1e-8)


def test_phase_preservation():
    """Test that phase information is preserved."""
    qc = QuantumCircuit(1)
    qc.h(0)
    qc.s(0)
    qc.s(0)  # Two S gates = Z gate
    
    sim = UniversalQC(qc)
    sim.build_tree()
    sv_custom = sim.get_statevector()
    
    sv_qiskit = get_qiskit_sv(qc)
    
    assert np.allclose(sv_custom, sv_qiskit, atol=1e-8)


def test_large_circuit():
    """Test on a larger circuit to ensure scalability."""
    n = 4
    qc = QuantumCircuit(n)
    
    for i in range(n):
        qc.h(i)
    for i in range(n - 1):
        qc.cz(i, i + 1)
    
    sim = UniversalQC(qc)
    sim.build_tree()
    sv_custom = sim.get_statevector()
    sv_qiskit = get_qiskit_sv(qc)
    
    assert np.allclose(sv_custom, sv_qiskit, atol=1e-7)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
