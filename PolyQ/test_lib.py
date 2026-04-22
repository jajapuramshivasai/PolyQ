import unittest
import numpy as np
from PolyQ.lib import DicksonOp, QC, _HAS_QISKIT
#python -m unittest PolyQ/test_lib.py
class TestDicksonOp(unittest.TestCase):
    def test_swap_operation(self):
        op = DicksonOp("SWAP", a=1, b=2)
        self.assertEqual(op.type, "SWAP")
        self.assertEqual(op.a, 1)
        self.assertEqual(op.b, 2)

    def test_add_operation(self):
        op = DicksonOp("ADD", a=None, pivot=3, target=4)
        self.assertEqual(op.type, "ADD")
        self.assertEqual(op.a, 3)
        self.assertEqual(op.b, 4)

class TestQC(unittest.TestCase):
    def setUp(self):
        self.qc = QC(num_qubits=3)

    def test_gate_constructors(self):
        self.qc.h(0)
        self.qc.z(1)
        self.qc.s(2)
        self.qc.cz(0, 1)
        self.qc.rz(2, np.pi / 4)
        self.assertEqual(len(self.qc.gates), 5)

    def test_compile_clifford_only(self):
        self.qc.h(0)
        self.qc.cz(0, 1)
        self.qc.s(2)
        self.qc.compile()
        self.assertTrue(self.qc.compiled)

    def test_get_amplitude_clifford(self):
        self.qc.h(0)
        self.qc.cz(0, 1)
        self.qc.compile()
        amp = self.qc.get_amplitude(1)
        self.assertIsInstance(amp, complex)

    def test_get_statevector_0(self):
        self.qc.h(0)
        self.qc.cz(0, 1)
        self.qc.compile()
        sv = self.qc.get_statevector_0()
        self.assertEqual(len(sv), 2 ** self.qc.num_qubits)

    def test_branching_with_rz(self):
        self.qc.h(0)
        self.qc.rz(1, np.pi / 4)
        branches = self.qc._build_clifford_branches(threshold=0.1)
        self.assertGreater(len(branches), 0)

    def test_verify_against_qiskit(self):
        self.qc.h(0)
        self.qc.cz(0, 1)
        self.qc.rz(1, np.pi / 4)
        self.qc.verify_against_qiskit()

# if __name__ == "__main__":
#     unittest.main()