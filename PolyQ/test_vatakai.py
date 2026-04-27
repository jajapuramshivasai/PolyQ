from qiskit import QuantumCircuit
from lib import UniversalQC

def test_analytic_decomposition():
    qc = QuantumCircuit(2)
    qc.h(0)
    qc.t(0)  # Causes branching in the middle
    qc.h(1)
    qc.cz(0, 1)
    qc.h(1)
    
    qc.h(0)
    qc.s(0)
    qc.h(0)
    
    qc.h(0)
    qc.z(0)
    qc.h(0)

    uqc = UniversalQC(qc)
    uqc.build_tree()
    
    # This will print the Weight, ID/Z path, and the XOR formulas for each branch
    uqc.print_full_analytic_decomposition(transition_mode=False)

if __name__ == "__main__":
    test_analytic_decomposition()