import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector
from sim import UniversalQC # Assuming the class is in universal.py

def run_universal_test():
    # 1. Create a "Magic" Circuit (Clifford + T + Rz)
    num_qubits = 2
    qc = QuantumCircuit(num_qubits)
    
    # Clifford part
    qc.h(0)
    qc.h(1)
    qc.cz(0, 1)
    
    # Non-Clifford part (Magic gates)
    qc.rz(np.pi/8, 0)                # Arbitrary Rz rotation
    qc.rz(np.pi/8, 1)
    
    qc.h(0)
    
    print("--- Circuit Summary ---")
    print(qc)

    # 2. Initialize Universal Simulator
    sim = UniversalQC(qc)
    
    # 3. Build the Execution Tree with Pruning
    # threshold=0.01 will prune very small branches
    print("\nBuilding execution tree...")
    sim.build_tree(threshold=0.000001) 
    
    # 4. Calculate Amplitude for state |11> (y=3, x=0)
    target_state = 3 # Binary 11
    print(f"\nCalculating amplitude for state |{bin(target_state)[2:]:>02}>...")
    
    # Test Full Summation
    amp_full = sim.get_amplitude(y_val=target_state)
    
    # Test Top-N Approximation (using only the most significant path)
    amp_top1 = sim.get_amplitude(y_val=target_state, top_n=1)
    
    # 5. Get Ground Truth from Qiskit
    qiskit_sv = Statevector(qc).data
    expected_amp = qiskit_sv[target_state]

    # 6. Results Validation
    print(f"{'='*40}")
    print(f"Universal Sim (Full):  {amp_full:.6f}")
    print(f"Universal Sim (Top-1): {amp_top1:.6f}")
    print(f"Qiskit Ground Truth:   {expected_amp:.6f}")
    print(f"{'='*40}")

    is_correct = np.allclose(amp_full, expected_amp, atol=1e-8)
    print(f"Test Passed: {is_correct}")

    # 7. Visualization
    print("\nGenerating Tree Visualization...")
    sim.visualize()

if __name__ == "__main__":
    run_universal_test()