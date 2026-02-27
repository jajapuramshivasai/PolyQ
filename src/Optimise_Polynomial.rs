//todo: seperate poly into Z8+z4+z2 which can be sent to their respective engine


mod Optimise_Polynomial {
    use crate::quantum_circuit::QuantumCircuit;
    use crate::Z2::engine::phase_polynomial;
    use crate::Z4::engine::phase_polynomial_z4;
    use crate::Z8::engine::phase_polynomial_z8;

    /// Optimise a quantum circuit by simplifying its phase polynomial.
    pub fn optimise_circuit(circuit: &QuantumCircuit) -> QuantumCircuit {
        // For simplicity, we will just return the original circuit.
        // In a real implementation, we would extract the phase polynomial,
        // simplify it, and then reconstruct an optimized circuit.
        circuit.clone()
    }
}