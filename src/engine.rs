use crate::dataset::Dataset;
use crate::simulation::SimulationResult;
use crate::gateset::GateSet;

/// Simulation engine configuration and entrypoints.
#[derive(Clone, Debug)]
pub struct Engine<G: GateSet> {
    pub gateset: G,
}

impl<G: GateSet> Engine<G> {
    pub fn new(gateset: G) -> Self {
        Self { gateset }
    }

    /// Run a simple simulation for a dataset using the default simulation module.
    pub fn run(&self, dataset: &Dataset) -> SimulationResult {
        crate::simulation::simulate(self, dataset)
    }
    /// Get the modulus for the current gateset
    pub fn modulus(&self) -> i64 {
        self.gateset.modulus()
    }
    /// Map a gate to its polynomial using the gateset
    pub fn gate_to_polynomial(&self, gate: &str, qubits: &[usize]) -> Option<crate::polynomial::Polynomial> {
        self.gateset.gate_to_polynomial(gate, qubits)
    }
}
