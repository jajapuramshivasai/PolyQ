use crate::engine::Engine;
use crate::gateset::GateSet;
use crate::dataset::Dataset;

/// Straightforward simulation result produced by the simulation module.
#[derive(Clone, Debug)]
pub struct SimulationResult {
    pub dataset_name: String,
    pub value: i64,
}

/// A minimal simulation: evaluate the dataset polynomial at x=2 modulo engine.modulus.
// ...existing code...
pub fn simulate<G: GateSet>(engine: &Engine<G>, dataset: &Dataset) -> SimulationResult {
    let value = dataset.polynomial.eval(2, Some(engine.modulus()));
    SimulationResult {
        dataset_name: dataset.name.clone(),
        value,
    }
}
