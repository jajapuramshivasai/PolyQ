use crate::dataset::Dataset;
use crate::simulation::SimulationResult;

/// Simulation engine configuration and entrypoints.
#[derive(Clone, Debug)]
pub struct Engine {
    pub modulus: i64,
}

impl Engine {
    pub fn new(modulus: i64) -> Self {
        Self { modulus }
    }

    /// Run a simple simulation for a dataset using the default simulation module.
    pub fn run(&self, dataset: &Dataset) -> SimulationResult {
        crate::simulation::simulate(self, dataset)
    }
}
