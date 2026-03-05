use crate::sim::PhasePolynomial;

/// Public representation of a gate in a quantum circuit.
///
/// This enum is re-exported from the library root so downstream
/// users can inspect or construct gates if they wish (e.g. for
/// testing or serialization).  Internally we also use it as the
/// "raw" gate type for the optimizer and circuit builder.
#[derive(Clone, Debug, PartialEq)]
pub enum Gate {
    H(usize),
    S(usize),
    Z(usize),
    CZ(usize, usize),
    T(usize),
}

// Convenient alias for the internal optimizer; using the same
// enum keeps conversions trivial.

struct CircuitOptimizer {
    raw_gates: Vec<Gate>,
}

impl CircuitOptimizer {
    fn new() -> Self {
        Self { raw_gates: Vec::new() }
    }

    fn add_gate(&mut self, gate: Gate) {
        self.raw_gates.push(gate);
    }

    /// Performs a pass to eliminate H-H pairs on the same qubit
    fn optimize_hadamards(&mut self) {
        let mut optimized = Vec::new();
        let mut last_h_on_qubit: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();

        for gate in self.raw_gates.drain(..) {
            match gate {
                Gate::H(q) => {
                    if let Some(&idx) = last_h_on_qubit.get(&q) {
                        // Check if any gates in between 'idx' and now commute/blocked
                        // For simplicity, this pass removes immediate adjacent H-H
                        if idx == optimized.len() - 1 {
                            optimized.remove(idx);
                            last_h_on_qubit.remove(&q);
                        } else {
                            optimized.push(gate);
                            last_h_on_qubit.insert(q, optimized.len() - 1);
                        }
                    } else {
                        optimized.push(gate);
                        last_h_on_qubit.insert(q, optimized.len() - 1);
                    }
                },
                _ => {
                    // Non-H gate blocks the H-H cancellation for that qubit in this simple pass
                    optimized.push(gate);
                    if let Gate::CZ(q1, q2) = &optimized[optimized.len()-1] {
                        last_h_on_qubit.remove(q1);
                        last_h_on_qubit.remove(q2);
                    } else if let Gate::S(q) | Gate::T(q) | Gate::Z(q) = &optimized[optimized.len()-1] {
                        // In Clifford+T, S and T don't commute with H
                        last_h_on_qubit.remove(q);
                    }
                }
            }
        }
        self.raw_gates = optimized;
    }

    /// Converted optimized gates to the PhasePolynomial structure
    fn build_poly_simulator(&self, n: usize) -> PhasePolynomial {
        let mut sim = PhasePolynomial::new(n);
        for gate in &self.raw_gates {
            match *gate {
                Gate::H(q) => sim.apply_clifford("h", &[q]),
                Gate::S(q) => sim.apply_clifford("s", &[q]),
                Gate::Z(q) => sim.apply_clifford("z", &[q]),
                Gate::CZ(q1, q2) => sim.apply_clifford("cz", &[q1, q2]),
                Gate::T(q) => sim.t_gate_vars.push(sim.output_mapping[q]),
            }
        }
        sim
    }
}

/// Simple builder for a quantum circuit that consumers of the library
/// will normally interact with.  It records gates in order and provides
/// a convenient method to turn the finished circuit into the polynomial
/// simulator defined in `sim.rs`.
#[derive(Clone, Debug)]
pub struct Circuit {
    num_qubits: usize,
    gates: Vec<Gate>,
}

impl Circuit {
    /// Create a fresh circuit acting on `n` qubits.
    pub fn new(n: usize) -> Self {
        Self { num_qubits: n, gates: Vec::new() }
    }

    /// Add a Hadamard to the given wire.
    pub fn h(&mut self, q: usize) -> &mut Self {
        self.gates.push(Gate::H(q));
        self
    }

    /// Add an S gate.
    pub fn s(&mut self, q: usize) -> &mut Self {
        self.gates.push(Gate::S(q));
        self
    }

    /// Add a Z gate.
    pub fn z(&mut self, q: usize) -> &mut Self {
        self.gates.push(Gate::Z(q));
        self
    }

    /// Add a controlled-Z between `q1` and `q2`.
    pub fn cz(&mut self, q1: usize, q2: usize) -> &mut Self {
        self.gates.push(Gate::CZ(q1, q2));
        self
    }

    /// Add a controlled-X gate.  We decompose using the standard
    /// Clifford identity `CNOT = (I ⊗ H) CZ (I ⊗ H)`.
    pub fn cx(&mut self, control: usize, target: usize) -> &mut Self {
        // Note: the extra `.h(target)` at the end is critical!
        self.h(target);
        self.cz(control, target);
        self.h(target)
    }

    /// Add a T gate.
    pub fn t(&mut self, q: usize) -> &mut Self {
        self.gates.push(Gate::T(q));
        self
    }

    /// Add an X gate (H Z H decomposition).
    pub fn x(&mut self, q: usize) -> &mut Self {
        self.h(q).z(q).h(q)
    }

    /// Consumes the circuit and produces a simulator instance ready
    /// to compute amplitudes.
    pub fn into_simulator(self) -> PhasePolynomial {
        let mut optimizer = CircuitOptimizer::new();
        for gate in &self.gates {
            optimizer.add_gate(gate.clone());
        }
        // perform the hadamard cancellation pass before generating the
        // polynomial representation.  this ensures the internal optimizer
        // API remains exercised and avoids dead-code warnings.
        optimizer.optimize_hadamards();

        let mut sim = optimizer.build_poly_simulator(self.num_qubits);
        // remember the unoptimized list so `simulate` can fall back
        sim.original_gates = Some(self.gates);
        sim
    }
}
