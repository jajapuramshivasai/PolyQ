//Clifford Gate accepted by poly q {H,Z,CZ}
//Clifford+S Gate accepted by poly q {H,Z,CZ,S}
//Clifford+T Gate accepted by poly q {H,Z,CZ,S,T}

use crate::quantum_circuit::{Gate, QuantumCircuit, QuantumCircuitClass};

/// Transpile a circuit to the accepted gateset for the given class.
pub fn transpile_to_gateset(
    circuit: &QuantumCircuit,
    class: QuantumCircuitClass,
) -> QuantumCircuit {
    let accepted_gates = match class {
        QuantumCircuitClass::Clifford => vec!["h", "z", "cz"],
        QuantumCircuitClass::CliffordPlusS => vec!["h", "z", "cz", "s"],
        QuantumCircuitClass::CliffordPlusT => vec!["h", "z", "cz", "t"],
    };
    let mut new_gates = Vec::new();
    for gate in &circuit.gates {
        let g = gate.name.to_lowercase();
        if accepted_gates.contains(&g.as_str()) {
            new_gates.push(gate.clone());
        } else if g == "x" && gate.qubits.len() == 1 {
            // X(q) = H(q) Z(q) H(q)
            let q = gate.qubits[0];
            new_gates.push(Gate {
                name: "h".to_string(),
                qubits: vec![q],
                params: vec![],
            });
            new_gates.push(Gate {
                name: "z".to_string(),
                qubits: vec![q],
                params: vec![],
            });
            new_gates.push(Gate {
                name: "h".to_string(),
                qubits: vec![q],
                params: vec![],
            });
        } else if g == "cx" && gate.qubits.len() == 2 {
            // CX(control, target) = H(target) CZ(control, target) H(target)
            let control = gate.qubits[0];
            let target = gate.qubits[1];
            new_gates.push(Gate {
                name: "h".to_string(),
                qubits: vec![target],
                params: vec![],
            });
            new_gates.push(Gate {
                name: "cz".to_string(),
                qubits: vec![control, target],
                params: vec![],
            });
            new_gates.push(Gate {
                name: "h".to_string(),
                qubits: vec![target],
                params: vec![],
            });
        } else if g == "y" && gate.qubits.len() == 1 {
            // Y(q) = S(q) H(q) X(q) H(q) S†(q)
            let q = gate.qubits[0];
            new_gates.push(Gate {
                name: "s".to_string(),
                qubits: vec![q],
                params: vec![],
            });
            new_gates.push(Gate {
                name: "h".to_string(),
                qubits: vec![q],
                params: vec![],
            });
            // X(q) = H(q) Z(q) H(q)
            new_gates.push(Gate {
                name: "h".to_string(),
                qubits: vec![q],
                params: vec![],
            });
            new_gates.push(Gate {
                name: "z".to_string(),
                qubits: vec![q],
                params: vec![],
            });
            new_gates.push(Gate {
                name: "h".to_string(),
                qubits: vec![q],
                params: vec![],
            });
            new_gates.push(Gate {
                name: "h".to_string(),
                qubits: vec![q],
                params: vec![],
            });
            new_gates.push(Gate {
                name: "sdag".to_string(),
                qubits: vec![q],
                params: vec![],
            });
        } else {
        }
    }
    QuantumCircuit {
        num_qubits: circuit.num_qubits,
        gates: new_gates,
    }
}

/*
FUNCTION Gridsynth(Target_Unitary U, Precision epsilon):

    1. CONVERSION TO Z[i, 1/sqrt(2)]
       Map the target rotation U to a candidate point (u, v)
       where |u|^2 + |v|^2 = 1.

    2. DETERMINE SEARCH REGION
       Define a region (an ellipse) in the complex plane that
       contains all points within epsilon-distance of U.

    3. GRID SEARCH (The Core "Solve" Step)
       # We look for a "Gaussian Integer" approximation
       FOR each candidate solution z in the search region:
           IF z satisfies the "unit requirement":
               # Check if z can be represented as a sum of squares
               # in the specific ring Z[sqrt(2), i]
               Valid_Points.append(z)

    4. SELECT OPTIMAL POINT
       Pick the point z from Valid_Points that minimizes T-count
       (the point with the smallest "denominator exponent").

    5. EXACT SYNTHESIS (Kliuchnikov's Algorithm)
       # Once we have the algebraic point, we turn it into gates
       Gate_Sequence = []
       WHILE Point != Identity:
           Find a Clifford+T gate G that reduces the "complexity"
           of the Point.
           Point = G * Point
           Gate_Sequence.prepend(Inverse(G))

    RETURN Gate_Sequence
*/
