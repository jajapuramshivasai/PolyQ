//Clifford Gate accepted by poly q {H,Z,CZ}
//Clifford+S Gate accepted by poly q {H,Z,CZ,S}
//Clifford+T Gate accepted by poly q {H,Z,CZ,T}

use crate::quantum_circuit::{QuantumCircuit, Gate, QuantumCircuitClass};

/// Transpile a circuit to the accepted gateset for the given class.
pub fn transpile_to_gateset(circuit: &QuantumCircuit, class: QuantumCircuitClass) -> QuantumCircuit {
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
			new_gates.push(Gate { name: "h".to_string(), qubits: vec![q], params: vec![] });
			new_gates.push(Gate { name: "z".to_string(), qubits: vec![q], params: vec![] });
			new_gates.push(Gate { name: "h".to_string(), qubits: vec![q], params: vec![] });
		} else if g == "cx" && gate.qubits.len() == 2 {
			// CX(control, target) = H(target) CZ(control, target) H(target)
			let control = gate.qubits[0];
			let target = gate.qubits[1];
			new_gates.push(Gate { name: "h".to_string(), qubits: vec![target], params: vec![] });
			new_gates.push(Gate { name: "cz".to_string(), qubits: vec![control, target], params: vec![] });
			new_gates.push(Gate { name: "h".to_string(), qubits: vec![target], params: vec![] });
		} else if g == "y" && gate.qubits.len() == 1 {
			// Y(q) = S(q) H(q) X(q) H(q) S†(q)
			let q = gate.qubits[0];
			new_gates.push(Gate { name: "s".to_string(), qubits: vec![q], params: vec![] });
			new_gates.push(Gate { name: "h".to_string(), qubits: vec![q], params: vec![] });
			// X(q) = H(q) Z(q) H(q)
			new_gates.push(Gate { name: "h".to_string(), qubits: vec![q], params: vec![] });
			new_gates.push(Gate { name: "z".to_string(), qubits: vec![q], params: vec![] });
			new_gates.push(Gate { name: "h".to_string(), qubits: vec![q], params: vec![] });
			new_gates.push(Gate { name: "h".to_string(), qubits: vec![q], params: vec![] });
			new_gates.push(Gate { name: "sdag".to_string(), qubits: vec![q], params: vec![] });
		} else {
		}
	}
	QuantumCircuit {
		num_qubits: circuit.num_qubits,
		gates: new_gates,
	}
}