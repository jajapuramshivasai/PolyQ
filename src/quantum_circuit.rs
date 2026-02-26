use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum QuantumCircuitClass {
    Clifford,
    CliffordPlusS,
    CliffordPlusT,
}

#[derive(Debug, Clone)]
pub struct QuantumCircuit {
    pub num_qubits: usize,
    pub gates: Vec<Gate>,
}

#[derive(Debug, Clone)]
pub struct Gate {
    pub name: String,
    pub qubits: Vec<usize>,
    pub params: Vec<f64>,
}

impl QuantumCircuit {
    /// Initialize a QuantumCircuit from a QASM file
    pub fn from_qasm_file<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let mut num_qubits = 0;
        let mut gates = Vec::new();

        for line in reader.lines() {
            let line = line?;
            let line = line.trim();
            if line.is_empty() || line.starts_with("//") {
                continue;
            }
            if line.starts_with("qreg") {
                // Example: qreg q[5];
                if let Some(start) = line.find('[') {
                    if let Some(end) = line.find(']') {
                        let n = &line[start + 1..end];
                        if let Ok(n) = n.parse::<usize>() {
                            num_qubits = n;
                        }
                    }
                }
            } else if line.ends_with(';') {
                // Parse gate line
                let line = &line[..line.len() - 1];
                let mut parts = line.split_whitespace();
                if let Some(name) = parts.next() {
                    let rest: Vec<&str> = parts.collect();
                    let mut qubits = Vec::new();
                    let mut params = Vec::new();
                    // Handle parameterized gates: e.g., rx(pi/2) q[0];
                    let (gate_name, param_str) = if let Some(idx) = name.find('(') {
                        let end_idx = name.find(')').unwrap_or(name.len());
                        (&name[..idx], Some(&name[idx + 1..end_idx]))
                    } else {
                        (name, None)
                    };
                    if let Some(param_str) = param_str {
                        // For simplicity, only handle numeric params
                        for p in param_str.split(',') {
                            if let Ok(val) = p.trim().parse::<f64>() {
                                params.push(val);
                            }
                        }
                    }
                    for arg in rest {
                        if let Some(start) = arg.find('[') {
                            if let Some(end) = arg.find(']') {
                                let idx = &arg[start + 1..end];
                                if let Ok(q) = idx.parse::<usize>() {
                                    qubits.push(q);
                                }
                            }
                        }
                    }
                    gates.push(Gate {
                        name: gate_name.to_string(),
                        qubits,
                        params,
                    });
                }
            }
        }
        Ok(QuantumCircuit { num_qubits, gates })
    }
    /// Classify the circuit as Clifford, Clifford+S, or Clifford+T
    pub fn classify(&self) -> QuantumCircuitClass {
        let mut has_s = false;
        let mut has_t = false;
        for gate in &self.gates {
            let g = gate.name.to_lowercase();
            match g.as_str() {
                // Clifford gates: h, cx, cz, x, y, z, swap, s, sdag, etc.
                "h" | "cx" | "cz" | "x" | "y" | "z" | "swap" | "sdag" => {}
                "s" => {
                    has_s = true;
                }
                "t" | "tdag" => {
                    has_t = true;
                }
                _ => {} // For now, ignore other gates
            }
        }
        if has_t {
            QuantumCircuitClass::CliffordPlusT
        } else if has_s {
            QuantumCircuitClass::CliffordPlusS
        } else {
            QuantumCircuitClass::Clifford
        }
    }
}
