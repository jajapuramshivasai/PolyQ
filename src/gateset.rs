// GateSet trait and implementations for Clifford and Clifford+T
use crate::polynomial::Polynomial;

pub trait GateSet {
    fn modulus(&self) -> i64;
    fn gate_to_polynomial(&self, gate: &str, qubits: &[usize]) -> Option<Polynomial>;
}

pub struct CliffordGateSet;

impl GateSet for CliffordGateSet {
    fn modulus(&self) -> i64 {
        4 // Z4
    }
    fn gate_to_polynomial(&self, gate: &str, qubits: &[usize]) -> Option<Polynomial> {
        match gate {
            // Example: S gate (phase gate) on qubit q: S = diag(1, i) => polynomial x_q^2
            "S" => Some(Polynomial::new(vec![0, 0, 1])), // x^2
            // H and CZ are Clifford but do not have a nontrivial phase polynomial
            "H" | "CZ" => Some(Polynomial::new(vec![0])),
            _ => None,
        }
    }
}

pub struct CliffordTGateSet;

impl GateSet for CliffordTGateSet {
    fn modulus(&self) -> i64 {
        8 // Z8
    }
    fn gate_to_polynomial(&self, gate: &str, qubits: &[usize]) -> Option<Polynomial> {
        match gate {
            // T gate: phase pi/4, polynomial x_q^3
            "T" => Some(Polynomial::new(vec![0, 0, 0, 1])), // x^3
            // S gate: phase pi/2, polynomial x_q^2
            "S" => Some(Polynomial::new(vec![0, 0, 1])), // x^2
            // H and CZ are Clifford but do not have a nontrivial phase polynomial
            "H" | "CZ" => Some(Polynomial::new(vec![0])),
            _ => None,
        }
    }
}
