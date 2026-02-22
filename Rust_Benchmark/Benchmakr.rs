
//! # Quantum Circuit Simulation Benchmark
//!
//! This program benchmarks the performance of a simple quantum circuit simulator in Rust.
//! It applies random single-qubit and CNOT gates to a statevector and measures execution time.
//!
//! ## Usage
//!
//! ```sh
//! cargo run --release --bin Benchmakr -- <N> <G> <SEED>
//! ```
//! - `<N>`: Number of qubits (default: 8)
//! - `<G>`: Number of gates (default: 1000)
//! - `<SEED>`: Random seed (default: 42)
//!
//! The output includes timing and statevector norm for correctness.

use num_complex::Complex64;
use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;
use std::env;
use std::time::Instant;

/// Applies a single-qubit gate to the statevector.
///
/// # Arguments
/// * `state` - Mutable reference to the statevector.
/// * `n` - Number of qubits.
/// * `target` - Target qubit index.
/// * `gate` - 2x2 unitary matrix representing the gate.
fn apply_single_qubit_gate(state: &mut [Complex64], n: usize, target: usize, gate: [[Complex64; 2]; 2]) {
	let size = state.len();
	for i in 0..size {
		if ((i >> target) & 1) == 0 {
			let j = i | (1 << target);
			let a = state[i];
			let b = state[j];
			state[i] = gate[0][0] * a + gate[0][1] * b;
			state[j] = gate[1][0] * a + gate[1][1] * b;
		}
	}
}

/// Applies a CNOT gate to the statevector.
///
/// # Arguments
/// * `state` - Mutable reference to the statevector.
/// * `n` - Number of qubits.
/// * `control` - Control qubit index.
/// * `target` - Target qubit index.
fn apply_cnot(state: &mut [Complex64], n: usize, control: usize, target: usize) {
	let size = state.len();
	if control == target { return; }
	if control < target {
		for i in 0..size {
			if ((i >> control) & 1) == 1 && ((i >> target) & 1) == 0 {
				let j = i | (1 << target);
				state.swap(i, j);
			}
		}
	} else {
		for i in 0..size {
			if ((i >> control) & 1) == 1 && ((i >> target) & 1) == 0 {
				let j = i | (1 << target);
				state.swap(i, j);
			}
		}
	}
}

/// Entry point for the benchmark program.
///
/// Parses command-line arguments, initializes the statevector, applies random gates,
/// measures execution time, and prints results.
fn main() {
	let args: Vec<String> = env::args().collect();
	let n: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(8);
	let g: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(1000);
	let seed: u64 = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(42);

	let dim = 1usize << n;
	let mut state = vec![Complex64::new(0.0, 0.0); dim];
	state[0] = Complex64::new(1.0, 0.0);

	let mut rng = StdRng::seed_from_u64(seed);

	let sqrt2_inv = 1.0 / 2.0f64.sqrt();
	let h = [[Complex64::new(sqrt2_inv, 0.0), Complex64::new(sqrt2_inv, 0.0)],
			 [Complex64::new(sqrt2_inv, 0.0), Complex64::new(-sqrt2_inv, 0.0)]];
	let s = [[Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)],
			 [Complex64::new(0.0, 0.0), Complex64::new(0.0, 1.0)]];
	let x = [[Complex64::new(0.0, 0.0), Complex64::new(1.0, 0.0)],
			 [Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)]];
	let y = [[Complex64::new(0.0, 0.0), Complex64::new(0.0, -1.0)],
			 [Complex64::new(0.0, 1.0), Complex64::new(0.0, 0.0)]];
	let z = [[Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)],
			 [Complex64::new(0.0, 0.0), Complex64::new(-1.0, 0.0)]];

	let start = Instant::now();
	for _ in 0..g {
		let gate_choice = rng.gen_range(0..6);
		match gate_choice {
			0 => { // H
				let target = rng.gen_range(0..n);
				apply_single_qubit_gate(&mut state, n, target, h);
			}
			1 => { // S
				let target = rng.gen_range(0..n);
				apply_single_qubit_gate(&mut state, n, target, s);
			}
			2 => { // X
				let target = rng.gen_range(0..n);
				apply_single_qubit_gate(&mut state, n, target, x);
			}
			3 => { // Y
				let target = rng.gen_range(0..n);
				apply_single_qubit_gate(&mut state, n, target, y);
			}
			4 => { // Z
				let target = rng.gen_range(0..n);
				apply_single_qubit_gate(&mut state, n, target, z);
			}
			_ => { // CNOT
				let control = rng.gen_range(0..n);
				let mut target = rng.gen_range(0..n);
				while target == control { target = rng.gen_range(0..n); }
				apply_cnot(&mut state, n, control, target);
			}
		}
	}
	let elapsed = start.elapsed();
	let secs = elapsed.as_secs_f64();
	println!("N={} qubits, G={} gates, time: {:.6} s, gates/s: {:.3}", n, g, secs, (g as f64) / secs);

	// compute simple fidelity (norm) to ensure correctness
	let mut norm = 0.0f64;
	for amp in &state {
		norm += amp.norm_sqr();
	}
	println!("statevector norm = {:.6}", norm);
}

