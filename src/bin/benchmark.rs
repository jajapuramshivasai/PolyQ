use std::time::Instant;
use PolyQ::circuit::Circuit;

// A simple benchmark that builds a random circuit via the public `Circuit`
// API and then runs the polynomial simulator on it.

fn run_benchmark(n_qubits: usize, t_count: usize, h_count: usize) {
    let mut circ = Circuit::new(n_qubits);
    // build a random sequence of H and T gates, plus one CZ
    for _ in 0..h_count {
        circ.h(rand::random::<usize>() % n_qubits);
    }
    for _ in 0..t_count {
        circ.t(rand::random::<usize>() % n_qubits);
    }
    circ.cz(0, 1);

    println!("--- Starting Benchmark (n={}, t={}, h={}) ---", n_qubits, t_count, h_count);

    let mut sim = circ.into_simulator();
    let start_sim = Instant::now();
    let amplitude = sim.simulate();
    let sim_time = start_sim.elapsed();

    println!("Simulation Result Amplitude: {:.4}", amplitude);
    println!("Polynomial Simulation Time: {:?}", sim_time);
}

fn main() {
    // Benchmark a high-qubit, T-sparse circuit (where this simulator excels)
    run_benchmark(22, 6, 20);
}