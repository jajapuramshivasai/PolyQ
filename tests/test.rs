use std::f64::consts::PI;
use PolyQ::circuit::Circuit;

#[test]
fn test_ghz_with_non_trivial_phase() {
    let n = 3;
    let mut circ = Circuit::new(n);
    
    // 1. Prepare the GHZ state: |000> -> 1/sqrt(2)(|000> + |111>)
    circ.h(0)
        .cx(0, 1)
        .cx(1, 2);
    
    // 2. Apply a non-trivial phase to the |111> component.
    // A T-gate on q2 applies a phase of exp(i * pi / 4) to the |1> state.
    // Since q2 is only 1 in the |111> branch, this isolates the phase.
    circ.t(2);
    
    // 3. Uncompute the GHZ state.
    // We reverse the preparation operations in reverse order.
    circ.cx(1, 2)
        .cx(0, 1)
        .h(0);
        
    // Mathematical Expectation:
    // The state before uncomputing is 1/sqrt(2)(|000> + exp(i * pi / 4)|111>)
    // After uncomputing, the amplitude of |000> becomes:
    // 1/2 * (1 + exp(i * pi / 4))
    
    let mut sim = circ.into_simulator();
    let amp = sim.simulate();
    
    let expected_re = 0.5 * (1.0 + (PI / 4.0).cos());
    let expected_im = 0.5 * (PI / 4.0).sin();
    
    // Check if the simulator matches the expected amplitude within a tight tolerance
    assert!((amp.re - expected_re).abs() < 1e-6, "Real part mismatch: {} != {}", amp.re, expected_re);
    assert!((amp.im - expected_im).abs() < 1e-6, "Imaginary part mismatch: {} != {}", amp.im, expected_im);
}

#[test]
fn test_bernstein_vazirani() {
    // We are trying to find the hidden string s = 1011
    let secret_string = vec![1, 0, 1, 1]; 
    let n = secret_string.len();
    let total_qubits = n + 1; // +1 for the |-> ancilla
    
    let mut circ = Circuit::new(total_qubits);
    let ancilla = n;
    
    // 1. Initialize the ancilla to |1> (which transpiles to H-Z-H)
    circ.x(ancilla);
    
    // 2. Apply Hadamard to all qubits to create superposition
    for i in 0..total_qubits {
        circ.h(i);
    }
    
    // 3. Apply the BV Oracle: f(x) = s \cdot x
    for (i, &bit) in secret_string.iter().enumerate() {
        if bit == 1 {
            circ.cx(i, ancilla);
        }
    }
    
    // 4. Apply Hadamard to input qubits
    // At this point, qubits 0..n are guaranteed to be in the state |s>
    for i in 0..n {
        circ.h(i);
    }
    
    // 5. Uncompute the result to check the simulator's all-zero amplitude
    // We flip the bits where s = 1 back to 0.
    for (i, &bit) in secret_string.iter().enumerate() {
        if bit == 1 {
            circ.x(i);
        }
    }
    
    // Uncompute the ancilla: |-> needs H then X to return to |0>
    circ.h(ancilla)
        .x(ancilla);
    
    // Because we completely uncomputed everything back to the ground state,
    // the simulator's computed amplitude for <00..0 | U | 00..0> MUST be 1.0
    
    let mut sim = circ.into_simulator();
    let amp = sim.simulate();
    
    assert!((amp.re.abs() - 1.0).abs() < 1e-6, "BV Amplitude magnitude should be 1.0, got {}", amp.re);
    assert!(amp.im.abs() < 1e-6, "BV Imaginary part should be 0.0, got {}", amp.im);
}