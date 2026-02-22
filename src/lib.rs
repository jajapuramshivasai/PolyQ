pub mod polynomial;
pub mod ring_math;
pub mod dataset;
pub mod engine;
pub mod simulation;
pub mod poly_opt;

pub use polynomial::Polynomial;
pub use dataset::Dataset;
pub use engine::Engine;

// PolyQ modular Rust library root.
//
// This crate exposes the following modules:
// - `polynomial` — polynomial representation and helpers
// - `ring_math` — small modular arithmetic helpers
// - `dataset` — dataset descriptors binding a name to a polynomial
// - `engine` — simulation `Engine` and configuration
// - `simulation` — simple simulation routines
// - `template` — optimization helpers for polynomials
