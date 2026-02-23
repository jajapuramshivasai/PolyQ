# Rust Benchmark Instructions

This directory contains a Rust benchmark for quantum circuit simulation.

## Prerequisites
- [Rust toolchain](https://www.rust-lang.org/tools/install) (install via `rustup`)

## Building the Rust Library

From the root of the repository or inside the `Rust_Benchmark` directory, run:

```
cargo build --release
```

This will build all Rust code in release mode for best performance. The compiled binaries will be located in `../target/release/`.

## Running the Benchmark

To run the benchmark executable:

```
cargo run --release --bin Benchmark -- <N> <G> <SEED>
```
- `<N>`: Number of qubits (default: 8)
- `<G>`: Number of gates (default: 1000)
- `<SEED>`: Random seed (default: 42)

Example:

```
cargo run --release --bin benchmarker -- 10 2000 123
```

Alternatively, you can run the binary directly after building:

```
../target/release/Benchmakr 10 2000 123
```

## Output
The benchmark prints the number of qubits, number of gates, elapsed time, gates per second, and the final statevector norm for correctness checking.

## Building and Viewing Rust Documentation

You can generate and view documentation for the Rust code using Cargo's built-in doc tool.

To build the documentation, run:

```
cargo doc --no-deps --open
```

This will build the documentation for the project and open it in your default web browser. You can also find the generated docs in the `target/doc/` directory.

For more details on documenting Rust code, see: https://doc.rust-lang.org/rustdoc/how-to-write-documentation.html



## PolyQ Architecture

The following flowchart illustrates the architecture of the PolyQ backend and its benchmark integration:

```mermaid
flowchart TD
A[PolyQ Rust Library]
A --> B[polynomial.rs]
A --> C[ring_math.rs]
A --> D[dataset.rs]
A --> E[engine.rs]
A --> F[simulation.rs]
A --> G[poly_opt.rs]
B --> H[Polynomial struct]
D --> I[Dataset struct]
E --> J[Engine struct]
F --> K[Simulation routines]
G --> L[Optimization helpers]
C --> M[Modular arithmetic]

%% Descriptions as tooltips
click B "#" "Handles polynomial arithmetic and representations"
click C "#" "Implements modular arithmetic for ring operations"
click D "#" "Manages datasets and data loading"
click E "#" "Core engine logic and orchestration"
click F "#" "Simulation routines for quantum/classical systems"
click G "#" "Helpers for polynomial optimization"
click H "#" "Main struct for polynomial objects"
click I "#" "Main struct for dataset objects"
click J "#" "Main struct for engine logic"
click K "#" "Functions for running simulations"
click L "#" "Utility functions for optimization"
click M "#" "Implements modular arithmetic functions"

%% Color styling with dark text
classDef main fill:#f9f,stroke:#333,stroke-width:2px,color:#111;
classDef file fill:#bbf,stroke:#222,stroke-width:1.5px,color:#111;
classDef struct fill:#bfb,stroke:#222,stroke-width:1.5px,color:#111;
class A main;
class B file;
class C file;
class D file;
class E file;
class F file;
class G file;
class H struct;
class I struct;
class J struct;
class K struct;
class L struct;
class M struct;
```