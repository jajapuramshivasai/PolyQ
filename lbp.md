## PolyQ: High-Performance Phase Polynomial Quantum Simulator

**Project Report**

---

### Abstract

PolyQ is a next-generation quantum circuit simulator that leverages phase polynomials over $\mathbb{Z}_4$ and Dickson’s theorem to provide high-performance simulation. By mapping quantum evolution to Boolean polynomials, the framework enables polynomial-time simulation for Clifford circuits and efficient universal simulation via a "Sum-over-Cliffords" branching approach. This methodology allows for the analytical representation of quantum states and transition amplitudes, significantly reducing memory requirements compared to traditional statevector simulators.

---

### 1. Introduction

The simulation of quantum circuits on classical hardware is traditionally limited by the exponential growth of the Hilbert space, requiring massive memory for statevector storage. PolyQ addresses this by representing quantum circuits as low-degree polynomials over finite fields. This approach shifts the computational burden from matrix-vector multiplication to algebraic evaluation, providing exponential advantages for specific circuit classes and offering deep theoretical insights into destructive interference through kernel analysis.

---

### 2. Core Workings: $\mathbb{F}_2$-Based Simulation

The foundation of the simulator rests on mapping a quantum circuit with $n$ qubits and $h$ Hadamard gates to a characteristic polynomial $f(x)$ with $n+h$ variables and a maximum degree of 3.

* **Polynomial Mapping Rules**:
* 
**Z Gate**: Adds a linear term $x_i$ to $f(x)$.


* 
**CZ Gate**: Adds a quadratic term $x_i x_j$.


* 
**CCZ Gate**: Adds a cubic term $x_i x_j x_k$.


* 
**H Gate**: Introduces a new variable $x_{new}$ and adds a quadratic term $x_{old}x_{new}$.




* **Amplitude Relation**: The final statevector is derived by fixing input variables and summing over internal variables:


$$|\psi\rangle = \frac{1}{2^{h/2}} \sum_{x_{internal}} (-1)^{f(x)} |w_0 w_1 \dots w_{l-1}\rangle$$


.



---

### 3. Reduction and Extension to Clifford/$\mathbb{Z}_4$

#### 3.1 Dickson's Theorem and Hamming Weight

To solve the summation in polynomial time for quadratic forms, PolyQ utilizes **Dickson’s Theorem** to reduce the quadratic form matrix $B$ to a canonical block-diagonal form.

* **Rank and Kernel**: The reduction identifies the rank $2k$ of the bilinear form. Variables outside this rank are "kernel variables".


* 
**Efficiency**: If a kernel variable persists in the linear part of the polynomial, the amplitude is zero (balanced function). This structural analysis reduces complexity from $O(2^h)$ to $O(h^3)$.



#### 3.2 $\mathbb{Z}_4$ Universal Extension

PolyQ extends this to universal sets (Clifford+T) by utilizing $\mathbb{Z}_4$-valued quadratic forms.

* 
**Sum-over-Cliffords**: Non-Clifford gates (e.g., $T, R_z$) are handled via Pauli decomposition, branching the simulation into Clifford sub-circuits with associated complex weights.


* **Canonical $\mathbb{Z}_4$ Sums**: Amplitudes are calculated using a closed-form result:

$$A = 2^{m-r} (1+i)^{N_0} (1-i)^{N_1}$$


where $N_0$ and $N_1$ are derived from the transformed linear vector.



---

### 4. Implementation and Optimizations

PolyQ incorporates several high-performance features to ensure scalability:

| Feature | Description |
| --- | --- |
| **Circuit Splitting** | Circuits are split into $E_L$ (left fringe), Core $U$, and $E_R$ (right fringe) to handle affine transforms and phases analytically.

 |
| **DicksonEngine Caching** | The structural reduction of the core quadratic matrix is performed once and reused across all simulation branches.

 |
| **Gray-Code Traversal** | High-speed statevector generation utilizes Gray-code incremental updates, allowing $O(1)$ updates per amplitude.

 |
| **Parallelization** | Operations in the $\mathbb{F}_2$ domain leverage SIMD and bitmap operations to accelerate linear algebra.

 |

---

### 5. Analytical Transition Amplitudes

A key differentiator for PolyQ is the ability to output **closed-form analytical expressions** for amplitudes.

* 
**Functionality**: Instead of a raw numerical value, the simulator provides XOR expressions for exponential sum variables and explicit kernel conditions.


* 
**Utility**: This allows researchers to identify exactly which gate interactions lead to destructive interference without full statevector extraction.



---

### 6. Scaling Analysis and Benchmarks

Benchmarking against industry standards (Qiskit-Aer and DDSIM) reveals:

* 
**Qubit Independence**: Unlike traditional simulators that scale exponentially with the number of qubits ($n$), PolyQ's wall time is primarily dependent on the Hadamard count ($h$).


* 
**Regime of Advantage**: PolyQ achieves significant speedups for high-qubit circuits with low-to-moderate Hadamard depth.


* 
**Memory Efficiency**: By avoiding full statevector storage, PolyQ can simulate $n > 50$ qubit circuits that crash standard simulators on high-memory systems.



---

### 7. Future Work and Applications

* 
**Match Circuits**: Future integration of matchgate simulation using similar algebraic frameworks.


* 
**Algorithm Verification**: High-precision verification of sub-circuits for tensor network contraction.


* 
**Symbolic Manipulation**: Development of symbolic amplitude engines for advanced circuit optimization and error characterization.



Would you like me to refine the technical description of the Gray-code traversal logic used in the statevector generation?