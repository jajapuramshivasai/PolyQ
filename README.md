# PolyQ: High-Performance Phase Polynomial Quantum Simulator

**PolyQ** is a next-generation quantum circuit simulator based on **phase polynomials over ℤ₄** and **Dickson’s theorem**. It achieves polynomial-time simulation for Clifford circuits and efficient universal simulation via weighted branching over non-Clifford gates.

---

## Abstract

PolyQ enables **polynomial-time evaluation of transition amplitudes** and provides **closed-form analytical representations** of quantum states and transition matrices. By leveraging phase polynomial representations and symplectic reduction (Dickson’s theorem), it offers exponential advantages over naive Feynman-path or full statevector simulators for many circuit classes. The framework supports both exact simulation and efficient observable computation with dramatically reduced memory requirements.

---

## Core Theory

### 1. Phase Polynomial Representation

Any Clifford circuit acting on an input basis state $|x\rangle$ can be represented algebraically as:

$$
|\psi\rangle = \frac{1}{\sqrt{2^h}} \sum_{v \in \{0,1\}^h} i^{Q(v, x)} |y(v, x)\rangle
$$

Where:
- $h$ = number of Hadamard gates (determines the number of summation variables)
- $Q(v, x)$ is a **quadratic phase polynomial over ℤ₄**
- $y(v, x)$ is an **affine transformation** of the input and internal variables

The phase polynomial takes the general form:

$$
Q(u) = 2 \sum_{j < k} B_{jk} u_j u_k + \sum_j L_j u_j + \epsilon \pmod{4}
$$

- $B \in \mathbb{F}_2^{m \times m}$: symmetric adjacency matrix of quadratic interactions
- $L \in \mathbb{Z}_4^m$: linear phase coefficients
- $\epsilon \in \mathbb{Z}_4$: global phase

### 2. Action of Clifford Gates

PolyQ maintains the polynomial structure by applying simple update rules:

| Gate     | Effect on Polynomial |
|----------|----------------------|
| **S**    | $L_q \gets L_q + 1 \pmod{4}$ |
| **Z**    | $L_q \gets L_q + 2 \pmod{4}$ |
| **CZ**   | $B_{ab} \gets B_{ab} \oplus 1$ |
| **H**    | Introduces new variable $u_{new}$ and quadratic edge $2 \cdot x_{old} \cdot u_{new}$ |

### 3. Amplitude Evaluation via Dickson’s Theorem

For a fixed output $y$, the amplitude is:

$$
\langle y | C | x \rangle = \frac{i^{\epsilon(y)}}{\sqrt{2^h}} \sum_{u} i^{u^T B u + L(y)^T u}
$$

**Dickson’s Theorem** allows reduction of the quadratic form $B$ to canonical block-diagonal form:

$$
B_{\text{canonical}} = \bigoplus_{i=1}^{k} \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix} \oplus \mathbf{0}
$$

This reduces the exponential sum to a **product of independent closed-form sub-sums**:

- **Paired variables** (rank $2k$): Each pair contributes a factor of $(1+i)$ or $(1-i)$ (scaled).
- **Kernel variables** (nullity $m-2k$): Contribute $1+i^{L_j}$. If any $L_j \equiv 2 \pmod{4}$, amplitude = **0** (destructive interference).

Final closed form:

$$
A(y) = \frac{2^{m-k}}{\sqrt{2^h}} \, i^{\epsilon(y)} \, (1+i)^{N_0(y)} \, (1-i)^{N_1(y)}
$$

---

## Universal Simulation (Sum-over-Cliffords)

Non-Clifford gates ($T$, $R_z(\theta)$, etc.) are handled via **Pauli decomposition**:

$$
R_z(\theta) = \cos(\theta/2) \, I - i \sin(\theta/2) \, Z
$$

Each non-Clifford gate **branches** the simulation into two Clifford sub-circuits with associated complex weights. The quadratic structure $B$ remains identical across branches — only the linear vector $L$ changes. This enables **single Dickson reduction** reused across all branches.

---

## Architecture & Workings

### 1. Circuit Splitting (E_L, Core, E_R)

PolyQ splits any circuit into three parts:

- **E_L (Left Fringe)**: Pre-core gates (affine transform + phases on input)
- **Core U**: Clifford + non-Clifford gates (phase polynomial engine)
- **E_R (Right Fringe)**: Post-core gates (permutations and diagonal phases on output)

This separation allows analytical handling of input/output transformations without inflating the core polynomial.

### 2. DicksonEngine

- Builds adjacency matrix $B_4$ from circuit gates
- Performs optimized **Dickson reduction** (with SWAP/ADD operation tracking)
- Supports Gray-code traversal for full statevector generation
- Provides fast single-amplitude and transition-matrix evaluation

### 3. UniversalQC

- Manages branching tree for non-Clifford gates
- Pre-caches v4 linear phase states per branch
- Applies E_L / E_R transformations analytically
- Combines branch contributions with proper weighting and global phase

---

## Key Features

### Analytical Capabilities
- **Closed-form amplitude formulas** with explicit XOR expressions for exponential sum variables and kernel conditions
- **Full transition matrix** in analytical form
- **Branch decomposition visualization** showing weights and formulas

### Performance
- **Clifford circuits**: Near-constant time per amplitude after initial reduction
- **Universal circuits**: Exponential in number of T-gates, but highly pruneable via weight thresholding
- Dramatic speedup vs. statevector simulators for medium-to-large qubit counts with moderate non-Clifford depth

### Memory Efficiency
- No need to store full $2^n$ statevector when only amplitudes or observables are required
- Gray-code incremental updates for full statevector when needed

---

## Applications

### 1. Quantum Algorithm Verification
- Exact or high-precision verification of small-to-medium quantum algorithms
- Analytical insight into **why** certain amplitudes are zero (kernel conditions)
- Debugging of phase kickback and interference patterns

### 2. Quantum Simulation Research
- Fast exploration of circuit families (e.g. stabilizer circuits , Match Circuits)
- Study of **stabilizer rank** and non-Clifford resource requirements
- Benchmarking of quantum advantage thresholds

### 3. Observable Computation
- Efficient computation of expectation values without full statevector
- Marginal probabilities and correlation functions via analytical formulas

### 4. Hybrid Classical-Quantum Workflows
- Pre-computation of sub-circuit amplitudes for tensor network contraction
- Error mitigation and characterization using exact branch weights

---

## Implementation Highlights

- **Dickson reduction** with operation replay for linear vector transformation
- **Gray code** statevector generation for efficient full simulation
- **Branch-and-bound style pruning** potential via small weight filtering
- **Qiskit integration** for easy circuit input
- **NumPy-accelerated** linear algebra over ℤ₄ and ℤ₂

---

## Future Directions

- Parallel branch evaluation (GPU/CPU)
- Adaptive pruning and importance sampling of branches
- Support for more gate families (MCRZ, etc.)
- Integration with tensor networks for deeper circuits
- Symbolic amplitude manipulation for circuit optimization

---

**PolyQ** represents a powerful bridge between algebraic quantum information theory and practical high-performance simulation, enabling both deep theoretical insight and scalable numerical computation.

---

*Developed as part of advanced quantum simulation research (LBP-2026 framework).*