The simulation method implemented in the engine is based on the algebraic representation of quantum states as **phase polynomials** over Boolean variables. This framework maps the exponential complexity of quantum statevectors into the evaluation of exponential sums over $\mathbb{Z}_4$ (integers modulo 4), allowing Clifford circuits to be simulated in polynomial time. Universal simulation is then achieved by decomposing non-Clifford gates into weighted sums of Clifford operations.

---

### 1. The Phase Polynomial Representation

In this formalism, the action of a quantum circuit on a standard basis state $|x\rangle$ is tracked algebraically rather than via matrix multiplication. For a Clifford circuit containing $n$ qubits and $h$ Hadamard gates, the resulting state is represented up to a global phase by:

$$|\psi\rangle = \frac{1}{\sqrt{2^h}} \sum_{v \in \{0,1\}^h} i^{Q(v, x)} |y(v, x)\rangle$$

Where:
* $v \in \{0,1\}^h$ represents the internal summation variables (paths) introduced by Hadamard gates.
* $Q(v, x)$ is a quadratic polynomial evaluated over $\mathbb{Z}_4$ (where $+1$ corresponds to a phase of $i$, $+2$ to $-1$, $+3$ to $-i$, and $+0$ to $1$).
* $|y(v, x)\rangle$ is a deterministic affine function mapping the inputs and internal variables to the output basis state.



The phase polynomial $Q(u)$ takes the general form:
$$Q(u) = 2 \sum_{j < k} B_{jk} u_j u_k + \sum_j L_j u_j + \epsilon \pmod 4$$
Here, $B_{jk} \in \mathbb{F}_2$ acts as an adjacency matrix of quadratic interactions, $L_j \in \mathbb{Z}_4$ is a vector of linear phase shifts, and $\epsilon \in \mathbb{Z}_4$ is the global phase. 

### 2. Action of Clifford Gates (The Transpilation Rules)

The circuit is compiled by sweeping through the gates and updating the variables $B$, $L$, and the connectivity graph. The rules are geometrically straightforward:

* **Phase Gate ($S_q$):** Adds a $\pi/2$ phase. $L_q \to L_q + 1 \pmod 4$.
* **Pauli-Z Gate ($Z_q$):** Adds a $\pi$ phase. $L_q \to L_q + 2 \pmod 4$.
* **Controlled-Z Gate ($CZ_{a, b}$):** Entangles variables $a$ and $b$, creating a quadratic edge. $B_{ab} \to B_{ab} \oplus 1$.
* **Hadamard Gate ($H_q$):** This is the only gate that expands the Hilbert space dimensions in the polynomial. It replaces the variable on wire $q$ with a new internal summation variable $u_{new}$, and adds a quadratic interaction between the old variable and the new one: $Q \to Q + 2x_{old}u_{new} \pmod 4$.

### 3. Amplitude Evaluation as an Exponential Sum

To calculate the specific transition amplitude from an input state $|0\dots0\rangle$ to a target output bitstring $|y\rangle$, we fix the boundary variables. This collapses the polynomial to depend *only* on the $m$ internal variables $u$.

The amplitude is given by the exponential sum:
$$\langle y | C | 0 \rangle = \frac{1}{\sqrt{2^h}} i^{\epsilon(y)} \sum_{u \in \{0,1\}^m} i^{u^T B u + L(y)^T u}$$

Evaluating this sum naively requires $O(2^m)$ operations. However, because $B$ is a symmetric matrix over $\mathbb{F}_2$, we can evaluate it in polynomial time using Dickson's Theorem.

### 4. Dickson's Theorem (Symplectic Reduction)

**Theorem:** Any symmetric quadratic form over $\mathbb{F}_2$ can be reduced to a canonical, block-diagonal form via an invertible linear transformation $T$. 

We apply operations (row/column swaps and additions) to $B$ to diagonalize it into $k$ disjoint pairs. The rank of $B$ is $r = 2k$. 

$$B_{canonical} = T^T B T = \bigoplus_{i=1}^{k} \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix} \oplus \mathbf{0}$$

Under this change of basis ($u = T w$), the quadratic polynomial becomes purely a sum of independent pairs and isolated kernel variables:
$$Q_{canonical}(w) = 2\sum_{i=1}^{k} w_{2i-1} w_{2i} + \sum_{j} L'_j w_j \pmod 4$$
where the new linear vector is $L' = L T$.



### 5. The Closed-Form Analytical Solution

Because the canonical polynomial is completely uncoupled, the exponential sum over $2^m$ states factors into a product of independent, easily solvable sub-sums:

$$\sum_{w} i^{Q_{can}(w)} = \left( \prod_{i=1}^{k} \sum_{w_{2i-1}, w_{2i} \in \{0,1\}} i^{2 w_{2i-1} w_{2i} + L'_{2i-1}w_{2i-1} + L'_{2i}w_{2i}} \right) \times \left( \prod_{j=2k}^{m} \sum_{w_j \in \{0,1\}} i^{L'_j w_j} \right)$$

**Proof of the Kernel Condition (Balanced Function):**
For any variable $j \ge 2k$ (a kernel variable with no quadratic terms), its sum is $\sum_{w_j} i^{L'_j w_j} = 1 + i^{L'_j}$. 
If the linear coefficient $L'_j = 2 \pmod 4$, the sum evaluates to $1 + i^2 = 1 - 1 = 0$. Therefore, if the output bitstring $y$ shifts any kernel variable's linear term to $2$, the entire amplitude collapses to exactly $0$ due to destructive interference.

**Proof of the Gaussian Pair Sum:**
For the paired rank variables ($j < 2k$), there are 4 combinations of $(w_{2i-1}, w_{2i})$. Evaluating the inner sum yields either a phase rotated by $(1+i)$ or $(1-i)$, scaled by a factor of 2. 

Let $N_1$ be the number of canonical pairs evaluating to $1-i$, and $N_0$ be the pairs evaluating to $1+i$. The exact closed-form solution for the amplitude becomes:

$$\text{Amplitude}(y) = \frac{1}{\sqrt{2^h}} i^{\epsilon(y)} 2^{m - k} (1+i)^{N_0(y)} (1-i)^{N_1(y)}$$

This reduces the complexity of finding any amplitude to $O(m^3)$ for the initial Dickson reduction, and effectively $O(1)$ (using bitwise F2 shifts) for each subsequent target bitstring $y$.

### 6. Universal Simulation (Sum-over-Cliffords)

Clifford gates alone cannot achieve universal quantum computation; they require the addition of non-Clifford operations like the $T$ gate or arbitrary $R_z(\theta)$ rotations. These gates do not map cleanly to $\mathbb{Z}_4$ arithmetic.

To bypass this, the simulator uses the **Pauli Decomposition** (Stabilizer Rank approach). Any single-qubit unitary can be expressed as a linear combination of Pauli matrices. For a Z-rotation:
$$R_z(\theta) = \cos(\theta/2) I - i \sin(\theta/2) Z$$

**The Branching Execution:**
When the simulator encounters an $R_z(\theta)$ gate, it does not calculate a new polynomial. Instead, it forks the simulation into two parallel Clifford universes:
1. **The $I$ Branch:** The gate is treated as the Identity. The resulting amplitude is scaled by $\cos(\theta/2)$.
2. **The $Z$ Branch:** The gate is treated as a Pauli-Z (adding $+2 \pmod 4$ to the linear vector). The resulting amplitude is scaled by $-i \sin(\theta/2)$.

note in all branches only the linear term dependence on x,y is changes the quadratic part stays the same so we have to dickson reduction only once and use it to evaluate all branches efficiently.




For a circuit with $t$ non-Clifford rotations, the total amplitude is the weighted sum of $2^t$ independent Clifford circuits. Because the $\\sin(\\theta/2)$ weight drops exponentially for small rotation angles, branches falling below a certain magnitude threshold are pruned. This prevents an exponential memory explosion and allows for fast, highly parallelized universal simulation. Alternatively, if the weights are mostly uniform, we can evaluate only the first $n$ branches we can quantify the error induces with operator norm distance.