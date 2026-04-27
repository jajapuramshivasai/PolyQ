import numpy as np
from math import comb


"""

# Z2 Closed-Form Simulation Theory and Complexity Analysis

This simulation framework leverages the **Sum-over-Paths** formalism to represent quantum circuits as Boolean polynomials over the finite field $\mathbb{F}_2$. For the gate set $\{H, Z, CZ\}$, every quantum circuit can be mapped to a characteristic polynomial of at most degree 2, allowing for efficient amplitude calculation via structural analysis.

---

### Core Concepts and Algorithms

#### 1. Sum-over-Paths Foundation
Any quantum circuit using the gate set $\{H, Z, CZ, CCZ\}$ with $n$ qubits and $h$ Hadamard gates can be mapped to a characteristic Boolean polynomial $f(x)$ with $n+h$ variables. The final statevector is obtained by fixing the initial input variables and summing over the internal variables introduced by the $H$ gates:

$$|\psi\rangle = \frac{1}{2^{h/2}} \sum_{x_{internal}} (-1)^{f(x)} |w_0 w_1 \dots w_{n-1}\rangle$$



#### 2. Symplectic and Alternating Forms
A bilinear form $B$ is defined as **alternating** if $B(x, x) = 0$ for all $x$ in the vector space. In the context of $\mathbb{F}_2$, this means the matrix representation of the bilinear form has a zero diagonal. A **Symplectic Matrix** is a symmetric, K-valued matrix that maintains this alternating property. qubit interactions (created by $CZ$ and $H$ gates) form this adjacency structure.

#### 3. Dickson's Reduction Algorithm
**Dickson’s Theorem** establishes that any symplectic matrix of rank $2k$ can be transformed into a canonical block-diagonal form consisting of $k$ blocks of $\begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$. This transformation preserves the quadratic nature of the polynomial while making the Hamming weight calculation computationally feasible. The process uses row and column operations (similar to Gauss-Jordan elimination) to isolate the rank of the bilinear form.



#### 4. Closed-Loop Hamming Weight Calculation
The transition amplitude is calculated by evaluating the **Hamming Weight** ($wt$) of the Boolean polynomial—the count of inputs for which the function results in 1. The exponential sum is related to the weight by:

$$\sum_{v \in \mathbb{F}_2^m} (-1)^{f'(v)} = 2 \times wt(f'(v)) - 2^m$$

The weight of the quadratic part $U$ for a rank $2k$ system (with constant phase $\epsilon=0$) is:
$$wt(U) = \sum_{r=1}^{\lceil k/2 \rceil} \binom{k}{2r-1} 3^{k-(2r-1)}$$

* **Balanced Condition:** If any "kernel variable" (a variable not part of the quadratic rank) persists in the linear part after reduction, the function is balanced, and the amplitude is exactly **zero**.
* **Normalization:** The final amplitude is normalized by the Hadamard count ($h$) using the factor $2^{-h/2}$.



---

### Complexity Reduction Analysis

Traditional statevector simulators are limited by the exponential growth of the Hilbert space, scaling at $O(2^n)$. This polynomial framework provides a significant alternative by focusing on the Hadamard count ($h$) and structural rank.

| Simulation Step | Standard Matrix Method | Z2 Closed-Form Method |
| :--- | :--- | :--- |
| **Parsing** | N/A | $\Theta(g)$ (Linear in gate count) |
| **Evaluation** | $O(2^n)$ (Qubit dependent) | $O(h^3)$ (Hadamard dependent) |
| **Memory** | Exponential ($2^n$ complex numbers) | Polynomial (Matrix of size $h \times h$) |

#### Theoretical Implications
1.  **Elimination of Exponential Sums:** By using Dickson's Theorem, we replace the need to evaluate $2^h$ paths with a rank-finding algorithm that operates in $O(h^3)$ time.
2.  **Qubit Independence:** The complexity of finding a single amplitude is largely independent of the total number of qubits $n$, scaling instead with the number of variables $m$ (where $h-n \le m \le h$) that remain after fixing inputs and outputs.
3.  **Classical Hardness Thresholds:** This method demonstrates that quantum circuits with quadratic polynomials (Clifford-like) are classically efficient to simulate. The introduction of cubic terms ($CCZ$ or $T$ gates) breaks this $O(h^3)$ scaling, reintroducing exponential dependence and highlighting the transition to "quantum advantage."
    

  
"""


class Z2Sim:
    """
    A Quantum Circuit Simulator based on Boolean Polynomials over F2.
    Specifically designed for the {H, Z, CZ} gate set.
    """
    def __init__(self, n):
        """
        Initializes the simulator.
        
        Args:
            n (int): Number of qubits.
        """
        self.n = n
        self.gates = []
        self.compiled = False
        self.n_vars = n
        self.Q = np.zeros((256, 256), dtype=int)
        self.L = np.zeros(256, dtype=int)
        self.h_count = 0

    def h(self, q):
        """Adds a Hadamard gate to the circuit."""
        self.gates.append(('H', q))

    def z(self, q):
        """Adds a Pauli-Z gate to the circuit."""
        self.gates.append(('Z', q))

    def cz(self, q1, q2):
        """Adds a Controlled-Z gate to the circuit."""
        self.gates.append(('CZ', q1, q2))

    def compile(self):
        """
        Maps the sequence of gates to a characteristic quadratic Boolean polynomial.
        H gates introduce new variables, while Z and CZ add linear and quadratic terms.
        """
        wires = [[i] for i in range(self.n)]
        v_next = self.n
        for g in self.gates:
            if g[0] == 'H':
                q = g[1]
                prev, cur = wires[q][-1], v_next
                v_next += 1
                wires[q].append(cur)
                self.h_count += 1
                self.Q[prev, cur] = self.Q[cur, prev] = 1
            elif g[0] == 'CZ':
                v1, v2 = wires[g[1]][-1], wires[g[2]][-1]
                self.Q[v1, v2] = self.Q[v2, v1] = 1
            elif g[0] == 'Z':
                self.L[wires[g[1]][-1]] ^= 1
        self.n_vars, self.output_vars = v_next, [w[-1] for w in wires]
        self.compiled = True

    def print_analytic(self, transition_mode=True):
        """
        Prints the analytical XOR conditions for the amplitude summation based on 
        Dickson's Theorem and the Sum-over-Paths formalism.
        
        The amplitude is derived from the structural weight of the Boolean polynomial:
        Sum = 2^m - 2 * wt(f')
        
        Where:
        - m: Number of internal variables (nu).
        - wt(f'): The Hamming weight of the quadratic form.
        - Rank 2k: The number of variables involved in quadratic pairs.
        - Kernel: Variables not part of the rank; if they have linear terms, 
          the function is "balanced" and the amplitude is zero.
        """
        if not self.compiled:
            print("Compile the circuit first.")
            return

        # 1. Identify boundary vs internal variables
        fixed_indices = set(range(self.n)) | set(self.output_vars)
        u_vars = [i for i in range(self.n_vars) if i not in fixed_indices]
        nu = len(u_vars)
        
        if nu == 0:
            print("\nNo internal variables to sum over (Trivial Circuit).")
            return

        # 2. Initialize Coefficient Matrix: [x_bits | y_bits | constant]
        # Tracks how boundary bits and constants contribute to each internal u_i
        n_cols = 2 * self.n + 1 if transition_mode else self.n + 1
        coeff_matrix = np.zeros((nu, n_cols), dtype=int)
        
        for i, u_idx in enumerate(u_vars):
            # Constant linear contribution
            coeff_matrix[i, -1] = self.L[u_idx]
            
            # Linear contributions from boundary bits via quadratic interactions (Q)
            if transition_mode:
                for xi in range(self.n):
                    if self.Q[u_idx, xi]: coeff_matrix[i, xi] ^= 1
            
            for yi, o_idx in enumerate(self.output_vars):
                col_idx = self.n + yi if transition_mode else yi
                if self.Q[u_idx, o_idx]: coeff_matrix[i, col_idx] ^= 1

        # 3. Mirror the Dickson Diagonalization on the coefficient matrix
        B = self.Q[np.ix_(u_vars, u_vars)].copy()
        rank, p = 0, 0
        while p < nu - 1:
            pivot = next(((i, j) for i in range(p, nu) for j in range(i+1, nu) if B[i, j]), None)
            if not pivot: break
            
            i, j = pivot
            # Mirror Swaps
            B[[p, i]] = B[[i, p]]; B[:, [p, i]] = B[:, [i, p]]
            coeff_matrix[[p, i]] = coeff_matrix[[i, p]]
            
            j_act = i if j == p else j
            B[[p+1, j_act]] = B[[j_act, p+1]]; B[:, [p+1, j_act]] = B[:, [j_act, p+1]]
            coeff_matrix[[p+1, j_act]] = coeff_matrix[[j_act, p+1]]
            
            # Mirror XOR Additions
            for k in range(p + 2, nu):
                if B[k, p]:
                    B[k] ^= B[p+1]; B[:, k] ^= B[:, p+1]
                    coeff_matrix[k] ^= coeff_matrix[p+1]
                if B[k, p+1]:
                    B[k] ^= B[p]; B[:, k] ^= B[:, p]
                    coeff_matrix[k] ^= coeff_matrix[p]
            p += 2; rank += 2

        # 4. Helper to format XOR expressions
        def build_xor_expr(idx):
            terms = []
            if transition_mode:
                terms += [f"x_{xi}" for xi in range(self.n) if coeff_matrix[idx, xi]]
                terms += [f"y_{yi}" for yi in range(self.n) if coeff_matrix[idx, self.n + yi]]
            else:
                terms += [f"y_{yi}" for yi in range(self.n) if coeff_matrix[idx, yi]]
            
            if coeff_matrix[idx, -1]: terms.append("1")
            return " ⊕ ".join(terms) if terms else "0"

        # 5. Final Output
        title = "<y|U|x>" if transition_mode else "Statevector <y|U|0>"
        print(f"\n{'='*60}\n   ANALYTICAL FUNCTIONS FOR {title} (Z2)\n{'='*60}")
        print(f"Hadamard Count (h)         : {self.h_count}")
        print(f"Dickson Rank (2k)          : {rank}")
        
        print("\n1. CANONICAL SUM VARIABLES (Pairs contribute to weight):")
        if rank > 0:
            for i in range(0, rank, 2):
                print(f"   Pair_{i//2}: ( {build_xor_expr(i)} ) * ( {build_xor_expr(i+1)} )")
        else:
            print("   No quadratic pairs found.")

        print("\n2. ZERO AMPLITUDE CONDITIONS (Kernel Variables):")
        kernel_exists = False
        for k in range(rank, nu):
            kernel_exists = True
            print(f"   Condition_{k-rank}: ({build_xor_expr(k)}) == 1  => Amplitude = 0")
        if not kernel_exists:
            print("   No kernel variables exist.")
            
        print("\n3. WEIGHT FORMULA:")
        print("   Sum = 2^(nu-2k) * (2^2k - 2 * wt(Pairs))")
        print("="*60 + "\n")

    def get_amplitude(self, y, x=0):
        """
        Calculates the transition amplitude <y|U|x> using Dickson's Theorem.
        
        Args:
            y (int): Output bitstring as an integer.
            x (int): Input bitstring as an integer.
            
        Returns:
            complex: The calculated amplitude.
        """
        if not self.compiled: self.compile()
        
        fixed = {i: (x >> i) & 1 for i in range(self.n)}
        for i, v_idx in enumerate(self.output_vars):
            fixed[v_idx] = (y >> i) & 1
            
        u_vars = [v for v in range(self.n_vars) if v not in fixed]
        m = len(u_vars)
        
        red_L, epsilon = np.zeros(m, dtype=int), 0
        for i in range(self.n_vars):
            for j in range(i + 1, self.n_vars):
                if self.Q[i, j]:
                    if i in fixed and j in fixed: epsilon ^= (fixed[i] & fixed[j])
                    elif i in fixed: 
                        if fixed[i]: red_L[u_vars.index(j)] ^= 1
                    elif j in fixed: 
                        if fixed[j]: red_L[u_vars.index(i)] ^= 1
        for i in range(self.n_vars):
            if self.L[i]:
                if i in fixed: epsilon ^= fixed[i]
                else: red_L[u_vars.index(i)] ^= 1

        if m == 0: return ((-1)**epsilon) / (2**(self.h_count / 2.0))

        B = self.Q[np.ix_(u_vars, u_vars)]
        rank, p = 0, 0
        while p < m - 1:
            pivot = next(((i, j) for i in range(p, m) for j in range(i+1, m) if B[i, j]), None)
            if not pivot: break
            i, j = pivot
            B[[p, i]] = B[[i, p]]; B[:, [p, i]] = B[:, [i, p]]
            red_L[p], red_L[i] = red_L[i], red_L[p]
            j_act = i if j == p else j
            B[[p+1, j_act]] = B[[j_act, p+1]]; B[:, [p+1, j_act]] = B[:, [j_act, p+1]]
            red_L[p+1], red_L[j_act] = red_L[j_act], red_L[p+1]
            for k in range(p+2, m):
                if B[k, p]: B[k] ^= B[p+1]; B[:, k] ^= B[:, p+1]; red_L[k] ^= red_L[p+1]
                if B[k, p+1]: B[k] ^= B[p]; B[:, k] ^= B[:, p]; red_L[k] ^= red_L[p]
            p += 2; rank += 2

        if any(red_L[rank:]): return 0.0

        for i in range(rank // 2):
            epsilon ^= (red_L[2*i] & red_L[2*i+1])

        k = rank // 2
        wt_u = sum(comb(k, 2*r-1) * (3**(k-(2*r-1))) for r in range(1, k//2+2) if 2*r-1 <= k)
        wt_final = (2**(m-rank)) * ((2**rank - wt_u) if epsilon else wt_u)
        return (2**m - 2*wt_final) / (2**(self.h_count / 2.0))
    
    
    
    
    
    
def test_print_analytic():
    # Test with a Bell State preparation circuit [cite: 85]
    sim = Z2Sim(2)
    sim.h(0)
    
    sim.h(1)
    sim.cz(0, 1)
    sim.h(1)
    sim.z(0)
    
    sim.h(0)
    sim.z(0)
    sim.h(0)
    
    sim.compile()
    for y in range(4): #print statevector amplitudes for all output states given input |00>
        amp = sim.get_amplitude(y, 0)
        print(f"<{y:02b}|U|00> = {amp:.4f}")
    # sim.get_amplitude(0, 0) # <00|U|00>
    sim.print_analytic()
    
def test_z2_bell_state():
    """Verifies Bell state prep: |00> -> (|00> + |11>) / sqrt(2)."""
    sim = Z2Sim(2)
    sim.h(0)
    sim.h(1)
    sim.cz(0, 1)
    sim.h(1)
    sim.compile()
    
    # Input x=0 (|00>), Output y=0 (|00>)
    amp_00 = sim.get_amplitude(0, 0)
    print(f"Amplitude <00|U|00>: {amp_00:.4f}")
    assert np.isclose(abs(amp_00), 1/np.sqrt(2))

    # Input x=0 (|00>), Output y=3 (|11>)
    amp_11 = sim.get_amplitude(3, 0)
    print(f"Amplitude <11|U|00>: {amp_11:.4f}")
    assert np.isclose(abs(amp_11), 1/np.sqrt(2))

def test_x_gate():
    """Verifies X-gate behavior: |0> -> |1>."""
    sim = Z2Sim(1)
    sim.h(0); sim.z(0); sim.h(0) # Equivalent to X
    sim.compile()
    
    amp_0 = sim.get_amplitude(0, 0) # <0|U|0>
    amp_1 = sim.get_amplitude(1, 0) # <1|U|0>
    print(f"X-gate: <0|U|0> = {amp_0}, <1|U|0> = {amp_1}")
    assert np.isclose(amp_0, 0.0)
    assert np.isclose(amp_1, 1.0)

def test_z_gate():
    """Verifies Z-gate behavior: |0> -> |0>, |1> -> -|1>."""
    sim = Z2Sim(1)
    sim.z(0)
    sim.compile()
    
    amp_0 = sim.get_amplitude(0, 0) # <0|U|0>
    amp_1 = sim.get_amplitude(1, 0) # <1|U|0>
    print(f"Z-gate: <0|U|0> = {amp_0}, <1|U|0> = {amp_1}")
    assert np.isclose(amp_0, 1.0)
    assert np.isclose(amp_1, -1.0)