import numpy as np
# from .QInfo import fidelity



class Z4Sim:
    def __init__(self, n):
        self.n = n
        self.gates = []
        self.compiled = False
        self.v4 = np.zeros(256, dtype=int)  
        self.b4 = np.zeros((256, 256), dtype=int)
        self.h_count = 0

    def h(self, q): self.gates.append(('H', q))
    def s(self, q): self.gates.append(('S', q))
    def z(self, q): self.gates.append(('Z', q))
    def cz(self, q1, q2): self.gates.append(('CZ', q1, q2))

    def compile(self):
        wires = [[i] for i in range(self.n)]
        v_next = self.n
        for g in self.gates:
            if g[0] == 'H':
                q = g[1]; prev, cur = wires[q][-1], v_next
                v_next += 1; wires[q].append(cur); self.h_count += 1
                self.b4[prev, cur] = self.b4[cur, prev] = 1 # Quadratic weight 4 
            elif g[0] == 'CZ':
                v1, v2 = wires[g[1]][-1], wires[g[2]][-1]
                self.b4[v1, v2] = self.b4[v2, v1] = 1
            elif g[0] == 'S': self.v4[wires[g[1]][-1]] = (self.v4[wires[g[1]][-1]] + 1) % 4
            elif g[0] == 'Z': self.v4[wires[g[1]][-1]] = (self.v4[wires[g[1]][-1]] + 2) % 4
        
        self.n_vars, self.output_vars = v_next, [w[-1] for w in wires]
        self.compiled = True

    def get_amplitude(self, y, x=0):
        """Calculates <y|U|x> using Z4 Schmidt Summation."""
        fixed = {}
        for i in range(self.n): fixed[i] = (x >> i) & 1
        for i, v in enumerate(self.output_vars): fixed[v] = (y >> i) & 1
            
        u_vars = [i for i in range(self.n_vars) if i not in fixed]
        nu = len(u_vars)
        
        # eps_base tracks constant phase; vu_base tracks linear terms for sum variables
        eps, vu = 0, np.zeros(nu, dtype=int)
        f_list = [v for v, b in fixed.items() if b == 1]
        
        for f in f_list:
            eps = (eps + self.v4[f]) % 4
            for f2 in f_list:
                if f < f2 and self.b4[f, f2]: eps = (eps + 2) % 4
        
        for i, u in enumerate(u_vars):
            vu[i] = self.v4[u]
            for f in f_list:
                if self.b4[u, f]: vu[i] = (vu[i] + 2) % 4

        if nu == 0: return {0:1, 1:1j, 2:-1, 3:-1j}[eps]

        # Diagonalize bilinear form of internal variables
        B_u = self.b4[np.ix_(u_vars, u_vars)]
        rank, p = 0, 0
        while p < nu:
            pivot = next((i for i in range(p, nu) if B_u[i, p]), None) # Non-alternating pivot 
            if pivot is None: break
            B_u[[p, pivot]] = B_u[[pivot, p]]; B_u[:, [p, pivot]] = B_u[:, [pivot, p]]
            vu[p], vu[pivot] = vu[pivot], vu[p]
            for k in range(p+1, nu):
                if B_u[k, p]: B_u[k] ^= B_u[p]; B_u[:, k] ^= B_u[:, p]; vu[k] = (vu[k] + vu[p]) % 4
            rank += 1; p += 1

        if any(vu[rank:] == 2): return 0j # Balanced condition 
        
        n1 = np.count_nonzero(vu[:rank] == 1)
        n0 = rank - n1
        s = (2**(nu-rank)) * ((1+1j)**n0) * ((1-1j)**n1)
        return {0:1, 1:1j, 2:-1, 3:-1j}[eps] * s * (2**(-self.h_count/2))

    def print_analytic(self, transition_mode=True):
        """
        Prints the analytical XOR conditions for the amplitude summation based on 
        Schmidt's Z4 exponential sum calculation.
        
        The amplitude is derived from the closed-form expression:
        Sum = 2^(m-r) * (1+i)^N0 * (1-i)^N1 
        
        Where:
        - m: Number of internal variables (nu).
        - r: Rank of the non-alternating bilinear form[.
        - N1: Count of sum variables u_j = 1 in the canonical form.
        - N0: Count of sum variables u_j = 0 in the canonical form (r - N1).
        """
        if not self.compiled:
            print("Circuit must be compiled first.")
            return

        # Identify internal variables (not fixed by input or output)
        # Note: In Z4Sim, we assume standard input/output mapping for the analytic view
        fixed_indices = set(range(self.n)) | set(self.output_vars)
        u_vars = [i for i in range(self.n_vars) if i not in fixed_indices]
        nu = len(u_vars)
        
        # Coefficient matrix: [x_bits | y_bits | constant]
        n_cols = 2 * self.n + 1 if transition_mode else self.n + 1
        coeff_matrix = np.zeros((nu, n_cols), dtype=int)
        
        # Initialize linear contributions based on b4 (quadratic weights) and v4 (linear weights)
        for ui, orig_u in enumerate(u_vars):
            coeff_matrix[ui, -1] = self.v4[orig_u] % 4
            if transition_mode:
                for xi in range(self.n):
                    if self.b4[orig_u, xi]: 
                        coeff_matrix[ui, xi] = 2 # Quadratic interaction adds 2 in Z4 
            
            for yi in range(self.n):
                if self.b4[orig_u, self.output_vars[yi]]:
                    col_idx = self.n + yi if transition_mode else yi
                    coeff_matrix[ui, col_idx] = (coeff_matrix[ui, col_idx] + 2) % 4

        # Perform the Z4 diagonalization (Schmidt's Basis Change) on coefficients 
        B_u = self.b4[np.ix_(u_vars, u_vars)].copy()
        rank, p = 0, 0
        while p < nu:
            pivot = next((i for i in range(p, nu) if B_u[i, p]), None) 
            if pivot is None: break
            
            # Apply swaps to both Bilinear matrix and Analytic coefficients
            B_u[[p, pivot]] = B_u[[pivot, p]]
            B_u[:, [p, pivot]] = B_u[:, [pivot, p]]
            coeff_matrix[[p, pivot]] = coeff_matrix[[pivot, p]]
            
            # Eliminate entries
            for k in range(p+1, nu):
                if B_u[k, p]:
                    B_u[k] ^= B_u[p]
                    B_u[:, k] ^= B_u[:, p]
                    coeff_matrix[k] = (coeff_matrix[k] + coeff_matrix[p]) % 4
            rank += 1; p += 1

        def build_xor_expr(j):
            # In Z4, linear bits (XOR) correspond to weight 2 [cite: 652]
            terms = ["1"] if coeff_matrix[j, -1] in (1, 3) else []
            if transition_mode:
                terms += [f"x_{xi}" for xi in range(self.n) if coeff_matrix[j, xi] == 2]
                terms += [f"y_{yi}" for yi in range(self.n) if coeff_matrix[j, self.n + yi] == 2]
            else:
                terms += [f"y_{yi}" for yi in range(self.n) if coeff_matrix[j, yi] == 2]
            return " XOR ".join(terms) if terms else "0"

        title = "<y|U|x>" if transition_mode else "Statevector <y|U|0>"
        print(f"\n{'='*60}\n   ANALYTICAL FUNCTIONS FOR {title} (Z4)\n{'='*60}")
        print(f"Hadamard Count (h)         : {self.h_count}")
        print(f"Norm Factor                : 2^(-{self.h_count}/2)")
        print(f"Schmidt Rank (r)           : {rank}")
        
        n_sum_terms = [f"v_{j} = ({build_xor_expr(j)})" for j in range(rank)]
        print(f"\n1. EXPONENTIAL SUM VARIABLES (v_j)")
        if n_sum_terms:
            for term in n_sum_terms: print(f"   {term}")
            print(f"\n   N1 = count(v_j == 1), N0 = {rank} - N1")
        else: print("   No sum variables (Rank 0).")
        
        zero_conditions = [build_xor_expr(k) for k in range(rank, nu)]
        print("\n2. ZERO AMPLITUDE CONDITIONS (Balanced Function) ")
        if zero_conditions:
            for idx, cond in enumerate(zero_conditions): 
                print(f"   Kernel_{idx}: ({cond}) == 1")
        else: print("   No kernel variables exist.")
        print("="*60 + "\n")

def test_z4_s_gate():
    sim = Z4Sim(1)
    sim.s(0); sim.compile()
    res = sim.get_amplitude(y=1, x=1) # <1|S|1>
    print(f"S-gate Result: {res}")
    # Fix: res should be exactly 1j
    assert np.isclose(res.real, 0.0) and np.isclose(res.imag, 1.0)
    print("✅ test_z4_s_gate passed!")
    



def test_bell_state_00():
    bell = Z4Sim(2)
    
    bell.h(0)
    
    bell.h(1)
    bell.cz(0,1)
    bell.h(1)
    bell.compile()
    
    amp_00 = bell.get_amplitude(0, 0) # <00|Bell|00>
    amp_11 = bell.get_amplitude(3, 0) # <11|Bell|00>
    
    print(f"Bell State Amplitudes: <00|Bell|00> = {amp_00}, <11|Bell|00> = {amp_11}")
    
    assert np.isclose(amp_00, 1/np.sqrt(2))
    assert np.isclose(amp_11, 1/np.sqrt(2))
    bell.print_analytic()
    
    
def test_bell_state_01():
    bell = Z4Sim(2)
    
    bell.h(0)
    
    bell.h(1)
    bell.cz(0,1)
    bell.h(1)
    
    bell.z(0) # Apply Z to the first qubit to flip the phase of  |11>
    
    bell.compile()
    sv = []
    expected_sv = [ 1/np.sqrt(2),0, 0, -1/np.sqrt(2)] # |01> + (-|11>) / sqrt(2)
    for i in range(4):
        sv.append(bell.get_amplitude(i, 0))
        assert np.isclose(fidelity(sv[-1], expected_sv[i]), 1.0), f"Expected {expected_sv[i]}, got {sv[-1]}"
        print(f"<{i:02b}|Bell|00> = {bell.get_amplitude(i, 0)}")
    
    bell.print_analytic()
    
def test_bell_state_11():
    bell = Z4Sim(2)
    
    bell.h(0)
    
    bell.h(1)
    bell.cz(0,1)
    bell.h(1)
    
    bell.z(1)
    
    bell.h(0)
    bell.z(0) 
    bell.h(0)
    
    bell.compile()
    sv = []
    expected_sv = [0, 1/np.sqrt(2), -1/np.sqrt(2),0] 
    
    # for i in range(4):
    #     sv.append(bell.get_amplitude(i, 0))
    
    # fiedility = fidelity(sv, expected_sv)
    # assert np.isclose(fiedility, 1.0)
        
    
    bell.print_analytic()
    
def fidelity(state1, state2):
    """Computes the fidelity between two quantum states."""
    
    return np.abs(np.dot(state1.conj(), state2))**2

def test_print_analytic():
    test_z4_s_gate()
    
    # Bell State Example with Analytic Representation
    bell = Z4Sim(2)
    bell.h(0)
    
    bell.h(1)
    bell.cz(0,1)
    bell.h(1)
    bell.compile()
    print(f"Bell <00|U|00>: {bell.get_amplitude(0, 0)}")
    bell.print_analytic()
    
    assert np.isclose(bell.get_amplitude(0, 0), 1/np.sqrt(2))
    assert np.isclose(bell.get_amplitude(3, 0), 1/np.sqrt(2))
    


# test_bell_state_11()

def test_ghz_state():
    ghz = Z4Sim(10)
    ghz.h(0)
    for i in range(1, 10):
        ghz.h(i)
        ghz.cz(0, i)
        ghz.h(i)
    ghz.compile()
    
    
    # for i in range(4):
    #     amp = ghz.get_amplitude(i, 0)
    #     print(f"<{i:010b}|GHZ|0000000000> = {amp}")
    
    amp_000 = ghz.get_amplitude(0, 0) # <000|GHZ|000>
    amp_111 = ghz.get_amplitude(1023, 0) # <111|GHZ|000>
    
    print(f"GHZ State Amplitudes: <000|GHZ|000> = {amp_000}, <111|GHZ|000> = {amp_111}")
    
    assert np.isclose(amp_000, 1/np.sqrt(2))
    assert np.isclose(amp_111, 1/np.sqrt(2))
    
test_ghz_state()