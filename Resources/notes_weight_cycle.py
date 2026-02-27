"""
̈-Theta Reduction sub Routine

A(O) = sigma { Œ^P(I,M,O) }

I:input space
M:intermediate space
O:output space

Œ: exp(i.2pi/n)
n: base of ring polynomials

Typically I = 0
P(I,M,O) = f(I,M,O) + g(M)

Compute g(M) and cache it as hash map G: M -> g(M) #  Ø(K.2^|M|) // parallalizable

when required calculate f(I,M',O) + G(M) # Ø(K.2^|M'+O|)) + Ø(1) // parallalizable

net complexity reduction: Ø(K.2^|M|) + Ø(K.2^|M'+O|)) + Ø(1) << Ø(K.2^|I+M+O|))

Gray Code Iteration

FUNCTION Simulate_Optimized_Z8(Polynomial P, Input_Bitstring x_in):
    1. PARTITION VARIABLES
       Identify indices for x_in (n bits), x_int (m bits), and x_out (n bits).
       m = num_hadamards.

    2. PRE-PROCESS TERMS
       - Group terms that ONLY involve (x_in, x_int). Call this P_core(x_int).
       - Group terms that involve interactions between x_out and x_int. 
         These are typically linear in x_out: Σ (w_j * x_{out,i} * x_{int,k}).

    3. COMPUTE BASE DISTRIBUTION (The "Kernel")
       Initialize Base_Dist[2^n][8] to zeros # One Z8 distribution per possible output assignment
       
       # This is the heavy lifting: O(2^m)
       FOR each assignment of intermediate variables s ∈ {0, 1}^m:
           - Calculate the internal phase: φ_int = P_core(s) mod 8.
           - Determine the output bitstring x_out produced by this path.
           - Base_Dist[x_out][φ_int] += 1.

    4. APPLY CYCLICAL SHIFTS (Linear Adjustments)
       # Adjust for terms involving x_in and x_out interaction: O(2^n)
       FOR each output state k in {0...2^n-1}:
           - Calculate shift = Σ (weight * x_in * x_out) mod 8.
           - IF shift != 0:
               Base_Dist[k] = Cyclical_Shift(Base_Dist[k], shift).

    5. COMPUTE AMPLITUDES (Z8 FFT)
       Initialize StateVector[2^n]
       Norm = 2^(-m/2).
       
       FOR each k in {0...2^n-1}:
           Amplitude = 0
           FOR j from 0 to 7:
               Phase_Angle = j * π / 4
               Amplitude += Base_Dist[k][j] * exp(i * Phase_Angle).
           StateVector[k] = Amplitude * Norm.

    RETURN StateVector
    
# 1. Theoretical Basis: The Variable Partition
A phase polynomial P is defined over three distinct variable sets:
- Fixed Variables (X_in): Determined by the initial input state (e.g., |0...0>).
- Intermediate Variables (X_int): Variables 's' introduced by Hadamard gates that 
  must be summed over to calculate interference.
- Output Variables (X_out): Variables 'y' that define the computational basis 
  state being measured.

The amplitude for a specific output bitstring 'y' is:
Ψ(y) = (1/√2^m) * Σ_{s ∈ {0,1}^m} exp(i * (π/4) * P(x_in, s, y))

# 2. The "Shortcut" Approach: Distribution Kernels
Instead of evaluating the full polynomial for every output state (2^n * 2^m), 
we isolate the "Core" of the polynomial:

A. Partition the Polynomial:
   P(x_in, s, y) = P_core(x_in, s) + P_inter(s, y) + P_bias(x_in, y)
   
   - P_core: Terms containing only input and intermediate variables.
   - P_inter: Interaction terms (typically weight-2 or weight-4) between 
     intermediate and output variables.
   - P_bias: Direct shifts between input and output variables.

B. Build the Base Distribution (Z8 Histogram):
   We iterate through the 2^m intermediate paths ONCE. For each path, we 
   evaluate P_core and determine which output bitstring 'y' it maps to. 
   We store this in a "Frequency Table" (Histogram) of phases mod 8.

C. Cyclical Weight Shifting:
   The interaction terms (P_inter) do not create new interference paths; 
   they only shift the phases of existing paths. For a specific output 'y', 
   the phase shift is a linear offset: 
   Phase_final = (Phase_core + Shift_y) mod 8.

# 3. Complexity Reduction Analysis

| Parameter           | Baseline Algorithm       | Shortcut Algorithm         |
|---------------------|--------------------------|----------------------------|
| Summation Space     | 2^t (Total variables)    | 2^m (Hadamard variables)   |
| Evaluation Cost     | O(2^(n+m) * Terms)       | O(2^m * Terms_core)        |
| Basis Mapping       | Performed 2^(n+m) times  | Performed 2^m times        |
| Output Generation   | Implicit in loop         | O(2^n) Cyclical Shifts     |

# 4. Impact on Resource Scaling
- Time Complexity: The exponent is reduced from (n + m) to (m). 
  In circuits where n (qubits) is large but m (Hadamards) is localized 
  or sparse, this provides an exponential speedup of 2^n.
- Memory: Shifts the bottleneck from high-frequency path evaluation 
  to a fixed-size Z8 frequency table per output state.
"""