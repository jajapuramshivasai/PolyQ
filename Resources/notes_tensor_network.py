"""
    
   
    ┌───────────────────────────────────────┐
    │Quantum Circuit (Matrix Product State) │
    └───────────────────────────────────────┘

    Qubit 0: ─[H]─[X]─[Y]─[Z]─[S]─[T]─

    Qubit 1: ─[RX]─[RY]─[RZ]─[H]─[X]─

    Qubit 2: ─[Y]─[Z]─[S]─[T]─[H]─[X]─
    
    nor P(0, M, O) = F(I_0, M_0, O_0) + F(I_1, M_1, O_1) + F(I_2, M_2, O_2)
    where V_0, V_1, V_2 are the variable spaces of qubits 0, 1, 2 respectively

    note: all the {V}_i are decoupled ie suitable for MPS style Simulation
    
    Fact: While 1D Tensor Network (MPS) contraction is efficient (P-time), 
    the exact contraction of 2D Tensor Networks (PEPS) and general 
    graphs is #P-complete.

   
    
"""