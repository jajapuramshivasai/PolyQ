import numpy as np

def fidelity(state1, state2):
    """Computes the fidelity between two quantum states."""
    
    return np.abs(np.dot(state1.conj(), state2))**2

