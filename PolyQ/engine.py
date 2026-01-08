# Core Logic File

"""
Core logic for simulating Clifford+S+T quantum circuits.

Provides functions for constructing polynomial representations of circuits,
evaluating polynomials, generating truth tables, and computing statevectors.
"""

import numpy as np

from polyq_backend import (
    create_poly as _create_poly_backend,
    get_statevector as _get_statevector_backend,
    get_statevector_file as _get_statevector_file_backend,
    get_truthtable as _get_truthtable_backend,
    get_truthtable_no_ivs as _get_truthtable_no_ivs_backend,
)

def create_poly(qc, n: int):
    """
    Create the polynomial array representation of a quantum circuit.

    Args:
        qc: QuantumCircuit object.
        n (int): Number of qubits.

    Returns:
        tuple: (terms, wire_array, t) where terms is the polynomial, wire_array tracks variable names, and t is the total number of variables.
    """
    instructions = [(instruction.operation.name,
                    [qc.find_bit(q).index for q in instruction.qubits]) 
                    for index, instruction in enumerate(qc.data)]

    terms, wire_array, t = _create_poly_backend(instructions, n)
    return terms, wire_array, t


def get_truthtable(terms, n, t, initial_state):
    """
    Generate the truth table for the polynomial given an initial state.

    Args:
        terms (list): Polynomial terms.
        n (int): Number of qubits.
        t (int): Total number of variables.
        initial_state (list): Initial state of the qubits.

    Returns:
        np.ndarray: Truth table values.
    """
    return np.asarray(_get_truthtable_backend(terms, n, t, list(initial_state)), dtype=np.uint8)

def get_truthtable_no_ivs(terms, n, t, initial_state):
    """
    Generate the truth table, removing input variables from the polynomial.

    Args:
        terms (list): Polynomial terms.
        n (int): Number of qubits.
        t (int): Total number of variables.
        initial_state (list): Initial state of the qubits.

    Returns:
        np.ndarray: Truth table values.
    """
    return np.asarray(_get_truthtable_no_ivs_backend(terms, n, t, list(initial_state)), dtype=np.uint8)

def get_statevector(ttb, n, t, ovs, starting_index=0):
    """
    Compute the statevector from the truth table.

    Args:
        ttb (np.ndarray): Truth table values.
        n (int): Number of qubits.
        t (int): Total number of variables.
        ovs (list): Output variable indices.
        starting_index (int, optional): Starting index. Defaults to 0.

    Returns:
        np.ndarray: Statevector as a complex array.
    """
    return np.asarray(_get_statevector_backend(list(ttb), n, t, list(ovs), starting_index), dtype=complex)

def get_statevector_file(ttb, n, t, ovs, starting_index=0):
    """
    Compute the statevector and write amplitudes to a file.

    Args:
        ttb (np.ndarray): Truth table values.
        n (int): Number of qubits.
        t (int): Total number of variables.
        ovs (list): Output variable indices.
        starting_index (int, optional): Starting index. Defaults to 0.

    Returns:
        str: Filename where the statevector is written.
    """
    return _get_statevector_file_backend(list(ttb), n, t, list(ovs), None, starting_index)

