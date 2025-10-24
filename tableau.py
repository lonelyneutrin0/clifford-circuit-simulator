import cirq 
import numpy as np 

from dataclasses import dataclass 

class Tableau: 
    """ 
    Clifford tableau for representing Pauli generators.
    It is a 2n x (2n + 1) binary matrix.

    Structure: 
    - Rows 1 to n: Destabilizers
    - Rows n+1 to 2n: Stabilizers
    - Column 2n+1: Phase bits, 1 if negative phase and 0 if positive.

    The ith row tells you the Pauli string for the ith generator, that is, if 
    R_j = ±P_1 ⊗ P_2 ⊗ ... ⊗ P_n, then:
    x_ij, z_ij = 0,0 -> P_j = I
    x_ij, z_ij = 1,0 -> P_j = X
    x_ij, z_ij = 0,1 -> P_j = Z
    x_ij, z_ij = 1,1 -> P_j = Y
    """

    n_qubits: int
    """ Number of qubits in the tableau."""

    matrix: np.ndarray
    """ The tableau matrix itself."""

    def __init__(self, n_qubits: int): 
        self.n_qubits = n_qubits

        self.matrix = np.eye((2 * n_qubits), (2 * n_qubits + 1), dtype=int)

    def _rowsum(self, h: int, i: int): 
        """Rowsum subroutine for tableau update."""
        # Zero indexing
        h -= 1
        i -= 1

        if h > 2 * self.n_qubits or i > 2 * self.n_qubits: 
            raise ValueError("Row index out of bounds.")
        
        # x_hj = x_hj ⊕ x_ij
        self.matrix[h, :self.n_qubits] ^= self.matrix[i, :self.n_qubits]

        # z_hj = z_hj ⊕ z_ij
        self.matrix[h, self.n_qubits:2*self.n_qubits] ^= self.matrix[i, self.n_qubits:2*self.n_qubits]

        def _g(x1, x2, z1, z2) -> np.ndarray:
            result = np.zeros_like(x1, dtype=int)
            
            # Case: x1 == 0, z1 == 0
            mask = (x1 == 0) & (z1 == 0)
            result[mask] = 0
            
            # Case: x1 == 1, z1 == 1
            mask = (x1 == 1) & (z1 == 1)
            result[mask] = z2[mask] - x2[mask]
            
            # Case: x1 == 1, z1 == 0
            mask = (x1 == 1) & (z1 == 0)
            result[mask] = z2[mask] * (2 * x2[mask] - 1)
            
            # Case: x1 == 0, z1 == 1
            mask = (x1 == 0) & (z1 == 1)
            result[mask] = x2[mask] * (1 - 2 * z2[mask])
            
            return result

        g = _g(
            self.matrix[h, :self.n_qubits], 
            self.matrix[i, :self.n_qubits], 
            self.matrix[h, self.n_qubits:2*self.n_qubits], 
            self.matrix[i, self.n_qubits:2*self.n_qubits]
            ).sum()

        p = (2 * self.matrix[h, -1] + 2 * self.matrix[i, -1] + g) % 4
        if p == 0:
            self.matrix[h, -1] = 0
        elif p == 2:
            self.matrix[h, -1] = 1