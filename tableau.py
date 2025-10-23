import cirq 
import numpy as np 

from dataclasses import dataclass 

class Tableau: 
    """ 
    Clifford tableau for representing Pauli generators.
    It is a 2n x (2n + 1) binary matrix.

    Structure: 
    - Rows 0 to n-1: X stabilizers
    - Rows n to 2n-1: Z stabilizers
    - Columns 0 to n-1: X components
    - Columns n to 2n-1: Z components
    - Last column: Phase bits
    """

    n_qubits: int
    """ Number of qubits in the tableau."""

    matrix: np.ndarray
    """ The tableau matrix itself."""

    def __init__(self, n_qubits: int): 
        self.n_qubits = n_qubits
        self.matrix = np.zeros((2 * n_qubits, 2 * n_qubits + 1), dtype=int)
        self._initialize_tableau() 

    def _initialize_tableau(self):
        """Initialize the tableau to |00000....> state."""

        self.matrix[:self.n_qubits, :self.n_qubits] = np.eye(self.n_qubits, dtype=int)  # X stabilizers
        self.matrix[self.n_qubits:, self.n_qubits:2*self.n_qubits] = np.eye(self.n_qubits, dtype=int)  # Z stabilizers
    
    def cnot(self, control: int, target: int): 
        """Apply CNOT gate with given control and target qubits."""

        # If the target qubit has an X generator, add the control's X to it
        for row in range(2 * self.n_qubits): 
            if self.matrix[row, target] == 1:
                self.matrix[row, control] ^= 1
        
        # If the control qubit 
        for row in range(2 * self.n_qubits):
            if self.matrix[row, self.n_qubits + control] == 1:
                self.matrix[row, self.n_qubits + target] ^= 1
    
    def hadamard(self, qubit: int): 
        """Apply Hadamard gate on the given qubit."""

        for row in range(2 * self.n_qubits):
            x = self.matrix[row, qubit]
            z = self.matrix[row, self.n_qubits + qubit]
            self.matrix[row, qubit] = z
            self.matrix[row, self.n_qubits + qubit] = x

            if x == 1 and z == 1:
                self.matrix[row, 2 * self.n_qubits] ^= 1 
            
    def phase(self, qubit: int):
        """Apply Phase (S) to the given qubit."""
        for row in range(2 * self.n_qubits):
            if self.matrix[row, qubit] == 1:
                self.matrix[row, self.n_qubits + qubit] ^= 1

                self.matrix[row, 2 * self.n_qubits] ^= 1
    
    def apply_gate(self, gate, *qubits): 
        """Apply a gate to the tableau"""

        if gate.lower() == 'cnot' or gate.lower() == 'cx': 
            self.cnot(qubits[0], qubits[1])
        elif gate.lower() == 'h' or gate.lower() == 'hadamard':
            self.hadamard(qubits[0])
        elif gate.lower() == 's' or gate.lower() == 'phase':
            self.phase(qubits[0])
        else:
            raise ValueError(f"Unsupported gate: {gate}")
        
    def __str__(self): 
        """Pretty print the tableau."""
        result = "Tableau:\n"
        result += "X block | Z block | Phase\n"
        result += "-" * (4 * self.n_qubits + 10) + "\n"
        
        for i in range(2 * self.n_qubits):

            x_block = " ".join(str(self.matrix[i, j]) for j in range(self.n_qubits))

            z_block = " ".join(str(self.matrix[i, j]) for j in range(self.n_qubits, 2 * self.n_qubits))

            phase = str(self.matrix[i, 2 * self.n_qubits])

            generator_type = "X" + str(i) if i < self.n_qubits else "Z" + str(i - self.n_qubits)
            result += f"{generator_type:2}: {x_block:>8} | {z_block:>8} | {phase:>3}\n"
        
        return result
    

