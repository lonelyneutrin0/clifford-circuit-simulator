# clifford-circuit-simulator
Quantum circuit simulator for Clifford gate circuits in Cirq.

## Background 
Group theory provides a very useful to classically simulate Clifford gates within polynomial time. Consider the Pauli subgroup   $$\mathcal U(2^n) \geq \mathcal P = \langle X_i, Z_i \rangle_{i=1}^n$$ where $$X_i = I \otimes \cdots X \otimes \cdots I$$ In other words, the $X$ gate acting on the ith qubit. Also consider its normalizer group, 
$$\mathcal U(2^n) \geq \mathcal C = \langle CNOT_{ij}, S_i, H_i \rangle$$
Then, by definition, $\forall P \in \mathcal P, \forall C \in \mathcal C$, 
$$CPC^{\dagger} \in \mathcal P $$
Therefore, instead of tracking all Pauli matrices, one can simply track the 2n generators of the Pauli subgroup (along with bits for the phases) to evolve Clifford circuits classically.
