from tableau import Tableau

# Create a 3-qubit tableau
tab = Tableau(3)
print("Initial state:")
print(tab)

# Apply some gates
tab.apply_gate('h', 0)  # Hadamard on qubit 0
print("\nAfter H(0):")
print(tab)

tab.apply_gate('cnot', 0, 1)  # CNOT(0,1)
print("\nAfter CNOT(0,1):")
print(tab)

tab.apply_gate('s', 2)  # Phase gate on qubit 2
print("\nAfter S(2):")
print(tab)