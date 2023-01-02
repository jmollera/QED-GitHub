from scipy import linalg
from qed.eng_elec import Network, VoltageSource, Impedance

Z_M = [[2, -1, 0],
       [-1, 4, -1],
       [0, -1, 2]]

E_M = [5, 0, -2]

I_M = linalg.solve(Z_M, E_M)
I = [I_M[0], -I_M[1], I_M[0] - I_M[1], I_M[1], I_M[1] - I_M[2], -I_M[2]]
print('Resolució algebraica')
print('-'*20)
for br, i in enumerate(I, 1):
    print(f'I{br} = {i:.2f} A')

circuit = Network()
circuit.add(VoltageSource(from_node=0, to_node=1, E=5, Z=1), branch=1)
circuit.add(Impedance(from_node=2, to_node=1, Z=1), branch=2)
circuit.add(Impedance(from_node=1, to_node=0, Z=1), branch=3)
circuit.add(Impedance(from_node=3, to_node=0, Z=1), branch=4)
circuit.add(Impedance(from_node=2, to_node=3, Z=1), branch=5)
circuit.add(VoltageSource(from_node=3, to_node=2, E=2, Z=1), branch=6)
circuit.solve()
print('\nMètode dels nusos')
print('-'*17)
for br in range(1, circuit.num_branches + 1):
    print(f'I{br} = {circuit.branch_current(br):.2f} A')
