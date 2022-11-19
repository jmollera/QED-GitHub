from qed.eng_elec import Network, VoltageSource, CurrentSource, Impedance

circuit = Network()

circuit.add(VoltageSource(from_node=0, to_node=1, E=200, Z=10), branch=1)
circuit.add(VoltageSource(from_node=1, to_node=2, E=-50, Z=20j), branch=2)
circuit.add(Impedance(from_node=2, to_node=0, Z=5j), branch=3)
circuit.add(Impedance(from_node=1, to_node=2, Z=20), branch=4)
circuit.add(CurrentSource(from_node=0, to_node=2, J=4, Y=1/10), branch=5)

circuit.solve()

for branch in (2, 5):
    Ibr = circuit.branch_current(branch)
    print(f'I{branch} = {Ibr:.4f} A')