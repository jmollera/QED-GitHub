from qed.eng_elec import Network, Branch, VoltageSource, CurrentSource, Impedance

circuit = Network()

circuit.add(Branch(VoltageSource(E=200, Z=10), from_to=(0, 1), branch=1))
circuit.add(Branch(VoltageSource(E=-50, Z=20j), from_to=(1, 2), branch=2))
circuit.add(Branch(Impedance(Z=5j), from_to=(2, 0), branch=3))
circuit.add(Branch(Impedance(Z=20), from_to=(1, 2), branch=4))
circuit.add(Branch(CurrentSource(J=4, Y=1 / 10), from_to=(0, 2), branch=5))

circuit.solve()

for b in (2, 5):
    Ib = circuit.current(branch=b)
    print(f'I{b} = {Ib:.4f} A')
