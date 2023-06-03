from qed.eng_elec import Network, VoltageSource, CurrentSource, Impedance

circuit = Network()

circuit.add(VoltageSource(from_to=(0, 1), E=200, Z=10), branch=1)
circuit.add(VoltageSource(from_to=(1, 2), E=-50, Z=20j), branch=2)
circuit.add(Impedance(from_to=(2, 0), Z=5j), branch=3)
circuit.add(Impedance(from_to=(1, 2), Z=20), branch=4)
circuit.add(CurrentSource(from_to=(0, 2), J=4, Y=1/10), branch=5)

circuit.solve()

for b in (2, 5):
    Ib = circuit.current(branch=b)
    print(f'I{b} = {Ib:.4f} A')
