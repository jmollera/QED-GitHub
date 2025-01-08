from qed.eng_elec import Network, VoltageSource, CurrentSource, Impedance, MutualCoupling, ShortCircuit

circuit = Network()

circuit.add(VoltageSource(E=200, Z=10), from_to=(0, 1), branch=1)
circuit.add(VoltageSource(E=-50, Z=20j), from_to=(1, 2), branch=2)
circuit.add(Impedance(Z=5j), from_to=(2, 0), branch=3)
circuit.add(Impedance(Z=20), from_to=(1, 2), branch=4)
circuit.add(CurrentSource(J=4, Y=1 / 10), from_to=(0, 2), branch=5)
circuit.add(MutualCoupling(X=5j), coupled_branches=(2, 3))
circuit.add(ShortCircuit(), from_to=(1, 2), branch=6)

circuit.solve()

print(circuit.results(u_fmt='z.4f', u_unit='V', i_fmt='z.4f', i_unit='A'))
