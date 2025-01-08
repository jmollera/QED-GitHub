from qed.eng_elec import Network, VoltageSource, Impedance, ShortCircuit

circuit = Network()
circuit.add(VoltageSource(E=5, Z=1), from_to=(0, 1), branch=1)
circuit.add(Impedance(Z=1), from_to=(2, 1), branch=2)
circuit.add(Impedance(Z=1), from_to=(1, 0), branch=3)
circuit.add(Impedance(Z=1), from_to=(3, 0), branch=4)
circuit.add(Impedance(Z=1), from_to=(2, 3), branch=5)
circuit.add(VoltageSource(E=2, Z=1), from_to=(3, 2), branch=6)
circuit.add(ShortCircuit(), from_to=(2, 0), branch=7)

circuit.solve()

print(circuit.results(u_fmt='z.3f', u_unit='V', i_fmt='z.3f', i_unit='A'))
