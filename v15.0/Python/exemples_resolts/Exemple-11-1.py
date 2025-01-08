from qed.eng_elec import Network, VoltageSource, CurrentSource, Impedance, MutualCoupling

circuit = Network()

circuit.add(VoltageSource(E=1.10, Z=0.25j), from_to=(0, 1), branch=1)
circuit.add(VoltageSource(E=1.05 + 0.10j, Z=0.20j), from_to=(0, 2), branch=2)
circuit.add(VoltageSource(E=1.08 + 0.12j, Z=0.25j), from_to=(0, 3), branch=3)
circuit.add(Impedance(Z=0.10j), from_to=(3, 4), branch=4)
circuit.add(Impedance(Z=0.405j), from_to=(2, 4), branch=5)
circuit.add(Impedance(Z=0.50j), from_to=(2, 3), branch=6)
circuit.add(Impedance(Z=0.16j), from_to=(1, 2), branch=7)
circuit.add(CurrentSource(J=2 - 0.9j, Y=1 / (-25j)), from_to=(4, 0), branch=8)
circuit.add(MutualCoupling(X=0.05j), coupled_branches=(5, 6))

circuit.solve()

print(circuit.results(polar=True, u_unit='pu', u_fmt='.4f', i_unit='pu', i_fmt='.4f'))

s_G1 = 1.10*circuit.current(branch=1).conjugate()
s_G2 = (1.05+0.10j)*circuit.current(branch=2).conjugate()
s_G3 = (1.08+0.12j)*circuit.current(branch=3).conjugate()
s_Q8 = circuit.voltage(branch=8)*circuit.current(branch=8).conjugate()
print(f'\ns_G1 = {s_G1:.4f}')
print(f's_G2 = {s_G2:.4f}')
print(f's_G3 = {s_G3:.4f}')
print(f's_Q8 = {s_Q8:.4f}')
