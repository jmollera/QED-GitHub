from qed.eng_elec import Network, MutualCoupling, VoltageSource, CurrentSource, Impedance

circuit = Network()

circuit.add(VoltageSource(from_node=0, to_node=1, E=1.10, Z=0.25j), branch=1)
circuit.add(VoltageSource(from_node=0, to_node=2, E=1.05+0.10j, Z=0.20j), branch=2)
circuit.add(VoltageSource(from_node=0, to_node=3, E=1.08+0.12j, Z=0.25j), branch=3)
circuit.add(Impedance(from_node=3, to_node=4, Z=0.10j), branch=4)
circuit.add(Impedance(from_node=2, to_node=4, Z=0.405j), branch=5)
circuit.add(Impedance(from_node=2, to_node=3, Z=0.50j), branch=6)
circuit.add(Impedance(from_node=1, to_node=2, Z=0.16j), branch=7)
circuit.add(CurrentSource(from_node=4, to_node=0, J=2-0.9j, Y=1/(-25j)), branch=8)
circuit.add(MutualCoupling(XM=0.05j), branch=5, with_branch=6)

circuit.solve()

print(circuit.results(polar=True, u_unit='pu', u_fmt='.4f', i_unit='pu', i_fmt='.4f'))

s_G1 = 1.10*circuit.branch_current(1).conjugate()
s_G2 = (1.05+0.10j)*circuit.branch_current(2).conjugate()
s_G3 = (1.08+0.12j)*circuit.branch_current(3).conjugate()
s_Q8 = circuit.branch_voltage(8)*circuit.branch_current(8).conjugate()
print(f'\ns_G1 = {s_G1:.4f}')
print(f's_G2 = {s_G2:.4f}')
print(f's_G3 = {s_G3:.4f}')
print(f's_Q8 = {s_Q8:.4f}')
