from qed.eng_elec import Network, MutualCoupling, VoltageSource, CurrentSource, Impedance

circuit = Network()

circuit.add(VoltageSource(from_to=(0, 1), E=200, Z=10), branch=1)
circuit.add(VoltageSource(from_to=(1, 2), E=-50, Z=20j), branch=2)
circuit.add(Impedance(from_to=(2, 0), Z=5j), branch=3)
circuit.add(Impedance(from_to=(1, 2), Z=20), branch=4)
circuit.add(CurrentSource(from_to=(0, 2), J=4, Y=1/10), branch=5)
circuit.add(MutualCoupling(XM=5j), coupled_branches=(2, 3))

circuit.solve()

eth, zth = circuit.thevenin(1, 2)
jno, yno = circuit.norton(1, 2)
print(f'E_Th(1,2) = {eth:.4f} V')
print(f'Z_Th(1,2) = {zth:.4f} Ω')
print(f'J_No(1,2) = {jno:.4f} A')
print(f'Y_No(1,2) = {yno:.4f} S')
