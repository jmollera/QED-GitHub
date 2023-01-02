import qed.eng_elec as ee

U_1 = 125.1  # V
U_2 = 124.8  # V
U_3 = 125.2  # V

R_1 = 34e-3  # Ω
R_2 = 41e-3  # Ω
R_3 = 29e-3  # Ω
R_Q = 50  # Ω

Z_Th = ee.z_parallel([R_1, R_2, R_3])
E_Th = -ee.millman([-U_1, -U_2, -U_3], [R_1, R_2, R_3])
Y_No = 1/Z_Th
J_No = E_Th/Z_Th
I_Q = E_Th/(Z_Th + R_Q)
U_Q = R_Q*I_Q
print('Teorema de Millman')
print('-'*18)
print(f'Z_Th = {Z_Th*1000:.4f} mΩ')
print(f'E_Th = {E_Th:.4f} V')
print(f'Y_No = {Y_No:.2f} S')
print(f'J_No = {J_No/1000:.2f} kA')
print(f'I_Q  = {I_Q:.4f} A')
print(f'U_Q  = {U_Q:.2f} V\n')

circuit = ee.Network()
circuit.add(ee.VoltageSource(from_node=0, to_node=1, E=U_1, Z=R_1), branch=1)
circuit.add(ee.VoltageSource(from_node=0, to_node=1, E=U_2, Z=R_2), branch=2)
circuit.add(ee.VoltageSource(from_node=0, to_node=1, E=U_3, Z=R_3), branch=3)
circuit.solve()
E_Th, Z_Th = circuit.thevenin(1, 0)
J_No, Y_No = circuit.norton(1, 0)
print('Mètode dels nusos')
print('-'*18)
print(f'Z_Th = {Z_Th*1000:.4f} mΩ')
print(f'E_Th = {E_Th:.4f} V')
print(f'Y_No = {Y_No:.2f} S')
print(f'J_No = {J_No/1000:.2f} kA')
circuit.add(ee.Impedance(from_node=1, to_node=0, Z=R_Q), branch=4)
circuit.solve()
I_Q = circuit.branch_current(branch=4)
U_Q = circuit.branch_voltage(branch=4)
print(f'I_Q  = {I_Q:.4f} A')
print(f'U_Q  = {U_Q:.2f} V\n')
