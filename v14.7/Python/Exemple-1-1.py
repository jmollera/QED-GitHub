import qed.eng_elec as ee

U_1 = 125.1  # V
U_2 = 124.8  # V
U_3 = 125.2  # V

R_1 = 34e-3  # Ω
R_2 = 41e-3  # Ω
R_3 = 29e-3  # Ω
R = 50  # Ω

Z_Th = ee.z_parallel([R_1, R_2, R_3])
E_Th = -ee.millman([-U_1, -U_2, -U_3], [R_1, R_2, R_3])
Y_No = 1/Z_Th
J_No = E_Th/Z_Th
I = E_Th / (Z_Th + R)
U = R * I
print('Teorema de Millman')
print('-'*18)
print(f'Z_Th = {Z_Th*1000:.4f} mΩ')
print(f'E_Th = {E_Th:.4f} V')
print(f'Y_No = {Y_No:.2f} S')
print(f'J_No = {J_No/1000:.2f} kA')
print(f'I = {I:.4f} A')
print(f'U = {U:.2f} V\n')

circuit = ee.Network()
circuit.add(ee.Branch(ee.VoltageSource(E=U_1, Z=R_1), from_to=(0, 1), branch=1))
circuit.add(ee.Branch(ee.VoltageSource(E=U_2, Z=R_2), from_to=(0, 1), branch=2))
circuit.add(ee.Branch(ee.VoltageSource(E=U_3, Z=R_3), from_to=(0, 1), branch=3))
circuit.solve()
E_Th, Z_Th = circuit.thevenin(1, 0)
J_No, Y_No = circuit.norton(1, 0)
print('Mètode dels nusos')
print('-'*18)
print(f'Z_Th = {Z_Th*1000:.4f} mΩ')
print(f'E_Th = {E_Th:.4f} V')
print(f'Y_No = {Y_No:.2f} S')
print(f'J_No = {J_No/1000:.2f} kA')
circuit.add(ee.Branch(ee.Impedance(Z=R), from_to=(1, 0), branch=4))
circuit.solve()
I = circuit.current(branch=4)
U = circuit.voltage(branch=4)
print(f'I = {I:.4f} A')
print(f'U = {U:.2f} V\n')
