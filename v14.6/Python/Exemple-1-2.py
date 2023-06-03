from qed.utils import ComplexD
import qed.eng_elec as ee

U_AN = 230  # V
U_BN = ComplexD(230, 240)  # V
U_CN = ComplexD(230, 120)  # V

R_A = 50e-3  # Ω
R_B = 80e-3  # Ω
R_C = 70e-3  # Ω
R_N = 46  # Ω

U_GN = ee.millman([U_AN, U_BN, U_CN, 0], [R_A, R_B, R_C, R_N])
I_A = (U_AN - U_GN)/R_A
I_B = (U_BN - U_GN)/R_B
I_C = (U_CN - U_GN)/R_C
I_N = U_GN/R_N
print('Teorema de Millman')
print('-'*30)
print(f'U_GN = {ComplexD(U_GN):.4f/.2f} V')
print(f'I_A  = {ComplexD(I_A):.4f/.2f} A')
print(f'I_B  = {ComplexD(I_B):.4f/.2f} A')
print(f'I_C  = {ComplexD(I_C):.4f/.2f} A')
print(f'I_N  = {ComplexD(I_N):.4f/.2f} A\n')

circuit = ee.Network()
circuit.add(ee.VoltageSource(from_to=(1, 0), E=U_AN, Z=R_A), branch=1)
circuit.add(ee.VoltageSource(from_to=(1, 0), E=U_BN, Z=R_B), branch=2)
circuit.add(ee.VoltageSource(from_to=(1, 0), E=U_CN, Z=R_C), branch=3)
circuit.add(ee.Impedance(from_to=(0, 1), Z=R_N), branch=4)
circuit.solve()
U_GN = circuit.voltage(branch=4)
I_A = circuit.current(branch=1)
I_B = circuit.current(branch=2)
I_C = circuit.current(branch=3)
I_N = circuit.current(branch=4)
print('Mètode dels nusos')
print('-'*30)
print(f'U_GN = {ComplexD(U_GN):.4f/.2f} V')
print(f'I_A  = {ComplexD(I_A):.4f/.2f} A')
print(f'I_B  = {ComplexD(I_B):.4f/.2f} A')
print(f'I_C  = {ComplexD(I_C):.4f/.2f} A')
print(f'I_N  = {ComplexD(I_N):.4f/.2f} A')
