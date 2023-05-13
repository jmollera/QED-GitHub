from qed.utils import ComplexD
import qed.eng_elec as ee

R = 10  # Ω

U_AB, U_BC, U_CA = ee.triangle_to_phasors(2760, 1840, 2300)  # V
print(f'U_AB = {ComplexD(U_AB):.2f} V')
print(f'U_BC = {ComplexD(U_BC):.2f} V')
print(f'U_CA = {ComplexD(U_CA):.2f} V\n')

U_AB_0, U_AB_1, U_AB_2 = ee.ABC_to_A012(U_AB, U_BC, U_CA)
print(f'U_AB_1 = {ComplexD(U_AB_1):.2f} V')
print(f'U_AB_2 = {ComplexD(U_AB_2):.2f} V')
print(f'U_AB_0 = {ComplexD(U_AB_0):.2f} V\n')

U_AG_1, U_AG_2 = ee.AB12_to_AN12(U_AB_1, U_AB_2)
U_AG_0 = 0
print(f'U_AG_1 = {ComplexD(U_AG_1):.2f} V')
print(f'U_AG_2 = {ComplexD(U_AG_2):.2f} V')
print(f'U_AG_0 = {ComplexD(U_AG_0):.2f} V\n')

I_A_1 = U_AG_1/R
I_A_2 = U_AG_2/R
I_A_0 = U_AG_0/R
print(f'I_A_1 = {ComplexD(I_A_1):.2f} A')
print(f'I_A_2 = {ComplexD(I_A_2):.2f} A')
print(f'I_A_0 = {ComplexD(I_A_0):.2f} A\n')

S_3F = 3*(U_AG_0*I_A_0.conjugate() + U_AG_1*I_A_1.conjugate() + U_AG_2*I_A_2.conjugate())
print(f'S_3F = {ComplexD(S_3F)/1000:.2f} kVA\n')

U_AG, U_BG, U_CG = ee.A012_to_ABC(U_AG_0, U_AG_1, U_AG_2)
print(f'U_AG = {ComplexD(U_AG):.2f} V')
print(f'U_BG = {ComplexD(U_BG):.2f} V')
print(f'U_CG = {ComplexD(U_CG):.2f} V\n')

U_AG, U_BG, U_CG = ee.LL_to_LG(U_AB, U_BC, U_CA)
print(f'U_AG = {ComplexD(U_AG):.2f} V')
print(f'U_BG = {ComplexD(U_BG):.2f} V')
print(f'U_CG = {ComplexD(U_CG):.2f} V')
