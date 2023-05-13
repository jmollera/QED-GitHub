from qed.utils import ComplexD
import qed.eng_elec as ee

R_1 = 0.5  # Ω/fase
X_1 = 1.5  # Ω/fase
R_2 = 0.625  # Ω/fase
X_2 = 1.25  # Ω/fase
R_Fe = 360  # Ω/fase
X_m = 40  # Ω/fase
s_n = 0.05  # Lliscament nominal

U_AB, U_BC, U_CA = ee.triangle_to_phasors(399, 370, 370)
print(f'U_AB = {ComplexD(U_AB):.2f} V')
print(f'U_BC = {ComplexD(U_BC):.2f} V')
print(f'U_CA = {ComplexD(U_CA):.2f} V\n')

U_AB_0, U_AB_1, U_AB_2 = ee.ABC_to_A012(U_AB, U_BC, U_CA)
print(f'U_AB,1 = {ComplexD(U_AB_1):.2f} V')
print(f'U_AB,2 = {ComplexD(U_AB_2):.2f} V')
print(f'U_AB,0 = {ComplexD(U_AB_0):.2f} V\n')

U_AN_1, U_AN_2 = ee.AB12_to_AN12(U_AB_1, U_AB_2)
U_rel = abs(U_AN_2/U_AN_1)
print(f'U_AN,1 = {ComplexD(U_AN_1):.2f} V')
print(f'U_AN,2 = {ComplexD(U_AN_2):.2f} V')
print(f'U_AN,2/U_AN,1 = {U_rel:.3f}\n')

Z_0 = R_Fe*1j*X_m/(R_Fe + 1j*X_m)

def Z_mot_1(s):
    return R_1 + 1j*X_1 + Z_0*(R_2/s + 1j*X_2)/(Z_0 + R_2/s + 1j*X_2)

def I_1_1(s):
    return U_AN_1/Z_mot_1(s)

def Z_mot_2(s):
    return R_1 + 1j*X_1 + Z_0*(R_2/(2-s) + 1j*X_2)/(Z_0 + R_2/(2-s) + 1j*X_2)

def I_1_2(s):
    return U_AN_2/Z_mot_2(s)

Z_mot_n_1 = Z_mot_1(s=s_n)
Z_mot_n_2 = Z_mot_2(s=s_n)
Z_mot_arr_1 = Z_mot_1(s=1)
print(f'Z_mot,n,1 = {Z_mot_n_1:.3f} Ω')
print(f'Z_mot,n,2 = {Z_mot_n_2:.3f} Ω')
print(f'Z_mot,arr,1 = {Z_mot_arr_1:.3f} Ω\n')

I_1_n_1 = I_1_1(s=s_n)
I_1_n_2 = I_1_2(s=s_n)
I_rel = abs(I_1_n_2/I_1_n_1)
print(f'I_1,n,1 = {ComplexD(I_1_n_1):.2f} A')
print(f'I_1,n,2 = {ComplexD(I_1_n_2):.2f} A')
print(f'I_1,n,2/I_1,n,1 = {I_rel:.3f}\n')

I_1_n_A, I_1_n_B, I_1_n_C = ee.A012_to_ABC(0, I_1_n_1, I_1_n_2)
print(f'I_1,n,A = {ComplexD(I_1_n_A):.2f} A')
print(f'I_1,n,B = {ComplexD(I_1_n_B):.2f} A')
print(f'I_1,n,C = {ComplexD(I_1_n_C):.2f} A')
