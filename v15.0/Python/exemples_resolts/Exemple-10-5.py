from qed.utils import Complex
from qed.eng_elec import triangle_to_phasors, ABC_to_A012, AB12_to_AN12, A012_to_ABC

R_1 = 0.5  # Ω/fase
X_1 = 1.5  # Ω/fase
R_2 = 0.625  # Ω/fase
X_2 = 1.25  # Ω/fase
R_Fe = 360  # Ω/fase
X_m = 40  # Ω/fase
s_n = 0.05  # Lliscament nominal

U_AB, U_BC, U_CA = triangle_to_phasors(399, 370, 370)
print(f'U_AB = {Complex(U_AB) :z.2f} V')
print(f'U_BC = {Complex(U_BC) :z.2f} V')
print(f'U_CA = {Complex(U_CA) :z.2f} V\n')

U_AB_0, U_AB_1, U_AB_2 = ABC_to_A012(U_AB, U_BC, U_CA)
print(f'U_AB,1 = {Complex(U_AB_1) :z.2f} V')
print(f'U_AB,2 = {Complex(U_AB_2) :z.2f} V')
print(f'U_AB,0 = {Complex(U_AB_0) :z.2f} V\n')

U_AN_1, U_AN_2 = AB12_to_AN12(U_AB_1, U_AB_2)
U_rel = abs(U_AN_2/U_AN_1)
print(f'U_AN,1 = {Complex(U_AN_1) :z.2f} V')
print(f'U_AN,2 = {Complex(U_AN_2) :z.2f} V')
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
print(f'Z_mot,n,1 = {Z_mot_n_1:z.3f} Ω')
print(f'Z_mot,n,2 = {Z_mot_n_2:z.3f} Ω')
print(f'Z_mot,arr,1 = {Z_mot_arr_1:z.3f} Ω\n')

I_1_n_1 = I_1_1(s=s_n)
I_1_n_2 = I_1_2(s=s_n)
I_rel = abs(I_1_n_2/I_1_n_1)
print(f'I_1,n,1 = {Complex(I_1_n_1) :z.2f} A')
print(f'I_1,n,2 = {Complex(I_1_n_2) :z.2f} A')
print(f'I_1,n,2/I_1,n,1 = {I_rel:.3f}\n')

I_1_n_A, I_1_n_B, I_1_n_C = A012_to_ABC(0, I_1_n_1, I_1_n_2)
print(f'I_1,n,A = {Complex(I_1_n_A) :z.2f} A')
print(f'I_1,n,B = {Complex(I_1_n_B) :z.2f} A')
print(f'I_1,n,C = {Complex(I_1_n_C) :z.2f} A')
