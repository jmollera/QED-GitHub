from qed.utils import Complex
from qed.eng_elec import triangle_to_phasors, ABC_to_A012, LL_to_LN

R_AN = 5   # Ω
R_BN = 10  # Ω
R_CN = 15  # Ω

U_AB, U_BC, U_CA = triangle_to_phasors(2760, 1840, 2300)
print(f'U_AB = {Complex(U_AB) :.2f} V')
print(f'U_BC = {Complex(U_BC) :.2f} V')
print(f'U_CA = {Complex(U_CA) :.2f} V\n')

U_AN, U_BN, U_CN = LL_to_LN(U_AB, U_BC, U_CA, R_AN, R_BN, R_CN)
print(f'U_AN = {Complex(U_AN) :.2f} V')
print(f'U_BN = {Complex(U_BN) :.2f} V')
print(f'U_CN = {Complex(U_CN) :.2f} V\n')

U_AN_0, U_AN_1, U_AN_2 = ABC_to_A012(U_AN, U_BN, U_CN)
print(f'U_AN_1 = {Complex(U_AN_1) :.2f} V')
print(f'U_AN_2 = {Complex(U_AN_2) :.2f} V')
print(f'U_AN_0 = {Complex(U_AN_0) :.2f} V\n')

U_AG = (U_AB - U_CA)/3
U_GN = -U_AG + U_AN
print(f'U_GN = {Complex(U_GN) :.2f} V')
