import numpy as np
from qed.utils import Complex

S_n = 400_000  # VA
U_n1 = 25_000  # V
U_n2 = 400  # V
ε_cc = 0.04
W_cc = 4000  # W
i_0 = 0.02
W_0 = 2000  # W
P_2 = 200_000  # W
U_2 = 380  # V
cos_𝜑 = 0.8

S_B = S_n
U_B1 = U_n1
U_B2 = U_n2

g_Fe = w_0 = W_0/S_B
b_m = np.sqrt(i_0**2 - w_0**2)
r = w_cc = W_cc/S_B
x = np.sqrt(ε_cc**2 - w_cc**2)
print(f'g_Fe = {g_Fe:.4f}')
print(f'b_m = {b_m:.4f}')
print(f'r = {r:.4f}')
print(f'x = {x:.4f}\n')
p_2 = P_2/S_B
s_2 = p_2 + 1j * p_2 * np.sqrt(1 - cos_𝜑 ** 2) / cos_𝜑
u_2 = U_2/U_B2
print(f's_2 = {s_2:.4f}')
print(f'u_2 = {u_2:.4f}\n')
i_2 = s_2.conjugate()/u_2.conjugate()
i_0 = u_2*(g_Fe - 1j*b_m)
i_1 = i_2 + i_0
print(f'i_2 = {Complex(i_2) :.4f/.2f}')
print(f'i_0 = {Complex(i_0) :.4f/.2f}')
print(f'i_1 = {Complex(i_1) :.4f/.2f}\n')
u_1 = (r + 1j*x)*i_1 + u_2
U_1 = u_1*U_B1
print(f'u_1 = {Complex(u_1) :.4f/.2f}')
print(f'U_1 = {Complex(U_1) :.1f/.2f} V\n')
Δ_u = abs(u_1) - abs(u_2)
Δ_U2 = Δ_u*U_B2
print(f'Δ_u = {Δ_u:.4f}')
print(f'Δ_U2 = {Δ_U2:.4f} V\n')
p_Cu = r*abs(i_1)**2
p_Fe = g_Fe*abs(u_2)**2
η = p_2/(p_2 + p_Cu + p_Fe)
print(f'p_Cu = {p_Cu:.6f}')
print(f'p_Fe = {p_Fe:.6f}')
print(f'η = {η:.2f}')
