import numpy as np
from qed.eng_elec import Motor3ph, A012_to_ABC
from qed.utils import Complex

R_1 = 0.5  # Ω/fase
X_1 = 1.5  # Ω/fase
R_2 = 0.625  # Ω/fase
X_2 = 1.25  # Ω/fase
R_Fe = 360  # Ω/fase
X_m = 40  # Ω/fase
p = 4  # Nombre de pols
s_n = 0.05  # Lliscament nominal
U = 380  # V (fase-fase)
f = 50  # Hz

U_1 = U / np.sqrt(3)

motor = Motor3ph(p, f, R_1, X_1, R_2, X_2, R_Fe, X_m)
mot_n = motor.run(s_n, U_1)

print(f's_n = {s_n}')
print(f'U_1 = {U_1:.3f} V')
print(f'ω_m,sinc = {motor.ω_m_sync:.3f} rad/s\n')

I_1_1 = U_1 / (2*(R_1 + R_2/(s_n*(2-s_n)) + 1j*(X_1+X_2)))
I_1_2 = -I_1_1
I_1_0 = 0
print(f'I_1_1 = {Complex(I_1_1) :.1f/.2f}')
print(f'I_1_2 = {Complex(I_1_2) :.1f/.2f}')
print(f'I_1_0 = {Complex(I_1_0) :.1f/.2f}\n')

I_1_A, I_1_B, I_1_C = A012_to_ABC(I_1_0, I_1_1, I_1_2)
print(f'I_1_A = {Complex(I_1_A) :.1f/.2f}')
print(f'I_1_B = {Complex(I_1_B) :.1f/.2f}')
print(f'I_1_C = {Complex(I_1_C) :.1f/.2f}\n')

P_m = 3 * R_2 * ((1-s_n)/s_n - (1-s_n)/(2-s_n)) * abs(I_1_1)**2
T_m = P_m / ((1-s_n) * motor.ω_m_sync)
print(f'P_m = {P_m:.1f} W')
print(f'T_m = {T_m:.1f} N·m\n')

print(f'I_1_B/I_1_n = {abs(I_1_B/mot_n.I):.2f}')
print(f'P_m/P_m_n = {P_m/mot_n.P_m:.2f}')
print(f'T_m/T_m_n = {T_m/mot_n.T_m:.2f}')
