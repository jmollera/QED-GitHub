import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
from qed.eng_elec import Motor3ph, SQRT3

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


def find_s(motor: Motor3ph, U_1: float, P_m: float, s_guess: float):
    s_found = optimize.root_scalar(lambda s: motor.run(s, U_1).P_m - P_m, x0=s_guess)
    return s_found.root


U_1 = U / SQRT3

motor = Motor3ph(p, f, R_1, X_1, R_2, X_2, R_Fe, X_m)
mot_n = motor.run(s_n, U_1)
P_m_n = mot_n.P_m

print(f'P_m,n = {P_m_n:.1f} W')

U_1_80 = 0.8*U_1
print(f'U_1,80% = {U_1_80:.3f} V')

s_U_1_80 = find_s(motor, U_1_80, P_m_n, s_n)
mot_U_1_80 = motor.run(s_U_1_80, U_1_80)

print(f's,80% = {s_U_1_80:.4f}\n')
print(f'Z_mot,80% = {mot_U_1_80.Z:.3f} Ω\n')
print(f'I_1,80% = {abs(mot_U_1_80.I):.1f} A')
print(f'I_1,80%/I_1,n = {abs(mot_U_1_80.I / mot_n.I):.1f}\n')
print(f'cos 𝜑,80% = {mot_U_1_80.PF:.4f}')
print(f'η,80% = {mot_U_1_80.Eff:.4f}\n')

pu = np.linspace(0.75, 1, 500)
s_U_1_x = np.vectorize(find_s)(motor, pu*U_1, P_m_n, s_n)
mot_s_U_1_x = motor.run(s_U_1_x , pu*U_1)

fig, ax1 = plt.subplots(figsize=(6, 4))
ax1.set_xlim(0.75, 1.0)
xticks = np.arange(0.75, 1.01, 0.05)
ax1.set_xticks(xticks, [f'{tick:.2f}'.replace('.', ',') for tick in xticks])
ax1.set_ylim(0, 0.14)
yticks1 = np.arange(0, 0.15, 0.02)
ax1.set_yticks(yticks1, [f'{tick:.2f}'.replace('.', ',') for tick in yticks1])
ax1.grid()
ax1.set_ylabel(r'$s$')
ax1.set_xlabel(r'$U_1/U_{\mathrm{n}}$')
ax1.plot(pu, s_U_1_x, 'red', label=r'$s$')
ax2 = ax1.twinx()
ax2.set_ylim(16, 30)
ax2.set_yticks(np.arange(16, 31, 2))
ax2.set_ylabel(r'$I_1$ / A')
ax2.plot(pu, abs(mot_s_U_1_x.I), 'green', label=r'$I_1$')
fig.legend(loc='upper right', bbox_to_anchor=(1, 1), bbox_transform=ax1.transAxes)
plt.show()

fig, ax = plt.subplots(figsize=(6, 4))
ax.set_xlim(0.75, 1.0)
xticks = np.arange(0.75, 1.01, 0.05)
ax.set_xticks(xticks, [f'{tick:.2f}'.replace('.', ',') for tick in xticks])
ax.set_ylim(0.76, 0.90)
yticks = np.arange(0.76, 0.91, 0.02)
ax.set_yticks(yticks, [f'{tick:.2f}'.replace('.', ',') for tick in yticks])
ax.grid()
ax.set_ylabel(r'$\cos\varphi$,   $\eta$')
ax.set_xlabel(r'$U_1/U_{\mathrm{n}}$')
ax.plot(pu, mot_s_U_1_x.PF, 'red', label=r'$\cos\varphi$')
ax.plot(pu, mot_s_U_1_x.Eff, 'green', label=r'$\eta$')
ax.legend(loc='lower right')
plt.show()
