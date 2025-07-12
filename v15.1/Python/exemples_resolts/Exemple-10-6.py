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
print(f'U_1 = {U_1:.3f} V\n')

motor = Motor3ph(p, f, R_1, X_1, R_2, X_2, R_Fe, X_m)
mot_n = motor.run(s_n, U_1)

P_m_n = mot_n.P_m
P_m_50 = 0.5*P_m_n
s_P_m_50 = find_s(motor, U_1, P_m_50, 0.5*s_n)
mot_s_P_m_50 = motor.run(s_P_m_50, U_1)

print(f'P_m,50% = {P_m_50:.1f} W')
print(f's,50% = {s_P_m_50:.4f}\n')
print(f'Z_mot,50% = {mot_s_P_m_50.Z:.3f} Ω\n')
print(f'I_1,50% = {abs(mot_s_P_m_50.I):.1f} A')
print(f'I_1,50%/I_1,n = {abs(mot_s_P_m_50.I / mot_n.I):.1f}\n')
print(f'cos 𝜑,50% = {mot_s_P_m_50.PF:.4f}')
print(f'η,50% = {mot_s_P_m_50.Eff:.4f}\n')

pu = np.linspace(0, 1, 500)
s_P_m_x = np.vectorize(find_s)(motor, U_1, pu*P_m_n, pu*s_n)
mot_s_P_m_x = motor.run(s_P_m_x, U_1)

fig, ax1 = plt.subplots(figsize=(6, 4))
ax1.set_xlim(0, 1)
xticks = np.arange(0, 1.01, 0.1)
ax1.set_xticks(xticks, [f'{tick:.1f}'.replace('.', ',') for tick in xticks])
ax1.set_ylim(0, 0.05)
yticks1 = np.arange(0, 0.051, 0.005)
ax1.set_yticks(yticks1, [f'{tick:.3f}'.replace('.', ',') for tick in yticks1])
ax1.grid()
ax1.set_ylabel(r'$s$')
ax1.set_xlabel(r'$P_{\mathrm{m}}/P_{\mathrm{m,n}}$')
ax1.plot(pu, s_P_m_x, 'red', label=r'$s$')
ax2 = ax1.twinx()
ax2.set_ylim(0, 20)
yticks2 = np.arange(0, 21, 2)
ax2.set_yticks(yticks2, [f'{tick:2.0f}' for tick in yticks2])
ax2.set_ylabel(r'$I_1$ / A')
ax2.plot(pu, abs(mot_s_P_m_x.I), 'green', label=r'$I_1$')
fig.legend(loc='lower right', bbox_to_anchor=(1, 0), bbox_transform=ax1.transAxes)
plt.show()

fig, ax = plt.subplots(figsize=(6, 4))
ax.set_xlim(0, 1)
xticks = np.arange(0, 1.01, 0.1)
ax.set_xticks(xticks, [f'{tick:.1f}'.replace('.', ',') for tick in xticks])
ax.set_ylim(0, 1)
yticks = np.arange(0, 1.01, 0.1)
ax.set_yticks(yticks, [f'{tick:.1f}'.replace('.', ',') for tick in xticks])
ax.grid()
ax.set_ylabel(r'$\cos\varphi$,   $\eta$')
ax.set_xlabel(r'$P_{\mathrm{m}}/P_{\mathrm{m,n}}$')
ax.plot(pu, mot_s_P_m_x.PF, 'red', label=r'$\cos\varphi$')
ax.plot(pu, mot_s_P_m_x.Eff, 'green', label=r'$\eta$')
ax.legend(loc='lower right')
plt.show()
