import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
from qed.eng_elec import Motor3ph

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

U_1 = U/np.sqrt(3)
print(f'U_1 = {U_1:.3f} V\n')

mot = Motor3ph(p, f, R_1, X_1, R_2, X_2, R_Fe, X_m)
mot_n = mot.run(s_n, U_1)

P_m_n = mot_n.P_m
P_m_50 = 0.5 * P_m_n

sol = optimize.root(lambda s: mot.run(s, U_1).P_m - P_m_50, x0=0.5*s_n)
s_50 = sol.x[0]
mot_50 = mot.run(s_50, U_1)
print(f'P_m,50% = {P_m_50:.1f} W')
print(f's,50% = {s_50:.4f}\n')
print(f'Z_mot,50% = {mot_50.Z:.3f} Ω\n')
print(f'I_1,50% = {abs(mot_50.I):.1f} A')
print(f'I_1,50%/I_1,n = {abs(mot_50.I / mot_n.I):.1f}\n')
print(f'cos 𝜑,50% = {mot_50.PF:.4f}')
print(f'η,50% = {mot_50.Eff:.4f}\n')

Pm_pu = np.linspace(0, 1, 500)  # (0 % a 100 %) Pm,n
s = np.empty_like(Pm_pu)
for i, pu in enumerate(Pm_pu):
    sol = optimize.root(lambda s: mot.run(s, U_1).P_m - pu*P_m_n, x0=pu*s_n)
    s[i] = sol.x[0]
mot_s = mot.run(s, U_1)

fig, ax1 = plt.subplots(figsize=(6, 4))
ax1.set_xlim(0, 1)
ax1.set_xticks(np.arange(0, 1.01, 0.1))
ax1.set_ylim(0, 0.05)
ax1.set_yticks(np.arange(0, 0.051, 0.005))
ax1.grid()
ax1.set_ylabel(r'$s$')
ax1.set_xlabel(r'$P_{\mathrm{m}}/P_{\mathrm{m,n}}$')
ax1.plot(Pm_pu, s, 'red', label=r'$s$')
ax2 = ax1.twinx()
ax2.set_ylim(0, 20)
ax2.set_yticks(np.arange(0, 21, 2))
ax2.set_ylabel(r'$I_1$ / A')
ax2.plot(Pm_pu, abs(mot_s.I), 'green', label=r'$I_1$')
fig.legend(loc='lower right', bbox_to_anchor=(1, 0), bbox_transform=ax1.transAxes)
plt.show()

fig, ax = plt.subplots(figsize=(6, 4))
ax.set_xlim(0, 1)
ax.set_xticks(np.arange(0, 1.01, 0.1))
ax.set_ylim(0, 1)
ax.set_yticks(np.arange(0, 1.01, 0.1))
ax.grid()
ax.set_ylabel(r'$\cos\varphi$,   $\eta$')
ax.set_xlabel(r'$P_{\mathrm{m}}/P_{\mathrm{m,n}}$')
ax.plot(Pm_pu, mot_s.PF, 'red', label=r'$\cos\varphi$')
ax.plot(Pm_pu, mot_s.Eff, 'green', label=r'$\eta$')
ax.legend(loc='lower right')
plt.show()
