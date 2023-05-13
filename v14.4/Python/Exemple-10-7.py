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

mot = Motor3ph(p, f, R_1, X_1, R_2, X_2, R_Fe, X_m)
mot_n = mot.run(s_n, U_1)
P_m_n = mot_n.P_m

print(f'P_m,n = {P_m_n:.1f} W')

U_1_80 = 0.8*U_1
print(f'U_1,80% = {U_1_80:.3f} V')

sol = optimize.root(lambda s: mot.run(s, U_1_80).P_m - P_m_n, x0=s_n)
s_80 = sol.x[0]
mot_80 = mot.run(s_80, U_1_80)

print(f's,80% = {s_80:.4f}\n')
print(f'Z_mot,80% = {mot_80.Z:.3f} Ω\n')
print(f'I_1,80% = {abs(mot_80.I):.1f} A')
print(f'I_1,80%/I_1,n = {abs(mot_80.I / mot_n.I):.1f}\n')
print(f'cos 𝜑,80% = {mot_80.PF:.4f}')
print(f'η,80% = {mot_80.Eff:.4f}\n')

U1_pu = np.linspace(0.75, 1, 500)  # (75 % a 100 %) U_1,n
s = np.empty_like(U1_pu)
for i, pu in enumerate(U1_pu):
    sol = optimize.root(lambda s: mot.run(s, pu*U_1).P_m - P_m_n, x0=s_n)
    s[i] = sol.x[0]
mot_s = mot.run(s, U1_pu * U_1)

fig, ax1 = plt.subplots(figsize=(6, 4))
ax1.set_xlim(0.75, 1.0)
ax1.set_xticks(np.arange(0.75, 1.01, 0.05))
ax1.set_ylim(0, 0.14)
ax1.set_yticks(np.arange(0, 0.15, 0.02))
ax1.grid()
ax1.set_ylabel(r'$s$')
ax1.set_xlabel(r'$U_1/U_{\mathrm{n}}$')
ax1.plot(U1_pu, s, 'red', label=r'$s$')
ax2 = ax1.twinx()
ax2.set_ylim(16, 30)
ax2.set_yticks(np.arange(16, 31, 2))
ax2.set_ylabel(r'$I_1$ / A')
ax2.plot(U1_pu, abs(mot_s.I), 'green', label=r'$I_1$')
fig.legend(loc='upper right', bbox_to_anchor=(1, 1), bbox_transform=ax1.transAxes)
plt.show()

fig, ax = plt.subplots(figsize=(6, 4))
ax.set_xlim(0.75, 1.0)
ax.set_xticks(np.arange(0.75, 1.01, 0.05))
ax.set_ylim(0.76, 0.90)
ax.set_yticks(np.arange(0.76, 0.91, 0.02))
ax.grid()
ax.set_ylabel(r'$\cos\varphi$,   $\eta$')
ax.set_xlabel(r'$U_1/U_{\mathrm{n}}$')
ax.plot(U1_pu, mot_s.PF, 'red', label=r'$\cos\varphi$')
ax.plot(U1_pu, mot_s.Eff, 'green', label=r'$\eta$')
ax.legend(loc='lower right')
plt.show()
