import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt

R_1 = 0.5  # Ω/fase
X_1 = 1.5  # Ω/fase
R_2 = 0.625  # Ω/fase
X_2 = 1.25  # Ω/fase
R_Fe = 360  # Ω/fase
X_m = 40  # Ω/fase
s_n = 0.05  # Lliscament nominal
U = 380  # V (fase-fase)

U_1 = U/np.sqrt(3)

Z_0 = R_Fe*1j*X_m/(R_Fe + 1j*X_m)
Z_Th = (R_1 + 1j*X_1)*Z_0/(R_1 + 1j*X_1 + Z_0)
print(f'Z_0 = {Z_0:.3f} Ω')
print(f'Z_Th = {Z_Th:.3f} Ω\n')

def Z_mot(s):
    return R_1 + 1j*X_1 + Z_0*(R_2/s + 1j*X_2)/(Z_0 + R_2/s + 1j*X_2)

def I_1(s, U1):
    return U1/Z_mot(s)

def I_2(s, U1):
    ETh = U1*Z_0/(R_1 + 1j*X_1 + Z_0)
    return ETh/(Z_Th + R_2/s + 1j*X_2)

def P_m(s, U1):
    return 3*(1-s)*R_2*abs(I_2(s, U1))**2/s

def equació(s, Pm, U1):
    return P_m(s, U1) - Pm

P_m_n = P_m(s_n, U_1)
print(f'P_m,n = {P_m_n:.1f} W\n')

U_1_80 = 0.8*U_1
print(f'U_1,80% = {U_1_80:.3f} V')

sol = optimize.root(equació, x0=s_n, args=(P_m_n, U_1_80))
s_80 = sol.x[0]
print(f's,80% = {s_80:.6f}\n')

Z_mot_80 = Z_mot(s_80)
I_1_80 = abs(I_1(s_80, U_1_80))
I_1_n = abs(I_1(s_n, U_1))
cos_φ_80= np.cos(np.angle(Z_mot(s_80)))
η_80 = (1-s_80)/s_80*R_2*abs(I_2(s_80,U_1_80))**2/U_1_80/I_1_80/cos_φ_80
print(f'Z_mot,80% = {Z_mot_80:.3f} Ω\n')
print(f'I_1,80% = {I_1_80:.3f} A')
print(f'I_1,80%/I_1,n = {I_1_80/I_1_n:.1f}\n')
print(f'cos φ,80% = {cos_φ_80:.6f}')
print(f'η,80% = {η_80:.6f}\n')

U1_pu = np.linspace(0.75, 1, 500)  # 75% a 100% U_1,n
s = np.empty_like(U1_pu)
for i, pu in enumerate(U1_pu):
    sol = optimize.root(equació, x0=s_n, args=(P_m_n, pu*U_1))
    s[i] = sol.x[0]
I1 = abs(I_1(s, U1_pu*U_1))
cos_φ = np.cos(np.angle(Z_mot(s)))
η = (1-s)/s*R_2*abs(I_2(s, U1_pu*U_1))**2/(U1_pu*U_1)/I1/cos_φ

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
ax2.plot(U1_pu, I1, 'green', label=r'$I_1$')
fig.legend(loc='upper right', bbox_to_anchor=(1,1), bbox_transform=ax1.transAxes)
plt.show()

fig, ax = plt.subplots(figsize=(6, 4))
ax.set_xlim(0.75, 1.0)
ax.set_xticks(np.arange(0.75, 1.01, 0.05))
ax.set_ylim(0.76, 0.90)
ax.set_yticks(np.arange(0.76, 0.91, 0.02))
ax.grid()
ax.set_ylabel(r'$\cos\varphi$,   $\eta$')
ax.set_xlabel(r'$U_1/U_{\mathrm{n}}$')
ax.plot(U1_pu, cos_φ, 'red', label=r'$\cos\varphi$')
ax.plot(U1_pu, η, 'green', label=r'$\eta$')
ax.legend(loc='lower right')
plt.show()