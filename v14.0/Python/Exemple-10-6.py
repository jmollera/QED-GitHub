import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt

R_1 = 0.5  # Ω/fase
X_1 = 1.5  # Ω/fase
R_2 = 0.625  # Ω/fase
X_2 = 1.25 # Ω/fase
R_Fe = 360  # Ω/fase
X_m = 40  # Ω/fase
s_n = 0.05  # Lliscament nominal
U = 380  # V (fase-fase)

U_1 = U/np.sqrt(3)
print(f'U_1 = {U_1:.3f} V\n')

Z_0 = R_Fe*1j*X_m/(R_Fe + 1j*X_m)
Z_Th = (R_1 + 1j*X_1)*Z_0/(R_1 + 1j*X_1 + Z_0)
print(f'Z_0 = {Z_0:.3f} Ω')
print(f'Z_Th = {Z_Th:.3f} Ω\n')

def Z_mot(s):
    return R_1 + 1j*X_1 + Z_0*(R_2/s + 1j*X_2)/(Z_0 + R_2/s + 1j*X_2)

def I_1(s, U1):
    return U1/Z_mot(s)

def I_2(s, U1):
    UTh = U1*Z_0/(R_1 + 1j*X_1 + Z_0)
    return UTh/(Z_Th + R_2/s + 1j*X_2)

def P_m(s, U1):
    return 3*(1-s)*R_2*abs(I_2(s, U1))**2/s

def equació(s, Pm, U1):
    return P_m(s, U1) - Pm

P_m_n = P_m(s_n, U_1)
P_m_50 = 0.5*P_m_n

sol = optimize.root(equació, x0=s_n, args=(P_m_50, U_1))
s_50 = sol.x[0]
print(f'P_m,50% = {P_m_50:.1f} W')
print(f's,50% = {s_50:.6f}\n')

Z_mot_50 = Z_mot(s_50)
I_1_50 = abs(I_1(s_50, U_1))
I_1_n = abs(I_1(s_n, U_1))
cos_φ_50= np.cos(np.angle(Z_mot(s_50)))
η_50 = (1-s_50)/s_50*R_2*abs(I_2(s_50,U_1))**2/U_1/I_1_50/cos_φ_50
print(f'Z_mot,50% = {Z_mot_50:.3f} Ω\n')
print(f'I_1,50% = {I_1_50:.3f} A')
print(f'I_1,50%/I_1,n = {I_1_50/I_1_n:.1f}\n')
print(f'cos φ,50% = {cos_φ_50:.6f}')
print(f'η,50% = {η_50:.6f}\n')

Pm_pu = np.linspace(0.001, 1, 500) # 0.1% a 100% Pm,n
s = np.empty_like(Pm_pu)
for i, pu in enumerate(Pm_pu):
    sol = optimize.root(equació, x0=s_n, args=(pu*P_m_n, U_1))
    s[i] = sol.x[0]
I1 = abs(I_1(s, U_1))
cos_φ = np.cos(np.angle(Z_mot(s)))
η = (1-s)/s*R_2*abs(I_2(s, U_1))**2/U_1/I1/cos_φ

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
ax2.plot(Pm_pu, I1, 'green', label=r'$I_1$')
fig.legend(loc='lower right', bbox_to_anchor=(1,0), bbox_transform=ax1.transAxes)
plt.show()

fig, ax = plt.subplots(figsize=(6, 4))
ax.set_xlim(0, 1)
ax.set_xticks(np.arange(0, 1.01, 0.1))
ax.set_ylim(0, 1)
ax.set_yticks(np.arange(0, 1.01, 0.1))
ax.grid()
ax.set_ylabel(r'$\cos\varphi$,   $\eta$')
ax.set_xlabel(r'$P_{\mathrm{m}}/P_{\mathrm{m,n}}$')
ax.plot(Pm_pu, cos_φ, 'red', label=r'$\cos\varphi$')
ax.plot(Pm_pu, η, 'green', label=r'$\eta$')
ax.legend(loc='lower right')
plt.show()