import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt

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
U_1_80 = 0.8*U_1
print(f'U_1,100% = {U_1:.3f} V')
print(f'U_1,80% = {U_1_80:.3f} V\n')

n_m_sinc = 120*f/p
ω_m_sinc = 4*np.pi*f/p
print(f'n_m,sinc = {n_m_sinc:.0f} r/min,  ω_m,sinc = {ω_m_sinc:.3f} rad/s\n')

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

def T_m(s, U1):
    return 3*R_2*abs(I_2(s, U1))**2/(s*ω_m_sinc)

def T_load(s):
    return 3.2*(1-s) + 58.9*(1-s)**2

def equació(s, U1):
    return T_m(s, U1) - T_load(s)

sol = optimize.root(equació, x0=s_n, args=(U_1))
s_100 = sol.x[0]
n_100 = (1 - s_100)*n_m_sinc
zmot_100 = Z_mot(s_100)
I_1_100 = abs(I_1(s_100, U_1))
I_1_arr_100 = abs(I_1(1, U_1))
Tm_100 = T_m(s_100, U_1)
Tm_arr_100 = T_m(1.0001, U_1)
cos_φ_100 = np.cos(np.angle(Z_mot(s_100)))
η_100 = (1-s_100)/s_100*R_2*abs(I_2(s_100, U_1))**2/U_1/I_1_100/cos_φ_100
print(f's,100% = {s_100:.3f}')
print(f'n,100% = {n_100:.1f} r/min')
print(f'Z_mot,100% = {zmot_100:.3f} Ω')
print(f'I_1,100% = {I_1_100:.1f} A')
print(f'I_1,arr,100% = {I_1_arr_100:.1f} A')
print(f'Tm,100% = {Tm_100:.1f} N·m')
print(f'Tm,arr,100% = {Tm_arr_100:.1f} N·m')
print(f'cos φ,100% = {cos_φ_100:.3f}')
print(f'η,100% = {η_100:.3f}\n')

sol = optimize.root(equació, x0=s_n, args=(U_1_80))
s_80 = sol.x[0]
n_80 = (1 - s_80)*n_m_sinc
zmot_80 = Z_mot(s_80)
I_1_80 = abs(I_1(s_80, U_1_80))
I_1_arr_80 = abs(I_1(1, U_1_80))
Tm_80 = T_m(s_80, U_1_80)
Tm_arr_80 = T_m(1.0001, U_1_80)
cos_φ_80 = np.cos(np.angle(Z_mot(s_80)))
η_80 = (1-s_80)/s_80*R_2*abs(I_2(s_80, U_1_80))**2/U_1_80/I_1_80/cos_φ_80
print(f's,80% = {s_80:.3f}')
print(f'n,80% = {n_80:.1f} r/min')
print(f'Z_mot,80% = {zmot_80:.3f} Ω')
print(f'I_1,80% = {I_1_80:.1f} A')
print(f'I_1,arr,80% = {I_1_arr_80:.1f} A')
print(f'Tm,80% = {Tm_80:.1f} N·m')
print(f'Tm,arr,80% = {Tm_arr_80:.1f} N·m')
print(f'cos φ,80% = {cos_φ_80:.3f}')
print(f'η,80% = {η_80:.3f}\n')

n = np.linspace(0.001, n_m_sinc-0.001, 500)
s = 1 - n/n_m_sinc
Tload = T_load(s)
Tm_100 = T_m(s, U_1)
Tm_80 = T_m(s, U_1_80)
I1_100 = abs(I_1(s, U_1))
I1_80 = abs(I_1(s, U_1_80))

fig, ax1 = plt.subplots(figsize=(6, 4))
plt.rc('xtick', labelsize=8)
ax1.set_xlim(0, 1500)
ax1.set_xticks(np.arange(0, 1501, 100))
ax1.set_ylim(0, 150)
ax1.set_yticks(np.arange(0, 151, 10))
ax1.grid()
ax1.set_ylabel(r'$T_{\mathrm{m}}$ / N m,   $T_{\mathrm{load}}$ / N m')
ax1.set_xlabel(r'$n_{\mathrm{m}}$ / r/min')
ax1.plot(n, Tload, 'black', label=r'$T_{\mathrm{load}}$')
ax1.plot(n, Tm_100, 'red', label=r'$T_{\mathrm{m},100 \%}$')
ax1.plot(n, Tm_80, 'green', label=r'$T_{\mathrm{m},80 \%}$')
ax2 = ax1.twinx()
ax2.set_ylim(0, 75)
ax2.set_yticks(np.arange(0, 76, 5))
ax2.set_ylabel(r'$I_1$ / A')
ax2.plot(n, I1_100, 'orange', label=r'$I_{1,100 \%}$')
ax2.plot(n, I1_80, 'lime', label=r'$I_{1,80 \%}$')
fig.legend(loc='lower left', bbox_to_anchor=(0,0), bbox_transform=ax1.transAxes)
plt.show()

fig, ax1 = plt.subplots(figsize=(6, 4))
plt.rc('xtick', labelsize=8)
ax1.set_xlim(1200, 1500)
ax1.set_xticks(np.arange(1200, 1501, 50))
ax1.set_ylim(0, 90)
ax1.set_yticks(np.arange(0, 91, 10))
ax1.grid()
ax1.set_ylabel(r'$T_{\mathrm{m}}$ / N m,   $T_{\mathrm{load}}$ / N m')
ax1.set_xlabel(r'$n_{\mathrm{m}}$ / r/min')
ax1.plot(n, Tload, 'black', label=r'$T_{\mathrm{load}}$')
ax1.plot(n, Tm_100, 'red', label=r'$T_{\mathrm{m},100 \%}$')
ax1.plot(n, Tm_80, 'green', label=r'$T_{\mathrm{m},80 \%}$')
ax2 = ax1.twinx()
ax2.set_ylim(0, 45)
ax2.set_yticks(np.arange(0, 46, 5))
ax2.set_ylabel(r'$I_1$ / A')
ax2.plot(n, I1_100, 'orange', label=r'$I_{1,100 \%}$')
ax2.plot(n, I1_80, 'lime', label=r'$I_{1,80 \%}$')
fig.legend(loc='lower left', bbox_to_anchor=(0,0), bbox_transform=ax1.transAxes)
plt.show()

U1_pu = np.linspace(0.75, 1, 500) # 75% a 100% U_1,n
s = np.empty_like(U1_pu)
for i, pu in enumerate(U1_pu):
    sol = optimize.root(equació, x0=s_n, args=(pu*U_1))
    s[i] = sol.x[0]
I1 = abs(I_1(s, U1_pu*U_1))
cos_φ = np.cos(np.angle(Z_mot(s)))
η = (1-s)/s*R_2*abs(I_2(s, U1_pu*U_1))**2/(U1_pu*U_1)/I1/cos_φ

fig, ax1 = plt.subplots(figsize=(6, 4))
ax1.set_xlim(0.75, 1.0)
ax1.set_xticks(np.arange(0.75, 1.01, 0.05))
ax1.set_ylim(0.04, 0.09)
ax1.set_yticks(np.arange(0.04, 0.091, 0.01))
ax1.grid()
ax1.set_ylabel(r'$s$')
ax1.set_xlabel(r'$U_1/U_{\mathrm{n}}$')
ax1.plot(U1_pu, s, 'red', label=r'$s$')
ax2 = ax1.twinx()
ax2.set_ylim(16, 21)
ax2.set_yticks(np.arange(16, 22, 1))
ax2.set_ylabel(r'$I_1$ / A')
ax2.plot(U1_pu, I1, 'green', label=r'$I_1$')
fig.legend(loc='upper right', bbox_to_anchor=(1,1), bbox_transform=ax1.transAxes)
plt.show()

fig, ax = plt.subplots(figsize=(6, 4))
ax.set_xlim(0.75, 1.0)
ax.set_xticks(np.arange(0.75, 1.01, 0.05))
ax.set_ylim(0.82, 0.9)
ax.set_yticks(np.arange(0.82, 0.91, 0.02))
ax.grid()
ax.set_ylabel(r'$\cos\varphi$,   $\eta$')
ax.set_xlabel(r'$U_1/U_{\mathrm{n}}$')
ax.plot(U1_pu, cos_φ, 'red', label=r'$\cos\varphi$')
ax.plot(U1_pu, η, 'green', label=r'$\eta$')
ax.legend(loc='lower right')
plt.show()