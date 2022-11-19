import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
from qed.utils import ComplexD

R_1 = 0.5  # Ω/fase
X_1 = 1.5  # Ω/fase
R_2 = 0.625  # Ω/fase
X_2 = 1.25  # Ω/fase
R_Fe = 360  # Ω/fase
X_m = 40  # Ω/fase
p = 4  # Nombre de pols
s_n = 0.05  # Lliscament nominal
f = 50  # Hz
U_sist = 380  # V (fase-fase)
Scc_sist = 5  # MVA
Rel_X_R = 9  # Relació X/R
L_cable = 80  # m
R_cable = 3.960  # Ω/km/fase
X_cable = 0.123  # Ω/km/fase

U_1_sist = U_sist/np.sqrt(3)
fase = np.arctan(Rel_X_R)
Z_sist = U_sist**2/(Scc_sist*1e6)*(np.cos(fase)+1j*np.sin(fase))
Z_cable = L_cable*1e-3*(R_cable+1j*X_cable)
print(f'U_1,sist = {U_1_sist:.3f} V')
print(f'Z_sist = {Z_sist:.3f} Ω')
print(f'Z_cable = {Z_cable:.3f} Ω\n')

n_m_sinc = 120*f/p
ω_m_sinc=4*np.pi*f/p
print(f'n_m,sinc = {n_m_sinc:.0f} r/min,  ω_m,sinc = {ω_m_sinc:.3f} rad/s\n')

Z_0 = R_Fe*1j*X_m/(R_Fe + 1j*X_m)
Z_Th = (Z_sist + Z_cable+R_1 + 1j*X_1)*Z_0/(Z_sist+Z_cable+R_1 + 1j*X_1 + Z_0)
E_Th = U_1_sist*Z_0/(Z_sist + Z_cable+R_1 + 1j*X_1 + Z_0)
print(f'Z_0 = {Z_0:.3f} Ω')
print(f'Z_Th = {Z_Th:.3f} Ω')
print(f'E_Th = {ComplexD(E_Th):.4f} V\n')

def Z_mot(s):
    return R_1 + 1j*X_1 + Z_0*(R_2/s + 1j*X_2)/(Z_0 + R_2/s + 1j*X_2)

def U_1(s):
    return U_1_sist*Z_mot(s)/(Z_sist+Z_cable+Z_mot(s))

def I_1(s):
    return U_1(s)/Z_mot(s)

def I_2(s):
    return E_Th/(Z_Th + R_2/s + 1j*X_2)

def T_m(s):
    return 3*R_2*abs(I_2(s))**2/(s*ω_m_sinc)

def T_load(s):
     return 3.2*(1-s) + 58.9*(1-s)**2

def equació(s):
     return T_m(s) - T_load(s)

sol = optimize.root(equació, x0=s_n)
s = sol.x[0]
U_1_n = U_1(s)
I_1_n = abs(I_1(s))
T_m_n = T_m(s)
print(f's,n = {s:.3f}')
print(f'U_1,n = {ComplexD(U_1_n):.4f} V')
print(f'I_1,n = {I_1_n:.2f} A')
print(f'T_m,n = {T_m_n:.1f} N·m\n')

U_1_arr = U_1(1)
I_1_arr = abs(I_1(1))
T_m_arr = T_m(1)
print(f'U_1,arr = {ComplexD(U_1_arr):.4f} V')
print(f'I_1,arr = {I_1_arr:.2f} A')
print(f'T_m,arr = {T_m_arr:.1f} N·m\n')

s_T_m_max = R_2/abs(Z_Th+1j*X_2)
n_T_m_max = (1-s_T_m_max)*n_m_sinc
ω_T_m_max = (1-s_T_m_max)*ω_m_sinc
T_m_max = T_m(s=s_T_m_max)
print(f's_T_m,max = {s_T_m_max:.4f}')
print(f'n_T_m,max = {n_T_m_max:.2f} r/min,  ω_T_m,max = {ω_T_m_max:.2f} rad/s')
print(f'T_m,max = {T_m_max:.1f} N·m')

n = np.linspace(0, n_m_sinc-0.01, 500)  # r/min
Tload = T_load(s=1-n/n_m_sinc)  # N·m
Tm = T_m(s=1-n/n_m_sinc)  # N·m
I1 = abs(I_1(s=1-n/n_m_sinc))  # A
U1 = abs(U_1(s=1-n/n_m_sinc))  # V

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
ax1.plot(n, Tm, 'red', label=r'$T_{\mathrm{m}}$')
ax2 = ax1.twinx()
ax2.set_ylim(0, 75)
ax2.set_yticks(np.arange(0, 76, 5))
ax2.set_ylabel(r'$I_1$ / A')
ax2.plot(n, I1, 'orange', label=r'$I_1$')
fig.legend(loc='lower left', bbox_to_anchor=(0,0), bbox_transform=ax1.transAxes)
plt.show()

fig, ax = plt.subplots(figsize=(6, 4))
ax.set_xlim(0, 1500)
ax.set_xticks(np.arange(0, 1502, 100))
ax.set_ylim(204, 220)
ax.set_yticks(np.arange(204, 221, 2))
ax.grid()
ax.set_ylabel(r'$U_1$ / V,   $U_{\mathrm{1,sist}}$ / V')
ax.set_xlabel(r'$n_{\mathrm{m}}$ / r/min')
ax.plot(n, U1, 'red', label=r'$U_1$')
ax.plot((0,1500), (U_1_sist,U_1_sist), 'green', label=r'$U_{\mathrm{1,sist}}$')
ax.legend(loc='lower right')
plt.show()