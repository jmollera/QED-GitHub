import numpy as np
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
U = 380  # V (fase-fase)
f = 50  # Hz

U_1 = U/np.sqrt(3)
n_m_sinc = 120*f/p
ω_m_sinc = 4*np.pi*f/p
n_m_n = (1-s_n)*n_m_sinc
ω_m_n = (1-s_n)*ω_m_sinc
print(f'U_1 = {U_1:.3f} V')
print(f'n_m,sinc = {n_m_sinc:.0f} r/min,  ω_m,sinc = {ω_m_sinc:.3f} rad/s')
print(f'n_m,n = {n_m_n:.0f} r/min,  ω_m,n = {ω_m_n:.3f} rad/s\n')
Z_0 = R_Fe*1j*X_m/(R_Fe + 1j*X_m)
Z_Th = (R_1 + 1j*X_1)*Z_0/(R_1 + 1j*X_1 + Z_0)
E_Th = U_1*Z_0/(R_1 + 1j*X_1 + Z_0)
print(f'Z_0 = {Z_0:.3f} Ω')
print(f'Z_Th = {Z_Th:.3f} Ω')
print(f'E_Th = {ComplexD(E_Th):.4f} V\n')

def Z_mot(s):
    return R_1 + 1j*X_1 + Z_0*(R_2/s + 1j*X_2)/(Z_0 + R_2/s + 1j*X_2)

def I_1(s):
    return U_1/Z_mot(s)

def I_2(s):
    return E_Th/(Z_Th + R_2/s + 1j*X_2)

def P_m(s):
    return 3*(1-s)*R_2*abs(I_2(s))**2/s

def T_m(s):
    return 3*R_2*abs(I_2(s))**2/(s*ω_m_sinc)

I_1_n = abs(I_1(s=s_n))
Z_arr = Z_mot(s=1)
I_1_arr = abs(I_1(s=1))
κ = 1.02 + 0.98*np.exp(-3*Z_arr.real/Z_arr.imag)
I_1_arr_pic_asim = κ*np.sqrt(2)*I_1_arr
print(f'I_1,n = {I_1_n:.1f} A')
print(f'I_1,arr = {I_1_arr:.1f} A,  I_1,arr/I_1,n = {I_1_arr/I_1_n:.1f}')
print(f'I_1,arr,pic,asim = {I_1_arr_pic_asim:.2f} A,  I_1,arr,pic,asim/I_1,arr = {I_1_arr_pic_asim/I_1_arr:.2f}\n')
P_m_n = P_m(s=s_n)
P_m_arr = P_m(s=1)
print(f'P_m,n = {P_m_n:.1f} W')
print(f'P_m,arr = {P_m_arr:.1f} W\n')
T_m_n = T_m(s=s_n)
T_m_arr = T_m(s=1)
print(f'T_m,n = {T_m_n:.1f} N·m')
print(f'T_m,arr = {T_m_arr:.1f} N·m,  T_m,arr/T_m,n = {T_m_arr/T_m_n:.2f}\n')
cos_φ_n = np.cos(np.angle(Z_mot(s=s_n)))
cos_φ_arr = np.cos(np.angle(Z_mot(s=1)))
print(f'cos φ,n = {cos_φ_n:.3f}')
print(f'cos φ,arr = {cos_φ_arr:.3f}\n')
P_n = 3*U_1*I_1_n*cos_φ_n
η_n = P_m_n/P_n
print(f'η,n = {η_n:.3f}\n')
s_T_m_max = R_2/abs(Z_Th+1j*X_2)
n_T_m_max = (1-s_T_m_max)*n_m_sinc
ω_T_m_max = (1-s_T_m_max)*ω_m_sinc
T_m_max = T_m(s=s_T_m_max)
print(f's,T_m,max = {s_T_m_max:.4f},  n,T_m,max = {n_T_m_max:.2f} r/min,  ω,T_m,max = {ω_T_m_max:.2f} rad/s')
print(f'T_m,max = {T_m_max:.1f} N·m,  T_m,max/T_m,n = {T_m_max/T_m_n:.2f}')

n = np.linspace(0, n_m_sinc-0.01, 500) # r/min
I1 = abs(I_1(s=1-n/n_m_sinc)) # A
Tm = T_m(s=1-n/n_m_sinc) # N·m
Pm = P_m(s=1-n/n_m_sinc)/1000 # kW
cos_φ = np.cos(np.angle(Z_mot(s=1-n/n_m_sinc)))
P = 3*U_1*I1*cos_φ/1000 # kW
Q = 3*U_1*I1*np.sqrt(1-cos_φ**2)/1000 # kvar
S = np.hypot(P,Q) # kVA
η = Pm/P

fig, ax1 = plt.subplots(figsize=(6, 4))
plt.rc('xtick', labelsize=8)
ax1.set_xlim(0, 1500)
ax1.set_xticks(np.arange(0, 1501, 100))
ax1.set_ylim(0, 150)
ax1.set_yticks(np.arange(0, 151, 10))
ax1.grid()
ax1.set_ylabel(r'$T_{\mathrm{m}}$ / N m')
ax1.set_xlabel(r'$n_{\mathrm{m}}$ / r/min')
ax1.plot(n, Tm, 'red', label=r'$T_{\mathrm{m}}$')
ax2 = ax1.twinx()
ax2.set_ylim(0, 75)
ax2.set_yticks(np.arange(0, 76, 5))
ax2.set_ylabel(r'$I_1$ / A')
ax2.plot(n, I1, 'green', label=r'$I_1$')
fig.legend(loc='lower left', bbox_to_anchor=(0,0), bbox_transform=ax1.transAxes)
plt.show()

fig, ax = plt.subplots(figsize=(6, 4))
plt.rc('xtick', labelsize=8)
ax.set_xlim(0, 1500)
ax.set_xticks(np.arange(0, 1501, 100))
ax.set_ylim(0, 50)
ax.set_yticks(np.arange(0, 51, 5))
ax.grid()
ax.set_ylabel(r'$P_{\mathrm{m}}$ / kW,   $P$ / kW,   $Q$ / kvar,   $S$ / kVA')
ax.set_xlabel(r'$n_{\mathrm{m}}$ / r/min')
ax.plot(n, Pm, 'red', label=r'$P_{\mathrm{m}}$')
ax.plot(n, P, 'green', label=r'$P$')
ax.plot(n, Q, 'orange', label=r'$Q$')
ax.plot(n, S, 'black', label=r'$S$')
ax.legend(loc='upper right')
plt.show()

fig, ax = plt.subplots(figsize=(6, 4))
plt.rc('xtick', labelsize=8)
ax.set_xlim(0, 1500)
ax.set_xticks(np.arange(0, 1501, 100))
ax.set_ylim(0, 1)
ax.set_yticks(np.arange(0, 1.01, 0.1))
ax.grid()
ax.set_ylabel(r'$\cos\varphi$,   $\eta$')
ax.set_xlabel(r'$n_{\mathrm{m}}$ / r/min')
ax.plot(n, cos_φ, 'red', label=r'$\cos\varphi$')
ax.plot(n, η, 'green', label=r'$\eta$')
ax.legend(loc='lower right')
plt.show()