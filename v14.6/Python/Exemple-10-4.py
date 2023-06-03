import numpy as np
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

U_1 = U / np.sqrt(3)

mot = Motor3ph(p, f, R_1, X_1, R_2, X_2, R_Fe, X_m)

print(f'U_1 = {U_1:.3f} V')
print(f'n_m,sinc = {mot.n_m_sync:.0f} r/min,  ω_m,sinc = {mot.ω_m_sync:.3f} rad/s')
print(f'n_m,n = {mot.to_n_m(s=s_n):.1f} r/min,  ω_m,n = {mot.to_ω_m(s=s_n):.3f} rad/s\n')

mot_n = mot.run(s_n, U_1)
mot_arr = mot.run(1, U_1)

I_1_n = abs(mot_n.I)
I_1_arr = abs(mot_arr.I)
Z_arr = mot_arr.Z
κ = 1.02 + 0.98*np.exp(-3*Z_arr.real/Z_arr.imag)
I_1_arr_pic_asim = κ*np.sqrt(2)*I_1_arr
print(f'I_1,n = {I_1_n:.1f} A')
print(f'I_1,arr = {I_1_arr:.1f} A,  I_1,arr/I_1,n = {I_1_arr/I_1_n:.1f}')
print(f'I_1,arr,pic,asim = {I_1_arr_pic_asim:.1f} A,  I_1,arr,pic,asim/I_1,arr = {I_1_arr_pic_asim/I_1_arr:.2f}\n')

print(f'P_m,n = {mot_n.P_m:.1f} W')
print(f'P_m,arr = {mot_arr.P_m:.1f} W\n')

T_m_n = mot_n.T_m
T_m_arr = mot_arr.T_m
print(f'T_m,n = {T_m_n:.1f} N·m')
print(f'T_m,arr = {T_m_arr:.1f} N·m,  T_m,arr/T_m,n = {T_m_arr/T_m_n:.2f}\n')

print(f'cos 𝜑,n = {mot_n.PF:.4f}')
print(f'cos 𝜑,arr = {mot_arr.PF:.4f}\n')

print(f'η,n = {mot_n.Eff:.4f}\n')

s_T_m_max = mot.s_T_m_max()
n_T_m_max = mot.to_n_m(s=s_T_m_max)
ω_T_m_max = mot.to_ω_m(s=s_T_m_max)
T_m_max = mot.run(s_T_m_max, U_1).T_m
print(f's,T_m,max = {s_T_m_max:.4f},  n,T_m,max = {n_T_m_max:.1f} r/min,  ω,T_m,max = {ω_T_m_max:.3f} rad/s')
print(f'T_m,max = {T_m_max:.1f} N·m,  T_m,max/T_m,n = {T_m_max/T_m_n:.2f}')

n_m = np.linspace(0, mot.n_m_sync, 500)  # r/min
mot_n_m = mot.run(mot.to_s(n_m=n_m), U_1)

fig, ax1 = plt.subplots(figsize=(6, 4))
plt.rc('xtick', labelsize=8)
ax1.set_xlim(0, 1500)
ax1.set_xticks(np.arange(0, 1501, 100))
ax1.set_ylim(0, 150)
ax1.set_yticks(np.arange(0, 151, 10))
ax1.grid()
ax1.set_ylabel(r'$T_{\mathrm{m}}$ / N m')
ax1.set_xlabel(r'$n_{\mathrm{m}}$ / r/min')
ax1.plot(n_m, mot_n_m.T_m, 'red', label=r'$T_{\mathrm{m}}$')
ax2 = ax1.twinx()
ax2.set_ylim(0, 75)
ax2.set_yticks(np.arange(0, 76, 5))
ax2.set_ylabel(r'$I_1$ / A')
ax2.plot(n_m, np.abs(mot_n_m.I), 'green', label=r'$I_1$')
fig.legend(loc='lower left', bbox_to_anchor=(0, 0), bbox_transform=ax1.transAxes)
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
S = mot_n_m.S / 1000
ax.plot(n_m, mot_n_m.P_m/1000, 'red', label=r'$P_{\mathrm{m}}$')
ax.plot(n_m, np.real(S), 'green', label=r'$P$')
ax.plot(n_m, np.imag(S), 'orange', label=r'$Q$')
ax.plot(n_m, np.abs(S), 'black', label=r'$S$')
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
ax.plot(n_m, mot_n_m.PF, 'red', label=r'$\cos\varphi$')
ax.plot(n_m, mot_n_m.Eff, 'green', label=r'$\eta$')
ax.legend(loc='lower right')
plt.show()
