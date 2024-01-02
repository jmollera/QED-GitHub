import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from qed.utils import Complex
from qed.eng_elec import Motor3ph

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
J = 2.1  # kg·m²

U_1_ext = U_sist / np.sqrt(3)
fase = np.arctan(Rel_X_R)
Z_sist = U_sist ** 2 / (Scc_sist * 1e6) * (np.cos(fase) + 1j * np.sin(fase))
Z_cable = L_cable * 1e-3 * (R_cable + 1j * X_cable)
print(f'U_1,sist = {U_1_ext:.3f} V')
print(f'Z_sist = {Z_sist:.3f} Ω')
print(f'Z_cable = {Z_cable:.3f} Ω\n')

motor = Motor3ph(p, f, R_1, X_1, R_2, X_2, R_Fe, X_m)
font_ext = (U_1_ext, Z_sist + Z_cable)

print(f'n_m,sinc = {motor.n_m_sync:.0f} r/min,  ω_m,sinc = {motor.ω_m_sync:.3f} rad/s\n')

def T_load(s):
    return 3.2 * (1 - s) + 58.9 * (1 - s) ** 2

s_fun = motor.s_match(T_load, *font_ext)
mot_fun = motor.run(s_fun, *font_ext)
mot_arr = motor.run(1, *font_ext)
print(f's,fun = {s_fun:.4f},   n_m_fun = {motor.to_n_m(s=s_fun):.1f} r/min')
print(f'U_1,fun = {Complex(mot_fun.U) :.3f} V')
print(f'I_1,fun = {abs(mot_fun.I):.1f} A')
print(f'T_m,fun = {mot_fun.T_m:.1f} N·m\n')
print(f'U_1,arr = {Complex(mot_arr.U) :.3f} V')
print(f'I_1,arr = {abs(mot_arr.I):.1f} A')
print(f'T_m,arr = {mot_arr.T_m:.1f} N·m\n')

s_T_m_max = motor.s_T_m_max(Z_sist + Z_cable)
n_T_m_max = motor.to_n_m(s=s_T_m_max)
ω_T_m_max = motor.to_ω_m(s=s_T_m_max)
T_m_max = motor.run(s_T_m_max, *font_ext).T_m
print(f's_T_m,max = {s_T_m_max:.4f}')
print(f'n_T_m,max = {n_T_m_max:.1f} r/min')
print(f'T_m,max = {T_m_max:.1f} N·m\n')

mot_arrencada = motor.start_up(T_load, J, *font_ext)
print(f't_arr = {mot_arrencada.t_start:.1f} s')

n_m = np.linspace(0, motor.n_m_sync, 500)  # r/min
s = motor.to_s(n_m=n_m)
mot = motor.run(s, *font_ext)

fig, ax1 = plt.subplots(figsize=(6, 4))
plt.rc('xtick', labelsize=8)
ax1.set_xlim(0, 1500)
ax1.set_xticks(np.arange(0, 1501, 100))
ax1.set_ylim(0, 150)
ax1.set_yticks(np.arange(0, 151, 10))
ax1.grid()
ax1.set_ylabel(r'$T_{\mathrm{m}}$ / N m,   $T_{\mathrm{load}}$ / N m')
ax1.set_xlabel(r'$n_{\mathrm{m}}$ / r/min')
ax1.plot(n_m, T_load(s), 'black', label=r'$T_{\mathrm{load}}$')
ax1.plot(n_m, mot.T_m, 'red', label=r'$T_{\mathrm{m}}$')
ax2 = ax1.twinx()
ax2.set_ylim(0, 75)
ax2.set_yticks(np.arange(0, 76, 5))
ax2.set_ylabel(r'$I_1$ / A')
ax2.plot(n_m, abs(mot.I), 'orange', label=r'$I_1$')
fig.legend(loc='lower left', bbox_to_anchor=(0, 0), bbox_transform=ax1.transAxes)
plt.show()

fig, ax = plt.subplots(figsize=(6, 4))
ax.set_xlim(0, 1500)
ax.set_xticks(np.arange(0, 1502, 100))
ax.set_ylim(205, 220)
ax.set_yticks(np.arange(205, 221, 1))
ax.grid()
ax.set_ylabel(r'$U_1$ / V,   $U_{\mathrm{1,sist}}$ / V')
ax.set_xlabel(r'$n_{\mathrm{m}}$ / r/min')
ax.plot(n_m, abs(mot.U), 'red', label=r'$U_1$')
ax.plot((0, 1500), (U_1_ext, U_1_ext), 'green', label=r'$U_{\mathrm{1,sist}}$')
ax.legend(loc='best')
plt.show()

def log_format(x, pos):
    decimal_places = int(np.maximum(-np.log10(x), 0))
    return f'{x:.{decimal_places:1d}f}'

t_arr = np.linspace(0, 10, 500)  # temps: de 0 s a 10 s
s_t_arr = mot_arrencada.s_start_up(t_arr)
mot_t_arr = motor.run(s_t_arr, *font_ext)
T_load_t_arr = T_load(s_t_arr)

fig, ax = plt.subplots(figsize=(6, 4))
ax.set(yscale='log')
ax.yaxis.set_major_formatter(ticker.FuncFormatter(log_format))
ax.set_xlim(15, 75)
ax.set_xticks(np.arange(15, 76, 5))
ax.set_ylim(0.01, 10)
ax.grid(which='minor', alpha=0.4)
ax.grid(which='major', alpha=1)
ax.set_xlabel(r'$I_1$ / A')
ax.set_ylabel(r'$t$ / s')
ax.plot(np.abs(mot_t_arr.I), t_arr, 'orange')
plt.show()

fig, ax = plt.subplots(figsize=(6, 4))
ax.set_xlim(0, 7)
ax.set_xticks(np.arange(0, 7.1, 0.5))
ax.set_ylim(205, 215)
ax.set_yticks(np.arange(205, 216, 1))
ax.grid()
ax.set_xlabel(r'$t$ / s')
ax.set_ylabel(r'$U_1$ / V')
ax.plot(t_arr, np.abs(mot_t_arr.U), 'red')
plt.show()

fig, ax = plt.subplots(figsize=(6, 4))
ax.set_xlim(0, 7)
ax.set_xticks(np.arange(0, 7.1, 0.5))
ax.set_ylim(0, 120)
ax.set_yticks(np.arange(0, 121, 10))
ax.grid()
ax.set_xlabel(r'$t$ / s')
ax.set_ylabel(r'$T_{\mathrm{m}}$ / N m,   $T_{\mathrm{load}}$ / N m')
ax.plot(t_arr, T_load_t_arr, 'black', label=r'$T_{\mathrm{load}}$')
ax.plot(t_arr, mot_t_arr.T_m, 'red', label=r'$T_{\mathrm{m}}$')
ax.legend(loc='best')
plt.show()
