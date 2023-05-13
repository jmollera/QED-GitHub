import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
from matplotlib import ticker
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
J = 2.1  # kg·m²

U_1 = U/np.sqrt(3)
U_1_80 = 0.8*U_1
print(f'U_1,100% = {U_1:.3f} V')
print(f'U_1,80% = {U_1_80:.3f} V\n')

mot = Motor3ph(p, f, R_1, X_1, R_2, X_2, R_Fe, X_m)

print(f'n_m,sinc = {mot.n_m_sync:.0f} r/min,  ω_m,sinc = {mot.ω_m_sync:.3f} rad/s\n')

def T_load(s):
    return 3.2*(1-s) + 58.9*(1-s)**2

s_100 = mot.s_match(T_load, U_1)
mot_100 = mot.run(s_100, U_1)
mot_100_arr = mot.run(1, U_1)
mot_100_arrencada = mot.start_up(T_load, J, U_1)
print(f's,100% = {s_100:.4f}')
print(f'n_m,100% = {mot.to_n_m(s=s_100):.1f} r/min')
print(f'Z_mot,100% = {mot_100.Z:.3f} Ω')
print(f'I_1,100% = {abs(mot_100.I):.1f} A')
print(f'I_1,arr,100% = {abs(mot_100_arr.I):.1f} A')
print(f'Tm,100% = {mot_100.T_m:.1f} N·m')
print(f'Tm,arr,100% = {mot_100_arr.T_m:.1f} N·m')
print(f'cos 𝜑,100% = {mot_100.PF:.4f}')
print(f'η,100% = {mot_100.Eff:.4f}')
print(f't_arr,100% = {mot_100_arrencada.t_start:.1f} s\n')

s_80 = mot.s_match(T_load, U_1_80)
mot_80 = mot.run(s_80, U_1_80)
mot_80_arr = mot.run(1, U_1_80)
mot_80_arrencada = mot.start_up(T_load, J, U_1_80)
print(f's,80% = {s_80:.4f}')
print(f'n_m,80% = {mot.to_n_m(s=s_80):.1f} r/min')
print(f'Z_mot,80% = {mot_80.Z:.3f} Ω')
print(f'I_1,80% = {abs(mot_80.I):.1f} A')
print(f'I_1,arr,80% = {abs(mot_80_arr.I):.1f} A')
print(f'Tm,80% = {mot_80.T_m:.1f} N·m')
print(f'Tm,arr,80% = {mot_80_arr.T_m:.1f} N·m')
print(f'cos 𝜑,80% = {mot_80.PF:.4f}')
print(f'η,80% = {mot_80.Eff:.4f}')
print(f't_arr,80% = {mot_80_arrencada.t_start:.1f} s\n')

n_m = np.linspace(0, mot.n_m_sync, 500)  # r/min
s = mot.to_s(n_m=n_m)
mot_100 = mot.run(s, U_1)
mot_80 = mot.run(s, U_1_80)

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
ax1.plot(n_m, mot_100.T_m, 'red', label=r'$T_{\mathrm{m},100 \%}$')
ax1.plot(n_m, mot_80.T_m, 'green', label=r'$T_{\mathrm{m},80 \%}$')
ax2 = ax1.twinx()
ax2.set_ylim(0, 75)
ax2.set_yticks(np.arange(0, 76, 5))
ax2.set_ylabel(r'$I_1$ / A')
ax2.plot(n_m, abs(mot_100.I), 'orange', label=r'$I_{1,100 \%}$')
ax2.plot(n_m, abs(mot_80.I), 'lime', label=r'$I_{1,80 \%}$')
fig.legend(loc='lower left', bbox_to_anchor=(0, 0), bbox_transform=ax1.transAxes)
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
ax1.plot(n_m, T_load(s), 'black', label=r'$T_{\mathrm{load}}$')
ax1.plot(n_m, mot_100.T_m, 'red', label=r'$T_{\mathrm{m},100 \%}$')
ax1.plot(n_m, mot_80.T_m, 'green', label=r'$T_{\mathrm{m},80 \%}$')
ax2 = ax1.twinx()
ax2.set_ylim(0, 45)
ax2.set_yticks(np.arange(0, 46, 5))
ax2.set_ylabel(r'$I_1$ / A')
ax2.plot(n_m, abs(mot_100.I), 'orange', label=r'$I_{1,100 \%}$')
ax2.plot(n_m, abs(mot_80.I), 'lime', label=r'$I_{1,80 \%}$')
fig.legend(loc='lower left', bbox_to_anchor=(0, 0), bbox_transform=ax1.transAxes)
plt.show()

U1_pu = np.linspace(0.75, 1, 500)  # (75 % a 100 %) U_1,n
s = np.empty_like(U1_pu)
for i, pu in enumerate(U1_pu):
    sol = optimize.root(lambda s: mot.run(s, pu*U_1).T_m - T_load(s), x0=s_n)
    s[i] = sol.x[0]
mot_sol = mot.run(s, U1_pu * U_1)

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
ax2.plot(U1_pu, abs(mot_sol.I), 'green', label=r'$I_1$')
fig.legend(loc='upper right', bbox_to_anchor=(1, 1), bbox_transform=ax1.transAxes)
plt.show()

fig, ax = plt.subplots(figsize=(6, 4))
ax.set_xlim(0.75, 1.0)
ax.set_xticks(np.arange(0.75, 1.01, 0.05))
ax.set_ylim(0.82, 0.9)
ax.set_yticks(np.arange(0.82, 0.91, 0.02))
ax.grid()
ax.set_ylabel(r'$\cos\varphi$,   $\eta$')
ax.set_xlabel(r'$U_1/U_{\mathrm{n}}$')
ax.plot(U1_pu, mot_sol.PF, 'red', label=r'$\cos\varphi$')
ax.plot(U1_pu, mot_sol.Eff, 'green', label=r'$\eta$')
ax.legend(loc='lower right')
plt.show()

def log_format(x, pos):
    decimal_places = int(np.maximum(-np.log10(x), 0))
    return f'{x:.{decimal_places:1d}f}'

t_arr = np.linspace(0, 10, 500)  # time form 0 s to 10 s
mot_100_arr = mot.run(mot_100_arrencada.s_start_up(t_arr), U_1)
mot_80_arr = mot.run(mot_80_arrencada.s_start_up(t_arr), U_1_80)

fig, ax = plt.subplots(figsize=(6, 4))
ax.set(yscale='log')
ax.yaxis.set_major_formatter(ticker.FuncFormatter(log_format))
ax.set_xlim(15, 80)
ax.set_xticks(np.arange(15, 81, 5))
ax.set_ylim(0.01, 10)
ax.grid(which='minor', alpha=0.4)
ax.grid(which='major', alpha=1)
ax.set_xlabel(r'$I_1$ / A')
ax.set_ylabel(r'$t$ / s')
ax.plot(np.abs(mot_100_arr.I), t_arr, 'red', label=r'$I_{1,100 \%}$')
ax.plot(np.abs(mot_80_arr.I), t_arr, 'green', label=r'$I_{1,80 \%}$')
ax.legend(loc='lower left')
plt.show()
