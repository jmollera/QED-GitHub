import numpy as np
import matplotlib.pyplot as plt

E_bat = 12  # V
R_bat = 2   # Ω
R = 4       # Ω
L = 60e-3   # H

τ = L / R

def i_u(t, I, E, R, L):
    u_L = (E - R*I) * np.exp(-t/τ)
    i_L = (E - u_L) / R
    return i_L, u_L

t = np.linspace(0, 0.1, 500)  # s
t_1 = 0.02  # s
t_2 = 0.04  # s

i_1, u_1 = i_u(t[t <= t_1], 0, E_bat, R_bat + R, L)
I_1 = i_1[-1]

i_2, u_2 = i_u(t[(t > t_1) & (t <= t_2)] - t_1, I_1, 0, R, L)
I_2 = i_2[-1]

i_3, u_3 = i_u(t[t > t_2] - t_2, I_2, E_bat, R_bat + R, L)

fig, ax1 = plt.subplots(figsize=(6, 4))
ax1.set_xlim(0, 100)
ax1.set_xticks(np.arange(0, 101, 10))
ax1.set_ylim(-12, 12)
yticks1 = np.arange(-12, 13, 3)
ax1.set_yticks(yticks1, [f'{tick: }' for tick in yticks1])
ax1.grid()
ax1.set_ylabel(r'$u(t)$ / V')
ax1.set_xlabel(r'$t$ / ms')
ax1.plot(t*1000, np.concatenate((u_1, u_2, u_3)), 'green', label=r'$u(t)$')
ax2 = ax1.twinx()
ax2.set_ylim(0, 2)
yticks2 = np.arange(0, 2.1, 0.25)
ax2.set_yticks(yticks2, [f'{tick:.2f}'.replace('.', ',') for tick in yticks2])
ax2.set_ylabel(r'$i(t)$ / A')
ax2.plot(t*1000, np.concatenate((i_1, i_2, i_3)), 'red', label=r'$i(t)$')
fig.legend(loc='lower right', bbox_to_anchor=(1,0), bbox_transform=ax1.transAxes)
plt.show()
