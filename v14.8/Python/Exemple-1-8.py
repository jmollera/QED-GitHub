import numpy as np
import matplotlib.pyplot as plt

E_bat = 12  # V
R_bat = 2  # Ω
R = 4  # Ω
L = 60e-3  # H
t = np.linspace(0, 0.1, 500)  # s
t_1 = 0.02  # s
t_2 = 0.04  # s

def i_L(t, I, E, R, L):
    τ = L/R
    return E/R - (E/R - I)*np.exp(-t/τ)

def u_L(t, I, E, R, L):
    τ = L/R
    return (E - R*I)*np.exp(-t/τ)

i_1 = i_L(t[t <= t_1], 0, E_bat, R_bat + R, L)
I_1 = i_L(t_1, 0, E_bat, R_bat + R, L)
u_1 = u_L(t[t <= t_1], 0, E_bat, R_bat + R, L)
i_2 = i_L(t[(t > t_1) & (t <= t_2)] - t_1, I_1, 0, R, L)
I_2 = i_L(t_2 - t_1, I_1, 0, R, L)
u_2 = u_L(t[(t > t_1) & (t <= t_2)] - t_1, I_1, 0, R, L)
i_3 = i_L(t[t > t_2] - t_2, I_2, E_bat, R_bat + R, L)
u_3 = u_L(t[t > t_2] - t_2, I_2, E_bat, R_bat + R, L)

fig, ax1 = plt.subplots(figsize=(6, 4))
ax1.set_xlim(0, 100)
ax1.set_xticks(np.arange(0, 101, 10))
ax1.set_ylim(-12, 12)
ax1.set_yticks(np.arange(-12, 13, 3))
ax1.grid()
ax1.set_ylabel(r'$u(t)$ / V')
ax1.set_xlabel(r'$t$ / ms')
ax1.plot(t*1000, np.concatenate((u_1, u_2, u_3)), 'green', label=r'$u(t)$')
ax2 = ax1.twinx()
ax2.set_ylim(0, 2)
ax2.set_yticks(np.arange(0, 2.1, 0.25))
ax2.set_ylabel(r'$i(t)$ / A')
ax2.plot(t*1000, np.concatenate((i_1, i_2, i_3)), 'red', label=r'$i(t)$')
fig.legend(loc='lower right', bbox_to_anchor=(1,0), bbox_transform=ax1.transAxes)
plt.show()
