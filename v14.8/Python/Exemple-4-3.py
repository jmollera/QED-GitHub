import numpy as np
import matplotlib.pyplot as plt


T = 10e-3  # s
U = 20  # V
R = 10  # Ω
L = 50e-3  # H

ω = 2*np.pi/T  # rad/s
k = np.arange(1, 60, 2)  # 1, 3, 5, ..., 59
U_k = 4*U/np.pi/k
Z_k = R + k*ω*L*1j
φ_k = np.angle(Z_k)
I_k = U_k/np.abs(Z_k)

I = np.sqrt(np.sum(I_k**2/2))
P = R*I**2
print(f'I = {I:.4f} A')
print(f'P = {P:.2f} W')
P = np.sum(U_k*I_k/2*np.cos(φ_k))
print(f'P = {P:.2f} W')

t = np.linspace(0, 16e-3, 500)  # s
u = np.vectorize(lambda time: np.sum(U_k * np.sin(k * ω * time)))(t)
i = np.vectorize(lambda time: np.sum(I_k * np.sin(k * ω * time - φ_k)))(t)

fig, ax1 = plt.subplots(figsize=(6, 4))
ax1.set_xlim(0, 16)
ax1.set_xticks(np.arange(0, 16.1, 1))
ax1.set_ylim(-25, 25)
ax1.set_yticks(np.arange(-25, 26, 5))
ax1.grid()
ax1.set_ylabel(r'$u(t)$ / V')
ax1.set_xlabel(r'$t$ / ms')
ax1.plot(t*1000, u, 'red', label=r'$u(t)$')
ax2 = ax1.twinx()
ax2.set_ylim(-1.25, 1.25)
ax2.set_yticks(np.arange(-1.25, 1.3, 0.25))
ax2.set_ylabel(r'$i(t)$ / A')
ax2.plot(t*1000, i, 'green', label=r'$i(t)$')
fig.legend(loc='upper center', bbox_to_anchor=(0.5,1), bbox_transform=ax1.transAxes)
plt.show()

t = np.linspace(0, 0.040, 1000)  # s
t_1 = 0.005  # s
t_2 = 0.010  # s
t_3 = 0.015  # s
t_4 = 0.020  # s
t_5 = 0.025  # s
t_6 = 0.030  # s
t_7 = 0.035  # s

def i_L(t, I, E, R, L):
    τ = L/R
    return E/R - (E/R - I)*np.exp(-t/τ)

def u(t, E):
    return E*np.ones_like(t)

i_1 = i_L(t[t <= t_1], 0, U, R, L)
I_1 = i_L(t_1, 0, U, R, L)
u_1 = u(t[t <= t_1], U)
i_2 = i_L(t[(t > t_1) & (t <= t_2)] - t_1, I_1, -U, R, L)
I_2 = i_L(t_2 - t_1, I_1, -U, R, L)
u_2 = u(t[(t > t_1) & (t <= t_2)] - t_1, -U)
i_3 = i_L(t[(t > t_2) & (t <= t_3)] - t_2, I_2, U, R, L)
I_3 = i_L(t_3 - t_2, I_2, U, R, L)
u_3 = u(t[(t > t_2) & (t <= t_3)] - t_2, U)
i_4 = i_L(t[(t > t_3) & (t <= t_4)] - t_3, I_3, -U, R, L)
I_4 = i_L(t_4 - t_3, I_3, -U, R, L)
u_4 = u(t[(t > t_3) & (t <= t_4)] - t_3, -U)
i_5 = i_L(t[(t > t_4) & (t <= t_5)] - t_4, I_4, U, R, L)
I_5 = i_L(t_5 - t_4, I_4, U, R, L)
u_5 = u(t[(t > t_4) & (t <= t_5)] - t_4, U)
i_6 = i_L(t[(t > t_5) & (t <= t_6)] - t_5, I_5, -U, R, L)
I_6 = i_L(t_6 - t_5, I_5, -U, R, L)
u_6 = u(t[(t > t_5) & (t <= t_6)] - t_5, -U)
i_7 = i_L(t[(t > t_6) & (t <= t_7)] - t_6, I_6, U, R, L)
I_7 = i_L(t_7 - t_6, I_6, U, R, L)
u_7 = u(t[(t > t_6) & (t <= t_7)] - t_6, U)
i_8 = i_L(t[t > t_7] - t_7, I_7, -U, R, L)
u_8 = u(t[t > t_7] - t_7, -U)
u_total = np.concatenate((u_1, u_2, u_3, u_4, u_5, u_6, u_7, u_8))
i_L_total = np.concatenate((i_1, i_2, i_3, i_4, i_5, i_6, i_7, i_8))

fig, ax1 = plt.subplots(figsize=(6, 4))
ax1.set_xlim(0, 40)
ax1.set_xticks(np.arange(0, 41, 5))
ax1.set_ylim(-25, 25)
ax1.set_yticks(np.arange(-25, 26, 5))
ax1.grid()
ax1.set_ylabel(r'$u(t)$ / V')
ax1.set_xlabel(r'$t$ / ms')
ax1.plot(t*1000, u_total, 'red', label=r'$u(t)$')
ax2 = ax1.twinx()
ax2.set_ylim(-1.25, 1.25)
ax2.set_yticks(np.arange(-1.25, 1.3, 0.25))
ax2.set_ylabel(r'$i(t)$ / A')
ax2.plot(t*1000, i_L_total, 'green', label=r'$i(t)$')
fig.legend(loc='lower left', bbox_to_anchor=(0,0), bbox_transform=ax1.transAxes)
plt.show()

t_sample = t[1] - t[0]
n = len(t)
freq = np.fft.rfftfreq(n, t_sample)[:n // 2]
F_u_total = np.fft.rfft(u_total)[:n // 2]
spectrum_u_total = 2 / n * np.abs(F_u_total)

fig, ax = plt.subplots(figsize=(6, 4))
ax.set_xlim(0, 1200)
ax.set_xticks(np.arange(0, 1201, 100))
ax.set_ylim(0, 26)
ax.set_yticks(np.arange(0, 27, 2))
ax.set_axisbelow(True)
ax.grid()
ax.set_ylabel(r'Amplitud / V')
ax.set_xlabel(r'Freqüència / Hz')
ax.bar(freq, spectrum_u_total, width=10)
plt.show()
