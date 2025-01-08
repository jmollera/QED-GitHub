import numpy as np
import matplotlib.pyplot as plt

T = 10e-3  # s
U = 20     # V
R = 10     # Ω
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
yticks1 = np.arange(-25, 26, 5)
ax1.set_yticks(yticks1, [f'{tick: }' for tick in yticks1])
ax1.grid()
ax1.set_ylabel(r'$u(t)$ / V')
ax1.set_xlabel(r'$t$ / ms')
ax1.plot(t*1000, u, 'red', label=r'$u(t)$')
ax2 = ax1.twinx()
ax2.set_ylim(-1.25, 1.25)
yticks2 = np.arange(-1.25, 1.3, 0.25)
ax2.set_yticks(yticks2, [f'{tick: .2f}'.replace('.', ',') for tick in yticks2])
ax2.set_ylabel(r'$i(t)$ / A')
ax2.plot(t*1000, i, 'green', label=r'$i(t)$')
fig.legend(loc='upper center', bbox_to_anchor=(0.5,1), bbox_transform=ax1.transAxes)
plt.show()

τ = L/R
def i_L(t, I, E, R):
     return E/R - (E/R - I)*np.exp(-t/τ)

t_0 = 0.000  # s
t_1 = 0.005  # s
t_2 = 0.010  # s
t_3 = 0.015  # s
t_4 = 0.020  # s
t_5 = 0.025  # s
t_6 = 0.030  # s
t_7 = 0.035  # s
t_8 = 0.040  # s
n_samples = 1000
t, t_step = np.linspace(t_0, t_8, n_samples, retstep=True)  # s
t_0_1 = t <= t_1
t_1_2 = (t > t_1) & (t <= t_2)
t_2_3 = (t > t_2) & (t <= t_3)
t_3_4 = (t > t_3) & (t <= t_4)
t_4_5 = (t > t_4) & (t <= t_5)
t_5_6 = (t > t_5) & (t <= t_6)
t_6_7 = (t > t_6) & (t <= t_7)
t_7_8 = t > t_7

i_1 = i_L(t[t_0_1], 0, U, R)
I_1 = i_1[-1]
i_2 = i_L(t[t_1_2] - t_1, I_1, -U, R)
I_2 = i_2[-1]
i_3 = i_L(t[t_2_3] - t_2, I_2, U, R)
I_3 = i_3[-1]
i_4 = i_L(t[t_3_4] - t_3, I_3, -U, R)
I_4 = i_4[-1]
i_5 = i_L(t[t_4_5] - t_4, I_4, U, R)
I_5 = i_5[-1]
i_6 = i_L(t[t_5_6] - t_5, I_5, -U, R)
I_6 = i_6[-1]
i_7 = i_L(t[t_6_7] - t_6, I_6, U, R)
I_7 = i_7[-1]
i_8 = i_L(t[t_7_8] - t_7, I_7, -U, R)

i_L_total = np.concatenate((i_1, i_2, i_3, i_4, i_5, i_6, i_7, i_8))

u_total = np.piecewise(t, [t_0_1, t_1_2, t_2_3, t_3_4, t_4_5, t_5_6, t_6_7, t_7_8],
                        [U, -U, U, -U, U, -U, U, -U])

fig, ax1 = plt.subplots(figsize=(6, 4))
ax1.set_xlim(0, 40)
ax1.set_xticks(np.arange(0, 41, 5))
ax1.set_ylim(-25, 25)
yticks1 = np.arange(-25, 26, 5)
ax1.set_yticks(yticks1, [f'{tick: }' for tick in yticks1])
ax1.grid()
ax1.set_ylabel(r'$u(t)$ / V')
ax1.set_xlabel(r'$t$ / ms')
ax1.plot(t*1000, u_total, 'red', label=r'$u(t)$')
ax2 = ax1.twinx()
ax2.set_ylim(-1.25, 1.25)
yticks2 = np.arange(-1.25, 1.3, 0.25)
ax2.set_yticks(yticks2, [f'{tick: .2f}'.replace('.', ',') for tick in yticks2])
ax2.set_ylabel(r'$i(t)$ / A')
ax2.plot(t*1000, i_L_total, 'green', label=r'$i(t)$')
fig.legend(loc='lower left', bbox_to_anchor=(0,0), bbox_transform=ax1.transAxes)
plt.show()

F_u_total = np.fft.rfft(u_total)
freq = np.fft.rfftfreq(n_samples, t_step)
spectrum_u_total = 2 / n_samples * np.abs(F_u_total)

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
