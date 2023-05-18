import numpy as np
import matplotlib.pyplot as plt

f = 50  # Hz
U = 400  # V
φ = 0  # rad
R = 9e-4  # Ω
L = 5e-5  # H

ω = 2*np.pi*f  # rad/s
α = φ - np.arctan2(ω*L, R)  # rad
I_sim = U/np.hypot(R, ω*L)  # A
I_asim_pic = 2*np.sqrt(2)*I_sim  # A
I_asim = np.sqrt(3)*I_sim  # A

print(f'Intensitat simètrica eficaç = {I_sim:.4f} A\n')
print(f'Intensitat asimètrica de pic = {I_asim_pic/1000:.1f} kA\n')
print(f'Intensitat asimètrica eficaç = {I_asim/1000:.1f} kA')

def i_RL_ac(t):
    return np.sqrt(2)*I_sim*np.sin(ω*t + α)

def i_RL_dc(t):
    return -np.sqrt(2)*I_sim*np.sin(α)*np.exp(-t*R/L)

t = np.linspace(0, 200, 500)  # ms
i_ac = i_RL_ac(t/1000)/1000  # kA
i_dc = i_RL_dc(t/1000)/1000  # kA

fig, ax = plt.subplots(figsize=(6, 4))
ax.set_xlim(0, 200)
ax.set_xticks(np.arange(0, 201, 20))
ax.set_ylim(-40, 80)
ax.set_yticks(np.arange(-40, 81, 20))
ax.grid()
ax.set_title('Corrent de curtcircuit')
ax.set_ylabel(r'$i_{\mathrm{dc}}(t)$ / kA,      $i_{\mathrm{ac}}(t)$ / kA,      $i(t)$ / kA')
ax.set_xlabel(r'$t$ / ms')
ax.plot(t, i_dc, 'green', label=r'$i_{\mathrm{dc}}(t)$')
ax.plot(t, i_ac, 'orange', label=r'$i_{\mathrm{ac}}(t)$')
ax.plot(t, i_ac + i_dc, 'red', label=r'$i(t)$')
ax.legend(loc='upper right')
plt.show()
