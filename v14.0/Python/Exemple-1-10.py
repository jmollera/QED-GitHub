import numpy as np
from scipy import optimize

f = 50  # Hz
U = 400  # V
φ = 0  # rad
R = 9e-4  # Ω
L = 5e-5  # H

ω = 2*np.pi*f  # rad/s
α = φ - np.arctan2(ω*L, R)  # rad
I_sim = U/np.hypot(R, ω*L)  # A
κ = 1.02 + 0.98*np.exp(-3*R/(ω*L))
I_asim_pic_CEI = κ*np.sqrt(2)*I_sim  # A

def i(t):
    return np.sqrt(2)*I_sim*np.sin(ω*t + α) - np.sqrt(2)*I_sim*np.sin(α)*np.exp(-t*R/L)

sol = optimize.minimize_scalar(lambda x: -i(x), bounds=[5e-3, 15e-3], method='bounded')
t_asim_pic, I_asim_pic = sol.x, -sol.fun

print(f'Factor κ (CEI 60909-1) = {κ:.4f}')
print(f'Intensitat asimètrica de pic (CEI 60909-1) = {I_asim_pic_CEI/1000:.1f} kA')
print(f'Intensitat asimètrica de pic (exacta) = {I_asim_pic:.1f} A (t = {t_asim_pic*1000:.4f} ms)')