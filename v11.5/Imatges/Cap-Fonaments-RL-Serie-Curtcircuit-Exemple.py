# -*- coding: utf-8 -*-

import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt

# Freqüència, Hz
f = 50
# Tensió, V
U = 400
# Angle inicial de la tensió, rad
φ = 0
# Resistència, Ω
R = 9e-4
# Inductància, H
L = 5e-5

# Velocitat angular, rad/s
ω = 2*np.pi*f

# Angle inicial del corrent AC, rad
α = φ - np.arctan2(ω*L, R)

# Intensitat simètrica (eficaç), A
I_sim = U / np.hypot(R, ω*L)

# Intensitat asimètrica (eficaç), A
I_asim = np.sqrt(3) * I_sim

# Intensitat asimètrica (pic), A
I_asim_pic = 2 * np.sqrt(2) * I_sim

# Factor κ segons CEI 60909-1
κ = 1.02 + 0.98 * np.exp(-3*R/(ω*L))

# Intensitat asimètrica (pic) segons CEI 60909-1, A
I_asim_pic_CEI = κ * np.sqrt(2) * I_sim

# Tensió en funció del temps, V
def u(t):
    return np.sqrt(2) * U * np.sin(ω*t + φ)

# Corrent AC en funció del temps, A
def i_ac(t):
    return np.sqrt(2) * I_sim * np.sin(ω*t + α)

# Corrent DC en funció del temps, A
def i_dc(t):
    return -np.sqrt(2) * I_sim * np.sin(α) * np.exp(-t*R/L)

# Corrent total en funció del temps, A
def i(t):
    return i_ac(t) + i_dc(t)

# Intensitat asimètrica (pic) exacte, A
sol = optimize.minimize_scalar(lambda x: -i(x),
                               bounds=[5e-3, 15e-3],
                               method='bounded')
t_max = sol.x
I_asim_pic_exact = -sol.fun

# Impressió dels resultats
print()
print("Intensitat simètrica (eficaç) = {:.4f} A".format(I_sim))
print("Intensitat asimètrica (eficaç) = {:.4f} A".format(I_asim))
print("Intensitat asimètrica (pic) = {:.4f} A".format(I_asim_pic))
print("Factor κ CEI 60909-1 = {:.4f}".format(κ))
print("Intensitat asimètrica (pic) CEI 60909-1 = {:.4f} A".format(I_asim_pic_CEI))
print("Intensitat asimètrica (pic) exacte = {0:.4f} A (t = {1:.4f} ms)".format(I_asim_pic_exact, t_max*1000))

# Rang de temps, ms
t = np.linspace(0, 200, 500)

# Rang de tensió, V
u_t = u(t/1000)

# Rang de corrent AC, kA
i_ac_t = i_ac(t/1000) / 1000

# Rang de corrent DC, kA
i_dc_t = i_dc(t/1000) / 1000

# Rang de corrent total, kA
i_t = i_ac_t + i_dc_t

# Gràfica de la tensió
plt.plot(t, u_t)
plt.xlim(0, 200)
plt.ylim(-600, 600)
plt.grid()
plt.title("Tensió")
plt.xlabel("t / ms")
plt.ylabel("u(t) / V")
plt.show()

# Gràfica del corrent AC, DC i total
plt.plot(t, i_ac_t, label="i ac")
plt.plot(t, i_dc_t, label="i dc")
plt.plot(t, i_t, label="i total")
plt.xlim(0, 200)
plt.ylim(-40, 80)
plt.grid()
plt.title("Corrent de curtcircuit")
plt.xlabel("t / ms")
plt.ylabel("i(t) / kA")
plt.legend(loc="upper right")
plt.show()
