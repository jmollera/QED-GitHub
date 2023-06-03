# Importem primer els mòduls necessaris
import numpy as np
import matplotlib.pyplot as plt
from qed.utils import ComplexD
import qed.eng_elec as ee

# Cable trifàsic de 25 mm² de secció i 100 m de longitud
z=ee.z_cable(25, 100)
print(f'Zcable = {z:.4f} Ω/fase = {ComplexD(z):.4f} Ω/fase')

# Cable monofàsic de 25 mm² de secció i 100 m de longitud
z=ee.z_cable(25, 2*100)
print(f'Zcable = {z:.4f} Ω = {ComplexD(z):.4f} Ω')

# Cable de corrent continu de 25 mm² de secció i 100 m de longitud
z=ee.z_cable(25, 2*100).real
print(f'Zcable = {z:.4f} Ω')

# Sistema trifàsic: U=380 V fase-fase, I=200 A, cos φ=0.87. Cable: S=240 mm², L=400 m
cdt = ee.voltage_drop(380/np.sqrt(3), 200, ee.z_cable(240, 400), 0.87)
print(f'Caiguda de tensió = {cdt:.2f} V fase-neutre')

# Sistema monofàsic: U=220 V, I=55 A, cos φ=0.85. Cable: S=35 mm², L=200 m
cdt = ee.voltage_drop(220, 55, ee.z_cable(35, 2*200), 0.85)
print(f'Caiguda de tensió = {cdt:.2f} V')

# Sistema de corrent continu: U=125 V, I=25 A. Cable: S=16 mm², L=50 m
cdt = ee.voltage_drop(125, 25, ee.z_cable(16, 2*50).real)
print(f'Caiguda de tensió = {cdt:.2f} V')

# Càlcul d'impedàncies sèrie i paral.lel
print(ee.z_series([1+1j, 1+1j, 1+1j]))
print(ee.z_parallel([3+3j ,3+3j, 3+3j]))

# Resistència unitària d'un fil de 0,5 mm² de secció, a 90
for material in ["Al", "Cu", "Ag", "Au"]:
    res = ee.r_cable(0.5, 1, 90, material)
    print(f'Material = {material}, Resistència = {res:.4f} Ω/m')

# Valors eficaç i mitjà d'una ona senoidal semirectificada
T = 20e-3 # 20 ms
def f(t):
    if t < T/2:
        return 10*np.sin(2*np.pi/T*t)
    else:
        return 0
print(ee.rms(f, 0, T)) # 10/2
print(ee.average(f, 0, T))  # 10/π

# Valors eficaç i mitjà d'una ona senoidal rectificada
T = 20e-3 # 20 ms
def f(t):
    if t < T/2:
        return 10*np.sin(2*np.pi/T*t)
    else:
        return -10*np.sin(2*np.pi/T*t)
print(ee.rms(f, 0, T)) # 10/√2
print(ee.average(f, 0, T)) # 2*10/π

# Comparació de les corbes 51 'very inverse' CEI i IEEE
I_p = 100 # A
I_min = 1.05*I_p # A
I_max = 5*I_p # A
I = np.linspace(I_min, I_max, 200) # A
t_CEI = ee.CEI_51_curve(I, 'VI', I_p) # s
t_IEEE = ee.IEEE_51_curve(I, 'VI', I_p) # s
fig, ax = plt.subplots(figsize=(6, 6))
ax.grid(which='minor', alpha=0.2)
ax.grid(which='major', alpha=0.9)
ax.set(yscale="log")
ax.set_ylim(1, 1000)
ax.set_xlim(I_p, I_max)
ax.plot(I, t_CEI, label='CEI')
ax.plot(I, t_IEEE, label='IEEE')
ax.set_title("Corbes 51 'very inverse'")
ax.set_ylabel("Temps / s")
ax.set_xlabel("Corrent / A")
ax.legend(loc="upper right")
plt.savefig("Exemples.pdf", dpi=300, bbox_inches="tight")
plt.show()
