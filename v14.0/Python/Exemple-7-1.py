import numpy as np
import qed.eng_elec as ee

U = 380  # V fase-faase
I = 630  # A
cos_φ = 0.87  # Factor de potència
L = 400  # m
n = 3  # cables en paraŀlel
R = 0.095  # Ω/km
X = 0.102  # Ω/km

U_FN = U/np.sqrt(3)
Z = (R + X*1j)*L/1000/n
Δ_U = ee.voltage_drop(U_FN, I, Z, cos_φ)
print(f'Caiguda de tensió: {Δ_U:.2f} V ({Δ_U/U_FN*100:.2f} %)')