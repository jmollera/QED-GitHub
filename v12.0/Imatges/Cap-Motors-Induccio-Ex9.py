# -*- coding: utf-8 -*-
"""
Created on Thu Jan 6 14:34:10 2020

@author: josep
"""

import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
from jmb.utils import ComplexD

# Impedàncies del motor en Ω/fase referits a l'estator
R_1 = 0.5
X_1 = 1.5
R_2 = 0.625
X_2 = 1.25
R_Fe = 360
X_m = 40
# pols
p = 4
# Lliscament nominal
s_N = 0.05
# Sistema d'alimentació
U_sist = 380   # Tensió fase-fase en V
Scc_sist = 5   # Potència de curtcicuit en MVA
Rel_X_R = 9    # Relació X/R
# Cable de 6 mm²
L_cable = 80     # m
R_cable = 3.960  # Ω/km
X_cable = 0.123  # Ω/km


# Freqüència en Hz
f = 50


U_1_sist = U_sist/np.sqrt(3)
fase = np.arctan(Rel_X_R)
Z_sist = U_sist**2/(Scc_sist*1e6)*(np.cos(fase)+1j*np.sin(fase))
print("U_1,sist = {:.3f} V".format(U_1_sist))
print("Z_sist = {:.3f} Ω".format(Z_sist))
print()

Z_cable = L_cable*1e-3*(R_cable+1j*X_cable)
print("Z_cable = {:.3f} Ω".format(Z_cable))
print()

n_m_sinc = 120*f/p
ω_m_sinc=4*np.pi*f/p
print("n_m,sinc = {:.0f} r/min,  ω_m,sinc = {:.3f} rad/s".format(n_m_sinc, ω_m_sinc))
print()

Z_0 = R_Fe*1j*X_m/(R_Fe + 1j*X_m)
Z_Th = (Z_sist + Z_cable+R_1 + 1j*X_1)*Z_0/(Z_sist+Z_cable+R_1 + 1j*X_1 + Z_0)
E_Th = U_1_sist*Z_0/(Z_sist + Z_cable+R_1 + 1j*X_1 + Z_0)

print("Z_0 = {:.3f} Ω".format(Z_0))
print("Z_Th = {:.3f} Ω".format(Z_Th))
print("E_Th = {:.4f} V".format(ComplexD(E_Th)))
print()
    

def Z_mot(s):
    return R_1 + 1j*X_1 + Z_0*(R_2/s + 1j*X_2)/(Z_0 + R_2/s + 1j*X_2)

def U_1(s):
    return U_1_sist*Z_mot(s)/(Z_sist+Z_cable+Z_mot(s))

def I_1(s):
    return U_1(s)/Z_mot(s)

def I_2(s):
    return E_Th/(Z_Th + R_2/s + 1j*X_2)

def P_m(s):
    return 3*(1-s)*R_2*abs(I_2(s))**2/s

def T_m(s):
    return P_m(s)/((1-s)*ω_m_sinc)
    
def T_load(s):
     return 3.2*(1-s) + 58.9*(1-s)**2

def equation(s):
     return T_m(s) - T_load(s)
 
s_T_m_max = R_2/abs(Z_Th+1j*X_2)
n_T_m_max = (1-s_T_m_max)*n_m_sinc
ω_T_m_max = (1-s_T_m_max)*ω_m_sinc
T_m_max = T_m(s=s_T_m_max)
print("s_T_m,max = {0:.4f},  n_T_m,max = {1:.2f} r/min,  ω_T_m,max = {2:.2f} rad/s".format(s_T_m_max, n_T_m_max, ω_T_m_max))
print("T_m,max = {0:.1f} N·m".format(T_m_max))
print()
    
sol = optimize.root(equation, x0=s_N)
s = sol.x[0]
zmot = Z_mot(s)
I = abs(I_1(s))
U1 = U_1(s)
Tm = T_m(s)
cos_φ = np.cos(np.angle(Z_mot(s)))
η = (1-s)/s*R_2*abs(I_2(s))**2/abs(U_1(s))/I/cos_φ
print ("s = {0:.3f}".format(s))
print ("Z_mot = {0:.3f} Ω".format(zmot))
print ("I = {0:.2f} A".format(I))
print("U_1 = {:.4f} V".format(ComplexD(U1)))
print ("Tm = {0:.1f} N m".format(Tm))
print ("cos φ = {0:.3f}".format(cos_φ))
print ("η = {0:.3f}".format(η))
print()

Iarr = abs(I_1(1))
Tarr = T_m(1.0001)
print ("Iarr = {0:.1f} A".format(Iarr))
print ("Tm,arr = {0:.1f} N m".format(Tarr))
print()


punts = 300
n_list = np.linspace(0.001, n_m_sinc-0.001, punts)  # r/min
T_load_list = T_load(s=1-n_list/n_m_sinc)  # N·m
T_m_list = T_m(s=1-n_list/n_m_sinc)  # N·m
I_list = abs(I_1(s=1-n_list/n_m_sinc))  # A
U1_list = abs(U_1(s=1-n_list/n_m_sinc))  # V

try:
    with open("Cap-Motors-Induccio-Ex9.txt", "w", encoding="utf8") as fh:
        for i in range(len(n_list)):
            fh.write("{0:8.3f} {1:7.3f} {2:7.3f} {3:6.3f}  {4:7.3f}\n".format(n_list[i], T_load_list[i], T_m_list[i], I_list[i], U1_list[i]))
except OSError as err:
    print(err)
    
# Gràfica dels parells
plt.plot(n_list, T_m_list, label="Tm")
plt.plot(n_list, T_load_list, label="Tload")
plt.xlim(0,1550)
plt.ylim(0,140)
plt.xticks(np.arange(0, 1600, 100))
plt.grid()
plt.title("Parell")
plt.xlabel("n / r/min")
plt.ylabel("T / N·m")
plt.legend(loc="upper left")
plt.show()

# Gràfica dels corrents
plt.plot(n_list, I_list, label="I")
plt.xlim(0,1550)
plt.ylim(0,80)
plt.xticks(np.arange(0, 1600, 100))
plt.grid()
plt.title("Corrent")
plt.xlabel("n / r/min")
plt.ylabel("I / A")
plt.legend(loc="lower left")
plt.show()

# Gràfica de la tensió
plt.plot(n_list, U1_list, label="U1")
plt.plot((0,1500), (U_1_sist,U_1_sist), label="U1_sist")
plt.xlim(0,1550)
plt.ylim(204,220)
plt.xticks(np.arange(0, 1600, 100))
plt.grid()
plt.title("Tensió")
plt.xlabel("n / r/min")
plt.ylabel("U1 / V")
plt.legend(loc="lower left")
plt.show()
