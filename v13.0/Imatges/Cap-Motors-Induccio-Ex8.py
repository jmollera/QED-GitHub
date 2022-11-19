# -*- coding: utf-8 -*-
"""
Created on Thu Jan 2 16:34:10 2020

@author: josep
"""

import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt

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
# Tensió fase-fase en V
U = 380
# Freqüència en Hz
f = 50

U_1 = U/np.sqrt(3)
print("100% U_1 = {:.3f} V".format(U_1))
print("80% U_1 = {:.3f} V".format(0.8*U_1))
print()

n_m_sinc = 120*f/p
ω_m_sinc=4*np.pi*f/p
print("n_m,sinc = {:.0f} r/min,  ω_m,sinc = {:.3f} rad/s".format(n_m_sinc, ω_m_sinc))
print()

Z_0 = R_Fe*1j*X_m/(R_Fe + 1j*X_m)
Z_Th = (R_1 + 1j*X_1)*Z_0/(R_1 + 1j*X_1 + Z_0)

print("Z_0 = {:.3f} Ω".format(Z_0))
print("Z_Th = {:.3f} Ω".format(Z_Th))
print()

def Z_mot(s):
    return R_1 + 1j*X_1 + Z_0*(R_2/s + 1j*X_2)/(Z_0 + R_2/s + 1j*X_2)

def I_1(s, U1):
    return U1/Z_mot(s)

def I_2(s, U1):
    UTh = U1*Z_0/(R_1 + 1j*X_1 + Z_0)
    return UTh/(Z_Th + R_2/s + 1j*X_2)

def P_m(s, U1):
    return 3*(1-s)*R_2*abs(I_2(s, U1))**2/s

def T_m(s, U1):
    return P_m(s, U1)/((1-s)*ω_m_sinc)
    
def T_load(s):
    return 3.2*(1-s) + 58.9*(1-s)**2

def equation(s, U1):
    return T_m(s, U1) - T_load(s)
 
    
sol = optimize.root(equation, x0=s_N, args=(U_1))
s = sol.x[0]
zmot = Z_mot(s)
I = abs(I_1(s, U_1))
Iarr = abs(I_1(1, U_1))
Tm = T_m(s, U_1)
Tarr = T_m(1.0001, U_1)
cos_φ = np.cos(np.angle(Z_mot(s)))
η = (1-s)/s*R_2*abs(I_2(s,U_1))**2/U_1/I/cos_φ
print ("U = 100 %")
print ("s = {0:.3f}".format(s))
print ("Z_mot = {0:.3f} Ω".format(zmot))
print ("I = {0:.1f} A".format(I))
print ("Iarr = {0:.1f} A".format(Iarr))
print ("Tm = {0:.1f} N m".format(Tm))
print ("Tm,arr = {0:.1f} N m".format(Tarr))
print ("cos φ = {0:.3f}".format(cos_φ))
print ("η = {0:.3f}".format(η))
print()

sol = optimize.root(equation, x0=s_N, args=(0.9*U_1))
s = sol.x[0]
zmot = Z_mot(s)
I = abs(I_1(s, 0.9*U_1))
Iarr = abs(I_1(1, 0.9*U_1))
Tm = T_m(s, 0.9*U_1)
Tarr = T_m(1.0001, 0.9*U_1)
cos_φ = np.cos(np.angle(Z_mot(s)))
η = (1-s)/s*R_2*abs(I_2(s,0.9*U_1))**2/(0.9*U_1)/I/cos_φ
print ("U = 90 %")
print ("s = {0:.3f}".format(s))
print ("Z_mot = {0:.3f} Ω".format(zmot))
print ("I = {0:.1f} A".format(I))
print ("Iarr = {0:.1f} A".format(Iarr))
print ("Tm = {0:.1f} N m".format(Tm))
print ("Tm,arr = {0:.1f} N m".format(Tarr))
print ("cos φ = {0:.3f}".format(cos_φ))
print ("η = {0:.3f}".format(η))
print()

sol = optimize.root(equation, x0=s_N, args=(0.8*U_1))
s = sol.x[0]
zmot = Z_mot(s)
I = abs(I_1(s, 0.8*U_1))
Iarr = abs(I_1(1, 0.8*U_1))
Tm = T_m(s, 0.8*U_1)
Tarr = T_m(1.0001, 0.8*U_1)
cos_φ = np.cos(np.angle(Z_mot(s)))
η = (1-s)/s*R_2*abs(I_2(s,0.8*U_1))**2/(0.8*U_1)/I/cos_φ
print ("U = 80 %")
print ("s = {0:.3f}".format(s))
print ("Z_mot = {0:.3f} Ω".format(zmot))
print ("I = {0:.1f} A".format(I))
print ("Iarr = {0:.1f} A".format(Iarr))
print ("Tm = {0:.1f} N m".format(Tm))
print ("Tm,arr = {0:.1f} N m".format(Tarr))
print ("cos φ = {0:.3f}".format(cos_φ))
print ("η = {0:.3f}".format(η))
print()

sol = optimize.root(equation, x0=s_N, args=(0.75*U_1))
s = sol.x[0]
zmot = Z_mot(s)
I = abs(I_1(s, 0.75*U_1))
Iarr = abs(I_1(1, 0.75*U_1))
Tm = T_m(s, 0.75*U_1)
Tarr = T_m(1.0001, 0.75*U_1)
cos_φ = np.cos(np.angle(Z_mot(s)))
η = (1-s)/s*R_2*abs(I_2(s,0.75*U_1))**2/(0.75*U_1)/I/cos_φ
print ("U = 75 %")
print ("s = {0:.3f}".format(s))
print ("Z_mot = {0:.3f} Ω".format(zmot))
print ("I = {0:.1f} A".format(I))
print ("Iarr = {0:.1f} A".format(Iarr))
print ("Tm = {0:.1f} N m".format(Tm))
print ("Tm,arr = {0:.1f} N m".format(Tarr))
print ("cos φ = {0:.3f}".format(cos_φ))
print ("η = {0:.3f}".format(η))
print()


punts = 300
n_list = np.linspace(0.001, n_m_sinc-0.001, punts)  # r/min
T_load_list = T_load(s=1-n_list/n_m_sinc)  # N·m
T_m_list_100 = T_m(s=1-n_list/n_m_sinc, U1=U_1)  # N·m
T_m_list_80 = T_m(s=1-n_list/n_m_sinc, U1=0.8*U_1)  # N·m
I_list_100 = abs(I_1(s=1-n_list/n_m_sinc, U1=U_1))  # A
I_list_80 = abs(I_1(s=1-n_list/n_m_sinc, U1=0.8*U_1))  # A

try:
    with open("Cap-Motors-Induccio-Ex8-1.txt", "w", encoding="utf8") as fh:
        for i in range(len(n_list)):
            fh.write("{0:8.3f} {1:7.3f} {2:7.3f} {3:7.3f} {4:6.3f} {5:6.3f}\n".format(n_list[i], T_load_list[i], T_m_list_100[i], T_m_list_80[i], I_list_100[i], I_list_80[i]))
except OSError as err:
    print(err)
    
# Gràfica dels parells
plt.plot(n_list, T_m_list_100, label="Tm 100%")
plt.plot(n_list, T_m_list_80, label="Tm 800%")
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
plt.plot(n_list, I_list_100, label="I 100 %")
plt.plot(n_list, I_list_80, label="I 80 %")
plt.xlim(0,1550)
plt.ylim(0,80)
plt.xticks(np.arange(0, 1600, 100))
plt.grid()
plt.title("Corrent")
plt.xlabel("n / r/min")
plt.ylabel("I / A")
plt.legend(loc="lower left")
plt.show()

punts = 200
U1_list = np.linspace(0.75, 1, punts) # 75% a 100%
s_list = np.empty(punts)
for i, pu in enumerate(U1_list):
    sol = optimize.root(equation, x0=s_N, args=(pu*U_1))
    s_list[i] = sol.x[0]
I_list = abs(I_1(s_list, U1_list*U_1))
cos_φ_list = np.cos(np.angle(Z_mot(s_list)))
η_list = (1-s_list)/s_list*R_2*abs(I_2(s_list,U1_list*U_1))**2/(U1_list*U_1)/I_list/cos_φ_list

try:
    with open("Cap-Motors-Induccio-Ex8-2.txt", "w", encoding="utf8") as fh:
        for i in range(punts):
            fh.write("{0:8.6f} {1:8.6f} {2:6.3f} {3:8.6f} {4:8.6f}\n".format(U1_list[i], s_list[i], I_list[i], cos_φ_list[i], η_list[i]))
except OSError as err:
    print(err)
    
# Gràfica del corrent
plt.plot(U1_list, I_list)
plt.xlim(0.75,1)
plt.xticks(np.arange(0.75, 1, 0.05))
plt.grid()
plt.title("Corrent")
plt.xlabel("U1/U1_N")
plt.ylabel("I / A")
plt.show()

# Gràfica del lliscament
plt.plot(U1_list, s_list)
plt.xlim(0.75,1)
plt.xticks(np.arange(0.75, 1, 0.05))
plt.grid()
plt.title("Lliscament")
plt.xlabel("U1/U1_N")
plt.ylabel("s")
plt.show()

# Gràfica del factor de potència i rendiement
plt.plot(U1_list, cos_φ_list, label="cos φ")
plt.plot(U1_list, η_list, label="η")
plt.xlim(0.75,1)
plt.xticks(np.arange(0.75, 1, 0.05))
plt.grid()
plt.title("Factor de potència i Rendiement")
plt.xlabel("U1/U1_N")
plt.ylabel("cos φ, η")
plt.legend(loc="lower right")
plt.show()