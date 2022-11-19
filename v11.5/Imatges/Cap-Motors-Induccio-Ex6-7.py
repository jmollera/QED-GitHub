# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 21:28:10 2019

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
print("U_1 = {:.3f} V".format(U_1))
print()

n_m_sinc = 120*f/p
ω_m_sinc = 4*np.pi*f/p
n_m_N = (1-s_N)*n_m_sinc
ω_m_N = (1-s_N)*ω_m_sinc
print("n_m,sinc = {:.0f} r/min,  ω_m,sinc = {:.3f} rad/s".format(n_m_sinc, ω_m_sinc))
print("n_m,N = {:.0f} r/min,  ω_m,N = {:.3f} rad/s".format(n_m_N, ω_m_N))
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


I_N = abs(I_1(s_N, U_1))
print("I_N = {0:.1f} A".format(I_N))
print()

P_m_N = P_m(s_N, U_1)
print("P_m,N = {0:.1f} W".format(P_m_N))
print()


def equation(s, Pm, U1):
    return 3*(1-s)/s*R_2*U1**2*(abs(Z_0)/(abs(R_1+1j*X_1+Z_0)*abs(Z_Th+R_2/s+1j*X_2)))**2-Pm
 
    
sol = optimize.root(equation, x0=s_N, args=(P_m_N, U_1))
s = sol.x[0]
zmot = Z_mot(s)
I = abs(I_1(s, U_1))
print ("U = 100 %")
print ("s = {0:.3f}".format(s))
print ("Z_mot = {0:.3f} Ω".format(zmot))
print ("I = {0:.1f} A".format(I))
print ("I/I_N = {0:.1f}".format(I/I_N))
print()

sol = optimize.root(equation, x0=s_N, args=(P_m_N, 0.9*U_1))
s = sol.x[0]
zmot = Z_mot(s)
I = abs(I_1(s, 0.9*U_1))
print ("U = 90 %")
print ("s = {0:.3f}".format(s))
print ("Z_mot = {0:.3f} Ω".format(zmot))
print ("I = {0:.1f} A".format(I))
print ("I/I_N = {0:.1f}".format(I/I_N))
print()

sol = optimize.root(equation, x0=s_N, args=(P_m_N, 0.8*U_1))
s = sol.x[0]
zmot = Z_mot(s)
I = abs(I_1(s, 0.8*U_1))
print ("U = 80 %")
print ("s = {0:.3f}".format(s))
print ("Z_mot = {0:.3f} Ω".format(zmot))
print ("I = {0:.1f} A".format(I))
print ("I/I_N = {0:.1f}".format(I/I_N))
print()

sol = optimize.root(equation, x0=s_N, args=(P_m_N, 0.75*U_1))
s = sol.x[0]
zmot = Z_mot(s)
I = abs(I_1(s, 0.75*U_1))
print ("U = 75 %")
print ("s = {0:.3f}".format(s))
print ("Z_mot = {0:.3f} Ω".format(zmot))
print ("I = {0:.1f} A".format(I))
print ("I/I_N = {0:.1f}".format(I/I_N))
print()

punts = 200
U1_list = np.linspace(75, 100, punts) # 75% a 100%
s_list = np.empty(punts)
for i, percent in enumerate(U1_list):
    sol = optimize.root(equation, x0=s_N, args=(P_m_N, percent/100*U_1))
    s_list[i] = sol.x[0]
I_list = abs(I_1(s_list, U1_list/100*U_1))
try:
    with open("Cap-Motors-Induccio-Ex6.txt", "w", encoding="utf8") as fh:
        for i in range(punts):
            fh.write("{0:8.3f} {1:8.6f} {2:6.3f}\n".format(U1_list[i], s_list[i], I_list[i]))
except Exception as err:
    print(err)
    
# Gràfica del corrent
plt.plot(U1_list, I_list)
plt.xlim(75,100)
plt.xticks(np.arange(75, 101, 5))
plt.grid()
plt.title("Corrent")
plt.xlabel("U1/U1_N / %")
plt.ylabel("I / A")
plt.show()

# Gràfica del lliscament
plt.plot(U1_list, s_list)
plt.xlim(75,100)
plt.xticks(np.arange(75, 101, 5))
plt.grid()
plt.title("Lliscament")
plt.xlabel("U1/U1_N / %")
plt.ylabel("s")
plt.show()


sol = optimize.root(equation, x0=s_N, args=(P_m_N, U_1))
s = sol.x[0]
zmot = Z_mot(s)
I = abs(I_1(s, U_1))
cos_φ = np.cos(np.angle(Z_mot(s)))
η = (1-s)/s*R_2*abs(I_2(s,U_1))**2/U_1/I/cos_φ
print ("P_m = 100 %")
print ("s = {0:.3f}".format(s))
print ("Z_mot = {0:.3f} Ω".format(zmot))
print ("I = {0:.1f} A".format(I))
print ("I/I_N = {0:.1f}".format(I/I_N))
print ("cos φ = {0:.3f}".format(cos_φ))
print ("η = {0:.3f}".format(η))
print()

sol = optimize.root(equation, x0=s_N, args=(0.75*P_m_N, U_1))
s = sol.x[0]
zmot = Z_mot(s)
I = abs(I_1(s, U_1))
cos_φ = np.cos(np.angle(Z_mot(s)))
η = (1-s)/s*R_2*abs(I_2(s,U_1))**2/U_1/I/cos_φ
print ("P_m = 75 %")
print ("s = {0:.3f}".format(s))
print ("Z_mot = {0:.3f} Ω".format(zmot))
print ("I = {0:.1f} A".format(I))
print ("I/I_N = {0:.1f}".format(I/I_N))
print ("cos φ = {0:.3f}".format(cos_φ))
print ("η = {0:.3f}".format(η))
print()

sol = optimize.root(equation, x0=s_N, args=(0.5*P_m_N, U_1))
s = sol.x[0]
zmot = Z_mot(s)
I = abs(I_1(s, U_1))
cos_φ = np.cos(np.angle(Z_mot(s)))
η = (1-s)/s*R_2*abs(I_2(s,U_1))**2/U_1/I/cos_φ
print ("P_m = 50 %")
print ("s = {0:.3f}".format(s))
print ("Z_mot = {0:.3f} Ω".format(zmot))
print ("I = {0:.1f} A".format(I))
print ("I/I_N = {0:.1f}".format(I/I_N))
print ("cos φ = {0:.3f}".format(cos_φ))
print ("η = {0:.3f}".format(η))
print()

sol = optimize.root(equation, x0=s_N, args=(0.25*P_m_N, U_1))
s = sol.x[0]
zmot = Z_mot(s)
I = abs(I_1(s, U_1))
cos_φ = np.cos(np.angle(Z_mot(s)))
η = (1-s)/s*R_2*abs(I_2(s,U_1))**2/U_1/I/cos_φ
print ("P_m = 25 %")
print ("s = {0:.3f}".format(s))
print ("Z_mot = {0:.3f} Ω".format(zmot))
print ("I = {0:.1f} A".format(I))
print ("I/I_N = {0:.1f}".format(I/I_N))
print ("cos φ = {0:.3f}".format(cos_φ))
print ("η = {0:.3f}".format(η))
print()

punts = 200
Pm_list = np.linspace(0.1, 100, punts) # 0.1% a 100%
s_list = np.empty(punts)
for i, percent in enumerate(Pm_list):
    sol = optimize.root(equation, x0=s_N, args=(percent/100*P_m_N, U_1))
    s_list[i] = sol.x[0]
I_list = abs(I_1(s_list, U_1))
cos_φ_list = np.cos(np.angle(Z_mot(s_list)))
η_list = (1-s_list)/s_list*R_2*abs(I_2(s_list,U_1))**2/U_1/I_list/cos_φ_list
try:
    with open("Cap-Motors-Induccio-Ex7.txt", "w", encoding="utf8") as fh:
        for i in range(punts):
            fh.write("{0:8.3f} {1:8.6f} {2:6.3f} {3:8.6f} {4:8.6f}\n".format(Pm_list[i], s_list[i], I_list[i], cos_φ_list[i], η_list[i]))
except Exception as err:
    print(err)
    
# Gràfica del corrent
plt.plot(Pm_list, I_list)
plt.xlim(0,105)
plt.xticks(np.arange(0, 105, 5))
plt.ylim(5,18)
plt.yticks(np.arange(5, 19, 1))
plt.grid()
plt.title("Corrent")
plt.xlabel("Pm/Pm_N / %")
plt.ylabel("I / A")
plt.show()

# Gràfica del lliscament
plt.plot(Pm_list, s_list)
plt.xlim(0,105)
plt.xticks(np.arange(0, 105, 5))
plt.grid()
plt.title("Lliscament")
plt.xlabel("Pm/Pm_N / %")
plt.ylabel("s")
plt.show()

# Gràfica del factor de potència i rendiement
plt.plot(Pm_list, cos_φ_list, label="cos φ")
plt.plot(Pm_list, η_list, label="η")
plt.xlim(0,105)
plt.xticks(np.arange(0, 105, 5))
plt.grid()
plt.title("Factor de potència i Rendiement")
plt.xlabel("Pm/Pm_N / %")
plt.ylabel("cos φ, η")
plt.legend(loc="upper left")
plt.show()