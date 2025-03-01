# -*- coding: utf-8 -*-

import numpy as np
from qed.utils import Complex

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
U_Th = U_1*Z_0/(R_1 + 1j*X_1 + Z_0)
print("Z_0 = {:.3f} Ω".format(Z_0))
print("Z_Th = {:.3f} Ω".format(Z_Th))
print("U_Th = {:.4f} V".format(Complex(U_Th)))
print()

def Z_mot(s):
    return R_1 + 1j*X_1 + Z_0*(R_2/s + 1j*X_2)/(Z_0 + R_2/s + 1j*X_2)

def I_1(s):
    return U_1/Z_mot(s)

def I_2(s):
    return U_Th/(Z_Th + R_2/s + 1j*X_2)

def P_m(s):
    return 3*(1-s)*R_2*abs(I_2(s))**2/s

def T_m(s):
    return 3*R_2*abs(I_2(s))**2/(s*ω_m_sinc)

I_N = abs(I_1(s=s_N))
Z_arr = Z_mot(s=1)
I_arr = abs(I_1(s=1))
kappa = 1.02 + 0.98*np.exp(-3*Z_arr.real/Z_arr.imag)
I_arr_pic_asim = kappa*np.sqrt(2)*I_arr
print("I_N = {0:.1f} A".format(I_N))
print("I_arr = {0:.1f} A,  I_arr/I_N = {1:.1f}".format(I_arr, I_arr/I_N))
print("I_arr,pic,asim = {0:.2f} A,  I_arr,pic,asim/I_arr = {1:.2f}".format(I_arr_pic_asim, I_arr_pic_asim/I_arr))
print()

P_m_N = P_m(s=s_N)
print("P_m,N = {0:.1f} W".format(P_m_N))
print()

T_m_N = T_m(s=s_N)
T_m_arr = T_m(s=1)
print("T_m,N = {0:.1f} N·m".format(T_m_N))
print("T_m,arr = {0:.1f} N·m,  T_m,arr/T_m,N = {1:.2f}".format(T_m_arr, T_m_arr/T_m_N))
print()

s_T_m_max = R_2/abs(Z_Th+1j*X_2)
n_T_m_max = (1-s_T_m_max)*n_m_sinc
ω_T_m_max = (1-s_T_m_max)*ω_m_sinc
T_m_max = T_m(s=s_T_m_max)
print("s_T_m,max = {0:.4f},  n_T_m,max = {1:.2f} r/min,  ω_T_m,max = {2:.2f} rad/s".format(s_T_m_max, n_T_m_max, ω_T_m_max))
print("T_m,max = {0:.1f} N·m,  T_m,max/T_m,N = {1:.2f}".format(T_m_max, T_m_max/T_m_N))
print()

cos_φ_N = np.cos(np.angle(Z_mot(s=s_N)))
cos_φ_arr = np.cos(np.angle(Z_mot(s=1)))
print("cos φ,N = {0:.2f}".format(cos_φ_N))
print("cos φ,arr = {0:.2f}".format(cos_φ_arr))
print()

Pm_N = P_m(s=s_N)/1000
Pm_arr = P_m(s=1)/1000
print("P_m,N = {0:.3f} kW".format(Pm_N))
print("P_m,arr = {0:.3f} kW".format(Pm_arr))
print()

P_N = 3*U_1*I_N*cos_φ_N/1000
P_arr = 3*U_1*I_arr*cos_φ_arr/1000
print("P,N = {0:.3f} kW".format(P_N))
print("P,arr = {0:.3f} kW".format(P_arr))
print()

Q_N = 3*U_1*I_N*np.sqrt(1-cos_φ_N**2)/1000
Q_arr = 3*U_1*I_arr*np.sqrt(1-cos_φ_arr**2)/1000
print("Q,N = {0:.3f} kvar".format(Q_N))
print("Q,arr = {0:.3f} kvar".format(Q_arr))
print()

η_N = Pm_N/P_N
η_arr = Pm_arr/P_arr
print("η,N = {0:.2f}".format(η_N))
print("η,arr = {0:.2f}".format(η_arr))
print()


n = np.linspace(0, n_m_sinc-0.0001, 300) # r/min
I = abs(I_1(s=1-n/n_m_sinc)) # A
Tm = T_m(s=1-n/n_m_sinc) # N·m
Pm = P_m(s=1-n/n_m_sinc)/1000 # kW
cos_φ = np.cos(np.angle(Z_mot(s=1-n/n_m_sinc)))
P = 3*U_1*I*cos_φ/1000 # kW
Q = 3*U_1*I*np.sqrt(1-cos_φ**2)/1000 # kvar
S = np.hypot(P,Q)
η = Pm/P
try:
    with open("Cap-Motors-Induccio-CaractMotor-Aux.txt", "w", encoding="utf8") as fh:
        for i in range(len(n)):
            fh.write("{0:8.3f} {1:7.3f} {2:6.3f} {3:6.3f} {4:6.3f} {5:6.3f} {6:6.3f} {7:8.6f} {8:8.6f}\n".format(n[i], Tm[i], I[i], Pm[i], P[i], Q[i], S[i], cos_φ[i], η[i]))
except OSError as err:
    print(err)
