# -*- coding: utf-8 -*-
"""
Created on Sat Nov  2 19:46:02 2019

@author: josep
"""

import numpy as np
from jmb.utils import ComplexD
import jmb.ee as EE

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
# Tensions fase-fase en V
U_AB, U_BC, U_CA = EE.triangle_to_phasors(399, 370, 370)
# Freqüència en Hz
f = 50

print("U_AB = {:.3f} V".format(ComplexD(U_AB)))
print("U_BC = {:.3f} V".format(ComplexD(U_BC)))
print("U_CA = {:.3f} V".format(ComplexD(U_CA)))
print()

U_AB_0, U_AB_1, U_AB_2 =  EE.ABC_to_A012(U_AB, U_BC, U_CA)
print("U_AB_0 = {:.3f} V".format(ComplexD(U_AB_0)))
print("U_AB_1 = {:.3f} V".format(ComplexD(U_AB_1)))
print("U_AB_2 = {:.3f} V".format(ComplexD(U_AB_2)))
print()

U_AN_1, U_AN_2 = EE.AB12_to_AN12(U_AB_1, U_AB_2)
print("U_AN_1 = {:.3f} V".format(ComplexD(U_AN_1)))
print("U_AN_2 = {:.3f} V".format(ComplexD(U_AN_2)))
print()

U_rel = ComplexD(U_AN_2).mod/ComplexD(U_AN_1).mod
print("U_AN_2/U_AN_1 = {:.3f}".format(U_rel))
print()

n_m_sinc = 120*f/p
ω_m_sinc = 4*np.pi*f/p
n_m_N = (1-s_N)*n_m_sinc
ω_m_N = (1-s_N)*ω_m_sinc
print("n_m,sinc = {:.0f} r/min,  ω_m,sinc = {:.3f} rad/s".format(n_m_sinc, ω_m_sinc))
print("n_m,N = {:.0f} r/min,  ω_m,N = {:.3f} rad/s".format(n_m_N, ω_m_N))
print()

Z_0 = R_Fe*1j*X_m/(R_Fe + 1j*X_m)
print("Z_0 = {:.3f} Ω".format(Z_0))
print()

def Z_mot_1(s):
    return R_1 + 1j*X_1 + Z_0*(R_2/s + 1j*X_2)/(Z_0 + R_2/s + 1j*X_2)

def I_1_1(s):
    return U_AN_1/Z_mot_1(s)

def Z_mot_2(s):
    return R_1 + 1j*X_1 + Z_0*(R_2/(2-s) + 1j*X_2)/(Z_0 + R_2/(2-s) + 1j*X_2)

def I_1_2(s):
    return U_AN_2/Z_mot_2(s)

Z_mot_N_1 = Z_mot_1(s=s_N)
Z_arr_1 = Z_mot_1(s=1)
Z_mot_N_2 = Z_mot_2(s=s_N)
print("Z_mot_N,1 = {0:.3f} Ω".format(Z_mot_N_1))
print("Z_arr,1 = {0:.3f} Ω".format(Z_arr_1))
print("Z_mot_N,2 = {0:.3f} Ω".format(Z_mot_N_2))
print()

I_N_1 = ComplexD(I_1_1(s=s_N))
I_N_2 = ComplexD(I_1_2(s=s_N))
print("I_N,1 = {0:.3f} A".format(I_N_1))
print("I_N,2 = {0:.3f} A".format(I_N_2))
print()

print("I_N,2/I_N,1 = {:.3f} %".format(100*I_N_2.mod/I_N_1.mod))
print()

I_N_A, I_N_B, I_N_C = EE.A012_to_ABC(0, I_N_1, I_N_2)
print("I_N,A = {0:.3f} A".format(ComplexD(I_N_A)))
print("I_N,B = {0:.3f} A".format(ComplexD(I_N_B)))
print("I_N,C = {0:.3f} A".format(ComplexD(I_N_C)))

print()