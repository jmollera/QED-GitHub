# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 17:18:33 2019

@author: josep
"""

from jmb.utils import ComplexD
from jmb.ee import millman

U_AN = 230
U_BN = ComplexD(230, 240)
U_CN = ComplexD(230, 120)

R_A = 50e-3
R_B = 80e-3
R_C = 70e-3
R_N = 46

U_GN = millman([U_AN, U_BN, U_CN, 0], [R_A, R_B, R_C, R_N])
print(f"U_GN = {ComplexD(U_GN):.4f} V")

I_N = U_GN / R_N
print(f"I_N  = {ComplexD(I_N):.4f} A")

I_A = (U_AN - U_GN) / R_A
print(f"I_A  = {ComplexD(I_A):.4f} A")

I_B = (U_BN - U_GN) / R_B
print(f"I_B  = {ComplexD(I_B):.4f} A")

I_C = (U_CN - U_GN) / R_C
print(f"I_C  = {ComplexD(I_C):.4f} A")