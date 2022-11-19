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

U_KN = millman([U_AN, U_BN, U_CN], [R_A, R_B, R_C])
print("U_KN = {:.4f} V".format(ComplexD(U_KN)))

I_A = (U_AN - U_KN) / R_A
print("I_A  = {:.4f} A".format(ComplexD(I_A)))

I_B = (U_BN - U_KN) / R_B
print("I_B  = {:.4f} A".format(ComplexD(I_B)))

I_C = (U_CN - U_KN) / R_C
print("I_C  = {:.4f} A".format(ComplexD(I_C)))
