# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 17:18:33 2019

@author: josep
"""

from jmb.utils import ComplexD
from jmb.ee import millman, Network


# 1r mètode de resolució

U_AN = 230
U_BN = ComplexD(230, 240)
U_CN = ComplexD(230, 120)

R_A = 50e-3
R_B = 80e-3
R_C = 70e-3

U_KN = millman([U_AN, U_BN, U_CN], [R_A, R_B, R_C])
print(f"U_KN = {ComplexD(U_KN):.4f} V")

I_A = (U_AN - U_KN) / R_A
print(f"I_A  = {ComplexD(I_A):.4f} A")

I_B = (U_BN - U_KN) / R_B
print(f"I_B  = {ComplexD(I_B):.4f} A")

I_C = (U_CN - U_KN) / R_C
print(f"I_C  = {ComplexD(I_C):.4f} A")

print()


# 2n mètode de resolució

net = Network()

net.add_voltage_source(1, 1, 0, 230, 50e-3)
net.add_voltage_source(2, 1, 0, ComplexD(230,240), 80e-3)
net.add_voltage_source(3, 1, 0, ComplexD(230,120), 70e-3)

net.solve()

branches = range(1, net.num_branches+1)
for branch in branches:
    print(f"Branch {branch} : Voltage = {ComplexD(net.branch_voltage(branch)):.4f} V")
print()
for branch in branches:
    print(f"Branch {branch} : Current = {ComplexD(net.branch_current(branch)):.4f} A")


