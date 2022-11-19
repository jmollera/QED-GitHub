# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 18:17:51 2021

@author: josep
"""

from jmb.utils import ComplexD
from jmb.ee import Network

net = Network()

net.add_voltage_source(1, 0, 1, 5, 4)
net.add_current_source(2, 0, 2, 2, 1e-6)
net.add_impedance(3, 1, 2, 5j)
net.add_impedance(4, 1, 0, -2j)

net.solve()

print(f"Corrent pel condensador: {ComplexD(net.branch_current(4)):.3f} A")


