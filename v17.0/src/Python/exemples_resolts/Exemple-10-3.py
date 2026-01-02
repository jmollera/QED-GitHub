import numpy as np
from scipy import integrate
from scipy import constants

J = 507  # lb·ft²

n_m = np.array([0, 90, 180, 270, 360, 450, 540, 630, 720, 810, 900,
                990, 1080, 1170, 1260, 1350, 1440, 1530, 1620, 1656,
                1692, 1728, 1764, 1780], dtype=np.double)  # r/min

T_m = np.array([266, 248, 234, 225, 221, 225, 234, 246, 258, 273,
                289, 304, 322, 338, 356, 373, 388, 396, 395, 388,
                370, 336, 248, 177], dtype=np.double)  # lbf·ft

T_load = np.array([0.0, 0.4, 1.7, 3.8, 6.8, 10.6, 15.2, 20.7, 27.0,
                   34.2, 42.2, 51.1, 60.8, 71.4, 82.8, 95.0, 108.1,
                   122.0, 136.8, 143.0, 149.2, 155.7, 162.2, 165.2], dtype=np.double)  # lbf·ft

J *= constants.lb * constants.foot**2  # kg·m²
ω_m = n_m*np.pi/30  # rad/s
T_acc_inv = 1/((T_m - T_load) * constants.lbf * constants.foot)  # 1/(N·m)

t_arr = J*integrate.trapezoid(y=T_acc_inv, x=ω_m)
print(f"Temps d'arrencada (trapezoid) = {t_arr:.1f} s")
t_arr = J*integrate.simpson(y=T_acc_inv, x=ω_m)
print(f"Temps d'arrencada (simpson)  = {t_arr:.1f} s")
