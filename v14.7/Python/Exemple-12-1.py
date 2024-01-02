import numpy as np
from scipy import optimize
import sympy
from sympy import sin as SIN
from sympy import cos as COS
import pandapower as pp

δ2, v3, δ3 = sympy.symbols('δ2, v3, δ3')

f_symb = [0.25 - 0.5 - 1.03*(1.05*(-2.8*COS(-δ2) - 9.6*SIN(-δ2)) + 1.03*7.8 + v3*(-5*COS(δ3-δ2)-15*SIN(δ3-δ2))),
          - 0.6 - v3*(1.05*(-10*COS(-δ3) - 30*SIN(-δ3)) + 1.03*(-5*COS(δ2-δ3)-15*SIN(δ2-δ3)) + v3*15),
          - 0.3 + v3*(1.05*(-10*SIN(-δ3) + 30*COS(-δ3)) + 1.03*(-5*SIN(δ2-δ3)+15*COS(δ2-δ3)) + v3*(-45))]

Jf_symb = sympy.Matrix(f_symb).jacobian(sympy.Matrix([δ2, v3, δ3]))

f_num = sympy.lambdify([(δ2, v3, δ3)], f_symb, 'numpy')
Jf_num = sympy.lambdify([(δ2, v3, δ3)], Jf_symb, 'numpy')
sol = optimize.root(f_num, x0=[0, 1.05, 0], jac=Jf_num)
δ_2, v_3, δ_3 = sol.x
print(f'δ_2 = {δ_2:.6f} rad = {np.degrees(δ_2):.6f}°')
print(f'v_3 = {v_3:.6f}')
print(f'δ_3 = {δ_3:.6f} rad = {np.degrees(δ_3):.6f}°\n')

s_1 = (1.05*((12.8-39.6j)*1.05 + (-2.8+9.6j)*1.03*np.exp(δ_2*1j) + (-10+30j)*v_3*np.exp(δ_3*1j))).conjugate()
s_2 = (1.03*np.exp(-δ_2*1j)*((-2.8+9.6j)*1.05 + (7.8-24.6j)*1.03*np.exp(δ_2*1j) + (-5+15j)*v_3*np.exp(δ_3*1j))).conjugate()
s_G1 = s_1
s_G2 = 0.5+0.25j + s_2
print(f's_1 = {s_1:.5f}')
print(f's_2 = {s_2:.5f}')
print(f's_G1 = {s_G1:.5f}')
print(f's_G2 = {s_G2:.5f}')

net = pp.create_empty_network(sn_mva=1)

b1 = pp.create_bus(net, vn_kv=1, index=1)
b2 = pp.create_bus(net, vn_kv=1, index=2)
b3 = pp.create_bus(net, vn_kv=1, index=3)

pp.create_ext_grid(net, b1, vm_pu=1.05, va_degree=0, index=1)
pp.create_gen(net, b2, p_mw=0.25, vm_pu=1.03, index=2) 
pp.create_line_from_parameters(net, b1, b2, length_km=1, r_ohm_per_km=0.028, x_ohm_per_km=0.096, c_nf_per_km=0, max_i_ka=1)
pp.create_line_from_parameters(net, b1, b3, length_km=1, r_ohm_per_km=0.01, x_ohm_per_km=0.03, c_nf_per_km=0, max_i_ka=1)
pp.create_line_from_parameters(net, b2, b3, length_km=1, r_ohm_per_km=0.02, x_ohm_per_km=0.06, c_nf_per_km=0, max_i_ka=1)
pp.create_load(net, b2, p_mw=0.5, q_mvar=0.25)
pp.create_load(net, b3, p_mw=0.6, q_mvar=0.3)

pp.runpp(net)

print('\nTensions i potències del nusos')
print(net.res_bus)
print('\nPotència generada en el nus 1')
print(net.res_ext_grid)
print('\nPotència generada en el nus 2')
print(net.res_gen)
