import numpy as np
from scipy import optimize
import sympy
from sympy import sin as SIN
from sympy import cos as COS
import pandapower as pp

v2, δ2 = sympy.symbols('v2, δ2')

print('Cas a)')
print('-'*45)

f_symb  = [-0.8 - v2*(1.05*(-2.80*COS(-δ2) - 9.60*SIN(-δ2)) + v2*2.80),
           -0.6 + v2*(1.05*(-2.80*SIN(-δ2) + 9.60*COS(-δ2)) + v2*(-9.55))]

Jf_symb = sympy.Matrix(f_symb).jacobian(sympy.Matrix([v2, δ2]))

f_num = sympy.lambdify([(v2, δ2)], f_symb, 'numpy')
Jf_num = sympy.lambdify([(v2, δ2)], Jf_symb, 'numpy')
sol = optimize.root(f_num, x0=[1.05, 0], jac=Jf_num)
v_2, δ_2 = sol.x
print(f'v_2 = {v_2:.6f}')
print(f'δ_2 = {δ_2:.6f} rad = {np.degrees(δ_2):.6f}°\n')

s_12 = 1.05*((0.1j/2*1.05 + (1.05 - v_2*np.exp(δ_2*1j))/(0.028+0.096j))).conjugate()
s_G1 = 1.2+0.3j + s_12
print(f's_12 = {s_12:.5f}')
print(f's_G1 = {s_G1:.5f}')

net = pp.create_empty_network(sn_mva=1)

b1 = pp.create_bus(net, vn_kv=1, index=1)
b2 = pp.create_bus(net, vn_kv=1, index=2)

pp.create_ext_grid(net, b1, vm_pu=1.05, va_degree=0, index=1)
cap_nf = 0.1/(2*np.pi*50)*1e9
pp.create_line_from_parameters(net, b1, b2, length_km=1, r_ohm_per_km=0.028, x_ohm_per_km=0.096, c_nf_per_km=cap_nf, max_i_ka=1)
pp.create_load(net, b1, p_mw=1.2, q_mvar=0.3)
pp.create_load(net, b2, p_mw=0.8, q_mvar=0.6)

pp.runpp(net)

print('\nTensions i potències del nusos')
print(net.res_bus)
print('\nPotència injectada al nus 1 per la xarxa externa')
print(net.res_ext_grid)

print('\nCas b)')
print('-'*45)

f_symb = -0.8 - 1.03*(1.05*(-2.80*COS(-δ2) - 9.60*SIN(-δ2)) + 1.03*2.80)

Df_symb = sympy.diff(f_symb, δ2)

f_num = sympy.lambdify(δ2, f_symb, 'numpy')
Df_num = sympy.lambdify(δ2, Df_symb, 'numpy')
sol = optimize.root_scalar(f_num, x0=0, fprime=Df_num)
δ_2 = sol.root
print(f'δ_2 = {δ_2:.6f} rad = {np.degrees(δ_2):.6f}°\n')

s_12 = 1.05*((0.1j/2*1.05 + (1.05 - 1.03*np.exp(δ_2*1j))/(0.028+0.096j))).conjugate()
s_21 = 1.03*np.exp(δ_2*1j)*((0.1j/2*1.03*np.exp(δ_2*1j) + (1.03*np.exp(δ_2*1j) - 1.05)/(0.028+0.096j))).conjugate()
s_G1 = 1.2+0.3j + s_12
s_C2 = 0.8+0.6j + s_21
print(f's_12 = {s_12:.5f}')
print(f's_21 = {s_21:.5f}')
print(f's_G1 = {s_G1:.5f}')
print(f's_C2 = {s_C2:.5f}')

pp.create_gen(net, b2, p_mw=0, vm_pu=1.03, index=2) 

pp.runpp(net)

print('\nTensions i potències del nusos')
print(net.res_bus)
print('\nPotència injectada al nus 1 per la xarxa externa')
print(net.res_ext_grid)
print('\nPotència injectada al nus 2 pel banc de condensadors')
print(net.res_gen)
