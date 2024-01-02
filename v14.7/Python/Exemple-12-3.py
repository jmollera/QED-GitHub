from scipy import optimize
import numpy as np
import sympy
import pandapower as pp
import pandapower.control as pc

m, δ2 = sympy.symbols('m, δ2')

f_symb = [-2 - (-50/m*sympy.sin(-δ2)),
          -1 + (50/m*sympy.cos(-δ2) - 50)]

Jf_symb = sympy.Matrix(f_symb).jacobian(sympy.Matrix([δ2, m]))

f_num = sympy.lambdify([(δ2, m)], f_symb, 'numpy')
Jf_num = sympy.lambdify([(δ2, m)], Jf_symb, 'numpy')
sol = optimize.root(f_num, x0=[0, 1], jac=Jf_num)
δ_2, m_ = sol.x
print(f'δ_2 = {δ_2:.6f} rad = {np.degrees(δ_2):.6f}°')
print(f'm = {m_:.6f}')

net = pp.create_empty_network(sn_mva=1)

b1 = pp.create_bus(net, vn_kv=1, index=1)
b2 = pp.create_bus(net, vn_kv=1, index=2)

pp.create_ext_grid(net, b1, vm_pu=1, va_degree=0, index=1)
tr1 = pp.create_transformer_from_parameters(net, b1, b2, sn_mva=1, vn_hv_kv=1, vn_lv_kv=1, vkr_percent=0, vk_percent=2, pfe_kw=0, i0_percent=0, tap_side='hv', tap_step_percent=1, tap_min=-10, tap_max=10, tap_neutral=0, tap_pos=0)
pp.create_load(net, b2, p_mw=2, q_mvar=1)

pc.ContinuousTapControl(net, tr1, vm_set_pu=1, tol=1e-6)

pp.runpp(net, run_control=True)

print('\nTensions i potències del nusos')
print(net.res_bus)
tap_pos = net.trafo.tap_pos[tr1]
print(f'\nPosició del "tap" = {tap_pos:.6f}')
rel = 1 + tap_pos/100
print(f'\nRelació de transformació = {rel:.6f}')
