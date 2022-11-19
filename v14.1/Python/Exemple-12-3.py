from scipy import optimize
import sympy

m, δ2 = sympy.symbols('m, δ2')

f_symb  = [-2 -(-50/m*sympy.sin(-δ2)),
           -1 +(50/m*sympy.cos(-δ2) -50)]

Jf_symb = sympy.Matrix(f_symb).jacobian(sympy.Matrix([δ2, m]))

f_num = sympy.lambdify([(δ2, m)], f_symb, 'numpy')
Jf_num = sympy.lambdify([(δ2, m)], Jf_symb, 'numpy')
sol = optimize.root(f_num, x0=[0, 1], jac=Jf_num)
δ_2, m_ = sol.x
print(f'δ_2 = {δ_2:.6f}')
print(f'm = {m_:.6f}')