import numpy as np
import sympy
import matplotlib.pyplot as plt

sympy.init_printing(use_unicode=True)

U_bat = 120  # V
R_bat = 5  # Ω
R_1 = 100  # Ω
R_2 = 75  # Ω
R_3 = 25  # Ω
C = sympy.Rational(10, 1000_000)  # F
L = sympy.Rational(250, 1000)  # H

iL_0 = sympy.Rational(U_bat, R_bat + R_3)
uC_0 = iL_0*R_3
print(f'iL_0 = {iL_0} A')
print(f'uC_0 = {uC_0} V')

s, iLs, uCs, ICs = sympy.symbols('s iLs uCs ICs')
t = sympy.symbols ('t', positive=True)

eq1 = sympy.Eq(R_1*iLs + s*L*iLs - L*iL_0 + R_2*(iLs + ICs), 0)
eq2 = sympy.Eq(uCs + R_3*ICs + R_2*(iLs + ICs), 0)
eq3 = sympy.Eq(uCs, ICs/(s*C) + uC_0/s)
solICs = sympy.solve(eq3, ICs)[0]
eq1 = eq1.subs(ICs, solICs)
eq2 = eq2.subs(ICs, solICs)
soluCs = sympy.solve(eq1, uCs)[0]
eq2 = eq2.subs(uCs, soluCs)
soliLs = sympy.simplify(sympy.solve(eq2, iLs)[0])
print(f'iL(s) = {soliLs}')

iLt = sympy.simplify(sympy.inverse_laplace_transform(soliLs, s, t, noconds=True))
print(f'iL(t) = {iLt}')
iL_0 = iLt.subs(t, 0)
print(f'iL(0) = {iL_0} A')

iLt_num  = sympy.lambdify(t, iLt, 'numpy')
t = np.linspace(0, 0.01, 500)

fig, ax = plt.subplots(figsize=(6, 4))
ax.set_xlim(0, 10)
ax.set_xticks(np.arange(0, 10.1, 1))
ax.set_ylim(-1, 4)
ax.set_yticks(np.arange(-1, 4.1, 1))
ax.grid()
ax.set_ylabel(r'$i_{\mathrm{L}}(t)$ / A')
ax.set_xlabel(r'$t$ / ms')
ax.plot(t*1000, iLt_num(t), 'red')
plt.show()