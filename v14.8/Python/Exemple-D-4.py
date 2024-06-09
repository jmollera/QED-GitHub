import numpy as np
from scipy import optimize
import sympy

def f(x):
    return 35953.6865*np.sin(314.1593*x - 1.5136) + 35894.8169*np.exp(-18*x)

def df(x):
    return 11295184.9833*np.cos(314.1593*x - 1.5136) - 646106.7042*np.exp(-18*x)

print('Resolució numèrica')
print('-'*64)
rs = optimize.root_scalar(f, method='newton', x0=0.015, fprime=df)
print('Newton: x =', rs.root, ', f(x) =', f(rs.root))
rs = optimize.root_scalar(f, method='secant', x0=0.015, x1=0.016)
print('Secant: x =', rs.root, ', f(x) =', f(rs.root))

x = sympy.symbols('x')
f_symb = 35953.6865*sympy.sin(314.1593*x - 1.5136) + 35894.8169*sympy.exp(-18*x)
df_symb = sympy.diff(f_symb, x)
f_num = sympy.lambdify(x, f_symb, 'numpy')
df_num = sympy.lambdify(x, df_symb, 'numpy')

print('\nResolució simbòlica + numèrica')
print('-'*64)
rs = optimize.root_scalar(f_num, method='newton', x0=0.015, fprime=df_num)
print('Newton: x =', rs.root, ', f(x) =', f(rs.root))
rs = optimize.root_scalar(f_num, method='secant', x0=0.015, x1=0.016)
print('Secant: x =', rs.root, ', f(x) =', f(rs.root))
