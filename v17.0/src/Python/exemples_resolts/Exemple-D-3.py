from scipy import integrate
import sympy

x_p = [1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
y_p = [1.0, 0.8333, 0.7143, 0.6250, 0.5556, 0.5]

int_trapz = integrate.trapezoid(y=y_p, x=x_p)
print(f'Integral (trapezoid) = {int_trapz:.4f}')

int_simps = integrate.simpson(y=y_p, x=x_p)
print(f'Integral (simpson) = {int_simps:.4f}')

int_quad, abserr = integrate.quad(lambda x: 1/x, a=1, b=2)
print(f'Integral (quad) = {int_quad}, Error = {abserr:.2e}')

z = sympy.symbols('z')
int_simb = sympy.integrate(1/z, (z, 1, 2))
int_num = sympy.N(int_simb, n=30)
print(f'Integral (simbòlica) = {int_simb} = {int_num}')
