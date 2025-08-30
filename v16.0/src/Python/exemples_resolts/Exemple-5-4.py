import numpy as np
import sympy
import matplotlib.pyplot as plt

sympy.init_printing(use_unicode=True)

s = sympy.symbols('s')
t = sympy.symbols('t', positive=True)
a = sympy.sqrt(2)*110_000_000
b = 100*sympy.pi

uCs = a*s/((s**2 + 1000*s + 500_000) * (s**2 + b**2))

uCt = sympy.inverse_laplace_transform(uCs, s, t, simplify=True)
print(f'uC(t) = {uCt}')

uCt_num = sympy.lambdify(t, uCt, 'numpy')
t = np.linspace(0, 0.05, 500)

fig, ax = plt.subplots(figsize=(6, 4))
ax.set_xlim(0, 50)
ax.set_xticks(np.arange(0, 50.1, 5))
ax.set_ylim(-400, 400)
yticks = np.arange(-400, 401, 100)
ax.set_yticks(yticks, [f'{tick: }' for tick in yticks])
ax.grid()
ax.set_ylabel(r'$u_{\mathrm{C}}(t)$ / V')
ax.set_xlabel(r'$t$ / ms')
ax.plot(t*1000, uCt_num(t), 'red')
plt.show()
