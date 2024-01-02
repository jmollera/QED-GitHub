import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

x = [1.2, 1.4, 1.6, 1.8]
y = [0.9329, 0.9854, 0.9996, 0.9738]

f_lin = interpolate.interp1d(x, y, kind='linear')
print(f'Interpolació lineal: {f_lin(np.pi/2):.4f}')

f_cub = interpolate.interp1d(x, y, kind='cubic')
print(f'Interpolació cúbica: {f_cub(np.pi/2):.4f}')

X = np.linspace(1.2, 1.8, 100)
fig, ax = plt.subplots(figsize=(6, 4))
ax.set_xlim(1.1, 1.9)
ax.set_xticks(np.arange(1.1, 1.91, 0.1))
ax.set_ylim(0.93, 1.01)
ax.set_yticks(np.arange(0.93, 1.011, 0.01))
ax.grid()
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')
ax.plot(X, np.sin(X), label=r'$y = \sin(x)$')
ax.plot(x, y,  'o', label='Punts d\'interpolació')
ax.plot(X, f_lin(X), label='Interpolació lineal')
ax.plot(X, f_cub(X), linestyle='--', label='Interpolació cúbica')
ax.legend(loc='best')
plt.show()
