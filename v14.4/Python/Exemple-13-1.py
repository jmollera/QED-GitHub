import numpy as np
from scipy import optimize, integrate
import matplotlib.pyplot as plt
from matplotlib import ticker

def Ia(t):
    """t en segons. Ia en ampere."""
    Id = 1.98*np.exp(-t/0.023) + 3.89*np.exp(-t/0.475) + 1.46 - 1.459*(1 - np.exp(-t/0.475))
    Iq = 3.159*np.exp(-t/0.023) + 0.781*np.exp(-t/0.106)
    return 4368*np.hypot(Id, Iq)

def t51(I):
    """I en ampere. t51 en segons. Cal definir el paràmetre T."""
    return T*(0.014/((I/6000)**0.022 - 1) + 0.0226)

def integral(tau):
    """ tau en segons. El valor de la integral és adimensional."""
    val, _ = integrate.quad(lambda t: 1/t51(Ia(t)), 0, tau)
    return val

T = 0.6

def log_format(x, pos):
    decimal_places = int(np.maximum(-np.log10(x), 0))
    return f'{x:.{decimal_places:1d}f}'

t = np.linspace(0, 3, 500)  # segons
I = np.linspace(6100, 100_000, 500)  # ampere
fig, ax = plt.subplots(figsize=(6, 4))
ax.set(xscale='log')
ax.set(yscale='log')
ax.xaxis.set_major_formatter(ticker.FuncFormatter(log_format))
ax.yaxis.set_major_formatter(ticker.FuncFormatter(log_format))
ax.set_xlim(100, 100000)
ax.set_ylim(0.01, 10)
ax.grid(which='minor', alpha=0.4)
ax.grid(which='major', alpha=1)
ax.set_xlabel(r'$I$ / A')
ax.set_ylabel(r'$t$ / s')
ax.plot(Ia(t), t, 'red', label=r'$I_{\mathrm{a}}(t)$')
ax.plot(I, t51(I), 'green', label=r'$t_{51}(I)$')
ax.legend(loc='upper right')
plt.show()

sol = optimize.root(lambda x: integral(x) - 1, x0=0.25)
print(f'T = {T}')
print('-'*35)
print(f"Temps calculat d'actuació: {sol.x[0]:.4f} s")
print(f'Valor de la integral: {integral(tau = sol.x[0]):.4f}\n')
print(sol)

T = 0.8

sol = optimize.root(lambda x: integral(x) - 1, x0=0.25)
print(f'\nT = {T}')
print('-'*35)
print(f"Temps calculat d'actuació: {sol.x[0]:.4f} s")
print(f'Valor de la integral: {integral(tau = sol.x[0]):.4f}\n')
print(sol)
