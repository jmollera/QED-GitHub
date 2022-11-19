import numpy as np
from scipy import interpolate

x = [1.2, 1.4, 1.6, 1.8]
y = [0.9329, 0.9854, 0.9996 ,0.9738]

π = np.pi

f_lin = interpolate.interp1d(x, y, kind='linear')
y_lin = f_lin(π/2)
print(f'Interpolació lineal: {y_lin:.4f}')

f_cub = interpolate.interp1d(x, y, kind='cubic')
y_cub = f_cub(π/2)
print(f'Interpolació cúbica: {y_cub:.4f}')