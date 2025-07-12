import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

temp = [60, 70, 80, 90]  # Temperatura en °F
conc = [40, 50, 60, 70]  # Concentració en %
visc = [[3.38, 4.55, 6.33, 8.97],
        [2.87, 3.81, 5.17, 7.22],
        [2.46, 3.23, 4.28, 5.88],
        [2.13, 2.76, 3.58, 4.58]]  # Viscositat dinàmica en cP = f(temp, conc)

punt = (76, 56.3)  # Punt d'interpolació (°F, %)

f_lin = interpolate.RegularGridInterpolator((temp, conc), visc, method='linear')
print(f'Interpolació lineal: {f_lin(punt):.2f} cP\n')

f_cub = interpolate.RegularGridInterpolator((temp, conc), visc, method='cubic')
print(f'Interpolació cúbica: {f_cub(punt):.2f} cP')

temp_ls, conc_ls = np.meshgrid(np.linspace(60, 90, 200), np.linspace(40, 70, 200))
visc_ls_cub = f_cub((temp_ls, conc_ls))
temp_scatter, conc_scatter = np.meshgrid(temp, conc)
visc_scatter = np.array(visc).T
fig = plt.figure(figsize=(7, 7))
ax = plt.axes(projection="3d")
ax.set_xlim(60, 90)
ax.set_xticks(np.arange(60, 91, 5))
ax.set_ylim(40, 70)
ax.set_yticks(np.arange(40, 71, 5))
ax.set_zlim(0, 10)
ax.set_xlabel('Temperatura / °F')
ax.set_ylabel('Concentració / %')
ax.set_zlabel('Viscositat / cP')
ax.plot_surface(temp_ls, conc_ls, visc_ls_cub, alpha=0.8)
ax.scatter(temp_scatter, conc_scatter, visc_scatter, color='red')
plt.show()
