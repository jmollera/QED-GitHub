from scipy import interpolate

x = [40, 50, 60, 70]  # Concentració en %
y = [60, 70, 80, 90]  # Temperatura en °F
z = [[3.38, 4.55, 6.33, 8.97],
     [2.87, 3.81, 5.17, 7.22],
     [2.46, 3.23, 4.28, 5.88],
     [2.13, 2.76, 3.58, 4.58]]  # Viscositat dinàmica en cP

p = (56.3, 76)  # Punt d'interpolació

f_lin = interpolate.RegularGridInterpolator((x, y), z, method='linear')
y_lin = f_lin(p)
print(f'Interpolació lineal: {y_lin:.2f} cP')

f_cub = interpolate.RegularGridInterpolator((x, y), z, method='cubic')
y_cub = f_cub(p)
print(f'Interpolació cúbica: {y_cub:.2f} cP')
