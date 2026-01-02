import locale
import numpy as np
from scipy import integrate

locale.setlocale(category=locale.LC_NUMERIC, locale="")

# moment d'inèrcia lb·ft²
J = 507
# velocitat r/min
ω = np.array([0, 90, 180, 270, 360, 450, 540, 630, 720, 810, 900,
              990, 1080, 1170, 1260, 1350, 1440, 1530, 1620, 1656,
              1692, 1728, 1764, 1780], dtype=float)
# parell motor lb·ft
Tm = np.array([266, 248, 234, 225, 221, 225, 234, 246, 258, 273,
               289, 304, 322, 338, 356, 373, 388, 396, 395, 388,
               370, 336, 248, 177], dtype=float)
# parell càrrega lb·ft
Tc = np.array([0.0, 0.4, 1.7, 3.8, 6.8, 10.6, 15.2, 20.7, 27.0,
               34.2, 42.2, 51.1, 60.8, 71.4, 82.8, 95.0, 108.1,
               122.0, 136.8, 143.0, 149.2, 155.7, 162.2, 165.2], dtype=float)
try:
    with open("Cap-Motors-Induccio-MotTempsArr-Aux-1.txt", "w", encoding="utf8") as fh:
        for i in range(len(ω)):
            fh.write(locale.format_string(
                    "%4.0f  &  %3.0f  &  %5.1f  \\\\ \n",
                    (ω[i], Tm[i], Tc[i])))
except Exception as err:
    print(err)

# moment d'inèrcia kg·m²
J = J*4.214011e-2
print(locale.format_string("moment d'inèrcia = %.2f kg·m²", J))
# velocitat rad/s
ω = ω*np.pi/30
# parell accelerador N·m
Tm = Tm*1.355818
Tc = Tc*1.355818
Tacc = Tm-Tc

try:
    with open("Cap-Motors-Induccio-MotTempsArr-Aux-2.txt", "w", encoding="utf8") as fh:
        for i in range(len(ω)):
            fh.write(locale.format_string(
                    "%7.3f  &  %7.3f  &  %7.3f  &  %8.6f  \\\\ \n",
                    (ω[i], Tm[i], Tc[i], 1/Tacc[i])))
except OSError as err:
    print(err)

# tarr s - integració pel mètode dels trapezis
tarr = J*integrate.trapezoid(1/Tacc, ω)
print (locale.format_string("temps d'arranc = %.1f s", tarr))
