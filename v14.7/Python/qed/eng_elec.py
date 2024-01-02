__author__ = 'Josep Mollera Barriga'
__version__ = '14.7'
__name__ = 'qed.eng_elc'
__all__ = ['EE', 'z_cable', 'r_cable', 'voltage_drop', 'kcmil_to_mm2',
           'mm2_to_kcmil', 'awg_to_mm2', 'mm2_to_awg', 'z_series',
           'z_parallel', 'millman', 'rms', 'average', 'apparent_power',
           'D_to_Y', 'Y_to_D', 'triangle_to_phasors', 'LN_to_LL',
           'LL_to_LG', 'LL_to_LN', 'ABC_to_A012', 'A012_to_ABC',
           'AN12_to_AB12', 'AB12_to_AN12', 'ezs_u', 'icc_tr', 'CEI_51_curve',
           'IEEE_51_curve', 'EEException', 'EEInvalidArguments',
           'EENetNotSolved', 'EENetMissingBranch', 'EENetMissingNode',
           'EENetMissingMutualCoupling', 'EENetInvalidElement',
           'EENetNotSolvable', 'Impedance', 'Admittance',
           'VoltageSource', 'VoltageSourceIdeal', 'CurrentSource',
           'CurrentSourceIdeal', 'ShortCircuit', 'Branch', 'CoupledBranches',
           'Network', 'Motor3phRun', 'Motor3phStartUp', 'Motor3ph']

import array
import tabulate
import enum
from dataclasses import dataclass
from typing import Callable, Optional
import pickle
import numpy as np
from numpy.typing import ArrayLike
from scipy import linalg, integrate, optimize, interpolate
from qed.utils import Complex

tabulate.PRESERVE_WHITESPACE = True

# Two dictionaries of impedances for copper cables.
#
# Index: section in mm²
# Value: impedance in Ω/km at 90 °C

_Z_CABLES_A = {4: 5.99 + 0.098j, 6: 3.96 + 0.123j, 10: 2.34 + 0.117j,
               16: 1.47 + 0.112j, 25: 0.936 + 0.107j, 35: 0.661 + 0.123j,
               50: 0.5 + 0.11j, 70: 0.34 + 0.114j, 120: 0.192 + 0.108j,
               150: 0.156 + 0.104j, 185: 0.125 + 0.103j, 240: 0.095 + 0.102j,
               300: 0.078 + 0.094j}

_Z_CABLES_V = {4: 5.993 + 0.168j, 6: 3.94 + 0.159j, 10: 2.33 + 0.151j,
               16: 1.46 + 0.142j, 25: 0.9361 + 0.0791j, 35: 0.6749 + 0.0765j,
               50: 0.4935 + 0.135j, 70: 0.3417 + 0.128j, 95: 0.249 + 0.130j,
               120: 0.1980 + 0.0825j}

# Dictionary of resistivity and coefficients of resistivity
# variation with temperature, for several materials.
#
# Index: chemical symbol of the material
# Value: resistivity in Ω·mm²/m at 20 °C, and coefficient of
#        resistivity variation with temperature in 1/°C at 20 °C

_R_MATERIAL = {'Al': (0.02825, 0.00391), 'Cu': (0.01723, 0.00393),
               'Ag': (0.01645, 0.00380), 'Au': (0.02440, 0.00340)}


class EE(enum.Enum):
    """Constants to be used in functions z_cable, CEI_51_curve and IEEE_51_curve"""
    CABLES_A = enum.auto()
    CABLES_V = enum.auto()
    CEI_INVERSE = enum.auto()
    CEI_VERY_INVERSE = enum.auto()
    CEI_EXTREMELY_INVERSE = enum.auto()
    IEEE_MODERATELY_INVERSE = enum.auto()
    IEEE_VERY_INVERSE = enum.auto()
    IEEE_EXTREMELY_INVERSE = enum.auto()


def z_cable(section: int, length: float, cable: EE = EE.CABLES_A) -> complex:
    """
    Find the impedance of a cable.

    Parameters
    ----------
    section: Cable section in mm².
    length:  Cable length in m.
    cable:   EE.CABLES_A → use dictionary _Z_CABLES_A,
             EE.CABLES_V → use dictionary _Z_CABLES_V,
             The default is EE.CABLES_A.

    Returns
    -------
    The cable impedance in Ω at 90 °C, or 0 if 'cable'
    is not EE.CABLES_A or EE.CABLES_V.

    Examples
    --------
    >>> z_cable(25, 100)
    (0.0936+0.0107j)
    >>> z_cable(25, 100, EE.CABLES_A)
    (0.0936+0.0107j)
    >>> z_cable(25,100, EE.CABLES_V)
    (0.09361+0.00791j)
    >>> z_cable(25,100, 'X')
    0
    """
    try:
        match cable:
            case EE.CABLES_A:
                return _Z_CABLES_A[section] * length / 1000
            case EE.CABLES_V:
                return _Z_CABLES_V[section] * length / 1000
            case _:
                raise KeyError
    except KeyError:
        return 0


def r_cable(section: float, length: float,
            temp: float = 20.0, material: str = 'Cu') -> float:
    """
    Find the resistance of a cable.

    Parameters
    ----------
    section:  Cable section in mm².
    length:   Cable length in m.
    temp:     Temperature in °C. The default is 20 °C.
    material: The chemical symbol of the cable material.
              The default is 'Cu'.

    Returns
    -------
    The cable resistance in Ω at the requested temperature,
    or 0 if 'material' is not valid.

    Examples
    --------
    >>> print(f'R = {r_cable(4, 100):.4f} Ω')
    R = 0.4307 Ω
    >>> print(f'R = {r_cable(4, 100, temp=90):.4f} Ω')
    R = 0.5492 Ω
    >>> print(f'R = {r_cable(4, 100, material="Al"):.4f} Ω')
    R = 0.7063 Ω
    >>> print(f'R = {r_cable(4, 100, temp=90, material="Al"):.4f} Ω')
    R = 0.8996 Ω
    >>> r_cable(4, 100, material='Fe')
    0
    """
    try:
        rho_20, alpha_20 = _R_MATERIAL[material]
        return rho_20 * (1 + alpha_20 * (temp - 20)) * length / section
    except KeyError:
        return 0


def voltage_drop(voltage: float, current: float,
                 impedance: complex, cos_phi: float = 1.0) -> float:
    """
    Find the voltage drop through a cable.

    Parameters
    ----------
    voltage:   Voltage at cable origin in V,
               (line-neutral voltage for a three-phase system).
    current:   Current through the cable in A.
    impedance: Cable impedance in Ω.
    cos_phi:   Load power factor. The default is 1.

    Returns
    -------
    Voltage drop through the cable in V,
    (line-neutral voltage for a three-phase system).

    Examples
    --------
    Three-phase: U=380 V (LL). I=200 A. Cable: 240 mm², 400 m. cos 𝜑=0.87
    >>> dv = voltage_drop(380/np.sqrt(3), 200, z_cable(240,400), 0.87) * np.sqrt(3)
    >>> print(f'{dv:.4f} V (L-L)')
    18.4652 V (L-L)

    One-phase: U=220 V. I=55 A. Cable: 35 mm², 200 m. cos 𝜑=0.85
    >>> dv = voltage_drop(220, 55, z_cable(35,2*200), 0.85)
    >>> print(f'{dv:.4f} V')
    13.8515 V

    Direct Current: U=125 V. I=25 A. Cable: 16 mm², 50 m
    >>> dv = voltage_drop(125, 25, z_cable(16,2*50).real)
    >>> print(f'{dv:.4f} V')
    3.6750 V
    """
    R = np.real(impedance)
    X = np.imag(impedance)
    U = np.abs(voltage)
    I = np.abs(current)
    sin_phi = np.sqrt(1 - cos_phi**2)
    return (R * cos_phi + X * sin_phi) * I + U - np.sqrt(U ** 2 - ((X * cos_phi - R * sin_phi) * I) ** 2)


def kcmil_to_mm2(kcmil: float) -> float:
    """
    Convert a circular section from kcmil to mm².

    Parameters
    ----------
    kcmil: Circular section in kcmil.

    Returns
    -------
    Equivalent section in mm².

    Examples
    --------
    >>> print(f'{kcmil_to_mm2(250):.4f}')
    126.6769
    """
    return kcmil * 0.506707479098


def mm2_to_kcmil(mm2: float) -> float:
    """
    Convert a circular section from mm² to kcmil.

    Parameters
    ----------
    mm2: Circular section in mm².

    Returns
    -------
    Equivalent section in kcmil.

    Examples
    --------
    >>> print(f'{mm2_to_kcmil(300):.4f}')
    592.0576
    """
    return mm2 / 0.506707479098


def awg_to_mm2(awg: int) -> float:
    """
    Convert a circular section from AWG to mm².

    Parameters
    ----------
    awg: Circular section in AWG; 0, -1, -2, and -3 must be
         used for AWG 1/0, 2/0, 3/0, and 4/0.

    Returns
    -------
    Equivalent section in mm².

    Examples
    --------
    >>> print(f'{awg_to_mm2(12):.4f}')
    3.3088
    """
    return 53.4751207321 / (92**(awg / 19.5))


def mm2_to_awg(mm2: float) -> int:
    """
    Convert a circular section from mm² to AWG.

    Parameters
    ----------
    mm2: Circular section in mm².

    Returns
    -------
    Equivalent section rounded to the nearest AWG; results 0, -1, -2,
    and -3 are equivalent to AWG 1/0, 2/0, 3/0, and 4/0.

    Examples
    --------
    >>> mm2_to_awg(2)
    14
    """
    return int(np.rint(4.31245284200 * np.log(53.4751207321 / mm2)))


def z_series(z: ArrayLike) -> complex:
    """
    Compute the series impedance of a list of impedances.

    Parameters
    ----------
    z: Impedances connected in series.

    Returns
    -------
    The series impedance.

    Examples
    --------
    >>> z_series([2, 3, 4])
    9
    >>> z_series([2+2j, 3+3j, 4+4j])
    (9+9j)
    """
    return np.sum(z)


def z_parallel(z: ArrayLike) -> complex:
    """
    Compute the parallel impedance of a list of impedances.

    Parameters
    ----------
    z: Impedances connected in parallel.

    Returns
    -------
    The parallel impedance.

    Examples
    --------
    >>> z_parallel([3, 3, 3])
    1.0
    >>> z_parallel([3+3j, 3+3j, 3+3j])
    (1+1j)
    """
    p = np.prod(z)
    return p / np.sum(p / z)


def millman(u: ArrayLike, z: ArrayLike) -> complex:
    """
    Apply the Millman theorem.

    Parameters
    ----------
    u: Voltages with respect to a reference point of the end points of
       a group of star connected impedances.
    z: Impedances of the star connected impedances.

    Returns
    -------
    The voltage with respect to the same reference point, of the
    neutral point of the start connected impedances.

    Examples
    --------
    >>> u = millman([125.1, 124.8, 125.2], [0.034, 0.041, 0.029])
    >>> print(f'{u:.4f}')
    125.0562
    >>> u = millman([24+10j, 26+9j, 27+7j], [1+1j, 2+2j ,3+3j])
    >>> print(f'{u:.4f}')
    25.0909+9.1818j
    """
    return np.sum(np.array(u) / z) * z_parallel(z)


def rms(f: Callable[[float], float], a: float, b: float, **kwargs: dict) -> float:
    """
    Compute the rms value of a periodic function.

    Parameters
    ----------
    f:      A periodic function of one variable.
    a:      Starting abscissa of the period.
    b:      Ending abscissa of the period.
    kwargs: Optional arguments accepted by scipy.integrate.quad

    Returns
    -------
    The rms value of 'f' with period from 'a' to 'b'.
    """
    val, _ = integrate.quad(lambda x: f(x)**2, a, b, **kwargs)
    return np.sqrt(val / (b - a))


def average(f: Callable[[float], float], a: float, b: float, **kwargs: dict) -> float:
    """
    Compute the average value of a periodic function.

    Parameters
    ----------
    f:      A periodic function of one variable.
    a:      Starting abscissa of the period.
    b:      Ending abscissa of the period.
    kwargs: Optional arguments accepted by scipy.integrate.quad

    Returns
    -------
    The average value of 'f' with period from 'a' to 'b'.
    """
    val, _ = integrate.quad(lambda x: f(x), a, b, **kwargs)
    return val / (b - a)


def apparent_power(u: ArrayLike, i: ArrayLike) -> complex:
    """
    Compute the apparent power drawn by a load.

    Parameters
    ----------
    u: Voltages with respect to a reference point of a set of
       branches feeding a load.
    i: Currents through the branches towards the load.

    Returns
    -------
    Te apparent power drawn by the load.

    Examples
    --------
    >>> apparent_power([24+10j, 26+9j, 27+7j], [1+1j, 2+2j ,3+3j])
    (206-108j)
    >>> apparent_power([110], [3])
    330
    """
    return np.array(u) @ np.array(i).conj()


def D_to_Y(z_ab: float | complex, z_bc: float | complex,
           z_ca: float | complex) -> tuple[float | complex, float | complex, float | complex]:
    """
    Convert three delta impedances into three wye impedances.

    Parameters
    ----------
    z_ab: Impedance between phases A and B.
    z_bc: Impedance between phases B and C.
    z_ca: Impedance between phases C and A.

    Returns
    -------
    Equivalent wye impedances: z_an, z_bn, and z_cn.

    Examples
    --------
    >>> D_to_Y(3, 3, 3)
    (1.0, 1.0, 1.0)
    >>> D_to_Y(10, -10j, -10j)
    ((4-2j), (4-2j), (-2-4j))
    """
    z = (z_ab * z_bc * z_ca) / (z_ab + z_bc + z_ca)
    return z / z_bc, z / z_ca, z / z_ab


def Y_to_D(z_an: float | complex, z_bn: float | complex,
           z_cn: float | complex) -> tuple[float | complex, float | complex, float | complex]:
    """
    Convert three wye impedances into three delta impedances.

    Parameters
    ----------
    z_an: Impedance between phase A and the neutral point.
    z_bn: Impedance between phase B and the neutral point.
    z_cn: Impedance between phase C and the neutral point.

    Returns
    -------
    Equivalent delta impedances: z_ab, z_bc, and z_ca.

    Examples
    --------
    >>> Y_to_D(1, 1, 1)
    (3.0, 3.0, 3.0)
    >>> Y_to_D(4-2j, 4-2j, -2-4j)
    ((10-0j), -10j, -10j)
    """
    z = z_an * (z_bn + z_cn) + z_bn * z_cn
    return z / z_cn, z / z_an, z / z_bn


def triangle_to_phasors(u_ab: float, u_bc: float, u_ca: float) -> tuple[complex, complex, complex]:
    """
    Convert the three sides of a voltage triangle into three voltage phasors.

    Parameters
    ----------
    u_ab: Side AB of the voltage triangle.
    u_bc: Side BC of the voltage triangle.
    u_ca: Side CA of the voltage triangle.

    Returns
    -------
    Voltage phasors: u_ab, u_bc, and u_ca.
    Phasor u_ab is taken as the reference phasor (angle equal to zero).

    Examples
    --------
    >>> print('  '.join(format(u, '.4f') for u in triangle_to_phasors(400, 500, 300)))
    400.0000+0.0000j  -400.0000-300.0000j  0.0000+300.0000j
    >>> print('  '.join(format(u, '.4f') for u in triangle_to_phasors(2760, 1840, 2300)))
    2760.0000+0.0000j  -1035.0000-1521.3070j  -1725.0000+1521.3070j
    """
    ang_BAC = np.arccos((u_ab**2 + u_ca**2 - u_bc**2) / (2 * u_ab * u_ca))
    ang_ABC = np.arccos((u_bc**2 + u_ab**2 - u_ca**2) / (2 * u_bc * u_ab))
    return (complex(u_ab, 0),
            u_bc * np.exp(1j * (np.pi + ang_ABC)),
            u_ca * np.exp(1j * (np.pi - ang_BAC)))


def LN_to_LL(u_an: complex, u_bn: complex, u_cn: complex) -> tuple[complex, complex, complex]:
    """
    Find the three line-line voltages of three line-neutral voltages.

    Parameters
    ----------
    u_an: Voltage AN phasor.
    u_bn: Voltage BN phasor.
    u_cn: Voltage CN phasor.

    Returns
    -------
    Voltage phasors: u_ab, u_bc, and u_ca.

    Examples
    --------
    >>> LN_to_LL(120j, 120, -120-120j)
    ((-120+120j), (240+120j), (-120-240j))
    """
    return u_an - u_bn, u_bn - u_cn, u_cn - u_an


def LL_to_LG(u_ab: complex, u_bc: complex, u_ca: complex) -> tuple[complex, complex, complex]:
    """
    Find the three line-G voltages of three line-line voltages.

    G is the barycentre of the triangle formed by the three line-line voltages.

    Parameters
    ----------
    u_ab: Voltage AB phasor.
    u_bc: Voltage BC phasor.
    u_ca: Voltage CA phasor.

    Returns
    -------
    Voltage phasors: u_ag, u_bg, and u_cg.

    Examples
    --------
    >>> LL_to_LG(-120+120j, 240+120j, -120-240j)
    (120j, (120+0j), (-120-120j))
    """
    return (u_ab - u_ca) / 3, (u_bc - u_ab) / 3, (u_ca - u_bc) / 3


def LL_to_LN(u_ab: complex, u_bc: complex, u_ca: complex,
             z_an: complex, z_bn: complex, z_cn: complex) -> tuple[complex, complex, complex]:
    """
    Find the three line-neutral voltages of three line-line voltages,
    with three star connected impedances.

    Parameters
    ----------
    u_ab: Voltage AB phasor.
    u_bc: Voltage BC phasor.
    u_ca: Voltage CA phasor.
    z_an: Impedance between phase A and neutral point.
    z_bn: Impedance between phase B and neutral point.
    z_cn: Impedance between phase C and neutral point.

    Returns
    -------
    Voltage phasors: u_an, u_bn, and u_cn.

    Examples
    --------
    >>> LL_to_LN(-120+120j, 240+120j, -120-240j, 10, 10, 10)
    (120j, (120+0j), (-120-120j))
    """
    z_p = z_parallel([z_an, z_bn, z_cn])
    return ((u_ab / z_bn - u_ca / z_cn) * z_p,
            (u_bc / z_cn - u_ab / z_an) * z_p,
            (u_ca / z_an - u_bc / z_bn) * z_p)


def ABC_to_A012(x_a: complex, x_b: complex, x_c: complex) -> tuple[complex, complex, complex]:
    """
    Find the sequence phasors of a three-phase unbalanced system.

    Parameters
    ----------
    x_a: Phasor A in a system of three unbalanced phasors.
    x_b: Phasor B in a system of three unbalanced phasors.
    x_c: Phasor C in a system of three unbalanced phasors.

    Returns
    -------
    Sequence phasors: x_a0, x_a1, and x_a2.

    Examples
    --------
    >>> print('  '.join(format(u, '.4f') for u in ABC_to_A012(2760, -1035-1521j, -1740+1535j)))
    -5.0000+4.6667j  2264.6912+201.1826j  500.3088-205.8493j
    """
    a = complex(-1/2, np.sqrt(3)/2)
    a2 = a.conjugate()
    return ((x_a + x_b + x_c) / 3,
            (x_a + a * x_b + a2 * x_c) / 3,
            (x_a + a2 * x_b + a * x_c) / 3)


def A012_to_ABC(x_a0: complex, x_a1: complex, x_a2: complex) -> tuple[complex, complex, complex]:
    """
    Find a three-phase unbalanced system of a sequence phasors.

    Parameters
    ----------
    x_a0: Homopolar sequence phasor.
    x_a1: Direct sequence phasor.
    x_a2: Inverse sequence phasor.

    Returns
    -------
    System of three unbalanced phasors: x_a, x_b and x_c.

    Examples
    --------
    >>> print('  '.join(format(u, '.4f') for u in A012_to_ABC(-5+4j, 2264+201j, 500-205j)))
    2759.0000+0.0000j  -1035.3937-1521.6688j  -1738.6063+1533.6688j
    """
    a = complex(-1/2, np.sqrt(3)/2)
    a2 = a.conjugate()
    return (x_a0 + x_a1 + x_a2,
            x_a0 + a2 * x_a1 + a * x_a2,
            x_a0 + a * x_a1 + a2 * x_a2)


def AN12_to_AB12(x_an_1: complex, x_an_2: complex) -> tuple[complex, complex]:
    """
    Find the line-line direct and inverse sequence phasors corresponding
    to the line-neutral direct and inverse sequence phasors.

    Parameters
    ----------
    x_an_1: Line-neutral direct sequence phasor.
    x_an_2: Line-neutral inverse sequence phasor.

    Returns
    -------
    Line-line sequence phasors: x_ab_1, and x_ab_2.

    Examples
    --------
    >>> print('  '.join(format(u, '.4f') for u in AN12_to_AB12(1187-552j, 308+44j)))
    2258.5460+199.9722j  500.1051-200.7358j
    """
    return (x_an_1 * np.sqrt(3) * np.exp(1j * np.pi/6),
            x_an_2 * np.sqrt(3) * np.exp(-1j * np.pi/6))


def AB12_to_AN12(x_ab_1: complex, x_ab_2: complex) -> tuple[complex, complex]:
    """
    Find the line-neutral direct and inverse sequence phasors
    corresponding to the line-line direct and inverse sequence phasors

    Parameters
    ----------
    x_ab_1: Line-line direct sequence phasor.
    x_ab_2: Line-line inverse sequence phasor.

    Returns
    -------
    Line-neutral sequence phasors: x_an_1, and x_an_2.

    Examples
    --------
    >>> print('  '.join(format(u, '.4f') for u in AB12_to_AN12(2258+200j, 500-200j)))
    1186.7350-551.8285j  307.7350+44.3376j
    """
    return (x_ab_1 / np.sqrt(3) / np.exp(1j * np.pi/6),
            x_ab_2 / np.sqrt(3) / np.exp(-1j * np.pi/6))


def ezs_u(e: float | complex, z: float | complex, s: float | complex) -> float | complex:
    """
    Find the voltage at a load, given the source
    voltage, the impedance between source and load, and
    the power drawn by the load.

    Parameters
    ----------
    e: Source voltage.
    z: Impedance between source and load.
    s: Power drawn by the load.

    Returns
    -------
    Voltage at the load.

    Examples
    --------
    >>> print(f'{ezs_u(0.4+0.3j, 0.1j, 0.6+0.45j):.4f}')
    0.3165+0.0874j
    >>> print(f'{ezs_u(125, 0.1, 100):.4f}')
    124.9199
    """
    ZcjS = z.conjugate() * s
    ImUe = ZcjS.imag / abs(e)
    E_abs = abs(e)
    ReUe = np.polynomial.Polynomial([ZcjS.real + ImUe ** 2, -E_abs, 1]).roots()

    if isinstance(ReUe[0], complex):
        print('No solution found')
        return 0

    U = complex(max(ReUe), ImUe) * e / E_abs
    is_complex = any(isinstance(var, complex) for var in [e, z, s])
    return U if is_complex else U.real


def icc_tr(up_tr: float, us_tr: float, s_tr: float, xcc_tr: float,
           u_ext: float, scc_ext: float = float('inf')) -> float:
    """
    Find the shortcircuit current at a transformer's secondary side.

    Parameters
    ----------
    up_tr:   Transformer's primary rated voltage, in V.
    us_tr:   Transformer's secondary rated voltage, in V.
    s_tr:    Transformer's rated power, in VA.
    xcc_tr:  Transformer's shortcircuit impedance, in %.
    u_ext:   Voltage of the external source connected to the
             transformer's primary side, in V.
    scc_ext: Shortcircuit power of the external source connected
             to the transformer's primary side, in VA.
             The default is 'inf' (infinite power source)

    Returns
    -------
    Shortcircuit current at the transformer's secondary side, in A

    Examples
    --------
    >>> print(f'Icc = {icc_tr(6900, 400, 850e3, 5, 6900, 200e6)/1000:.2f} kA')
    Icc = 22.62 kA
    >>> print(f'Icc = {icc_tr(6900, 400, 850e3, 5, 6900)/1000:.2f} kA')
    Icc = 24.54 kA
    """
    icc_tr_pu = u_ext / up_tr / (s_tr / scc_ext * (u_ext / up_tr) ** 2 + xcc_tr / 100)
    return icc_tr_pu * s_tr/(np.sqrt(3) * us_tr)


def CEI_51_curve(I: ArrayLike, curve: EE, I_p: float, T: float = 1.0,
                 C: float = 0.0, B: float = 0.0) -> ArrayLike:
    """
    Compute the actuation time 't' of a CEI 51 curve,
    for a given current 'I', according to:
        t = T*(k/((I/I_p)**a - 1) + C) + B,
    where 'k' and 'a' are set by the 'curve' parameter.

    Parameters
    ----------
    I: Current, in A.
    curve: The CEI 51 curve to use. Allowed values:
              EE.CEI_INVERSE
              EE.CEI_VERY_INVERSE
              EE.CEI_EXTREMELY_INVERSE
    I_p:   Pick up current of the CEI 51 curve, in A.
    T:     Parameter of the CEI 51 curve. The default is 1.
    C:     Parameter of the CEI 51 curve. The default is 0.
    B:     Parameter of the CEI 51 curve. The default is 0.

    Returns
    -------
    The actuation time of the CEI 51 curve, in s.

    Example
    -------
    >>> print(f'{CEI_51_curve(200, EE.CEI_VERY_INVERSE, 100):.4f}')
    13.5000
    """
    try:
        match curve:
            case EE.CEI_INVERSE:
                k = 0.14
                a = 0.02
            case EE.CEI_VERY_INVERSE:
                k = 13.5
                a = 1.0
            case EE.CEI_EXTREMELY_INVERSE:
                k = 80.0
                a = 2.0
            case _:
                raise AttributeError

        return T * (k / ((I / I_p) ** a - 1) + C) + B

    except AttributeError:
        print(f'Unknown curve "{curve}"')
        return 0


def IEEE_51_curve(I: ArrayLike, curve: EE, I_p: float,
                  T: float = 1.0, B: float = 0.0) -> ArrayLike:
    """
    Compute the actuation time 't' 'of an IEEE 51 curve,
    for a given current 'I', according to:
        t = T*(k/((I/I_p)**a - 1) + C) + B,
    where 'k', 'a' and 'C' are set by the 'curve' parameter.

    Parameters
    ----------
    I:     Current, in A.
    curve: The IEEE 51 curve to use. Allowed values:
              EE.IEEE_MODERATELY_INVERSE
              EE.IEEE_VERY_INVERSE
              EE.IEEE_EXTREMELY_INVERSE
    I_p:   Pick up current of the IEEE 51 curve, in A.
    T:     Parameter of the IEEE 51 curve. The default is 1.
    B:     Parameter of the IEEE 51 curve. The default is 0.

    Returns
    -------
    The actuation time of the IEEE 51 curve, in s.

    Example
    -------
    >>> print(f'{IEEE_51_curve(200, EE.IEEE_VERY_INVERSE, 100):.4f}')
    7.0277
    """
    try:
        match curve:
            case EE.IEEE_MODERATELY_INVERSE:
                k = 0.0515
                a = 0.02
                C = 0.114
            case EE.IEEE_VERY_INVERSE:
                k = 19.61
                a = 2.0
                C = 0.491
            case EE.IEEE_EXTREMELY_INVERSE:
                k = 28.2
                a = 2.0
                C = 0.1217
            case _:
                raise AttributeError

        return T * (k / ((I / I_p) ** a - 1) + C) + B

    except AttributeError:
        print(f'Unknown curve "{curve}"')
        return 0


class EEException(Exception):
    """Base exception class."""


class EEInvalidArguments(EEException):
    """
    Exception raised when a function is called with invalid arguments.
    """

    def __init__(self, invalid_args):
        if len(invalid_args) == 1:
            message = f'Invalid argument: {invalid_args}'
        else:
            message = f'Invalid arguments: {invalid_args}'
        super().__init__(message)
        self.invalid_args = invalid_args


class EENetNotSolved(EEException):
    """
    Exception raised in class Network when results want to be
    obtained before the network is solved.
    """

    def __init__(self):
        message = 'Network not solved yet!'
        super().__init__(message)


class EENetMissingBranch(EEException):
    """
    Exception raised in class Network when trying to access missing branches.
    """

    def __init__(self, missing_branches):
        if len(missing_branches) == 1:
            message = f'Branch {missing_branches} is invalid or has no information!'
        else:
            message = f'Branches {missing_branches} are invalid or have no information!'
        super().__init__(message)
        self.missing_branches = missing_branches


class EENetMissingNode(EEException):
    """
    Exception raised in class Network when trying to access missing nodes.
    """

    def __init__(self, missing_nodes):
        if len(missing_nodes) == 1:
            message = f'Node {missing_nodes} is invalid or has no connections!'
        else:
            message = f'Nodes {missing_nodes} are invalid or have no connections!'
        super().__init__(message)
        self.missing_nodes = missing_nodes


class EENetMissingMutualCoupling(EEException):
    """
    Exception raised in class Network when trying to access missing mutual couplings.
    """

    def __init__(self, missing_mutual_coupling):
        message = f'Mutual coupling between branches {missing_mutual_coupling} does not exist!'
        super().__init__(message)
        self.missing_mutual_coupling = missing_mutual_coupling


class EENetInvalidElement(EEException):
    """
    Exception raised in class Network when trying to access an invalid element
    """

    def __init__(self, invalid_element):
        message = f'Element "{invalid_element}" is not valid'
        super().__init__(message)
        self.invalid_element = invalid_element


class EENetNotSolvable(EEException):
    """
    Exception raised in class Network when trying to solve an invalid or empty network.
    Examples: Two ideal voltage sources in parallel, two ideal current sources
    in series, a short-circuited ideal voltage source, or an open current source.
    """

    def __init__(self, not_solvable):
        message = f'Non solvable network ({not_solvable})'
        super().__init__(message)
        self.not_solvable = not_solvable


@dataclass(frozen=True)
class Element:
    """Base class for the elements to be used with class Branch/CoupledBranches."""

    @property
    def is_complex(self) -> bool:
        return any(isinstance(value, complex) for value in self.__dict__.values())


@dataclass(frozen=True)
class Impedance(Element):
    """Create an impedance to be used with class Branch/CoupledBranches."""

    Z: complex  # Impedance


@dataclass(frozen=True)
class Admittance(Element):
    """Create an admittance to be used with class Branch."""

    Y: complex  # Admittance


@dataclass(frozen=True)
class VoltageSource(Element):
    """Create a voltage source to be used with class Branch."""

    E: complex  # Source voltage in an oriented graph
    Z: complex  # Source series impedance


@dataclass(frozen=True)
class VoltageSourceIdeal(Element):
    """Create an ideal voltage source to be used with class Branch."""

    E: complex  # Source voltage in an oriented graph


@dataclass(frozen=True)
class CurrentSource(Element):
    """Create a current source to be used with class Branch."""

    J: complex  # Source current in an oriented graph
    Y: complex  # Source parallel admittance


@dataclass(frozen=True)
class CurrentSourceIdeal(Element):
    """Create an ideal current source to be used with class Branch."""

    J: complex  # Source current in an oriented graph


@dataclass(frozen=True)
class ShortCircuit(Element):
    """Create a short circuit to be used with class Branch."""


@dataclass(frozen=True)
class Branch:
    """Create a branch to be used with class Network."""

    elem: Element  # Element of the branch
    from_to: tuple[int, int]  # From and to nodes of the branch in an oriented graph
    branch: int  # The branch number


@dataclass(frozen=True)
class CoupledBranches:
    """Create a coupling between two branches to be used with class Network."""

    elem: Impedance  # Impedance of the coupling (must be a pure imaginary number)
    coupled_branches: tuple[int, int]  # The numbers of the two coupled branches


class Network:
    """
    Solve an electrical networks using the nodes method.

    >>> net = Network('Test network')
    >>> net.add(Branch(VoltageSource(E=200, Z=10), from_to=(0, 1), branch=1))
    >>> net.add(Branch(VoltageSource(E=-50, Z=20j), from_to=(1, 2), branch=2))
    >>> net.add(Branch(Impedance(Z=5j), from_to=(2, 0), branch=3))
    >>> net.add(Branch(Impedance(Z=20), from_to=(1, 2), branch=4))
    >>> net.add(Branch(CurrentSource(J=4, Y=0.1), from_to=(0, 2), branch=5))
    >>> net.add(CoupledBranches(Impedance(Z=5j), coupled_branches=(2, 3)))
    >>> net.remove(branch=2)
    >>> net.remove(coupled_branches=(3, 2))
    >>> net.add(Branch(VoltageSource(E=-50, Z=20j), from_to=(1, 2), branch=2))
    >>> net.add(CoupledBranches(Impedance(Z=5j), coupled_branches=(2, 3)))
    >>> net.solve()
    >>> net.num_branches
    5
    >>> net.num_nodes
    3
    >>> net.is_solved
    True
    >>> print(net)
    Test network
    ------------------------------------------------------------------
    Branch(elem=VoltageSource(E=200, Z=10), from_to=(0, 1), branch=1)
    Branch(elem=VoltageSource(E=-50, Z=20j), from_to=(1, 2), branch=2)
    Branch(elem=Impedance(Z=5j), from_to=(2, 0), branch=3)
    Branch(elem=Impedance(Z=20), from_to=(1, 2), branch=4)
    Branch(elem=CurrentSource(J=4, Y=0.1), from_to=(0, 2), branch=5)
    CoupledBranches(elem=Impedance(Z=5j), coupled_branches=(2, 3))
    >>> print(f'{net.voltage(node=2):.4f}')
    33.5644+20.6436j
    >>> print(f'{net.voltage(branch=4):.4f}')
    119.2079+2.0792j
    >>> print(f'{net.current(branch=4) :.4f}')
    5.9604+0.1040j
    >>> print('  '.join(format(x, '.4f') for x in net.thevenin(node_1=1, node_2=2)))
    119.2079+2.0792j  4.2574+2.5743j
    >>> print('  '.join(format(x, '.4f') for x in net.norton(node_1=1, node_2=2)))
    20.7200-12.0400j  0.1720-0.1040j
    >>> print(net.results(polar=True, u_fmt='.2f/.1f', u_unit='V', i_fmt='.2f/.1f', i_unit='A'))
     Branch   Voltage / V     Current / A
    --------  --------------  -------------
       1      154.45∠-171.5°  5.24∠-25.7°
       2      119.23∠1.0°     2.68∠-117.5°
       3      39.40∠31.6°     6.90∠-38.9°
       4      119.23∠1.0°     5.96∠1.0°
       5      39.40∠-148.4°   2.16∠-72.7°
    """

    def __init__(self, name: str = '') -> None:
        """Create an AC or DC Network"""
        self._name = name
        self._net = {}
        self._is_solved = False
        self._is_complex = 0
        self._extra_nodes = 0
        self._extra_branches = 0
        self._list_of_nodes = array.array('Q')
        self._list_of_branches = array.array('Q')
        self._list_of_coupled_branches = array.array('Q')
        self._ZN = None
        self._VN = None
        self._IB = None

    def __str__(self) -> str:
        """Return a multiline string representing the class"""
        s = [[self._name]] if self._name else [['Network']]
        for idx in sorted(self._net, key=lambda x: str(len(x))+str(min(x))):
            s.append([self._net[idx]])
        return tabulate.tabulate(s, headers='firstrow')

    @property
    def name(self) -> str:
        """Return the name of the network"""
        return self._name

    @property
    def num_branches(self) -> int:
        """Return the umber of branches in the network."""
        return len(set(self._list_of_branches))

    @property
    def num_nodes(self) -> int:
        """Return the number of nodes in the network."""
        return len(set(self._list_of_nodes))

    @property
    def is_solved(self) -> bool:
        """Return True if the network has been solved, and False otherwise."""
        return self._is_solved

    def add(self, item_cls: Branch | CoupledBranches) -> None:
        """
        Add an element to the network (an instance of classes Branch or CoupledBranches).
        Exceptions EENetInvalidElement and EEInvalidArguments may be raised.
        """
        match (item_cls, item_cls.elem):
            case (Branch(), Element()):
                if (m := min(item_cls.from_to)) < 0 or m == max(item_cls.from_to):
                    raise EEInvalidArguments({'from_to': item_cls.from_to})
                if item_cls.branch < 1:
                    raise EEInvalidArguments({'branch': item_cls.branch})
                idx = frozenset([item_cls.branch])
                if idx in self._net:
                    self.remove(branch=item_cls.branch)
                self._net[idx] = item_cls
                self._is_solved = False
                self._is_complex += item_cls.elem.is_complex
                self._list_of_branches.append(item_cls.branch)
                self._list_of_nodes.extend(item_cls.from_to)
                match item_cls.elem:
                    case VoltageSourceIdeal() | ShortCircuit():
                        self._extra_branches += 1
                        self._extra_nodes += 1
                    case CurrentSourceIdeal():
                        self._extra_branches += 1
            case (CoupledBranches(), Impedance()):
                if item_cls.elem.Z.real != 0:
                    raise EEInvalidArguments({'Z': item_cls.elem.Z})
                if (m := min(item_cls.coupled_branches)) < 1 or m == max(item_cls.coupled_branches):
                    raise EEInvalidArguments({'coupled_branches': item_cls.coupled_branches})
                idx = frozenset(item_cls.coupled_branches)
                if idx in self._net:
                    self.remove(coupled_branches=item_cls.coupled_branches)
                self._net[idx] = item_cls
                self._is_solved = False
                self._is_complex += item_cls.elem.is_complex
                self._list_of_coupled_branches.extend(item_cls.coupled_branches)
            case _:
                raise EENetInvalidElement(item_cls)

    def remove(self, *, branch: Optional[int] = None,
               coupled_branches: Optional[tuple[int, int]] = None) -> None:
        """
        Remove a branch (int) or a coupling between two branches (tuple of int) from the network.

        Exceptions EENetMissingBranch, EENetMissingMutualCoupling and
        EEInvalidArguments may be raised.
        """
        match (branch, coupled_branches):
            case (None, None):
                raise EEInvalidArguments({'branch': None, 'coupled_branches': None})
            case (branch, None):
                try:
                    idx = frozenset([branch])
                    item_cls = self._net.pop(idx)
                    self._is_solved = False
                    self._is_complex -= item_cls.elem.is_complex
                    for n in item_cls.from_to:
                        self._list_of_nodes.remove(n)
                    self._list_of_branches.remove(item_cls.branch)
                    match item_cls.elem:
                        case VoltageSourceIdeal() | ShortCircuit():
                            self._extra_branches -= 1
                            self._extra_nodes -= 1
                        case CurrentSourceIdeal():
                            self._extra_branches -= 1
                except KeyError as err:
                    raise EENetMissingBranch({branch}) from err
            case (None, coupled_branches):
                try:
                    idx = frozenset(coupled_branches)
                    item_cls = self._net.pop(idx)
                    self._is_solved = False
                    self._is_complex -= item_cls.elem.is_complex
                    for branch in item_cls.coupled_branches:
                        self._list_of_coupled_branches.remove(branch)
                except KeyError as err:
                    raise EENetMissingMutualCoupling({coupled_branches}) from err
            case _:
                raise EEInvalidArguments({'branch': branch, 'coupled_branches': coupled_branches})

    def solve(self) -> None:
        """
        Solve the network.

        The voltage across every branch, and the current through every branch
        is calculated.
        Exception EENetMissingBranch is raised if not all branches
        1, 2, 3, ... max(self._branches) are present.
        Exception EENetMissingNode is raised if not all nodes
        0, 1, 2, ... max(self._nodes) are present.
        Exception EENetNotSolvable is raised if linalg.inv raises the
        linalg.LinAlgError exception.
        """

        # Check for missing nodes and branches, set the number of nodes and
        # branches, and check the data type (complex or real).
        if not self._net:
            raise EENetNotSolvable('empty network')

        actual_branches = set(self._list_of_branches)
        num_branches = max(actual_branches)
        expected_branches = set(b for b in range(1, num_branches + 1))
        if diff := actual_branches ^ expected_branches:
            raise EENetMissingBranch(diff)
        if diff := set(self._list_of_coupled_branches) - expected_branches:
            raise EENetMissingBranch(diff)

        actual_nodes = set(self._list_of_nodes)
        num_nodes = max(actual_nodes)
        expected_nodes = set(n for n in range(0, num_nodes + 1))
        if diff := actual_nodes ^ expected_nodes:
            raise EENetMissingNode(diff)

        # Set up the matrices needed to solve the network
        num_nodes_plus_extra = num_nodes + self._extra_nodes
        num_branches_plus_extra = num_branches + self._extra_branches
        A = np.zeros((num_nodes_plus_extra, num_branches_plus_extra), dtype=np.byte)
        data_type = np.cdouble if self._is_complex else np.double
        ZB = np.zeros((num_branches_plus_extra, num_branches_plus_extra), dtype=data_type)
        EB1 = np.zeros(num_branches_plus_extra, dtype=data_type)
        JB1 = np.zeros(num_branches_plus_extra, dtype=data_type)

        def fill_A(branch: int, from_to: tuple[int, int]) -> None:
            n_from, n_to = from_to
            if n_from > 0:
                A[n_from - 1, branch - 1] = 1
            if n_to > 0:
                A[n_to - 1, branch - 1] = -1

        for item_cls in self._net.values():
            match (item_cls, item_cls.elem):
                case (Branch(), Impedance()):
                    fill_A(item_cls.branch, item_cls.from_to)
                    ZB[item_cls.branch - 1, item_cls.branch - 1] = item_cls.elem.Z
                case (Branch(), Admittance()):
                    fill_A(item_cls.branch, item_cls.from_to)
                    ZB[item_cls.branch - 1, item_cls.branch - 1] = 1 / item_cls.elem.Y
                case (Branch(), VoltageSource()):
                    fill_A(item_cls.branch, item_cls.from_to)
                    EB1[item_cls.branch - 1] = item_cls.elem.E
                    ZB[item_cls.branch - 1, item_cls.branch - 1] = item_cls.elem.Z
                case (Branch(), VoltageSourceIdeal()):
                    fill_A(item_cls.branch, (item_cls.from_to[0], num_nodes_plus_extra))
                    fill_A(num_branches_plus_extra, (num_nodes_plus_extra, item_cls.from_to[1]))
                    EB1[item_cls.branch - 1] = item_cls.elem.E
                    ZB[item_cls.branch - 1, item_cls.branch - 1] = 1
                    ZB[num_branches_plus_extra - 1, num_branches_plus_extra - 1] = -1
                    num_branches_plus_extra -= 1
                    num_nodes_plus_extra -= 1
                case (Branch(), CurrentSource()):
                    fill_A(item_cls.branch, item_cls.from_to)
                    JB1[item_cls.branch - 1] = item_cls.elem.J
                    ZB[item_cls.branch - 1, item_cls.branch - 1] = 1 / item_cls.elem.Y
                case (Branch(), CurrentSourceIdeal()):
                    fill_A(item_cls.branch, item_cls.from_to)
                    fill_A(num_branches_plus_extra, item_cls.from_to)
                    JB1[item_cls.branch - 1] = item_cls.elem.J
                    ZB[item_cls.branch - 1, item_cls.branch - 1] = 1
                    ZB[num_branches_plus_extra - 1, num_branches_plus_extra - 1] = -1
                    num_branches_plus_extra -= 1
                case (Branch(), ShortCircuit()):
                    fill_A(item_cls.branch, (item_cls.from_to[0], num_nodes_plus_extra))
                    fill_A(num_branches_plus_extra, (num_nodes_plus_extra, item_cls.from_to[1]))
                    ZB[item_cls.branch - 1, item_cls.branch - 1] = 1
                    ZB[num_branches_plus_extra - 1, num_branches_plus_extra - 1] = -1
                    num_branches_plus_extra -= 1
                    num_nodes_plus_extra -= 1
                case (CoupledBranches(), Impedance()):
                    b1, b2 = item_cls.coupled_branches
                    ZB[b1 - 1, b2 - 1] = item_cls.elem.Z
                    ZB[b2 - 1, b1 - 1] = item_cls.elem.Z

        # Solve the network
        try:
            YB = linalg.inv(ZB, overwrite_a=True)
            JB = YB @ EB1 + JB1
            ZN = linalg.inv(A @ YB @ A.T, overwrite_a=True)
            VN = ZN @ (-A @ JB)
            IB = YB @ (A.T @ VN) + JB
            self._is_solved = True
            # Adjust result for branches with ideal current sources
            for item_cls in self._net.values():
                if isinstance(item_cls.elem, CurrentSourceIdeal):
                    IB[item_cls.branch - 1] = item_cls.elem.J
            # Remove extra nodes and branches
            self._ZN = np.delete(np.delete(ZN, np.s_[num_nodes:], axis=0), np.s_[num_nodes:], axis=1)
            self._VN = np.delete(VN, np.s_[num_nodes:])
            self._IB = np.delete(IB, np.s_[num_branches:])
        except linalg.LinAlgError as err:
            raise EENetNotSolvable(err) from err

    def thevenin(self, node_1: int, node_2: int) -> tuple[complex, complex]:
        """
        Return the Thevenin equivalent circuit between two nodes.

        Exceptions EENetNotSolved and EENetMissingNode may be raised.
        """
        self._check_is_solved()
        self._check_nodes(node_1, node_2)

        match (node_1, node_2):
            case (node_1, node_2) if node_1 == node_2:
                eth = zth = 0
            case (0, node_2):
                eth = -self._VN[node_2 - 1]
                zth = self._ZN[node_2 - 1, node_2 - 1]
            case (node_1, 0):
                eth = self._VN[node_1 - 1]
                zth = self._ZN[node_1 - 1, node_1 - 1]
            case (node_1, node_2):
                eth = self._VN[node_1 - 1] - self._VN[node_2 - 1]
                zth = (self._ZN[node_1 - 1, node_1 - 1] +
                       self._ZN[node_2 - 1, node_2 - 1] -
                       self._ZN[node_1 - 1, node_2 - 1] -
                       self._ZN[node_2 - 1, node_1 - 1])
        return eth, zth

    def norton(self, node_1: int, node_2: int) -> tuple[complex, complex]:
        """
        Return the Norton equivalent circuit between two nodes.

        May raise the same exceptions as the function 'thevenin'.
        """
        eth, zth = self.thevenin(node_1, node_2)
        match (eth, zth):
            case (0, 0):
                jno = 0
                yno = float('inf')
            case (eth, 0):
                jno = yno = float('inf')
            case _:
                jno = eth / zth
                yno = 1 / zth
        return jno, yno

    def voltage(self, *, branch: Optional[int] = None, node: Optional[int] = None) -> complex:
        """
        Return the voltage across a branch, or the voltage of a node (with respect to node 0).

        Exception EENetNotSolved, EENetMissingBranch, EENetMissingNode and
        EEInvalidArguments may be raised.
        """
        self._check_is_solved()

        match (branch, node):
            case(None, None):
                raise EEInvalidArguments({'branch': None, 'node': None})
            case (branch, None):
                self._check_branches(branch)
                idx = frozenset([branch])
                return self.thevenin(*self._net[idx].from_to)[0]
            case (None, node):
                self._check_nodes(node)
                return self.thevenin(node, 0)[0]
            case _:
                raise EEInvalidArguments({'branch': branch, 'node': node})

    def current(self, branch: int) -> complex:
        """
        Return the current through a branch.

        Exceptions EENetNotSolved and EENetMissingBranch may be raised.
        """
        self._check_is_solved()
        self._check_branches(branch)

        return self._IB[branch - 1]

    def results(self, polar: bool = False, u_fmt: str = '', u_unit: str = '', u_scale: float = 1,
                i_fmt: str = '', i_unit: str = '', i_scale: float = 1) -> str:
        """
        Return a multi-line string formatted as a table, with all the network's voltages and currents.

        A multi-line formatted string is returned with the voltages and
        currents of all the network branches.
        A different scale factor can be applied to voltages 'u_scale' and currents 'i_scale'.
        Polar display can be chosen 'polar', and different formats for voltages 'u_fmt'
        and currents 'i_fmt' are possible. Also, different units can be added
        to voltages 'u_unit' and currents 'u_unit'.
        Exception EENetNotSolved may be raised.
        """
        self._check_is_solved()

        if u_unit:
            u_unit = ' / ' + u_unit
        if i_unit:
            i_unit = ' / ' + i_unit
        s = [['Branch', 'Voltage' + u_unit, 'Current' + i_unit]]
        for branch in range(1, self.num_branches + 1):
            u = self.voltage(branch=branch) * u_scale
            i = self.current(branch) * i_scale
            if polar:
                u = Complex(u)
                i = Complex(i)
            s.append([branch, format(u, u_fmt), format(i, i_fmt)])
        return tabulate.tabulate(s, numalign='center', headers='firstrow')

    def save(self, filename: str) -> bool:
        """Save the contents of the network to a binary file."""
        try:
            with open(filename, 'wb') as fh:
                pickle.dump(self._name, fh, pickle.HIGHEST_PROTOCOL)
                pickle.dump(self._net, fh, pickle.HIGHEST_PROTOCOL)
                pickle.dump(self._is_solved, fh, pickle.HIGHEST_PROTOCOL)
                pickle.dump(self._is_complex, fh, pickle.HIGHEST_PROTOCOL)
                pickle.dump(self._extra_nodes, fh, pickle.HIGHEST_PROTOCOL)
                pickle.dump(self._extra_branches, fh, pickle.HIGHEST_PROTOCOL)
                pickle.dump(self._list_of_nodes, fh, pickle.HIGHEST_PROTOCOL)
                pickle.dump(self._list_of_branches, fh, pickle.HIGHEST_PROTOCOL)
                pickle.dump(self._list_of_coupled_branches, fh, pickle.HIGHEST_PROTOCOL)
                pickle.dump(self._ZN, fh, pickle.HIGHEST_PROTOCOL)
                pickle.dump(self._VN, fh, pickle.HIGHEST_PROTOCOL)
                pickle.dump(self._IB, fh, pickle.HIGHEST_PROTOCOL)
                return True
        except (OSError, pickle.PicklingError) as err:
            print(f'{filename}: save error: {err}')
            return False

    def load(self, filename: str) -> bool:
        """Load the network with the contents of a binary file."""
        try:
            with open(filename, 'rb') as fh:
                self._name = pickle.load(fh)
                self._net = pickle.load(fh)
                self._is_solved = pickle.load(fh)
                self._is_complex = pickle.load(fh)
                self._extra_nodes = pickle.load(fh)
                self._extra_branches = pickle.load(fh)
                self._list_of_nodes = pickle.load(fh)
                self._list_of_branches = pickle.load(fh)
                self._list_of_coupled_branches = pickle.load(fh)
                self._ZN = pickle.load(fh)
                self._VN = pickle.load(fh)
                self._IB = pickle.load(fh)
                return True
        except (OSError, pickle.UnpicklingError) as err:
            print(f'{filename}: load error: {err}')
            return False

    def _check_is_solved(self) -> None:
        """Checks if the network is solved."""
        if not self._is_solved:
            raise EENetNotSolved

    def _check_nodes(self, *nodes: int) -> None:
        """Checks if a list of nodes are valid."""
        for node in nodes:
            if node not in self._list_of_nodes:
                raise EENetMissingNode({node})

    def _check_branches(self, *branches: int) -> None:
        """Checks if a list of branches are valid."""
        for branch in branches:
            if branch not in self._list_of_branches:
                raise EENetMissingBranch({branch})


@dataclass(frozen=True)
class Motor3phRun:
    """Class used to store the result of method 'run' from class Motor3ph."""

    Z: np.ndarray    # Total motor impedance (Line-Neutral), Ω
    PF: np.ndarray   # Motor Power factor
    I: np.ndarray    # Current drawn by the motor, A
    U: np.ndarray    # Voltage at the motor terminals (Line-Neutral), V
    S: np.ndarray    # Apparent power drawn by the motor, VA
    T_m: np.ndarray  # Mechanical torque delivered by the motor, N·m
    P_m: np.ndarray  # Mechanical power delivered by the motor, W
    Eff: np.ndarray  # Motor Efficiency (between 0 and 1)


@dataclass(frozen=True)
class Motor3phStartUp:
    """Class used to store the result of method 'start_up' from class Motor3ph."""

    t_start: float  # The starting time of the motor, s
    s_start_up: interpolate.interp1d  # An interpolating function that returns the slip
                                      # at a given time in seconds
                                      # For t<0, the slip at 0 is returned
                                      # for t>t_start, the slip at t_start is returned


class Motor3ph:
    """
    Compute electrical and mechanical quantities of a 3-phase asynchronous motor.

    >>> R_1 = 0.095  # Ω/phase
    >>> X_1 = 0.680  # Ω/phase
    >>> R_2 = 0.3  # Ω/phase
    >>> X_2 = 0.672  # Ω/phase
    >>> R_Fe = 620  # Ω/phase
    >>> X_m = 18.7  # Ω/phase
    >>> p = 4  # number of poles
    >>> f = 50  # Hz
    >>> s_n = 0.04 # rated slip
    >>> U = 230  # V, phase-neutral
    >>> motor = Motor3ph(p, f, R_1, X_1, R_2, X_2, R_Fe, X_m)
    >>> mot_n = motor.run(s_n, U)
    >>> print(f'n_sync = {motor.n_m_sync:.0f} r/min')
    n_sync = 1500 r/min
    >>> print(f'I = {abs(mot_n.I):.1f} A')
    I = 32.3 A
    >>> print(f'T_m = {mot_n.T_m:.1f} N·m')
    T_m = 118.8 N·m
    """

    def __init__(self, p: int, f: float, r_1: float, x_1: float,
                 r_2: float, x_2: float, r_fe: float, x_m: float) -> None:
        """
        Create a motor.

        Parameters
        ----------
        p:    Number of poles
        f:    Frequency, Hz
        r_1:  Stator resistance, Ω/phase
        x_1   Stator reactance, Ω/phase
        r_2:  Rotor resistance, Ω/phase
        x_2:  Rotor reactance, Ω/phase
        r_fe: Core-loss resistance, Ω/phase
        x_m:  Magnetizing reactance, Ω/phase
        """
        self._p = p
        self._f = f
        self._z_1 = complex(r_1, x_1)
        self._r_2 = r_2
        self._jx_2 = 1j * x_2
        self._z_0 = z_parallel([r_fe, 1j * x_m])
        self._ω_m_sync = 4 * np.pi * f / p
        self._n_m_sync = 120 * f / p

    @property
    def ω_m_sync(self) -> float:
        """Return the motor synchronous speed in rad/s."""
        return self._ω_m_sync

    @property
    def n_m_sync(self) -> float:
        """Return the motor synchronous speed in r/min."""
        return self._n_m_sync

    def to_s(self, *, ω_m: Optional[float] = None,
             n_m: Optional[float] = None) -> float:
        """
        Return the slip corresponding to a mechanical speed
        given in rad/s (ω_m) or in r/min (n_m).
        """
        match (ω_m, n_m):
            case (None, None):
                raise EEInvalidArguments({'ω_m': None, 'n_m': None})
            case (ω_m, None):
                return 1 - ω_m / self._ω_m_sync
            case (None, n_m):
                return 1 - n_m / self._n_m_sync
            case _:
                raise EEInvalidArguments({'ω_m': ω_m, 'n_m': n_m})

    def to_ω_m(self, *, s: Optional[float] = None,
               n_m: Optional[float] = None) -> float:
        """
        Return the mechanical speed in rad/s corresponding to a
        given slip (s) or a mechanical speed in r/min (n_m).
        """
        match (s, n_m):
            case (None, None):
                raise EEInvalidArguments({'s': None, 'n_m': None})
            case (s, None):
                return (1 - s) * self._ω_m_sync
            case (None, n_m):
                return n_m * np.pi / 30
            case _:
                raise EEInvalidArguments({'s': s, 'n_m': n_m})

    def to_n_m(self, *, s: Optional[float] = None,
               ω_m: Optional[float] = None) -> float:
        """
        Return the mechanical speed in r/min corresponding to a
        given slip (s) or a mechanical speed in rad/s (ω_m).
        """
        match (s, ω_m):
            case (None, None):
                raise EEInvalidArguments({'s': None, 'ω_m': None})
            case (s, None):
                return (1 - s) * self._n_m_sync
            case (None, ω_m):
                return ω_m * 30 / np.pi
            case _:
                raise EEInvalidArguments({'s': s, 'ω_m': ω_m})

    def s_T_m_max(self, z_ext: complex = 0) -> float:
        """
        Return the slip at the maximum motor mechanical torque.

        Parameters
        ----------
        z_ext: Impedance of the external power system feeding the motor, Ω/phase.
               The default is 0.
        """
        z_th = z_parallel([self._z_0, self._z_1 + z_ext])
        return self._r_2 / abs(z_th + self._jx_2)

    def s_match(self, T_load: Callable[[float], float], e_ext: complex,
                z_ext: complex = 0, s_guess: float = 0.04) -> float:
        """
        Return the slip at which the motor torque matches the load torque.

        Parameters
        ----------
        T_load:  A function that returns the load torque for a given slip, N·m
        e_ext:   Voltage (phase-neutral) of the external power system feeding the motor, V
        z_ext:   Impedance of the external power system feeding the motor, Ω/phase.
                 The default is 0.
        s_guess: Matching slip initial guess. The default is 0.04.
        """
        sol = optimize.root(lambda s: self.run(s, e_ext, z_ext).T_m - T_load(s), x0=s_guess)
        return sol.x[0]

    def run(self, s: float, e_ext: complex, z_ext: complex = 0) -> Motor3phRun:
        """
        Return several motor quantities when running at a given slip, in a Motor3phRun class.

        Parameters
        ----------
        s:     Slip at which the motor quantities are to be computed
        e_ext: Voltage (line-neutral) of the external power system feeding the motor, V
        z_ext: Impedance of the external power system feeding the motor, Ω/phase.
               The default is 0.

        Returns
        -------
        A Motor3phRun class.
        """
        r_2_s = self._r_2 / np.where(np.isclose(s, 0.0, atol=1e-12), 1e-12, s)
        z_2_s = r_2_s + self._jx_2
        e_th = e_ext * self._z_0 / (self._z_0 + self._z_1 + z_ext)
        z_th = z_parallel([self._z_0, self._z_1 + z_ext])
        Z = self._z_1 + self._z_0 * z_2_s / (self._z_0 + z_2_s)
        PF = np.cos(np.angle(Z))
        I = e_ext / (z_ext + Z)
        U = np.where(z_ext == 0, e_ext, Z * I)
        S = 3 * U * np.conj(I)
        T_m = 3 * r_2_s * np.abs(e_th / (z_th + z_2_s)) ** 2 / self._ω_m_sync
        P_m = T_m * self.to_ω_m(s=s)
        Eff = P_m / np.real(S)
        return Motor3phRun(Z, PF, I, U, S, T_m, P_m, Eff)

    def start_up(self, T_load: Callable[[float], float], J: float,
                 e_ext: complex, z_ext: complex = 0) -> Motor3phStartUp:
        """
        Return the starting time of the motor, and a function of the time evolution of
        the slip during the starting up process, in a Motor3phStartUp class.

        Parameters
        ----------
        T_load: A function that returns the load torque for a given slip, N·m
        J:      Total moment of inertia of the motor plus the load, kg·m²
        e_ext:  Voltage (line-neutral) of the external power system feeding the motor, V
        z_ext:  Impedance of the external power system feeding the motor, Ω/phase.
                The default is 0.

        Returns
        -------
        A Motor3phStartUp class.
        """
        n_p1, n_p2 = 200, 800
        s_lin = 0.3  # The function to integrate is very smooth up to approximately this value
        s_match = self.s_match(T_load, e_ext, z_ext)
        s = np.concatenate([np.linspace(1, s_lin, n_p1, endpoint=False),
                            np.linspace(s_lin, s_match, n_p2, endpoint=False)])
        T_accel = self.run(s, e_ext, z_ext).T_m - T_load(s)
        t = -J * self._ω_m_sync * integrate.cumulative_trapezoid(1 / T_accel, s, initial=0)
        t_start = t[-1]
        s_start_up = interpolate.interp1d(t, s, kind='cubic', copy=False,
                                          assume_sorted=True, bounds_error=False,
                                          fill_value=(s[0], s[-1]))
        return Motor3phStartUp(t_start, s_start_up)


if __name__ == '__main__':
    import doctest

    doctest.testmod()
