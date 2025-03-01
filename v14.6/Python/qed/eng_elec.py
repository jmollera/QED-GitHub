__author__ = 'Josep Mollera Barriga'
__version__ = '14.6'
__name__ = 'qed.eng_elc'
__all__ = ['z_cable', 'r_cable', 'voltage_drop', 'kcmil_to_mm2',
           'mm2_to_kcmil', 'awg_to_mm2', 'mm2_to_awg', 'z_series',
           'z_parallel', 'millman', 'rms', 'average', 'apparent_power',
           'D_to_Y', 'Y_to_D', 'triangle_to_phasors', 'LN_to_LL', 'LL_to_LG',
           'LL_to_LN', 'ABC_to_A012', 'A012_to_ABC', 'AN12_to_AB12',
           'AB12_to_AN12', 'ezs_u', 'CEI_51_curve', 'IEEE_51_curve',
           'Impedance', 'Admittance', 'VoltageSource', 'VoltageSourceIdeal',
           'CurrentSource', 'CurrentSourceIdeal', 'ShortCircuit',
           'MutualCoupling', 'NetworkException', 'NetworkNotSolved',
           'NetworkMissingBranch', 'NetworkMissingNode',
           'NetworkMissingMutualCoupling', 'NetworkInvalidElement',
           'NetworkNotSolvable', 'Network', 'Motor3phRun', 'Motor3phStartUp',
           'Motor3ph']

import array
import tabulate
from dataclasses import dataclass
from abc import abstractmethod
import pickle
import numpy as np
from scipy import linalg, integrate, optimize, interpolate
from qed.utils import ComplexD

tabulate.PRESERVE_WHITESPACE = True

# Two dictionaries of impedances for copper cables.
#
# Index: section in mm²
# Value: impedance in Ω/km at 90 °C

Z_CABLE_A = {4: 5.99 + 0.098j, 6: 3.96 + 0.123j, 10: 2.34 + 0.117j,
             16: 1.47 + 0.112j, 25: 0.936 + 0.107j, 35: 0.661 + 0.123j,
             50: 0.5 + 0.11j, 70: 0.34 + 0.114j, 120: 0.192 + 0.108j,
             150: 0.156 + 0.104j, 185: 0.125 + 0.103j, 240: 0.095 + 0.102j,
             300: 0.078 + 0.094j}

Z_CABLE_V = {4: 5.993 + 0.168j, 6: 3.94 + 0.159j, 10: 2.33 + 0.151j,
             16: 1.46 + 0.142j, 25: 0.9361 + 0.0791j, 35: 0.6749 + 0.0765j,
             50: 0.4935 + 0.135j, 70: 0.3417 + 0.128j, 95: 0.249 + 0.130j,
             120: 0.1980 + 0.0825j}

# Dictionary of resistivity and coefficients of resistivity
# variation with temperature, for several materials.
#
# Index: chemical symbol of the material
# Value: resistivity in Ω·mm²/m at 20 °C, and coefficient of
#        resistivity variation with temperature in 1/°C at 20 °C

R_CABLE = {'Al': (0.02825, 0.00391), 'Cu': (0.01723, 0.00393),
           'Ag': (0.01645, 0.00380), 'Au': (0.02440, 0.00340)}


def z_cable(section, length, kind='A'):
    """
    Find the impedance of a cable.

    Parameters
    ----------
    section: int
        Cable section in mm².
    length: float
        Cable length in m.
    kind: str (optional)
        'A' or 'a' → use Z_CABLE_A,
        'V' or 'v' → use Z_CABLE_V,
        The default is 'A'.

    Returns
    -------
    float
        The cable impedance in Ω at 90 °C, or 0 if kind is not valid
        or section is no in Z_CABLE_A or Z_CABLE_V.

    Examples
    --------
    >>> z_cable(25, 100)
    (0.0936+0.0107j)
    >>> z_cable(25, 100, 'A')
    (0.0936+0.0107j)
    >>> z_cable(25, 100, 'V')
    (0.09361+0.00791j)
    >>> z_cable(25, 100, 'X')
    0
    """
    try:
        match kind.upper():
            case 'A':
                return Z_CABLE_A[section] * length / 1000
            case 'V':
                return Z_CABLE_V[section] * length / 1000
            case _:
                raise KeyError
    except KeyError:
        return 0


def r_cable(section, length, temp=20.0, material='Cu'):
    """
    Find the resistance of a cable.

    Parameters
    ----------
    section: float
        Cable section in mm².
    length: float
        Cable length in m.
    temp: float (optional)
        Temperature in °C.
        The default is 20 °C.
    material: str (optional)
        The chemical symbol of the cable material.
        The default is 'Cu'.

    Returns
    -------
    float
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
    >>> r_cable(4,100,material='Fe')
    0
    """
    try:
        ρ_20, α_20 = R_CABLE[material]
        return ρ_20 * (1 + α_20 * (temp - 20)) * length / section
    except KeyError:
        return 0


def voltage_drop(voltage, current, impedance, cos_phi=1.0):
    """
    Find the voltage drop through a cable.

    Parameters
    ----------
    voltage: float
        Voltage at cable origin (line-neutral voltage
        for a three-phase system) in V.
    current: float
        Current through the cable in A.
    impedance: float
        Cable impedance in Ω.
    cos_phi: float (optional)
        Load power factor. The default is 1.

    Returns
    -------
    float
        Voltage drop through the cable (line-neutral voltage
        for a three-phase system) in V.

    Examples
    --------
    Three-phase: U=380 V (LL). I=200 A. Cable: 240 mm², 400 m. cos 𝜑=0.87
    >>> dv = voltage_drop(380/np.sqrt(3),200,z_cable(240, 400),0.87) * np.sqrt(3)
    >>> print(f'{dv:.4f} V (L-L)')
    18.4652 V (L-L)

    One-phase: U=220 V. I=55 A. Cable: 35 mm², 200 m. cos 𝜑=0.85
    >>> dv = voltage_drop(220,55,z_cable(35, 2*200),0.85)
    >>> print(f'{dv:.4f} V')
    13.8515 V

    Direct Current: U=125 V. I=25 A. Cable: 16 mm², 50 m
    >>> dv = voltage_drop(125,25,z_cable(16, 2*50).real)
    >>> print(f'{dv:.4f} V')
    3.6750 V
    """
    R = np.real(impedance)
    X = np.imag(impedance)
    U = np.abs(voltage)
    I = np.abs(current)
    sin_phi = np.sqrt(1 - cos_phi ** 2)
    return (R * cos_phi + X * sin_phi) * I + U - np.sqrt(U ** 2 - ((X * cos_phi - R * sin_phi) * I) ** 2)


def kcmil_to_mm2(kcmil):
    """
    Convert a circular section from kcmil to mm².

    Parameters
    ----------
    kcmil: float
        Circular section in kcmil.

    Returns
    -------
    float
        Equivalent section in mm².

    Examples
    --------
    >>> print(f'{kcmil_to_mm2(250):.4f}')
    126.6769
    """
    return kcmil * 0.506707479098


def mm2_to_kcmil(mm2):
    """
    Convert a circular section from mm² to kcmil.

    Parameters
    ----------
    mm2: float
        Circular section in mm².

    Returns
    -------
    float
        Equivalent section in kcmil.

    Examples
    --------
    >>> print(f'{mm2_to_kcmil(300):.4f}')
    592.0576
    """
    return mm2 / 0.506707479098


def awg_to_mm2(awg):
    """
    Convert a circular section from AWG to mm².

    Parameters
    ----------
    awg: int
        Circular section in AWG; 0, -1, -2, and -3 must be
        used for AWG 1/0, 2/0, 3/0, and 4/0.

    Returns
    -------
    float
        Equivalent section in mm².

    Examples
    --------
    >>> print(f'{awg_to_mm2(12):.4f}')
    3.3088
    """
    return 53.4751207321 / (92 ** (awg / 19.5))


def mm2_to_awg(mm2):
    """
    Convert a circular section from mm² to AWG.

    Parameters
    ----------
    mm2: float
        Circular section in mm².

    Returns
    -------
    int
        Equivalent section rounded to the nearest AWG;
        results 0, -1, -2, and -3 are equivalent to
        AWG 1/0, 2/0, 3/0, and 4/0.

    Examples
    --------
    >>> mm2_to_awg(2)
    14
    """
    return int(np.rint(4.31245284200 * np.log(53.4751207321 / mm2)))


def z_series(z):
    """
    Compute the series impedance of a list of impedances.

    Parameters
    ----------
    z: float|complex array-like
        Impedances connected in series.

    Returns
    -------
    float|complex
        The series impedance.

    Examples
    --------
    >>> z_series([2, 3, 4])
    9
    >>> z_series([2+2j, 3+3j, 4+4j])
    (9+9j)
    """
    return np.sum(z)


def z_parallel(z):
    """
    Compute the parallel impedance of a list of impedances.

    Parameters
    ----------
    z: float|complex array-like
        Impedances connected in parallel.

    Returns
    -------
    float|complex
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


def millman(u, z):
    """
    Apply the Millman theorem.

    Parameters
    ----------
    u: float|complex array-like
        Voltages with respect to a reference point
        of a group of star connected impedances.
    z: float|complex array-like
        Impedances of the star connected impedances.

    Returns
    -------
    float|complex
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


def rms(f, a, b, **kwargs):
    """
    Compute the rms value of a periodic function.

    Parameters
    ----------
    f: function
        A periodic function of time.
    a: float
        Starting time of the period.
    b: float
        Ending time of the period.
    kwargs: keyword arguments
        Optional arguments accepted by scipy.integrate.quad

    Returns
    -------
    float
        The rms value of f with period b - a.
    """
    val, _ = integrate.quad(lambda x: f(x) ** 2, a, b, **kwargs)
    return np.sqrt(val / (b - a))


def average(f, a, b, **kwargs):
    """
    Compute the average value of a periodic function.

    Parameters
    ----------
    f: function
        A periodic function of time.
    a: float
        Starting time of the period.
    b: float
        Ending time of the period.
    kwargs: keyword arguments
        Optional arguments accepted by scipy.integrate.quad

    Returns
    -------
    float
        The average value of f with period b - a.
    """
    val, _ = integrate.quad(f, a, b, **kwargs)
    return val / (b - a)


def apparent_power(u, i):
    """
    Compute the apparent power drawn by a load.

    Parameters
    ----------
    u: float|complex array-like
        Voltages with respect to a reference point
        of a set of branches feeding a load.
    i: float|complex array-like
        Currents through the branches towards the load.

    Returns
    -------
    float|complex
        Te apparent power drawn by the load.

    Examples
    --------
    >>> apparent_power([24+10j, 26+9j, 27+7j], [1+1j, 2+2j ,3+3j])
    (206-108j)
    >>> apparent_power([110], [3])
    330
    """
    return np.array(u) @ np.array(i).conj()


def D_to_Y(z_ab, z_bc, z_ca):
    """
    Convert three delta impedances into three wye impedances.

    Parameters
    ----------
    z_ab: float|complex
        Impedance between phases A and B.
    z_bc: float|complex
        Impedance between phases B and C.
    z_ca: float|complex
        Impedance between phases C and A.

    Returns
    -------
    (float|complex, float|complex, float|complex)
        Equivalent wye impedances z_an, z_bn, and z_cn.

    Examples
    --------
    >>> D_to_Y(3, 3, 3)
    (1.0, 1.0, 1.0)
    >>> D_to_Y(10, -10j, -10j)
    ((4-2j), (4-2j), (-2-4j))
    """
    z = (z_ab * z_bc * z_ca) / (z_ab + z_bc + z_ca)
    return z / z_bc, z / z_ca, z / z_ab


def Y_to_D(z_an, z_bn, z_cn):
    """
    Convert three wye impedances into three delta impedances.

    Parameters
    ----------
    z_an: float|complex
        Impedance between phase A and neutral point.
    z_bn: float|complex
        Impedance between phase B and neutral point.
    z_cn: float|complex
        Impedance between phase C and neutral point.

    Returns
    -------
    (float|complex, float|complex, float|complex)
        Equivalent delta impedances z_ab, z_bc, and z_ca.

    Examples
    --------
    >>> Y_to_D(1, 1, 1)
    (3.0, 3.0, 3.0)
    >>> Y_to_D(4-2j, 4-2j, -2-4j)
    ((10-0j), -10j, -10j)
    """
    z = z_an * (z_bn + z_cn) + z_bn * z_cn
    return z / z_cn, z / z_an, z / z_bn


def triangle_to_phasors(u_ab, u_bc, u_ca):
    """
    Convert the three sides of a voltage triangle into three voltage phasors.

    Parameters
    ----------
    u_ab: float
        Side AB of the voltage triangle.
    u_bc: float
        Side BC of the voltage triangle.
    u_ca: float
        Side CA of the voltage triangle.

    Returns
    -------
    (complex, complex, complex)
        Voltage phasors u_ab, u_bc, and u_ca.
        Phasor u_ab is taken as the reference phasor (angle equal to zero).

    Examples
    --------
    >>> print('  '.join(format(u, '.4f') for u in triangle_to_phasors(400, 500, 300)))
    400.0000+0.0000j  -400.0000-300.0000j  0.0000+300.0000j
    >>> print('  '.join(format(u, '.4f') for u in triangle_to_phasors(2760, 1840, 2300)))
    2760.0000+0.0000j  -1035.0000-1521.3070j  -1725.0000+1521.3070j
    """
    φ_a = np.arccos((u_ab ** 2 + u_ca ** 2 - u_bc ** 2) / (2 * u_ab * u_ca))
    φ_b = np.arccos((u_bc ** 2 + u_ab ** 2 - u_ca ** 2) / (2 * u_bc * u_ab))
    return (complex(u_ab, 0),
            u_bc * np.exp(1j * (np.pi + φ_b)),
            u_ca * np.exp(1j * (np.pi - φ_a)))


def LN_to_LL(u_an, u_bn, u_cn):
    """
    Find the three line-line voltages of three line-neutral voltages.

    Parameters
    ----------
    u_an: complex
        Voltage AN phasor.
    u_bn: complex
        Voltage BN phasor.
    u_cn: complex
        Voltage CN phasor.

    Returns
    -------
    (complex, complex, complex)
        Voltage phasors u_ab, u_bc, and u_ca.

    Examples
    --------
    >>> LN_to_LL(120j, 120, -120-120j)
    ((-120+120j), (240+120j), (-120-240j))
    """
    return u_an - u_bn, u_bn - u_cn, u_cn - u_an


def LL_to_LG(u_ab, u_bc, u_ca):
    """
    Find the three line-G voltages of three line-line voltages.

    G is the barycentre of the line-line voltages.

    Parameters
    ----------
    u_ab: complex
        Voltage AB phasor.
    u_bc: complex
        Voltage BC phasor.
    u_ca: complex
        Voltage CA phasor.

    Returns
    -------
    (complex, complex, complex)
        Voltage phasors u_ag, u_bg, and u_cg.

    Examples
    --------
    >>> LL_to_LG(-120+120j, 240+120j, -120-240j)
    (120j, (120+0j), (-120-120j))
    """
    return (u_ab - u_ca) / 3, (u_bc - u_ab) / 3, (u_ca - u_bc) / 3


def LL_to_LN(u_ab, u_bc, u_ca, z_an, z_bn, z_cn):
    """
    Find the three line-neutral voltages of three line-line voltages.

    Parameters
    ----------
    u_ab: complex
        Voltage AB phasor.
    u_bc: complex
        Voltage BC phasor.
    u_ca: complex
        Voltage CA phasor.
    z_an: complex
        Impedance between phase A and neutral point.
    z_bn: complex
        Impedance between phase B and neutral point.
    z_cn: complex
        Impedance between phase C and neutral point.

    Returns
    -------
    (complex, complex, complex)
        Voltage phasors u_an, u_bn, and u_cn.

    Examples
    --------
    >>> LL_to_LN(-120+120j, 240+120j, -120-240j, 10, 10, 10)
    (120j, (120+0j), (-120-120j))
    """
    z_p = z_parallel([z_an, z_bn, z_cn])
    return ((u_ab / z_bn - u_ca / z_cn) * z_p,
            (u_bc / z_cn - u_ab / z_an) * z_p,
            (u_ca / z_an - u_bc / z_bn) * z_p)


def ABC_to_A012(x_a, x_b, x_c):
    """
    Find the sequence phasors of a three-phase unbalanced system.

    Parameters
    ----------
    x_a: complex
        Phasor A in a system of three unbalanced phasors.
    x_b: complex
        Phasor B in a system of three unbalanced phasors.
    x_c: complex
        Phasor C in a system of three unbalanced phasors.

    Returns
    -------
    (complex, complex, complex)
        Sequence phasors x_a0, x_a1, and x_a2.

    Examples
    --------
    >>> print('  '.join(format(u, '.4f') for u in ABC_to_A012(2760, -1035-1521j, -1740+1535j)))
    -5.0000+4.6667j  2264.6912+201.1826j  500.3088-205.8493j
    """
    a = complex(-1 / 2, np.sqrt(3) / 2)
    a2 = a.conjugate()
    return ((x_a + x_b + x_c) / 3,
            (x_a + a * x_b + a2 * x_c) / 3,
            (x_a + a2 * x_b + a * x_c) / 3)


def A012_to_ABC(x_a0, x_a1, x_a2):
    """
    Find a three-phase unbalanced system of a sequence phasors.

    Parameters
    ----------
    x_a0: complex
        Homopolar sequence phasor.
    x_a1: complex
        Direct sequence phasor.
    x_a2: complex
        Inverse sequence phasor.

    Returns
    -------
    (complex, complex, complex)
        System of three unbalanced phasors x_a, x_b and x_c.

    Examples
    --------
    >>> print('  '.join(format(u, '.4f') for u in A012_to_ABC(-5+4j, 2264+201j, 500-205j)))
    2759.0000+0.0000j  -1035.3937-1521.6688j  -1738.6063+1533.6688j
    """
    a = complex(-1 / 2, np.sqrt(3) / 2)
    a2 = a.conjugate()
    return (x_a0 + x_a1 + x_a2,
            x_a0 + a2 * x_a1 + a * x_a2,
            x_a0 + a * x_a1 + a2 * x_a2)


def AN12_to_AB12(x_an_1, x_an_2):
    """
    Find the line-line direct and inverse sequence phasors
    corresponding to the line-neutral direct and inverse
    sequence phasors.

    Parameters
    ----------
    x_an_1: complex
        Line-neutral direct sequence phasor.
    x_an_2: complex
        Line-neutral inverse sequence phasor.

    Returns
    -------
    (complex, complex)
        Line-line sequence phasors x_ab_1, and x_ab_2.

    Examples
    --------
    >>> print('  '.join(format(u, '.4f') for u in AN12_to_AB12(1187-552j, 308+44j)))
    2258.5460+199.9722j  500.1051-200.7358j
    """
    return (x_an_1 * np.sqrt(3) * np.exp(1j * np.pi / 6),
            x_an_2 * np.sqrt(3) * np.exp(-1j * np.pi / 6))


def AB12_to_AN12(x_ab_1, x_ab_2):
    """
    Find the line-neutral direct and inverse sequence phasors
    corresponding to the line-line direct and inverse
    sequence phasors

    Parameters
    ----------
    x_ab_1: complex
        Line-line direct sequence phasor.
    x_ab_2: complex
        Line-line inverse sequence phasor.

    Returns
    -------
    (complex, complex)
        Line-neutral sequence phasors x_an_1, and x_an_2.

    Examples
    --------
    >>> print('  '.join(format(u, '.4f') for u in AB12_to_AN12(2258+200j, 500-200j)))
    1186.7350-551.8285j  307.7350+44.3376j
    """
    return (x_ab_1 / np.sqrt(3) / np.exp(1j * np.pi / 6),
            x_ab_2 / np.sqrt(3) / np.exp(-1j * np.pi / 6))


def ezs_u(e, z, s):
    """
    Find the voltage at a load, given the source
    voltage, the impedance between source and load, and
    the power drawn by the load.

    Parameters
    ----------
    e: float|complex
        Source voltage.
    z: float|complex
        Impedance between source and load.
    s: float|complex
        Power drawn by the load.

    Returns
    -------
    float|complex
        Voltage at the load.

    Examples
    --------
    >>> print(f'{ezs_u(0.4+0.3j, 0.1j, 0.6+0.45j):.4f}')
    0.3165+0.0874j
    >>> print(f'{ezs_u(125, 0.1, 100):.4f}')
    124.9199
    """
    is_complex = any(isinstance(var, complex) for var in [e, z, s])
    ZcjS = z.conjugate() * s
    ImUe = ZcjS.imag / abs(e)
    E_abs = abs(e)
    ReUe = np.polynomial.Polynomial([ZcjS.real + ImUe ** 2, -E_abs, 1]).roots()

    if isinstance(ReUe[0], complex):
        print('No solution found')
        return None

    U = complex(max(ReUe), ImUe) * e / E_abs
    return U if is_complex else U.real


def CEI_51_curve(I, curve_type, I_p, T=1.0, C=0.0, B=0.0):
    """
    Compute the actuation time 't' of a CEI 51 curve,
    for a given current 'I', according to:
        t = T*(k/((I/I_p)**a - 1) + C) + B,
    where 'k' and 'a' are set by the 'curve_type' parameter.

    Parameters
    ----------
    I: float
        Current, in A.
    curve_type: str
        The CEI 51 curve to use. Allowed values are (uppercase or lowercase):
            'I'  → Inverse
            'VI' → Very inverse
            'EI' → Extremely inverse
    I_p: float
        Pick up current of the CEI 51 curve, in A.
    T: float (optional)
        Parameter of the CEI 51 curve. The default is 1.
    C: float (optional)
        Parameter of the CEI 51 curve. The default is 0.
    B: float (optional)
        Parameter of the CEI 51 curve. The default is 0.

    Returns
    -------
    float
        The actuation time of the CEI 51 curve, in s.

    Examples
    --------
    >>> print(f'{CEI_51_curve(200, "VI", 100):.4f}')
    13.5000
    >>> CEI_51_curve(200, 'Inverse', 100)
    Unknown curve "Inverse"
    0
    """
    try:
        match curve_type.upper():
            case 'I':
                k = 0.14
                a = 0.02
            case 'VI':
                k = 13.5
                a = 1.0
            case 'EI':
                k = 80.0
                a = 2.0
            case _:
                raise AttributeError

        return T * (k / ((I / I_p) ** a - 1) + C) + B

    except AttributeError:
        print(f'Unknown curve "{curve_type}"')
        return 0


def IEEE_51_curve(I, curve_type, I_p, T=1.0, B=0.0):
    """
    Compute the actuation time 't' 'of an IEEE 51 curve,
    for a given current 'I', according to:
        t = T*(k/((I/I_p)**a - 1) + C) + B,
    where 'k', 'a' and 'C' are set by the 'curve_type' parameter.

    Parameters
    ----------
    I: float
        Current, in A.
    curve_type: str
        The IEEE 51 curve to use. Allowed values are (uppercase or lowercase):
            'MI' → Moderately inverse
            'VI' → Very inverse
            'EI' → Extremely inverse
    I_p: float
        Pick up current of the IEEE 51 curve, in A.
    T: float (optional)
        Parameter of the IEEE 51 curve. The default is 1.
    B: float (optional)
        Parameter of the IEEE 51 curve. The default is 0.

    Returns
    -------
    float
        The actuation time of the IEEE 51 curve, in s.

    Examples
    --------
    >>> print(f'{IEEE_51_curve(200, "VI", 100):.4f}')
    7.0277
    >>> IEEE_51_curve(200, 'Inverse', 100)
    Unknown curve "Inverse"
    0
    """
    try:
        match curve_type.upper():
            case 'MI':
                k = 0.0515
                a = 0.02
                C = 0.114
            case 'VI':
                k = 19.61
                a = 2.0
                C = 0.491
            case 'EI':
                k = 28.2
                a = 2.0
                C = 0.1217
            case _:
                raise AttributeError

        return T * (k / ((I / I_p) ** a - 1) + C) + B

    except AttributeError:
        print(f'Unknown curve "{curve_type}"')
        return 0


class NetworkException(Exception):
    """Base exception class for all the class Network exceptions."""


class NetworkNotSolved(NetworkException):
    """
    Exception raised in class Network when results want to be
    obtained before the network is solved.
    """

    def __init__(self):
        message = 'Network not solved yet!'
        super().__init__(message)


class NetworkMissingBranch(NetworkException):
    """
    Exception raised in class Network when trying to access missing
    branches.
    """

    def __init__(self, missing_branches):
        if len(missing_branches) == 1:
            message = f'Branch {missing_branches} is invalid or has no information!'
        else:
            message = f'Branches {missing_branches} are invalid or have no information!'
        super().__init__(message)
        self.missing_branches = missing_branches


class NetworkMissingNode(NetworkException):
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


class NetworkMissingMutualCoupling(NetworkException):
    """
    Exception raised in class Network when trying to access missing mutual couplings.
    """

    def __init__(self, missing_mutual_coupling):
        message = f'Mutual coupling between branches {missing_mutual_coupling} does not exist!'
        super().__init__(message)
        self.missing_mutual_coupling = missing_mutual_coupling


class NetworkInvalidElement(NetworkException):
    """
    Exception raised in class Network when trying to access an invalid element
    """

    def __init__(self, invalid_element):
        message = f'Element "{invalid_element}" is not valid'
        super().__init__(message)
        self.invalid_element = invalid_element


class NetworkNotSolvable(NetworkException):
    """
    Exception raised in class Network when trying to solve an invalid or empty network.
    Examples: Two ideal voltage sources in parallel, two ideal current sources
    in series, a short-circuited ideal voltage source, or an open current source.
    """

    def __init__(self, not_solvable):
        message = f'Non solvable network ({not_solvable})'
        super().__init__(message)
        self.not_solvable = not_solvable


@dataclass
class BranchElement:
    """Base data class for branch elements."""

    from_to: tuple[int, int]  # Origin and destination nodes in an oriented graph

    @abstractmethod
    def is_complex(self):
        pass


@dataclass
class Impedance(BranchElement):
    """Create an impedance to be used with class Network."""

    Z: complex  # Impedance

    def is_complex(self):
        return isinstance(self.Z, complex)


@dataclass
class Admittance(BranchElement):
    """Create an admittance to be used with class Network."""

    Y: complex  # Admittance

    def is_complex(self):
        return isinstance(self.Y, complex)


@dataclass
class VoltageSource(BranchElement):
    """Create a voltage source to be used with class Network."""

    E: complex  # Source voltage in an oriented graph
    Z: complex  # Source series impedance

    def is_complex(self):
        return isinstance(self.E, complex) or isinstance(self.Z, complex)


@dataclass
class VoltageSourceIdeal(BranchElement):
    """Create an ideal voltage source to be used with class Network."""

    E: complex  # Source voltage in an oriented graph

    def is_complex(self):
        return isinstance(self.E, complex)


@dataclass
class CurrentSource(BranchElement):
    """Create a current source to be used with class Network."""

    J: complex  # Source current in an oriented graph
    Y: complex  # Source parallel admittance

    def is_complex(self):
        return isinstance(self.J, complex) or isinstance(self.Y, complex)


@dataclass
class CurrentSourceIdeal(BranchElement):
    """Create an ideal current source to be used with class Network."""

    J: complex  # Source current in an oriented graph

    def is_complex(self):
        return isinstance(self.J, complex)


@dataclass
class ShortCircuit(BranchElement):
    """Create a short circuit to be used with class Network."""

    def is_complex(self):
        return False


@dataclass
class MutualCoupling:
    """Create a mutual coupling to be used with class Network."""

    XM: complex  # Mutual coupling impedance between two branches in an oriented graph

    def is_complex(self):
        return True


class Network:
    """
    Solve an electrical networks using the nodes method.

    >>> net = Network()
    >>> net.add(VoltageSource(from_to=(0, 1), E=200, Z=10), branch=1)
    >>> net.add(VoltageSource(from_to=(1, 2), E=-50, Z=20j), branch=2)
    >>> net.add(Impedance(from_to=(2, 0), Z=5j), branch=3)
    >>> net.add(Impedance(from_to=(1, 2), Z=20), branch=4)
    >>> net.add(CurrentSource(from_to=(0, 2), J=4, Y=1/10), branch=5)
    >>> net.add(MutualCoupling(XM=5j), coupled_branches=(2, 3))
    >>> net.remove(branch=2)
    >>> net.remove(coupled_branches=(3, 2))
    >>> net.add(VoltageSource(from_to=(1, 2), E=-50, Z=20j), branch=2)
    >>> net.add(MutualCoupling(XM=5j), coupled_branches=(2, 3))
    >>> net.solve()
    >>> net.num_branches
    5
    >>> net.num_nodes
    3
    >>> net.is_solved
    True
    >>> print(net)
    Branch    Element
    --------  -------------------------------------------
     2, 3     MutualCoupling(XM=5j)
      1       VoltageSource(from_to=(0, 1), E=200, Z=10)
      2       VoltageSource(from_to=(1, 2), E=-50, Z=20j)
      3       Impedance(from_to=(2, 0), Z=5j)
      4       Impedance(from_to=(1, 2), Z=20)
      5       CurrentSource(from_to=(0, 2), J=4, Y=0.1)
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

    def __init__(self):
        self._net = {}
        self._is_solved = False
        self._is_complex = 0
        self._extra_nodes = 0
        self._extra_branches = 0
        self._list_nodes = array.array('Q')
        self._list_branches = array.array('Q')
        self._list_coupled_branches = array.array('Q')
        self._ZN = None
        self._VN = None
        self._IB = None

    def __str__(self):
        def sort_net(tup):
            k, _ = tup
            return k if isinstance(k, int) else -min(k)

        s = [['Branch', 'Element']]
        for branch, element in sorted(self._net.items(), key=sort_net):
            if isinstance(branch, int):
                s.append([f'{branch:^6d}', f'{element}'])
            else:
                b1, b2 = branch
                s.append([f'{b1}, {b2}'.center(6), f'{element}'])
        return tabulate.tabulate(s, headers='firstrow')

    @property
    def num_branches(self):
        """
        Return the umber of branches in the network.

        Exception NetworkNotSolved is raised if the network is not solved.
        """
        if not self._is_solved: raise NetworkNotSolved

        return len(set(self._list_branches))

    @property
    def num_nodes(self):
        """
        Return the number of nodes in the network.

        Exception NetworkNotSolved is raised if the network is not solved.
        """
        if not self._is_solved: raise NetworkNotSolved

        return len(set(self._list_nodes))

    @property
    def is_solved(self):
        """Return True if the network has been solved, and False otherwise."""
        return self._is_solved

    def add(self, element, *, branch=None, coupled_branches=None):
        """
        Add an element to the network.

        Allowed 'element' to be added to 'branch' are subclasses from
        class BranchElement.
        Allowed 'element' to be added between 'coupled_branches'
        is class MutualCoupling.
        Exception NetworkInvalidElement is raised if 'element' is not valid.
        """
        if coupled_branches is not None:
            if not isinstance(element, MutualCoupling):
                raise NetworkInvalidElement(element)
            branches = frozenset(coupled_branches)
            if branches in self._net:
                self.remove(coupled_branches=coupled_branches)
            self._net[branches] = element
            self._is_solved = False
            self._is_complex += int(element.is_complex())
            self._list_coupled_branches.extend(coupled_branches)
        elif branch is not None:
            if not isinstance(element, BranchElement):
                raise NetworkInvalidElement(element)
            if branch in self._net:
                self.remove(branch=branch)
            self._net[branch] = element
            self._is_solved = False
            self._is_complex += int(element.is_complex())
            self._list_branches.append(branch)
            self._list_nodes.extend(element.from_to)
            match element:
                case VoltageSourceIdeal() | ShortCircuit():
                    self._extra_branches += 1
                    self._extra_nodes += 1
                case CurrentSourceIdeal():
                    self._extra_branches += 1
        else:
            pass

    def remove(self, *, branch=None, coupled_branches=None):
        """
        Remove an element from the network.

        Exception NetworkMissingBranch is raised if 'branch' does not exist,
        when 'coupled_branches' is None.
        Exception NetworkMissingMutualCoupling is raised if 'coupled_branches'
        do not exit, when 'branch' is None.
        """
        if coupled_branches is not None:
            try:
                element = self._net.pop(frozenset(coupled_branches))
                self._is_solved = False
                self._is_complex -= int(element.is_complex())
                for b in coupled_branches:
                    self._list_coupled_branches.remove(b)
            except KeyError as err:
                raise NetworkMissingMutualCoupling({coupled_branches}) from err
        elif branch is not None:
            try:
                element = self._net.pop(branch)
                self._is_solved = False
                self._is_complex -= int(element.is_complex())
                for n in element.from_to:
                    self._list_nodes.remove(n)
                self._list_branches.remove(branch)
                match element:
                    case VoltageSourceIdeal() | ShortCircuit():
                        self._extra_branches -= 1
                        self._extra_nodes -= 1
                    case CurrentSourceIdeal():
                        self._extra_branches -= 1
            except KeyError as err:
                raise NetworkMissingBranch({branch}) from err
        else:
            pass

    def solve(self):
        """
        Solve the network.

        The voltage across every branch, and the current through every branch
        is calculated.
        Exception NetworkMissingBranch is raised if not all branches
        1, 2, 3, ... max(self._branches) are present.
        Exception NetworkMissingNode is raised if not all nodes
        0, 1, 2, ... max(self._nodes) are present.
        Exception NetworkNotSolvable is raised if linalg.inv raises the
        linalg.LinAlgError exception.
        """
        # Check for missing nodes and branches, set the number of nodes and
        # branches, and check the data type (complex or real). Computes also
        # the number of extra nodes anb branches needed for ideal elements.
        if not self._net:
            raise NetworkNotSolvable('empty network')

        actual_b = set(self._list_branches)
        num_b = max(actual_b)
        expected_b = set(b for b in range(1, num_b + 1))
        if diff := actual_b ^ expected_b:
            raise NetworkMissingBranch(diff)
        if diff := set(self._list_coupled_branches) - expected_b:
            raise NetworkMissingBranch(diff)

        actual_n = set(self._list_nodes)
        num_n = max(actual_n)
        expected_n = set(n for n in range(0, num_n + 1))
        if diff := actual_n ^ expected_n:
            raise NetworkMissingNode(diff)

        # Set up the matrices needed to solve the network
        total_n = num_n + self._extra_nodes
        total_b = num_b + self._extra_branches
        A = np.zeros((total_n, total_b), dtype=np.int_)
        data_type = np.complex_ if self._is_complex else np.float_
        ZB = np.zeros((total_b, total_b), dtype=data_type)
        EB1 = np.zeros(total_b, dtype=data_type)
        JB1 = np.zeros(total_b, dtype=data_type)

        def fill_A(branch, from_node, to_node):
            if from_node > 0:
                A[from_node - 1, branch - 1] = 1
            if to_node > 0:
                A[to_node - 1, branch - 1] = -1

        for branch, element in self._net.items():
            match element:
                case Impedance():
                    fill_A(branch, *element.from_to)
                    ZB[branch - 1, branch - 1] = element.Z
                case Admittance():
                    fill_A(branch, *element.from_to)
                    ZB[branch - 1, branch - 1] = 1 / element.Y
                case VoltageSource():
                    fill_A(branch, *element.from_to)
                    EB1[branch - 1] = element.E
                    ZB[branch - 1, branch - 1] = element.Z
                case VoltageSourceIdeal():
                    num_b += 1
                    num_n += 1
                    fill_A(branch, element.from_to[0], num_n)
                    fill_A(num_b, num_n, element.from_to[1])
                    EB1[branch - 1] = element.E
                    ZB[branch - 1, branch - 1] = 1
                    ZB[num_b - 1, num_b - 1] = -1
                case CurrentSource():
                    fill_A(branch, *element.from_to)
                    JB1[branch - 1] = element.J
                    ZB[branch - 1, branch - 1] = 1 / element.Y
                case CurrentSourceIdeal():
                    num_b += 1
                    fill_A(branch, *element.from_to)
                    fill_A(num_b, *element.from_to)
                    JB1[branch - 1] = element.J
                    ZB[branch - 1, branch - 1] = 1
                    ZB[num_b - 1, num_b - 1] = -1
                case ShortCircuit():
                    num_b += 1
                    num_n += 1
                    fill_A(branch, element.from_to[0], num_n)
                    fill_A(num_b, num_n, element.from_to[1])
                    ZB[branch - 1, branch - 1] = 1
                    ZB[num_b - 1, num_b - 1] = -1
                case MutualCoupling():
                    b1, b2 = branch
                    ZB[b1 - 1, b2 - 1] = element.XM
                    ZB[b2 - 1, b1 - 1] = element.XM

        # Solve the network
        try:
            YB = linalg.inv(ZB, overwrite_a=True)
            JB = YB @ EB1 + JB1
            self._ZN = linalg.inv(A @ YB @ A.T, overwrite_a=True)
            self._VN = self._ZN @ (-A @ JB)
            self._IB = YB @ (A.T @ self._VN) + JB
            # Adjust result for branches with ideal current sources
            for branch, element in self._net.items():
                if isinstance(element, CurrentSourceIdeal):
                    self._IB[branch - 1] = element.J
            self._is_solved = True
        except linalg.LinAlgError as err:
            raise NetworkNotSolvable(err) from err

    def thevenin(self, node_1, node_2):
        """
        Return the Thevenin equivalent circuit between two nodes.

        Exception NetworkNotSolved is raised if the network is not solved.
        Exception NetworkMissingNode is raised if 'node_1' or 'node_2'
        do not exist.
        """
        if not self._is_solved: raise NetworkNotSolved
        if node_1 not in self._list_nodes:
            raise NetworkMissingNode({node_1})
        if node_2 not in self._list_nodes:
            raise NetworkMissingNode({node_2})

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

    def norton(self, node_1, node_2):
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

    def voltage(self, *, branch=None, node=None):
        """
        Return the voltage across a branch, or the voltage of a node (with respect to node 0).

        Exception NetworkNotSolved is raised if the network is not solved.
        Exception NetworkMissingBranch is raised if 'branch' does not exist.
        Exception NetworkMissingNode is raised if 'node' does not exist.
        """
        if not self._is_solved:
            raise NetworkNotSolved

        if branch is not None:
            if branch not in self._list_branches:
                raise NetworkMissingBranch({branch})
            return self.thevenin(*self._net[branch].from_to)[0]
        elif node is not None:
            return self.thevenin(node, 0)[0]
        else:
            raise NetworkMissingBranch({})

    def current(self, branch):
        """
        Return the current through a branch.

        Exception NetworkNotSolved is raised if the network is not solved.
        Exception NetworkMissingBranch is raised if 'branch' does not exist.
        """
        if not self._is_solved:
            raise NetworkNotSolved
        if diff := {branch} - set(self._list_branches):
            raise NetworkMissingBranch(diff)

        return self._IB[branch - 1]

    def results(self, polar=False, u_fmt='', u_unit='', u_scale=1,
                i_fmt='', i_unit='', i_scale=1):
        """
        Return a multi-line string formatted as a table, with all the network's voltages and currents.

        A multi-line formatted string is returned with the voltages and
        currents of all the network branches.
        A different scale factor can be applied to voltages and currents.
        Polar display can be chosen, and different formats for voltages
        and currents are possible. Also, different units can be added
        to voltages and currents.
        Exception NetworkNotSolved is raised if the network is not solved.
        """
        if not self._is_solved:
            raise NetworkNotSolved

        if u_unit:
            u_unit = ' / ' + u_unit
        if i_unit:
            i_unit = ' / ' + i_unit
        s = [['Branch', 'Voltage' + u_unit, 'Current' + i_unit]]
        for branch in range(1, self.num_branches + 1):
            u = self.voltage(branch=branch) * u_scale
            i = self.current(branch) * i_scale
            if polar:
                u = ComplexD(u)
                i = ComplexD(i)
            s.append([branch, format(u, u_fmt), format(i, i_fmt)])
        return tabulate.tabulate(s, numalign='center', headers='firstrow')

    def save(self, filename):
        """Save the contents of the network to a binary file."""
        try:
            with open(filename, 'wb') as fh:
                pickle.dump(self._net, fh, pickle.HIGHEST_PROTOCOL)
                pickle.dump(self._is_solved, fh, pickle.HIGHEST_PROTOCOL)
                pickle.dump(self._is_complex, fh, pickle.HIGHEST_PROTOCOL)
                pickle.dump(self._extra_nodes, fh, pickle.HIGHEST_PROTOCOL)
                pickle.dump(self._extra_branches, fh, pickle.HIGHEST_PROTOCOL)
                pickle.dump(self._list_nodes, fh, pickle.HIGHEST_PROTOCOL)
                pickle.dump(self._list_branches, fh, pickle.HIGHEST_PROTOCOL)
                pickle.dump(self._list_coupled_branches, fh, pickle.HIGHEST_PROTOCOL)
                pickle.dump(self._ZN, fh, pickle.HIGHEST_PROTOCOL)
                pickle.dump(self._VN, fh, pickle.HIGHEST_PROTOCOL)
                pickle.dump(self._IB, fh, pickle.HIGHEST_PROTOCOL)
                return True
        except (OSError, pickle.PicklingError) as err:
            print(f'{filename}: save error: {err}')
            return False

    def load(self, filename):
        """Load the network with the contents of a binary file."""
        try:
            with open(filename, 'rb') as fh:
                self._net = pickle.load(fh)
                self._is_solved = pickle.load(fh)
                self._is_complex = pickle.load(fh)
                self._extra_nodes = pickle.load(fh)
                self._extra_branches = pickle.load(fh)
                self._list_nodes = pickle.load(fh)
                self._list_branches = pickle.load(fh)
                self._list_coupled_branches = pickle.load(fh)
                self._ZN = pickle.load(fh)
                self._VN = pickle.load(fh)
                self._IB = pickle.load(fh)
                return True
        except (OSError, pickle.UnpicklingError) as err:
            print(f'{filename}: load error: {err}')
            return False


@dataclass(frozen=True)
class Motor3phRun:
    """Class used to store the result of method 'run' from class Motor3ph."""

    Z: np.ndarray  # Total motor impedance (Line-Neutral), Ω
    PF: np.ndarray  # Motor Power factor
    I: np.ndarray  # Current drawn by the motor, A
    U: np.ndarray  # Voltage at the motor terminals (Line-Neutral), V
    S: np.ndarray  # Apparent power drawn by the motor, VA
    T_m: np.ndarray  # Mechanical torque delivered by the motor, N·m
    P_m: np.ndarray  # Mechanical power delivered by the motor, W
    Eff: np.ndarray  # Motor Efficiency (between 0 and 1)


@dataclass(frozen=True)
class Motor3phStartUp:
    """Class used to store the result of method 'start_up' from class Motor3ph."""

    t_start: float  # The starting time of the motor, s
    s_start_up: interpolate.interp1d  # An interpolating function returning the slip for a given time t in seconds
    # For t < 0, the slip at t = 0 is returned
    # for t > t_start, the slip at t = t_start is returned


class Motor3ph:
    """
    Compute electrical and mechanical quantities of 3-phase asynchronous motors.

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
    >>> mot = Motor3ph(p, f, R_1, X_1, R_2, X_2, R_Fe, X_m)
    >>> mot_n = mot.run(s_n, U)
    >>> print(f'n_sync = {mot.n_m_sync:.0f} r/min')
    n_sync = 1500 r/min
    >>> print(f'I = {abs(mot_n.I):.1f} A')
    I = 32.3 A
    >>> print(f'T_m = {mot_n.T_m:.1f} N·m')
    T_m = 118.8 N·m
    """

    def __init__(self, p, f, r_1, x_1, r_2, x_2, r_fe, x_m):
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
    def ω_m_sync(self):
        """Return the motor synchronous speed in rad/s."""
        return self._ω_m_sync

    @property
    def n_m_sync(self):
        """Return the motor synchronous speed in r/min."""
        return self._n_m_sync

    def to_s(self, ω_m=None, n_m=None):
        """
        Return the slip corresponding to a mechanical speed
        given in rad/s (ω_m) or in r/min (n_m).
        """
        if ω_m is not None:
            return 1 - ω_m / self._ω_m_sync
        elif n_m is not None:
            return 1 - n_m / self._n_m_sync
        else:
            return -1

    def to_ω_m(self, s=None, n_m=None):
        """
        Return the mechanical speed in rad/s corresponding to a
        given slip (s) or a mechanical speed in r/min (n_m).
        """
        if s is not None:
            return (1 - s) * self._ω_m_sync
        elif n_m is not None:
            return n_m * np.pi / 30
        else:
            return -1

    def to_n_m(self, s=None, ω_m=None):
        """
        Return the mechanical speed in r/min corresponding to a
        given slip (s) or a mechanical speed in rad/s (ω_m).
        """
        if s is not None:
            return (1 - s) * self._n_m_sync
        elif ω_m is not None:
            return ω_m * 30 / np.pi
        else:
            return -1

    def s_T_m_max(self, z_ext=0):
        """
        Return the slip at the maximum motor mechanical torque.

        Parameters
        ----------
        z_ext: Impedance of the external power system feeding the motor, Ω/phase
        """
        z_th = z_parallel([self._z_0, self._z_1 + z_ext])
        return self._r_2 / abs(z_th + self._jx_2)

    def s_match(self, T_load, e_ext, z_ext=0, s_guess=0.04):
        """
        Return the slip at which the motor torque matches the load torque.

        Parameters
        ----------
        T_load:  A function that returns the load torque for a given slip, N·m
        e_ext:   Voltage (phase-neutral) of the external power system feeding the motor, V
        z_ext:   Impedance of the external power system feeding the motor, Ω/phase
        s_guess: Slip initial guess
        """
        sol = optimize.root(lambda s: self.run(s, e_ext, z_ext).T_m - T_load(s), x0=s_guess)
        return sol.x[0]

    def run(self, s, e_ext, z_ext=0):
        """
        Return several motor quantities when running at a given slip, in a Motor3phRun class.

        Parameters
        ----------
        s:     Slip at which the motor quantities are to be computed
        e_ext: Voltage (line-neutral) of the external power system feeding the motor, V
        z_ext: Impedance of the external power system feeding the motor, Ω/phase

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

    def start_up(self, T_load, J, e_ext, z_ext=0):
        """
        Return the starting time of the motor, and a function of the time evolution of
        the slip during the starting up process, in a Motor3phStartUp class.

        Parameters
        ----------
        T_load: A function that returns the load torque for a given slip, N·m
        J:      Total moment of inertia of the motor plus the load, kg·m²
        e_ext:  Voltage (line-neutral) of the external power system feeding the motor, V
        z_ext:  Impedance of the external power system feeding the motor, Ω/phase

        Returns
        -------
        A Motor3phStartUp class.
        """
        n_p1, n_p2 = 200, 800
        s_lin = 0.3  # The function to integrate is very smooth up to approximately this value
        s_match = self.s_match(T_load, e_ext, z_ext)
        s = np.concatenate([np.linspace(1, s_lin, n_p1, endpoint=False),
                            np.linspace(s_lin, s_match, n_p2, endpoint=False)])
        mot_run = self.run(s, e_ext, z_ext)
        t = -J * self._ω_m_sync * integrate.cumulative_trapezoid(1 / (mot_run.T_m - T_load(s)), s, initial=0)
        t_start = t[-1]
        s_start_up = interpolate.interp1d(t, s, kind='cubic', copy=False,
                                          assume_sorted=True, bounds_error=False,
                                          fill_value=(s[0], s[-1]))
        return Motor3phStartUp(t_start, s_start_up)


if __name__ == '__main__':
    import doctest

    doctest.testmod()
