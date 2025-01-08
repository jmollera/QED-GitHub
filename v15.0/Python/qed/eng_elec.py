__author__ = 'Josep Mollera'
__version__ = '15.0'
__all__ = ['SQRT3', 'CableType', 'IECcurve', 'IEEEcurve', 'Speed', 'EEException',
           'EEInvalidArguments', 'EENetNotSolved', 'EENetMissingBranch',
           'EENetMissingNode', 'EENetMissingMutualCoupling',
           'EENetNotSolvable', 'z_cable', 'r_cable', 'voltage_drop',
           'kcmil_to_mm2', 'mm2_to_kcmil', 'awg_to_mm2', 'mm2_to_awg',
           'z_series', 'z_parallel', 'millman', 'rms', 'average',
           'apparent_power', 'D_to_Y', 'Y_to_D', 'triangle_to_phasors',
           'LN_to_LL', 'LL_to_LG', 'LL_to_LN', 'ABC_to_A012', 'A012_to_ABC',
           'AN12_to_AB12', 'AB12_to_AN12', 'ezs_u', 'icc_tr', 'curve_51',
           'Impedance', 'Admittance', 'VoltageSource', 'CurrentSource',
           'ShortCircuit', 'MutualCoupling', 'Network', 'Motor3phRun', 'Motor3ph']

import array
import tabulate
import enum
import math
import cmath
from dataclasses import dataclass
from collections.abc import Callable
import pickle
import numpy as np
from numpy.typing import ArrayLike
from scipy import linalg, integrate, optimize, interpolate
from qed.utils import Complex

tabulate.PRESERVE_WHITESPACE = True

SQRT3 = math.sqrt(3)  # The square root of three

# Two dictionaries of impedances for copper cables.
#
# Key: section in mm²
# Value: impedance in Ω/km at 90 °C
_Z_CABLE_A = {4: 5.99 + 0.098j, 6: 3.96 + 0.123j, 10: 2.34 + 0.117j,
              16: 1.47 + 0.112j, 25: 0.936 + 0.107j, 35: 0.661 + 0.123j,
              50: 0.5 + 0.11j, 70: 0.34 + 0.114j, 120: 0.192 + 0.108j,
              150: 0.156 + 0.104j, 185: 0.125 + 0.103j, 240: 0.095 + 0.102j,
              300: 0.078 + 0.094j}
_Z_CABLE_V = {4: 5.993 + 0.168j, 6: 3.94 + 0.159j, 10: 2.33 + 0.151j,
              16: 1.46 + 0.142j, 25: 0.9361 + 0.0791j, 35: 0.6749 + 0.0765j,
              50: 0.4935 + 0.135j, 70: 0.3417 + 0.128j, 95: 0.249 + 0.130j,
              120: 0.1980 + 0.0825j}

# Dictionary of resistivity and coefficients of resistivity
# variation with temperature, for several materials.
#
# Key: chemical symbol of the material
# Value: resistivity in Ω·mm²/m at 20 °C, and coefficient of
#        resistivity variation with temperature in 1/°C at 20 °C
_R_MATERIAL = {'Al': (0.02825, 0.00391), 'Cu': (0.01723, 0.00393),
               'Ag': (0.01645, 0.00380), 'Au': (0.02440, 0.00340)}


class CableType(enum.Enum):
    """Constants to be used in function z_cable"""
    A = enum.auto()
    V = enum.auto()


class IECcurve(enum.Enum):
    """Constants to be used in function curve_51"""
    Inverse = enum.auto()
    Very_Inverse = enum.auto()
    Extremely_Inverse = enum.auto()


class IEEEcurve(enum.Enum):
    """Constants to be used in function curve_51"""
    Moderately_Inverse = enum.auto()
    Very_Inverse = enum.auto()
    Extremely_Inverse = enum.auto()


class Speed(enum.Enum):
    """Constants to be used in method speed of class Motor3ph"""
    slip_to_rpm = enum.auto()        # slip to r/min
    slip_to_rad_per_s = enum.auto()  # slip to rad/s
    rpm_to_slip = enum.auto()        # r/min to slip
    rpm_to_rad_per_s = enum.auto()   # r/min to rad/s
    rad_per_s_to_slip = enum.auto()  # rad/s to slip
    rad_per_s_to_rpm = enum.auto()   # rad/s to r/min


class EEException(Exception):
    """Base exception class for eng_elec module."""


class EEInvalidArguments(EEException):
    """
    Exception raised when a function is called with invalid arguments.
    """

    def __init__(self, invalid_args: dict) -> None:
        invalid_args.pop('self', None)  # Drop "self" if present
        message = f'{invalid_args}'
        super().__init__(message)
        self.invalid_args = invalid_args


class EENetNotSolved(EEException):
    """
    Exception raised in class Network when results want to be
    obtained before the network is solved.
    """

    def __init__(self) -> None:
        message = 'Solve the network first!'
        super().__init__(message)


class EENetMissingBranch(EEException):
    """
    Exception raised in class Network when trying to access missing branches.
    """

    def __init__(self, missing_branches: set) -> None:
        pl1, pl2, pl3 = ('es', 'are', 'have') if len(missing_branches) > 1 else ('', 'is', 'has')
        message = f'Branch{pl1} {missing_branches} {pl2} invalid or {pl3} no information!'
        super().__init__(message)
        self.missing_branches = missing_branches


class EENetMissingNode(EEException):
    """
    Exception raised in class Network when trying to access missing nodes.
    """

    def __init__(self, missing_nodes: set) -> None:
        pl1, pl2, pl3 = ('s', 'are', 'have') if len(missing_nodes) > 1 else ('', 'is', 'has')
        message = f'Node{pl1} {missing_nodes} {pl2} invalid or {pl3} no connections!'
        super().__init__(message)
        self.missing_nodes = missing_nodes


class EENetMissingMutualCoupling(EEException):
    """
    Exception raised in class Network when trying to access missing mutual couplings.
    """

    def __init__(self, missing_mutual_coupling: set) -> None:
        message = f'Mutual coupling between branches {missing_mutual_coupling} does not exist!'
        super().__init__(message)
        self.missing_mutual_coupling = missing_mutual_coupling


class EENetNotSolvable(EEException):
    """
    Exception raised in class Network when trying to solve an invalid or empty network.
    Examples: Two ideal voltage sources in parallel, two ideal current sources
    in series, a short-circuited ideal voltage source, or an open current source.
    Raised also by function ezs_u when no solution in found.
    """

    def __init__(self, not_solvable: str) -> None:
        message = f'Non solvable network ({not_solvable})!'
        super().__init__(message)
        self.not_solvable = not_solvable


def z_cable(section: int, length: float, cable: CableType = CableType.A) -> complex:
    """
    Find the impedance of a cable.

    Parameters
    ----------
    section: Cable section, mm².
    length:  Cable length, m.
    cable:   CABLE.A → use dictionary _Z_CABLE_A,
             CABLE.V → use dictionary _Z_CABLE_V,
             The default is CABLE.A.

    Returns
    -------
    The cable impedance, Ω at 90 °C.
    Exception EEInvalidArguments may be raised.

    Examples
    --------
    >>> z_cable(25, 100)
    (0.0936+0.0107j)
    >>> z_cable(25, 100, CableType.A)
    (0.0936+0.0107j)
    >>> z_cable(25,100, CableType.V)
    (0.09361+0.00791j)
    """
    try:
        match cable:
            case CableType.A:
                return _Z_CABLE_A[section] * length / 1000
            case CableType.V:
                return _Z_CABLE_V[section] * length / 1000
            case _:
                raise EEInvalidArguments({'cable': cable})
    except KeyError:
        raise EEInvalidArguments({'section': section})


def r_cable(section: float, length: float,
            temp: float = 20.0, material: str = 'Cu') -> float:
    """
    Find the resistance of a cable.

    Parameters
    ----------
    section:  Cable section, mm².
    length:   Cable length, m.
    temp:     Temperature, °C. The default is 20 °C.
    material: The chemical symbol of the cable material: 'Al', 'Cu', 'Ag' or 'Au'.
              The default is 'Cu'.

    Returns
    -------
    The cable resistance at the requested temperature, Ω.
    Exception EEInvalidArguments may be raised.

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
    """
    try:
        rho_20, alpha_20 = _R_MATERIAL[material]
        return rho_20 * (1 + alpha_20 * (temp - 20)) * length / section
    except KeyError:
        raise EEInvalidArguments({'material': material})


def voltage_drop(voltage: float, current: float,
                 impedance: complex, PF: float = 1.0,
                 ret_approx: bool = False) -> float | tuple[float, float]:
    """
    Find the voltage drop through a cable.

    Parameters
    ----------
    voltage:    Voltage at cable origin.
                (line-neutral voltage for a three-phase system).
    current:    Current through the cable.
    impedance:  Cable impedance.
    PF:         Load power factor. The default is 1.
    ret_approx: If True, return a tuple with the exact and the approximate values
                of the voltage drop; otherwise return only the exact value.
                The approximate value is usually good enough when PF > 0.8.
                The default is False.

    Returns
    -------
    Voltage drop through the cable (one or two values, depending on ret_approx).
    (line-neutral voltage for a three-phase system).

    Examples
    --------
    Three-phase: U=380 V (LL), I=200 A, Cable: 240 mm², 400 m. cos 𝜑=0.87
    >>> dv = voltage_drop(380/SQRT3, 200, z_cable(240, 400), 0.87) * SQRT3
    >>> print(f'{dv:.4f} V (L-L)')
    18.4652 V (L-L)

    One-phase: U=220 V, I=55 A, Cable: 35 mm², 200 m. cos 𝜑=0.85
    >>> dv = voltage_drop(220, 55, z_cable(35, 2*200), 0.85)
    >>> print(f'{dv:.4f} V')
    13.8515 V

    Direct Current: U=125 V, I=25 A, Cable: 16 mm², 50 m
    >>> dv = voltage_drop(125, 25, z_cable(16, 2*50).real)
    >>> print(f'{dv:.4f} V')
    3.6750 V
    """
    R = impedance.real
    X = impedance.imag
    sin_phi = math.sqrt(1.0 - PF ** 2)
    vd = ((R * PF + X * sin_phi)*current + voltage -
          math.sqrt(voltage**2 - ((X * PF - R * sin_phi)*current)**2))
    if ret_approx:
        return vd, current*(R*PF + X*sin_phi)
    else:
        return vd


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
    return 53.4751207321 / (92 ** (awg / 19.5))


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
    return round(4.31245284200 * math.log(53.4751207321 / mm2))


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
    return np.sum(np.asarray(u) / z) * z_parallel(z)


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
    return math.sqrt(val / (b - a))


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


def apparent_power(u: ArrayLike | complex, i: ArrayLike | complex) -> complex:
    """
    Compute the apparent power drawn by a load.

    Parameters
    ----------
    u: Voltages with respect to a reference point of a set of
       nodes feeding a load.
    i: Currents from the nodes towards the load.

    Returns
    -------
    Te apparent power drawn by the load.

    Examples
    --------
    >>> apparent_power([24+10j, 26+9j, 27+7j], [1+1j, 2+2j ,3+3j])
    (206-108j)
    >>> apparent_power(110, 3)
    330
    """
    return np.vdot(i, u)


def D_to_Y(z_ab: complex, z_bc: complex, z_ca: complex) -> tuple[complex, complex, complex]:
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


def Y_to_D(z_an: complex, z_bn: complex, z_cn: complex) -> tuple[complex, complex, complex]:
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
    angle_bac = math.acos((u_ab ** 2 + u_ca ** 2 - u_bc ** 2) / (2 * u_ab * u_ca))
    angle_abc = math.acos((u_bc ** 2 + u_ab ** 2 - u_ca ** 2) / (2 * u_bc * u_ab))
    return (complex(u_ab, 0),
            cmath.rect(u_bc, math.pi + angle_abc),
            cmath.rect(u_ca, math.pi - angle_bac))


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
    a = complex(-1, SQRT3) / 2
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
    a = complex(-1, SQRT3) / 2
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
    m = complex(3, SQRT3) / 2
    return x_an_1 * m, x_an_2 * m.conjugate()


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
    m = complex(3, SQRT3) / 2
    return x_ab_1 / m, x_ab_2 / m.conjugate()


def ezs_u(e: complex, z: complex, s: complex) -> complex:
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
    Exception EENetNotSolvable is raised if no solution is found.

    Examples
    --------
    >>> print(f'{ezs_u(0.4+0.3j, 0.1j, 0.6+0.45j):.4f}')
    0.3165+0.0874j
    >>> print(f'{ezs_u(125, 0.1, 100):.4f}')
    124.9199
    """
    Zconj_S = z.conjugate() * s
    Im_Ue = Zconj_S.imag / abs(e)
    E_abs = abs(e)
    Re_Ue = np.polynomial.Polynomial([Zconj_S.real + Im_Ue ** 2, -E_abs, 1]).roots()

    if isinstance(Re_Ue[0], complex):
        raise EENetNotSolvable('no solution found')

    U = complex(max(Re_Ue), Im_Ue) * e / E_abs
    is_complex = any(isinstance(var, complex) for var in [e, z, s])
    return U if is_complex else U.real


def icc_tr(up_tr: float, us_tr: float, s_tr: float, xcc_tr: float,
           up_ext: float, scc_ext: float = float('inf')) -> float:
    """
    Find the shortcircuit current at a 3-phase transformer's secondary side.

    Parameters
    ----------
    up_tr:   Transformer's primary rated voltage, V.
    us_tr:   Transformer's secondary rated voltage, V.
    s_tr:    Transformer's rated power, VA.
    xcc_tr:  Transformer's shortcircuit impedance, %.
    up_ext:  Voltage of the external source connected to the
             transformer's primary side, V.
    scc_ext: Shortcircuit power of the external source connected
             to the transformer's primary side, VA.
             The default is 'inf' (infinite power source)

    Returns
    -------
    3-phase shortcircuit current at the transformer's secondary side, A

    Examples
    --------
    >>> print(f'Icc = {icc_tr(6900, 400, 850e3, 5, 6900, 200e6)/1000:.2f} kA')
    Icc = 22.62 kA
    >>> print(f'Icc = {icc_tr(6900, 400, 850e3, 5, 6900)/1000:.2f} kA')
    Icc = 24.54 kA
    """
    icc_tr_pu = up_ext / up_tr / (s_tr / scc_ext * (up_ext / up_tr) ** 2 + xcc_tr / 100)
    return icc_tr_pu * s_tr / (SQRT3 * us_tr)


def curve_51(I: ArrayLike | float, curve: IECcurve | IEEEcurve, I_p: float,
             T: float = 1.0, D: float = 0.0) -> ArrayLike | float:
    """
    Compute the actuation time t of an IEC or IEEE 51 curve.
    According to IEC 60255-151 and IEEE C37.112, for a given current I:
        t = T * (A / ((I/I_p)**P - 1) + B) + D,
    where A, P and B are set by the 'curve' parameter.

    Parameters
    ----------
    I:     Current, A.
    curve: The 51 curve to use. Allowed values:
              IECcurve.Inverse
              IECcurve.Very_Inverse
              IECcurve.Extremely_Inverse
              IEEEcurve.Moderately_Inverse
              IEEEcurve.Very_Inverse
              IEEEcurve.Extremely_Inverse
    I_p:   Pick up current, A.
    T:     Multiplicative parameter. The default is 1.
    D:     Additive parameter, s. The default is 0 s.

    Returns
    -------
    The actuation time of the 51 curve, s.
    Exception EEInvalidArguments may be raised.

    Example
    -------
    >>> print(f'{curve_51(200, IECcurve.Very_Inverse, 100):.4f}')
    13.5000
    >>> print(f'{curve_51(200, IEEEcurve.Very_Inverse, 100):.4f}')
    7.0277
    """
    match curve:
        case IECcurve.Inverse:
            A = 0.14
            P = 0.02
            B = 0.0
        case IECcurve.Very_Inverse:
            A = 13.5
            P = 1.0
            B = 0.0
        case IECcurve.Extremely_Inverse:
            A = 80.0
            P = 2.0
            B = 0.0
        case IEEEcurve.Moderately_Inverse:
            A = 0.0515
            P = 0.02
            B = 0.114
        case IEEEcurve.Very_Inverse:
            A = 19.61
            P = 2.0
            B = 0.491
        case IEEEcurve.Extremely_Inverse:
            A = 28.2
            P = 2.0
            B = 0.1217
        case _:
            raise EEInvalidArguments({'curve': curve})

    return T * (A / ((np.asarray(I) / I_p) ** P - 1) + B) + D


@dataclass(frozen=True)
class Element:
    """Base class for the elements to be used with class Network."""

    @property
    def is_complex(self) -> bool:
        return any(isinstance(value, complex) for value in self.__dict__.values())


@dataclass(frozen=True)
class Impedance(Element):
    """Create an impedance to be used with class Network."""

    Z: complex  # Impedance (must be different from 0)


@dataclass(frozen=True)
class Admittance(Element):
    """Create an admittance to be used with class Network."""

    Y: complex  # Admittance (must be different from 0)


@dataclass(frozen=True, repr=False)
class VoltageSource(Element):
    """Create a voltage source to be used with class Network."""

    E: complex                # Electromotive force in an oriented graph
    Z: complex | None = None  # Series impedance. None implies ideal source (Z=0)

    def __repr__(self):
        s = f'{type(self).__name__}(E={self.E!r}'
        if self.Z is not None:
            s = s + f', Z={self.Z!r}'
        return s + ')'


@dataclass(frozen=True, repr=False)
class CurrentSource(Element):
    """Create a current source to be used with class Network."""

    J: complex                # Intensivity in an oriented graph
    Y: complex | None = None  # Parallel admittance. None implies ideal source (Y=0)

    def __repr__(self):
        s = f'{type(self).__name__}(J={self.J!r}'
        if self.Y is not None:
            s = s + f', Y={self.Y!r}'
        return s + ')'


@dataclass(frozen=True)
class ShortCircuit(Element):
    """Create a short circuit to be used with class Network."""


@dataclass(frozen=True)
class MutualCoupling(Element):
    """Create a mutual coupling between two branches to be used with class Network."""

    X: complex  # Inductive reactance of the coupling (must be a pure imaginary number)

    def __post_init__(self):
        if self.X.real != 0:
            raise EEInvalidArguments({'X': self.X})

class Network:
    """
    Solve an electrical networks using the nodes method.

    >>> net = Network('Test network')
    >>> net.add(VoltageSource(E=200, Z=10), from_to=(0, 1), branch=1)
    >>> net.add(VoltageSource(E=-50, Z=20j), from_to=(1, 2), branch=2)
    >>> net.add(Impedance(Z=5j), from_to=(2, 0), branch=3)
    >>> net.add(Impedance(Z=20), from_to=(1, 2), branch=4)
    >>> net.add(CurrentSource(J=4, Y=0.1), from_to=(0, 2), branch=5)
    >>> net.add(MutualCoupling(X=5j), coupled_branches=(2, 3))
    >>> net.remove(branch=2)
    >>> net.remove(coupled_branches=(3, 2))
    >>> net.add(VoltageSource(E=-50, Z=20j), from_to=(1, 2), branch=2)
    >>> net.add(MutualCoupling(X=5j), coupled_branches=(2, 3))
    >>> net.solve()
    >>> net.num_branches
    5
    >>> net.num_nodes
    3
    >>> net.is_solved
    True
    >>> print(net)
    Element                      Connection
    ---------------------------  -----------------------
    VoltageSource(E=200, Z=10)   Branch: 1, Nodes: 0 ► 1
    VoltageSource(E=-50, Z=20j)  Branch: 2, Nodes: 1 ► 2
    Impedance(Z=5j)              Branch: 3, Nodes: 2 ► 0
    Impedance(Z=20)              Branch: 4, Nodes: 1 ► 2
    CurrentSource(J=4, Y=0.1)    Branch: 5, Nodes: 0 ► 2
    MutualCoupling(X=5j)         Coupled branches: 2 & 3
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
     Branch   Nodes    Voltage / V     Current / A
    --------  -------  --------------  -------------
       1      0 ► 1    154.45∠-171.5°  5.24∠-25.7°
       2      1 ► 2    119.23∠1.0°     2.68∠-117.5°
       3      2 ► 0    39.40∠31.6°     6.90∠-38.9°
       4      1 ► 2    119.23∠1.0°     5.96∠1.0°
       5      0 ► 2    39.40∠-148.4°   2.16∠-72.7°
    """

    __slots__ = ['_name', '_net_branches', '_net_couplings', '_is_solved', '_is_complex',
                 '_extra_nodes', '_extra_branches', '_list_of_nodes', '_list_of_branches',
                 '_list_of_coupled_branches', '_ZN', '_VN', '_IB']

    def __init__(self, name: str = '') -> None:
        """Create an AC or DC Network"""
        self._name = name
        self._net_branches = {}
        self._net_couplings = {}
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
        s = [['Element', 'Connection']]

        for branch, (element, (from_node, to_node)) in sorted(self._net_branches.items()):
            s.append([element, 'Branch: ' + str(branch) + ', Nodes: ' +
                      str(from_node) + ' ► ' + str(to_node)])
        for (branch_1, branch_2), element in sorted(self._net_couplings.items(),
                                                    key=lambda tup: min(tup[0])):
            s.append([element, 'Coupled branches: ' + str(branch_1) + r' & ' + str(branch_2)])
        return tabulate.tabulate(s, headers='firstrow')

    def __bool__(self) -> bool:
        """Return False for an empty network, and True otherwise"""
        return (len(self._list_of_branches) + len(self._list_of_coupled_branches)) > 0

    @property
    def name(self) -> str:
        """Return the name of the network"""
        return self._name

    @property
    def num_branches(self) -> int:
        """Return the number of branches in the network."""
        return len(set(self._list_of_branches))

    @property
    def num_nodes(self) -> int:
        """Return the number of nodes in the network."""
        return len(set(self._list_of_nodes))

    @property
    def is_solved(self) -> bool:
        """Return True if the network has been solved, and False otherwise."""
        return self._is_solved

    def add(self, element: Element, *,
            from_to: tuple[int, int] | None = None,
            branch: int | None = None,
            coupled_branches: tuple[int, int] | None = None) -> None:
        """
        Add an element to the network.
        'from_to' and 'branch' are mandatory for elements: Impedance, Admittance,
        VoltageSource, CurrentSource and ShortCircuit; 'coupled_branches' is ignored.
        'coupled_branches' is mandatory for element: MutualCoupling; 'from_to'
        and 'branch' are ignored.

        Exception EEInvalidArguments may be raised.
        """
        match element:
            case Impedance() | Admittance() | VoltageSource() | CurrentSource() | ShortCircuit():
                if from_to is None or branch is None:
                    raise EEInvalidArguments({'from_to': from_to, 'branch': branch})
                if len(from_to) != 2 or (m := min(from_to)) < 0 or m == max(from_to):
                    raise EEInvalidArguments({'from_to': from_to})
                if branch < 1:
                    raise EEInvalidArguments({'branch': branch})
                if branch in self._net_branches:
                    self.remove(branch=branch)
                self._net_branches[branch] = (element, from_to)
                self._is_solved = False
                self._is_complex += element.is_complex
                self._list_of_branches.append(branch)
                self._list_of_nodes.extend(from_to)
                match element:
                    case VoltageSource(_, None) | ShortCircuit():  # Ideal voltage source or shortcircuit
                        self._extra_branches += 1
                        self._extra_nodes += 1
                    case CurrentSource(_, None):  # Ideal current source
                        self._extra_branches += 1
            case MutualCoupling():
                if coupled_branches is None:
                    raise EEInvalidArguments({'coupled_branches': coupled_branches})
                if (len(coupled_branches) != 2 or (m := min(coupled_branches)) < 1 or
                        m == max(coupled_branches)):
                    raise EEInvalidArguments({'coupled_branches': coupled_branches})
                if (idx := frozenset(coupled_branches)) in self._net_branches:
                    self.remove(coupled_branches=coupled_branches)
                self._net_couplings[idx] = element
                self._is_solved = False
                self._is_complex += element.is_complex
                self._list_of_coupled_branches.extend(coupled_branches)
            case _:
                raise EEInvalidArguments({'element': element})

    def remove(self, *, branch: int | None = None,
               coupled_branches: tuple[int, int] | None = None) -> None:
        """
        Remove a branch (int) or a coupling between two branches (tuple of int) from the network.

        Exceptions EENetMissingBranch, EENetMissingMutualCoupling and
        EEInvalidArguments may be raised.
        """
        match (branch, coupled_branches):
            case (None, None):
                raise EEInvalidArguments(locals())
            case (branch, None):
                try:
                    element, from_to = self._net_branches.pop(branch)
                    self._is_solved = False
                    self._is_complex -= element.is_complex
                    for node in from_to:
                        self._list_of_nodes.remove(node)
                    self._list_of_branches.remove(branch)
                    match element:
                        case VoltageSource(_, None) | ShortCircuit():  # Ideal voltage source or shortcircuit
                            self._extra_branches -= 1
                            self._extra_nodes -= 1
                        case CurrentSource(_, None):  # Ideal current source
                            self._extra_branches -= 1
                except KeyError as err:
                    raise EENetMissingBranch({branch}) from err
            case (None, coupled_branches):
                try:
                    element = self._net_couplings.pop(frozenset(coupled_branches))
                    self._is_solved = False
                    self._is_complex -= element.is_complex
                    for branch in coupled_branches:
                        self._list_of_coupled_branches.remove(branch)
                except KeyError as err:
                    raise EENetMissingMutualCoupling({coupled_branches}) from err
            case _:
                raise EEInvalidArguments(locals())

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
        if not self:
            raise EENetNotSolvable('empty network')

        actual_branches = set(self._list_of_branches)
        if actual_branches == set():
            raise EENetMissingBranch({})
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
        EB0 = np.zeros(num_branches_plus_extra, dtype=data_type)
        JB0 = np.zeros(num_branches_plus_extra, dtype=data_type)

        def fill_matrices(branch: int, from_node: int, to_node: int, Zb: complex,
                          Eb: complex | None = None, Jb: complex | None = None) -> None:
            if from_node > 0:
                A[from_node - 1, branch - 1] = 1
            if to_node > 0:
                A[to_node - 1, branch - 1] = -1
            ZB[branch - 1, branch - 1] = Zb
            if Eb is not None:
                EB0[branch - 1] = Eb
            if Jb is not None:
                JB0[branch - 1] = Jb

        for branch, (element, (from_node, to_node)) in self._net_branches.items():
            match element:
                case Impedance(Z):
                    fill_matrices(branch, from_node, to_node, Zb=Z)
                case Admittance(Y):
                    fill_matrices(branch, from_node, to_node, Zb=1/Y)
                case VoltageSource(E, None):  # Ideal voltage source
                    fill_matrices(branch, from_node, num_nodes_plus_extra, Zb=1, Eb=E)
                    fill_matrices(num_branches_plus_extra, num_nodes_plus_extra, to_node, Zb=-1)
                    num_branches_plus_extra -= 1
                    num_nodes_plus_extra -= 1
                case VoltageSource(E, Z):
                    fill_matrices(branch, from_node, to_node, Zb=Z, Eb=E)
                case CurrentSource(J, None):  # Ideal current source
                    fill_matrices(branch, from_node, to_node, Zb=1, Jb=J)
                    fill_matrices(num_branches_plus_extra, from_node, to_node, Zb=-1)
                    num_branches_plus_extra -= 1
                case CurrentSource(J, Y):
                        fill_matrices(branch, from_node, to_node, Zb=1/Y, Jb=J)
                case ShortCircuit():
                    fill_matrices(branch, from_node, num_nodes_plus_extra, Zb=1)
                    fill_matrices(num_branches_plus_extra, num_nodes_plus_extra, to_node, Zb=-1)
                    num_branches_plus_extra -= 1
                    num_nodes_plus_extra -= 1
        for (branch_1, branch_2), element in self._net_couplings.items():
            ZB[branch_1 - 1, branch_2 - 1] = element.X
            ZB[branch_2 - 1, branch_1 - 1] = element.X

        # Solve the network
        try:
            YB = linalg.inv(ZB, overwrite_a=True)
            JB = JB0 + YB @ EB0
            self._ZN = linalg.inv(A @ YB @ A.T, overwrite_a=True)
            self._VN = self._ZN @ (-A @ JB)
            self._IB = YB @ (A.T @ self._VN) + JB
            # Adjust result for branches with ideal current sources
            for branch, (element, _) in self._net_branches.items():
                match element:
                    case CurrentSource(J, None):
                        self._IB[branch - 1] = J
            self._is_solved = True
        except linalg.LinAlgError as err:
            raise EENetNotSolvable(err) from err

    def thevenin(self, node_1: int, node_2: int) -> tuple[complex, complex]:
        """
        Return the Thevenin equivalent circuit between two nodes.

        Exceptions EENetNotSolved and EENetMissingNode may be raised.
        """
        if not self._is_solved:
            raise EENetNotSolved
        if (node_1 not in self._list_of_nodes) or (node_2 not in self._list_of_nodes):
            raise EENetMissingNode({node_1, node_2})

        match (node_1, node_2):
            case (node_1, node_2) if node_1 == node_2:
                eth = zth = 0
            case (node_1, 0):
                eth = self._VN[node_1 - 1]
                zth = self._ZN[node_1 - 1, node_1 - 1]
            case (0, node_2):
                eth = -self._VN[node_2 - 1]
                zth = self._ZN[node_2 - 1, node_2 - 1]
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

    def voltage(self, *, branch: int | None = None, node: int | None = None) -> complex:
        """
        Return the voltage across a branch, or the voltage of a node (with respect
        to node 0).

        Exception EENetNotSolved, EENetMissingBranch, EENetMissingNode and
        EEInvalidArguments may be raised.
        """
        if not self._is_solved:
            raise EENetNotSolved

        match (branch, node):
            case (None, None):
                raise EEInvalidArguments(locals())
            case (branch, None):
                if branch not in self._list_of_branches:
                    raise EENetMissingBranch({branch})
                from_node, to_node = self._net_branches[branch][1]
                return self.thevenin(from_node, to_node)[0]
            case (None, node):
                return self.thevenin(node, 0)[0]
            case _:
                raise EEInvalidArguments(locals())

    def current(self, branch: int) -> complex:
        """
        Return the current through a branch.

        Exceptions EENetNotSolved and EENetMissingBranch may be raised.
        """
        if not self._is_solved:
            raise EENetNotSolved
        if branch not in self._list_of_branches:
            raise EENetMissingBranch({branch})

        return self._IB[branch - 1]

    def results(self, polar: bool = False, u_fmt: str = '', u_unit: str = '',
                u_scale: float = 1, i_fmt: str = '', i_unit: str = '',
                i_scale: float = 1) -> str:
        """
        Return a multi-line string formatted as a table, with all the network's
        voltages and currents.

        A multi-line formatted string is returned with the voltages and
        currents of all the network branches.
        A different scale factor can be applied to voltages (u_scale)
        and currents (i_scale).
        Polar display can be chosen (polar = True/False), and different formats
        for voltages (u_fmt) and currents (i_fmt) are possible. Also, different
        units can be added to voltages (u_unit) and currents (u_unit).

        Exception EENetNotSolved may be raised.
        """
        if not self._is_solved:
            raise EENetNotSolved

        if u_unit:
            u_unit = ' / ' + u_unit
        if i_unit:
            i_unit = ' / ' + i_unit
        s = [['Branch', 'Nodes', 'Voltage' + u_unit, 'Current' + i_unit]]
        for branch in range(1, self.num_branches + 1):
            from_node, to_node = self._net_branches[branch][1]
            u = self.voltage(branch=branch) * u_scale
            i = self.current(branch) * i_scale
            if polar:
                u = Complex(u)
                i = Complex(i)
            s.append([branch, str(from_node) + ' ► ' + str(to_node), format(u, u_fmt), format(i, i_fmt)])
        return tabulate.tabulate(s, numalign='center', headers='firstrow')

    def save(self, filename: str) -> bool:
        """Save the contents of the network to a binary file."""
        try:
            with open(filename, 'wb') as fh:
                pickle.dump(self._name, fh, pickle.HIGHEST_PROTOCOL)
                pickle.dump(self._net_branches, fh, pickle.HIGHEST_PROTOCOL)
                pickle.dump(self._net_couplings, fh, pickle.HIGHEST_PROTOCOL)
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
                self._net_branches = pickle.load(fh)
                self._net_couplings = pickle.load(fh)
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


@dataclass(frozen=True)
class Motor3phRun:
    """Class used to store the result of method 'run' from class Motor3ph."""

    Z: np.ndarray    # Total motor impedance (Line-Neutral), Ω
    PF: np.ndarray   # Motor Power factor (between 0 and 1)
    I: np.ndarray    # Current drawn by the motor, A
    U: np.ndarray    # Voltage at the motor terminals (Line-Neutral), V
    S: np.ndarray    # Apparent power drawn by the motor, VA
    T_m: np.ndarray  # Mechanical torque delivered by the motor, N·m
    P_m: np.ndarray  # Mechanical power delivered by the motor, W
    Eff: np.ndarray  # Motor Efficiency (between 0 and 1)


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
    >>> print(motor)
    Motor parameter             Value
    --------------------------  -------------
    Number of poles, p          4
    Frequency, f                50 Hz
    Stator resistance, R_1      0.095 Ω/phase
    Stator reactance, X_1       0.68 Ω/phase
    Rotor resistance, R_2       0.3 Ω/phase
    Rotor reactance, X_2        0.672 Ω/phase
    Core-loss resistance, R_Fe  620.0 Ω/phase
    Magnetizing reactance, X_m  18.7 Ω/phase
    """

    __slots__ = ['_p', '_f', '_z_1', '_z_2', '_z_m', '_z_0']

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
        self._z_2 = complex(r_2, x_2)
        self._z_m = complex(r_fe, x_m)
        self._z_0 = z_parallel([r_fe, 1j*x_m])

    def __str__(self) -> str:
        """Return a multiline string representing the class"""
        s = [['Motor parameter', 'Value'],
             ['Number of poles, p', str(self._p)],
             ['Frequency, f', str(self._f) + ' Hz'],
             ['Stator resistance, R_1', str(self._z_1.real) + ' Ω/phase'],
             ['Stator reactance, X_1', str(self._z_1.imag) + ' Ω/phase'],
             ['Rotor resistance, R_2', str(self._z_2.real) + ' Ω/phase'],
             ['Rotor reactance, X_2', str(self._z_2.imag) + ' Ω/phase'],
             ['Core-loss resistance, R_Fe', str(self._z_m.real) + ' Ω/phase'],
             ['Magnetizing reactance, X_m', str(self._z_m.imag) + ' Ω/phase']]

        return tabulate.tabulate(s, headers='firstrow')

    @property
    def ω_m_sync(self) -> float:
        """Return the motor synchronous mechanical speed in rad/s."""
        return 4 * math.pi * self._f / self._p

    @property
    def n_m_sync(self) -> float:
        """Return the motor synchronous mechanical speed in r/min."""
        return 120 * self._f / self._p

    def convert(self, vel: ArrayLike | float, from_to: Speed) -> ArrayLike | float:
        """
        Convert vel from/to: slip (s), mechanical speed in rad/s (ω_m),
        or mechanical speed in r/min (n_m)
        """
        match from_to:
            case Speed.slip_to_rpm:
                return (1 - vel) * self.n_m_sync
            case Speed.slip_to_rad_per_s:
                return (1 - vel) * self.ω_m_sync
            case Speed.rpm_to_slip:
                return 1 - vel / self.n_m_sync
            case Speed.rpm_to_rad_per_s:
                return vel * np.pi / 30
            case Speed.rad_per_s_to_slip:
                return 1 - vel / self.ω_m_sync
            case Speed.rad_per_s_to_rpm:
                return vel * 30 / np.pi
            case _:
                raise EEInvalidArguments({'from_to': from_to})

    def s_T_m_max(self, z_ext: complex = 0) -> float:
        """
        Return the slip at the maximum motor mechanical torque.

        Parameters
        ----------
        z_ext: Impedance of the external power system feeding the motor, Ω/phase.
               The default is 0.
        """
        z_th = z_parallel([z_ext + self._z_1, self._z_0])
        return self._z_2.real / abs(z_th + 1j * self._z_2.imag)

    def s_match(self, T_load: Callable[[float], float], e_ext: complex,
                z_ext: complex = 0, s_guess: float = 0.04) -> float:
        """
        Return the slip at which the motor torque matches the load torque.

        Parameters
        ----------
        T_load:  A function that returns the load torque for a given slip, N·m
        e_ext:   Voltage (phase-neutral) of the external power system feeding
                 the motor, V
        z_ext:   Impedance of the external power system feeding the motor, Ω/phase.
                 The default is 0.
        s_guess: Matching slip initial guess. The default is 0.04.
        """
        sol = optimize.root_scalar(lambda s: self.run(s, e_ext, z_ext).T_m - T_load(s), x0=s_guess)
        return sol.root

    def run(self, s: ArrayLike | float, e_ext: complex, z_ext: complex = 0) -> Motor3phRun:
        """
        Return several motor quantities when running at a given slip,
        in a Motor3phRun class.

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
        z_th = z_parallel([z_ext + self._z_1, self._z_0])
        e_th = e_ext * self._z_0 / (self._z_0 + self._z_1 + z_ext)
        r_2_s = self._z_2.real / np.where(np.isclose(s, 0.0, atol=1e-12), 1e-12, s)
        z_2_s = r_2_s + 1j * self._z_2.imag
        Z = self._z_1 + self._z_0 * z_2_s / (self._z_0 + z_2_s)
        PF = np.cos(np.angle(Z))
        I = e_ext / (z_ext + Z)
        U = np.where(z_ext == 0, e_ext, Z * I)
        S = 3 * U * I.conj()
        T_m = 3 * r_2_s * np.abs(e_th / (z_th + z_2_s)) ** 2 / self.ω_m_sync
        P_m = T_m * self.convert(s, Speed.slip_to_rad_per_s)
        Eff = P_m / np.real(S)
        return Motor3phRun(Z, PF, I, U, S, T_m, P_m, Eff)

    def t_start(self, T_load: Callable[[float], float], J: float,
                 e_ext: complex, z_ext: complex = 0) -> float:
        """
        Return the starting time of the motor.

        Parameters
        ----------
        T_load: A function that returns the load torque for a given slip, N·m
        J:      Total moment of inertia of the motor plus the load, kg·m²
        e_ext:  Voltage (line-neutral) of the external power system feeding
                the motor, V
        z_ext:  Impedance of the external power system feeding the motor, Ω/phase.
                The default is 0.

        Returns
        -------
        The starting time of the motor, s.
        """
        n_p1, n_p2, factor = 200, 800, 1.25
        s_match = self.s_match(T_load, e_ext, z_ext)
        # The function to integrate is very smooth up to slightly before s_match
        s_lin = min(1.0, factor * s_match)
        s = np.concatenate([np.linspace(start=1, stop=s_lin, num=n_p1, endpoint=False),
                            np.linspace(start=s_lin, stop=s_match, num=n_p2, endpoint=False)])
        y_fun = 1 / (self.run(s, e_ext, z_ext).T_m - T_load(s))
        return -J * self.ω_m_sync * integrate.simpson(y=y_fun, x=s)

    def s_start_up(self, T_load: Callable[[float], float], J: float,
                 e_ext: complex, z_ext: complex = 0) -> interpolate.interp1d:
        """
        Return an interpolating function of the time evolution of
        the slip during the starting up process.

        Parameters
        ----------
        T_load: A function that returns the load torque for a given slip, N·m
        J:      Total moment of inertia of the motor plus the load, kg·m²
        e_ext:  Voltage (line-neutral) of the external power system feeding
                the motor, V
        z_ext:  Impedance of the external power system feeding the motor, Ω/phase.
                The default is 0.

        Returns
        -------
        An interpolating function (scipy.interpolate.interp1d) that returns the slip
        at a given time in seconds, during the starting up process.
        If time < 0, slip=1 (motor stopped) is returned.
        If time > starting time, the slip at the starting time is returned.
        """
        n_p1, n_p2, factor = 200, 800, 1.25
        s_match = self.s_match(T_load, e_ext, z_ext)
        # The function to integrate is very smooth up to slightly before s_match
        s_lin = min(1.0, factor * s_match)
        s = np.concatenate([np.linspace(start=1, stop=s_lin, num=n_p1, endpoint=False),
                            np.linspace(start=s_lin, stop=s_match, num=n_p2, endpoint=False)])
        y_fun = 1 / (self.run(s, e_ext, z_ext).T_m - T_load(s))
        # The variable x in cumulative_simpson must be strictly increasing,
        # so we need to change the sign of s, and then the sign of the integral
        t = J * self.ω_m_sync * integrate.cumulative_simpson(y=y_fun, x=-s, initial=0)
        return interpolate.interp1d(t, s, kind='cubic', copy=False,
                                    assume_sorted=True, bounds_error=False,
                                    fill_value=(s[0], s[-1]))


if __name__ == '__main__':
    import doctest

    doctest.testmod()
