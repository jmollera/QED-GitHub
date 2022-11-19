__author__ = 'Josep Mollera Barriga'
__version__ = '14.1'
__all__ = ['z_cable', 'r_cable', 'voltage_drop', 'kcmil_to_mm2', 'mm2_to_kcmil',
           'AWG_to_mm2', 'mm2_to_AWG', 'z_series', 'z_parallel',
           'millman', 'apparent_power', 'D_to_Y', 'Y_to_D',
           'triangle_to_phasors', 'LN_to_LL', 'LL_to_LG', 'LL_to_LN',
           'ABC_to_A012', 'A012_to_ABC', 'AN12_to_AB12', 'AB12_to_AN12',
           'ezs_u', 'CEI_51_curve', 'IEEE_51_curve',
           'Impedance', 'Admitance', 'VoltageSource', 'VoltageSourceIdeal',
           'CurrentSource', 'CurrentSourceIdeal', 'ShortCircuit', 'MutualCoupling',
           'NetworkException', 'NetworkNotSolved', 'NetworkMissingBranch', 
           'NetworkMissingNode', 'NetworkMissingMutualCoupling',  
           'NetworkInvalidElement', 'Network']

from dataclasses import dataclass, field
import numpy as np
from scipy import linalg
import pickle
from qed.utils import ComplexD

"""
Two dictionaries of copper cable impedances.

    Index: section in mm²

    Value: impedance at 90 °C in Ω/km
"""
Z_CABLE_A = {4   : 5.99+0.098j,    6   : 3.96+0.123j,
             10  : 2.34+0.117j,    16  : 1.47+0.112j,
             25  : 0.936+0.107j,   35  : 0.661+0.123j,
             50  : 0.5+0.11j,      70  : 0.34+0.114j,
             120 : 0.192+0.108j,   150 : 0.156+0.104j,
             185 : 0.125+0.103j,   240 : 0.095+0.102j,
             300 : 0.078+0.094j}

Z_CABLE_V = {4   : 5.993+0.168j,   6   : 3.94+0.159j,
             10  : 2.33+0.151j,    16  : 1.46+0.142j,
             25  : 0.9361+0.0791j, 35  : 0.6749+0.0765j,
             50  : 0.4935+0.135j,  70  : 0.3417+0.128j,
             95  : 0.249+0.130j,   120 : 0.1980+0.0825j}

"""
Dictionary of resistivity and coefficients of resistivity
variation with temperature, for several materials.

    Index: chemical symbol of the material

    Value: resistivity at 20 °C in Ω·mm²/m, and coefficient of
    resistivity variation with temperature at 20 °C in 1/°C
"""
R_CABLE = {'Al' : (0.02825, 0.00391), 'Cu' : (0.01723, 0.00393),
           'Ag' : (0.01645, 0.00380), 'Au' : (0.02440, 0.00340)}


def z_cable(section, length, kind='A'):
    """
    Finds the impedance of a cable.

    Parameters
    ----------
    section : float
        Cable section in mm².
    length : float
        Cable length in m.
    kind : str, optional
        'A' or 'a' → use Z_CABLE_A,
        'V' or 'v' → use Z_CABLE_V,
        The default is 'A'.

    Returns
    -------
    float
        The cable impedance at 90 °C in Ω, or 0 if kind is not valid
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
                return Z_CABLE_A[section] * length/1000
            case 'V':
                return Z_CABLE_V[section] * length/1000
            case _:
                raise KeyError
    except KeyError:
        return 0


def r_cable(section, length, temp=20, kind='Cu'):
    """
    Finds the resistance of a cable.

    Parameters
    ----------
    section : float
        Cable section in mm².
    length : float
        Cable length in m.
    temp : float, optional
        Temperature in °C.
        The default is 20 °C.
    kind : str, optional
        The chemical symbol of the cable material.
        The default is 'Cu'.

    Returns
    -------
    float
        The cable resistance at the requested temperature in Ω,
        or 0 if kind is not valid.

    Examples
    --------
    >>> print(f'{r_cable(4, 100):.4f}')
    0.4307
    >>> print(f'{r_cable(4, 100, temp=90):.4f}')
    0.5492
    >>> print(f'{r_cable(4, 100, kind="Al"):.4f}')
    0.7063
    >>> print(f'{r_cable(4, 100, temp=90, kind="Al"):.4f}')
    0.8996
    >>> r_cable(4, 100, kind='Fe')
    0
    """
    try:
        ρ_20, α_20 = R_CABLE[kind]
        return ρ_20*(1 + α_20*(temp - 20)) * length/section
    except KeyError:
        return 0


def voltage_drop(voltage, current, impedance, cos_φ=1):
    """
    Finds the voltage drop through a cable.

    Parameters
    ----------
    voltage : float
        Voltage at cable origin (line-neutral voltage
        for a three-phase system) in V.
    current : float
        Current through the cable in A.
    impedance : float
        Cable impedance in Ω.
    cos_φ : float, optional
        Load power factor. The default is 1.

    Returns
    -------
    float
        Voltage drop through the cable (line-neutral voltage
        for a three-phase system) in V.

    Examples
    --------
    Three-phase: U=380 V (LL). I=200 A. Cable: 240 mm², 400 m. cos φ=0.87
    >>> dv = voltage_drop(380/np.sqrt(3), 200, z_cable(240, 400), 0.87)
    >>> print(f'{dv:.4f}')
    10.6609

    One-phase: U=220 V. I=55 A. Cable: 35 mm², 200 m. cos φ=0.85
    >>> dv = voltage_drop(220, 55, z_cable(35, 2*200), 0.85)
    >>> print(f'{dv:.4f}')
    13.8515

    Direct Current: U=125 V. I=25 A. Cable: 16 mm², 50 m
    >>> dv = voltage_drop(125, 25, z_cable(16, 2*50).real)
    >>> print(f'{dv:.4f}')
    3.6750
    """
    R = np.real(impedance)
    X = np.imag(impedance)
    U = np.abs(voltage)
    I = np.abs(current)
    sin_φ = np.sqrt(1 - cos_φ**2)
    return (R*cos_φ + X*sin_φ)*I + U - np.sqrt(U**2 - ((X*cos_φ - R*sin_φ)*I)**2)


def kcmil_to_mm2(kcmil):
    """
    Converts a circular section from kcmil to mm².

    Parameters
    ----------
    kcmil : float
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
    return kcmil*0.506707479098


def mm2_to_kcmil(mm2):
    """
    Converts a circular section from mm² to kcmil.

    Parameters
    ----------
    mm2 : float
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
    return mm2/0.506707479098


def AWG_to_mm2(AWG):
    """
     Converts a circular section from AWG to mm².

    Parameters
    ----------
    AWG : int
        Circular section in AWG; 0, -1, -2, and -3 must be
        used for AWG 1/0, 2/0, 3/0, and 4/0.

    Returns
    -------
    float
        Equivalent section in mm².

    Examples
    --------
    >>> print(f'{AWG_to_mm2(12):.4f}')
    3.3088
    """
    return 53.4751207321/(92**(AWG/19.5))


def mm2_to_AWG(mm2):
    """
    Converts a circular section from mm² to AWG

    Parameters
    ----------
    mm2 : float
        Circular section in mm².

    Returns
    -------
    int
        Equivalent section rounded to the nearest AWG;
        results 0, -1, -2, and -3 are equivalent to
        AWG 1/0, 2/0, 3/0, and 4/0.

    Examples
    --------
    >>> mm2_to_AWG(2)
    14.0
    """
    return np.rint(4.31245284200*np.log(53.4751207321/mm2))


def z_series(z):
    """
    Computes the series impedance of a list of impedances

    Parameters
    ----------
    z : float/complex array-like
        Impedances connected in series.

    Returns
    -------
    float/complex
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
    Computes the parallel impedance of a list of impedances.

    Parameters
    ----------
    z : float/complex array-like
        Impedances connected in parallel.

    Returns
    -------
    float/complex
        The parallel impedance.

    Examples
    --------
    >>> z_parallel([3, 3, 3])
    1.0
    >>> z_parallel([3+3j, 3+3j, 3+3j])
    (1+1j)
    """
    p = np.prod(z)
    return p/np.sum(p/z)


def millman(u, z):
    """
    Applies the Millman theorem to a list of voltages
    and a list of star connected impedances.

    Parameters
    ----------
    u : float/complex array-like
        Voltages with respect to a reference point
        of a group of star connected impedances.
    z : float/complex array-like
        Impedances of the star connected impedances.

    Returns
    -------
    float/complex
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
    return np.sum(np.array(u)/z)*z_parallel(z)


def apparent_power(u, i):
    """
    Computes the apparent power for a list of voltages
    and a list of currents.

    Parameters
    ----------
    u : float/complex array-like
        Voltages with respect to a reference point
        of a set of branches feeding a load.
    i : float/complex array-like
        Currents through the branches towards the load.

    Returns
    -------
    float/complex
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
    Converts three delta connected impedances into three
    equivalent star connected impedances.

    Parameters
    ----------
    z_ab : float/complex
        Impedance between phases A and B.
    z_bc : float/complex
        Impedance between phases B and C.
    z_ca : float/complex
        Impedance between phases C and A.

    Returns
    -------
    float/complex
        Impedance between phase A and neutral point.
    float/complex
        Impedance between phase B and neutral point.
    float/complex
        Impedance between phase C and neutral point.

    Examples
    --------
    >>> D_to_Y(3, 3, 3)
    (1.0, 1.0, 1.0)
    >>> D_to_Y(10, -10j, -10j)
    ((4-2j), (4-2j), (-2-4j))
    """
    z = (z_ab*z_bc*z_ca)/(z_ab + z_bc + z_ca)
    return z/z_bc, z/z_ca, z/z_ab


def Y_to_D(z_an, z_bn, z_cn):
    """
    Converts three star connected impedances into three
    equivalent delta connected impedances.

    Parameters
    ----------
    z_an : float/complex
        Impedance between phase A and neutral point.
    z_bn : float/complex
        Impedance between phase B and neutral point.
    z_cn : float/complex
        Impedance between phase C and neutral point.

    Returns
    -------
    float/complex
        Impedance between phases A and B.
    float/complex
        Impedance between phases B and C.
    float/complex
        Impedance between phases C and A.

    Examples
    --------
    >>> Y_to_D(1, 1, 1)
    (3.0, 3.0, 3.0)
    >>> Y_to_D(4-2j, 4-2j, -2-4j)
    ((10-0j), -10j, -10j)
    """
    z = z_an*(z_bn + z_cn) + z_bn*z_cn
    return z/z_cn, z/z_an, z/z_bn


def triangle_to_phasors(u_ab, u_bc, u_ca):
    """
    Converts the three sides of a voltage triangle into
    three voltage phasors.

    Parameters
    ----------
    u_ab : float
        Side AB of the voltage triangle.
    u_bc : float
        Side BC of the voltage triangle.
    u_ca : float
        Side CA of the voltage triangle.

    Returns
    -------
    complex
        Voltage AB phasor. This voltage is taken as the reference
        phasor (phasor angle equal to zero).
    complex
        Voltage BC phasor.
    complex
        Voltage CA phasor.

    Examples
    --------
    >>> print('  '.join(format(u, '.4f') for u in triangle_to_phasors(400, 500, 300)))
    400.0000+0.0000j  -400.0000-300.0000j  0.0000+300.0000j
    >>> print('  '.join(format(u, '.4f') for u in triangle_to_phasors(2760, 1840, 2300)))
    2760.0000+0.0000j  -1035.0000-1521.3070j  -1725.0000+1521.3070j
    """
    φ_a = np.arccos((u_ab**2 + u_ca**2 - u_bc**2)/(2*u_ab*u_ca))
    φ_b = np.arccos((u_bc**2 + u_ab**2 - u_ca**2)/(2*u_bc*u_ab))
    return (complex(u_ab, 0),
            u_bc * np.exp(1j*(np.pi + φ_b)),
            u_ca * np.exp(1j*(np.pi - φ_a)))


def LN_to_LL(u_an, u_bn, u_cn):
    """
    Finds the three line-line voltages corresponding
    to three line-neutral voltages.

    Parameters
    ----------
    u_an : complex
        Voltage AN phasor.
    u_bn : complex
        Voltage BN phasor.
    u_cn : complex
        Voltage CN phasor.

    Returns
    -------
    complex
        Voltage AB phasor.
    complex
        Voltage BC phasor.
    complex
        Voltage CA phasor.

    Examples
    --------
    >>> LN_to_LL(120j, 120, -120-120j)
    ((-120+120j), (240+120j), (-120-240j))
    """
    return u_an-u_bn, u_bn-u_cn, u_cn-u_an


def LL_to_LG(u_ab, u_bc, u_ca):
    """
    Finds the three line-G voltages corresponding
    to three line-line voltages, where G is the
    barycentre of the line-line voltages

    Parameters
    ----------
    u_ab : complex
        Voltage AB phasor.
    u_bc : complex
        Voltage BC phasor.
    u_ca : complex
        Voltage CA phasor.

    Returns
    -------
    complex
        Voltage AG phasor.
    complex
        Voltage BG phasor.
    complex
        Voltage CG phasor.

    Examples
    --------
    >>> LL_to_LG(-120+120j, 240+120j, -120-240j)
    (120j, (120+0j), (-120-120j))
    """
    return (u_ab - u_ca)/3, (u_bc - u_ab)/3, (u_ca - u_bc)/3


def LL_to_LN(u_ab, u_bc, u_ca, z_an, z_bn, z_cn):
    """
    Finds the three line-neutral voltages corresponding to
    three line-line voltages, with three star connected impedances

    Parameters
    ----------
    u_ab : complex
        Voltage AB phasor.
    u_bc : complex
        Voltage BC phasor.
    u_ca : complex
        Voltage CA phasor.
    z_an : complex
        Impedance AN phasor.
    z_bn : complex
        Impedance BN phasor.
    z_cn : complex
        Impedance CN phasor.

    Returns
    -------
    complex
        Voltage AN phasor.
    complex
        Voltage BN phasor.
    complex
        Voltage CN phasor.

    Examples
    ---------
    >>> LL_to_LN(-120+120j, 240+120j, -120-240j, 10, 10, 10)
    (120j, (120+0j), (-120-120j))
    """
    z_p = z_parallel([z_an, z_bn, z_cn])
    return ((u_ab/z_bn - u_ca/z_cn) * z_p,
            (u_bc/z_cn - u_ab/z_an) * z_p,
            (u_ca/z_an - u_bc/z_bn) * z_p)


def ABC_to_A012(x_a, x_b, x_c):
    """
    Finds the three sequence phasors corresponding to
    a system of three unbalanced phasors.

    Parameters
    ----------
    x_a : complex
        Phasor A in a system of three unbalanced phasors.
    x_b : complex
        Phasor B in a system of three unbalanced phasors.
    x_c : complex
        Phasor C in a system of three unbalanced phasors.

    Returns
    -------
    complex
        Homopolar sequence phasor.
    complex
        Direct sequence phasor.
    complex
        Inverse sequence phasor.

    Examples
    --------
    >>> print('  '.join(format(u, '.4f') for u in ABC_to_A012(2760, -1035-1521j, -1740+1535j)))
    -5.0000+4.6667j  2264.6912+201.1826j  500.3088-205.8493j
    """
    a = complex(-1/2, np.sqrt(3)/2)
    a2 = a.conjugate()
    return ((x_a + x_b + x_c)/3,
            (x_a + a*x_b + a2*x_c)/3,
            (x_a + a2*x_b + a*x_c)/3)


def A012_to_ABC(x_a0, x_a1, x_a2):
    """
    Finds a system of three unbalanced phasors
    corresponding to three sequence phasors.

    Parameters
    ----------
    x_a0 : complex
        Homopolar sequence phasor.
    x_a1 : complex
        Direct sequence phasor.
    x_a2 : complex
        Inverse sequence phasor.

    Returns
    -------
    complex
        Phasor A in a system of three unbalanced phasors.
    complex
        Phasor B in a system of three unbalanced phasors.
    complex
        Phasor C in a system of three unbalanced phasors.

    Examples
    --------
    >>> print('  '.join(format(u, '.4f') for u in A012_to_ABC(-5+4j, 2264+201j, 500-205j)))
    2759.0000+0.0000j  -1035.3937-1521.6688j  -1738.6063+1533.6688j
    """
    a = complex(-1/2, np.sqrt(3)/2)
    a2 = a.conjugate()
    return (x_a0 + x_a1 + x_a2,
            x_a0 + a2*x_a1 + a*x_a2,
            x_a0 + a*x_a1 + a2*x_a2)


def AN12_to_AB12(x_an_1, x_an_2):
    """
    Finds the Line-Line direct and inverse sequence phasors
    corresponding to the Line-Neutral direct and inverse
    sequence phasors.

    Parameters
    ----------
    x_an_1 : complex
        Line-Neutral direct sequence phasor.
    x_an_2 : complex
        Line-Neutral inverse sequence phasor.

    Returns
    -------
    complex
        Line-Line direct sequence phasor.
    complex
        Line-Line inverse sequence phasor.

    Examples
    --------
    >>> print('  '.join(format(u, '.4f') for u in AN12_to_AB12(1187-552j, 308+44j)))
    2258.5460+199.9722j  500.1051-200.7358j
    """
    return (x_an_1*np.sqrt(3)*np.exp(1j*np.pi/6),
            x_an_2*np.sqrt(3)*np.exp(-1j*np.pi/6))


def AB12_to_AN12(x_ab_1, x_ab_2):
    """
    Finds the line-Neutral direct and inverse sequence phasors
    corresponding to the Line-Line direct and inverse
    sequence phasors

    Parameters
    ----------
    x_ab_1 : complex
        Line-Line direct sequence phasor.
    x_ab_2 : complex
        Line-Line inverse sequence phasor.

    Returns
    -------
    complex
        Line-Neutral direct sequence phasor.
    complex
        Line-Neutral inverse sequence phasor

    Examples
    --------
    >>> print('  '.join(format(u, '.4f') for u in AB12_to_AN12(2258+200j, 500-200j)))
    1186.7350-551.8285j  307.7350+44.3376j
    """
    return (x_ab_1/np.sqrt(3)/np.exp(1j*np.pi/6),
            x_ab_2/np.sqrt(3)/np.exp(-1j*np.pi/6))


def ezs_u(e, z, s):
    """
    Finds the voltage at a load, given the source
    voltage, the impedance between source and load, and
    the power drawn by the load.

    Parameters
    ----------
    e : float/complex
        Source voltage.
    z : float/complex
        Impedance between source and load.
    s : float/complex
        Power drawn by the load.

    Returns
    -------
    float/complex
        Voltage at the load.

    Examples
    --------
      >>> print(f'{ezs_u(0.4+0.3j, 0.1j, 0.6+0.45j):.4f}')
      0.3165+0.0874j
      >>> print(f'{ezs_u(125, 0.1, 100):.4f}')
      124.9199
    """
    ZcjS = z.conjugate()*s
    ImUe = ZcjS.imag/abs(e)
    Eabs = abs(e)
    ReUe = np.polynomial.Polynomial([ZcjS.real + ImUe**2, -Eabs, 1]).roots()
    if isinstance(ReUe[0], complex):
        print('No solution found')
        return None
    else:
        U = complex(max(ReUe), ImUe)*e/Eabs
        if any(isinstance(var, complex) for var in [e, z, s]):
            return U
        else:
            return U.real


def CEI_51_curve(I, curve_type, I_p, T=1, L=0, C=0):
    """
    Returns the actuation time 't' of a CEI 51 curve,
    for a given current 'I', according to:
        t = T*(k/(I/I_p)**a + L) + C
    where 'k' and 'a' are set by the 'curve_type' parameter.

    Parameters
    ----------
    I : float
        Current, in A.
    curve_type : str
        The CEI 51 curve to use. Allowed values are (uppercase or lowercase):
            'S' → Standard inverse
            'V' → Very inverse
            'E' → Extremely inverse
    I_p : float
        Pick up current of the CEI 51 curve, in A.
    T : float, optional.
        Parameter of the CEI 51 curve. The default is 1.
    L : float, optional.
        Parameter of the CEI 51 curve. The default is 0.
    C : float, optional.
        Parameter of the CEI 51 curve. The default is 0.

    Returns
    -------
    float
        The actuation time of the CEI 51 curve, in s.

    Examples
    --------
      >>> print(f'{CEI_51_curve(200, "V", 100):.4f}')
      13.5000
      >>> CEI_51_curve(200, 'Inverse', 100)
      Unknown curve "Inverse"
    """
    try:
        match curve_type.upper():
            case 'S':
                k = 0.14
                a = 0.02
            case 'V':
                k = 13.5
                a = 1.0
            case 'E':
                k = 80.0
                a = 2.0
            case _:
                raise AttributeError
        
        return T*(k/((I/I_p)**a - 1) + L) + C
    
    except AttributeError:
        print(f'Unknown curve "{curve_type}"')


def IEEE_51_curve(I, curve_type, I_p, T=1, C=0):
    """
    Returns the actuation time 't' of an IEEE 51 curve,
    for a given current 'I', according to:
        t = T*(k/(I/I_p)**a + L) + C
    where 'k', 'a' and 'L' are set by the 'curve_type' parameter.

    Parameters
    ----------
    I : float
        Current, in A.
    curve_type : str
        The IEEE 51 curve to use. Allowed values are (uppercase or lowercase):
            'M' → Moderately inverse
            'V' → Very inverse
            'E' → Extremely inverse
    I_p : float
        Pick up current of the IEEE 51 curve, in A.
    T : float, optional.
        Parameter of the IEEE 51 curve. The default is 1.
    C : float, optional.
        Parameter of the IEEE 51 curve. The default is 0.

    Returns
    -------
    float
        The actuation time of the IEEE 51 curve, in s.

    Examples
    --------
      >>> print(f'{IEEE_51_curve(200, "V", 100):.4f}')
      7.0277
      >>> IEEE_51_curve(200, 'Inverse', 100)
      Unknown curve "Inverse"
    """
    try:
        match curve_type.upper():
            case 'M':
                k = 0.0515
                a = 0.02
                L = 0.114
            case 'V':
                k = 19.61
                a = 2.0
                L = 0.491
            case 'E':
                k = 28.2
                a = 2.0
                L = 0.1217
            case _:
                raise AttributeError
       
        return T*(k/((I/I_p)**a - 1) + L) + C
    
    except AttributeError:
        print(f'Unknown curve "{curve_type}"')


class NetworkException(Exception):
    """
    Base exception class for all the class Network exceptions
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class NetworkNotSolved(NetworkException):
    """
    Exception raised in class Network when results want to be
    obtained before the network is solved.
    """
    def __init__(self):
        mssg = 'Netwok not solved yet!'
        super().__init__(mssg)
        self.mssg = mssg


class NetworkMissingBranch(NetworkException):
    """
    Exception raised in class Network when trying to access missing branches.
    """
    def __init__(self, missing_branches):
        if len(missing_branches) == 1:
            mssg = f'Branch {missing_branches} is invalid or has no information!'
        else:
            mssg = f'Branches {missing_branches} are invalid or have no information!'
        super().__init__(mssg)
        self.missing_branches = missing_branches
        self.mssg = mssg


class NetworkMissingNode(NetworkException):
    """
    Exception raised in class Network when trying to access missing nodes.
    """
    def __init__(self, missing_nodes):
        if len(missing_nodes) == 1:
            mssg = f'Node {missing_nodes} is invalid or has no connections!'
        else:
            mssg = f'Nodes {missing_nodes} are invalid or have no connections!'
        super().__init__(mssg)
        self.missing_nodes = missing_nodes
        self.mssg = mssg


class NetworkMissingMutualCoupling(NetworkException):
    """
    Exception raised in class Network when trying to access missing mutual couplings.
    """
    def __init__(self, missing_mutual_coupling):
        mssg = f'Mutual coupling between branches {missing_mutual_coupling} does not exist!'
        super().__init__(mssg)
        self.missing_mutual_coupling = missing_mutual_coupling
        self.mssg = mssg


class NetworkInvalidElement(NetworkException):
    """
    Exception raised in class Network when trying to access an invalid element
    (one of the named tuples to be used with class Network).
         
    """
    def __init__(self, invalid_element):
        mssg = f'Element "{invalid_element}" is not valid'
        super().__init__(mssg)
        self.invalid_element = invalid_element
        self.mssg = mssg


@dataclass
class Impedance:
    """ class to define an impedance, to be used with the class Network """
    from_node: int
    to_node: int
    Z: complex
    _is_complex: bool = field(init=False, repr=False)

    def __post_init__(self):
         self._is_complex = isinstance(self.Z, complex)


@dataclass
class Admitance:
    """ class to define an admitance, to be used with the class Network """

    from_node: int
    to_node: int
    Y: complex
    _is_complex: bool = field(init=False, repr=False)

    def __post_init__(self):
         self._is_complex = isinstance(self.Y, complex)
         
         
@dataclass
class VoltageSource:
    """ class to define a voltage source, to be used with the class Network """
    from_node: int
    to_node: int
    E: complex
    Z: complex
    _is_complex: bool = field(init=False, repr=False)

    def __post_init__(self):
         self._is_complex = isinstance(self.E, complex) or isinstance(self.Z, complex)       
         
         
@dataclass
class VoltageSourceIdeal:
    """ class to define an ideal voltage source, to be used with the class Network """
    from_node: int
    to_node: int
    E: complex
    _is_complex: bool = field(init=False, repr=False)

    def __post_init__(self):
         self._is_complex = isinstance(self.E, complex)


@dataclass
class CurrentSource:
    """ class to define a current source, to be used with the class Network """
    from_node: int
    to_node: int
    J: complex
    Y: complex
    _is_complex: bool = field(init=False, repr=False)

    def __post_init__(self):
         self._is_complex = isinstance(self.J, complex) or isinstance(self.Y, complex)       
         
         
@dataclass
class CurrentSourceIdeal:
    """ class to define an idealcurrent source, to be used with the class Network """
    from_node: int
    to_node: int
    J: complex
    _is_complex: bool = field(init=False, repr=False)

    def __post_init__(self):
         self._is_complex = isinstance(self.J, complex)


@dataclass
class ShortCircuit:
    """ class to define a short circuit, to be used with the class Network """
    from_node: int
    to_node: int
    _is_complex: bool = field(init=False, repr=False)

    def __post_init__(self):
         self._is_complex = False
         
         
@dataclass
class MutualCoupling:
    """ class to define a mutual coupling, to be used with the class Network """
    XM: complex
    _is_complex: bool = field(init=False, repr=False)

    def __post_init__(self):
         self._is_complex = True
                  
         
class Network:
    """
    Class to solve electrical networks using the nodes method
    
    >>> net = Network()
    >>> net.add(VoltageSource(from_node=0, to_node=1, E=200, Z=10), branch=1)
    >>> net.add(VoltageSource(from_node=1, to_node=2, E=-50, Z=20j), branch=2)
    >>> net.add(Impedance(from_node=2, to_node=0, Z=5j), branch=3)
    >>> net.add(Impedance(from_node=1, to_node=2, Z=20), branch=4)
    >>> net.add(CurrentSource(from_node=0, to_node=2, J=4, Y=1/10), branch=5)
    >>> net.add(MutualCoupling(XM=5j), branch=2, with_branch=3)
    >>> net.remove(branch=2)
    >>> net.remove(branch=2, with_branch=3)
    >>> net.add(VoltageSource(from_node=1, to_node=2, E=-50, Z=20j), branch=2)
    >>> net.add(MutualCoupling(XM=5j), branch=2, with_branch=3)
    >>> net.solve()
    >>> net.num_branches
    5
    >>> net.num_nodes
    3
    >>> net.is_solved
    True
    >>> print(net)
    Branch  Element
    ------  --------------------------------------------------
     2, 3   MutualCoupling(XM=5j)
      1     VoltageSource(from_node=0, to_node=1, E=200, Z=10)
      2     VoltageSource(from_node=1, to_node=2, E=-50, Z=20j)
      3     Impedance(from_node=2, to_node=0, Z=5j)
      4     Impedance(from_node=1, to_node=2, Z=20)
      5     CurrentSource(from_node=0, to_node=2, J=4, Y=0.1)
    >>> print(f'{net.node_voltage(2):.4f}')
    33.5644+20.6436j
    >>> print(f'{net.branch_voltage(4):.4f}')
    119.2079+2.0792j
    >>> print(f'{net.branch_current(4):.4f}')
    5.9604+0.1040j
    >>> print('  '.join(format(x, '.4f') for x in net.thevenin(1, 2)))
    119.2079+2.0792j  4.2574+2.5743j
    >>> print('  '.join(format(x, '.4f') for x in net.norton(1, 2)))
    20.7200-12.0400j  0.1720-0.1040j
    >>> print(net.results(polar=True, u_fmt='.2f/.1f', u_unit='V', i_fmt='.2f/.1f', i_unit='A'))
    U1 = 154.45∠-171.5° V,  I1 = 5.24∠-25.7° A
    U2 = 119.23∠1.0° V,  I2 = 2.68∠-117.5° A
    U3 = 39.40∠31.6° V,  I3 = 6.90∠-38.9° A
    U4 = 119.23∠1.0° V,  I4 = 5.96∠1.0° A
    U5 = 39.40∠-148.4° V,  I5 = 2.16∠-72.7° A
    """
    def __init__(self):
        self._net = {}
        self._is_solved = False
        self._is_complex = 0
        self._extra_nodes = 0
        self._extra_branches = 0
        self._list_nodes = []
        self._list_branches = []
        self._list_coupled_branches = []
    
    def __str__(self):
        def sort_net(tup):
            k, _ = tup
            if isinstance(k, int):
                return k
            else:
                return -min(k)
        s = ['Branch  Element\n------  ' + '-'*50]
        for branch, element in sorted(self._net.items(), key=sort_net):
            if isinstance(branch, int):
                s.append(f'{branch:^6d}  {element}')
            else:
                b1, b2 = branch
                s.append(f'{b1}, {b2}'.center(6) + f'  {element}')
        return '\n'.join(s)
        
    @property
    def num_branches(self):
        """ Returns the number of branches in the network.
            Exception NetworkNotSolved is raised if the network is not solved.
        """
        if not self._is_solved: raise NetworkNotSolved
        
        return len(set(self._list_branches))
        
    @property
    def num_nodes(self):
        """ Returns the number of nodes in the network,
            Exception NetworkNotSolved is raised if the network is not solved.
        """
        if not self._is_solved: raise NetworkNotSolved
        
        return len(set(self._list_nodes))

    @property
    def is_solved(self):
        """ Returns True if the network has been solved, and False otherwise """
        return self._is_solved

    def add(self, element, branch, with_branch=None):
        """ Adds 'element' (named tuple: Impedance, Admitance, VoltageSource,
            VoltageSourceIdeal, CurrentSource, CurrentSourceIdeal or
            ShortCircuit) in 'branch' to the network, or adds 'element' 
            (named tuple MutualCoupling) between 'branch' and 'with_branch'
            to the network.
            Exception NetworkInvalidElement is raised if 'element' is not 
            one of the above named tuples.
        """
        match element:
            case MutualCoupling():
                branches = tuple(sorted([branch, with_branch]))
                if branches in self._net: self.remove(branch, with_branch)
                self._net[branches] = element
                self._is_solved = False
                self._is_complex += int(element._is_complex)  
                self._list_coupled_branches.extend([branch, with_branch])
            case Impedance() | Admitance() | VoltageSource() | CurrentSource():
                if branch in self._net: self.remove(branch)
                self._net[branch] = element
                self._is_solved = False
                self._is_complex += int(element._is_complex)  
                self._list_branches.append(branch)
                self._list_nodes.extend([element.from_node, element.to_node])
            case VoltageSourceIdeal() | ShortCircuit():
                if branch in self._net: self.remove(branch)
                self._net[branch] = element
                self._is_solved = False
                self._is_complex += int(element._is_complex)  
                self._list_branches.append(branch)
                self._list_nodes.extend([element.from_node, element.to_node])
                self._extra_branches += 1
                self._extra_nodes += 1
            case CurrentSourceIdeal():
                if branch in self._net: self.remove(branch)
                self._net[branch] = element
                self._is_solved = False
                self._is_complex += int(element._is_complex)  
                self._list_branches.append(branch)
                self._list_nodes.extend([element.from_node, element.to_node])
                self._extra_branches += 1
            case _:
                raise NetworkInvalidElement(element)

    def remove(self, branch, with_branch=None):
        """ Removes an element in 'brancn' from the network, or a mutual coupling 
            between 'brancn' and 'with_branch' from the network.
            Exception NetworkMissingBranch is raised if 'branch' does not exist,
            Exception NetworkMissingMutualCoupling if raised if 'branch' or 
            'with_branch' do not exit.
        """
        if with_branch is not None:
            try:
                element = self._net.pop(tuple(sorted([branch, with_branch])))
                self._is_solved = False
                self._is_complex -= int(element._is_complex)
                self._list_coupled_branches.remove(branch)
                self._list_coupled_branches.remove(with_branch) 
            except KeyError:
                raise NetworkMissingMutualCoupling({branch, with_branch})
        else:
            try:
                element = self._net.pop(branch)
                self._is_solved = False
                self._is_complex -= int(element._is_complex)
                self._list_nodes.remove(element.from_node)
                self._list_nodes.remove(element.to_node)
                self._list_branches.remove(branch)
                match element:
                    case VoltageSourceIdeal() | ShortCircuit():
                        self._extra_branches -= 1
                        self._extra_nodes -= 1
                    case CurrentSourceIdeal():
                        self._extra_branches -= 1
            except KeyError:  
                raise NetworkMissingBranch({branch})
                
    def clear(self):
        """ Removes all the network's elements """
        self._net.clear()
        self._is_solved = False
        self._is_complex = 0
        self._extra_nodes = 0
        self._extra_branches = 0
        self._list_nodes.clear()
        self._list_branches.clear()
        self._list_coupled_branches.clear()
        
    def solve(self):
        """Solves the network, finding the voltage accros every branch,
           and the current through every branch.
           Exception NetworkMissingBranch is raised if not all branches
           1, 2, 3, ... max(self._branches) are present.
           Exception NetworkMissingNode is raised if not all nodes
           0, 1, 2, ... max(self._nodes) are present.
        """
            
        # Check for missing nodes and branches, set the number of
        # nodes and branches, and check the data type (complex or real).
        # Computes also the number of extra nodes anb branches needed for ideal elements.
        assert self._net, "Network is empty!"
        
        actual_branches = set(self._list_branches)
        num_branches = max(actual_branches)
        expected_branches = set(b for b in range(1, num_branches + 1))
        if diff := actual_branches ^ expected_branches: raise NetworkMissingBranch(diff)        
        if diff := set(self._list_coupled_branches) - expected_branches: raise NetworkMissingBranch(diff)
        
        actual_nodes = set(self._list_nodes)
        num_nodes = max(actual_nodes)
        expected_nodes = set(n for n in range(0, num_nodes + 1))
        if diff := actual_nodes ^ expected_nodes: raise NetworkMissingNode(diff)
        
        # Set up the matrices to solve the network
        A = np.zeros((num_nodes + self._extra_nodes, num_branches + self._extra_branches), dtype=np.int_)
        data_type = np.complex_ if self._is_complex else np.float_
        ZB = np.zeros((num_branches + self._extra_branches, num_branches + self._extra_branches), dtype=data_type)
        EB1 = np.zeros(num_branches + self._extra_branches, dtype=data_type)
        JB1 = np.zeros(num_branches + self._extra_branches, dtype=data_type)
            
        def set_A(branch, from_node, to_node):
            if from_node > 0:
                A[from_node - 1, branch - 1] = 1
            if to_node > 0:
                A[to_node - 1, branch - 1] = -1
     
        for branch, element in self._net.items():
            match element:
                case Impedance():
                    set_A(branch, element.from_node, element.to_node)
                    ZB[branch - 1, branch - 1] = element.Z
                case Admitance():
                    set_A(branch, element.from_node, element.to_node)
                    ZB[branch - 1, branch - 1] = 1/element.Y
                case VoltageSource():
                    set_A(branch, element.from_node, element.to_node)
                    EB1[branch - 1] = element.E
                    ZB[branch - 1, branch - 1] = element.Z
                case VoltageSourceIdeal():    
                    num_branches += 1
                    num_nodes += 1
                    set_A(branch, element.from_node, num_nodes)
                    set_A(num_branches, num_nodes, element.to_node)
                    EB1[branch - 1] = element.E
                    ZB[branch - 1, branch - 1] = 1
                    ZB[num_branches - 1, num_branches - 1] = -1
                case CurrentSource():
                    set_A(branch, element.from_node, element.to_node)
                    JB1[branch - 1] = element.J
                    ZB[branch - 1, branch - 1] = 1/element.Y
                case CurrentSourceIdeal():   
                    num_branches += 1
                    set_A(branch, element.from_node, element.to_node)
                    set_A(num_branches, element.from_node, element.to_node)
                    JB1[branch - 1] = element.J
                    ZB[branch - 1, branch - 1] = 1
                    ZB[num_branches - 1, num_branches - 1] = -1
                case ShortCircuit():
                    num_branches += 1
                    num_nodes += 1
                    set_A(branch, element.from_node, num_nodes)
                    set_A(num_branches, num_nodes, element.to_node)
                    ZB[branch - 1, branch - 1] = 1
                    ZB[num_branches - 1, num_branches - 1] = -1
                case MutualCoupling():
                    b1, b2 = branch
                    ZB[b1 - 1, b2 - 1] = element.XM
                    ZB[b2 - 1, b1 - 1] = element.XM
                    
        # Solve the network
        YB = linalg.inv(ZB, overwrite_a=True)
        JB = YB @ EB1 + JB1
        self._ZN = linalg.inv(A @ YB @ A.T, overwrite_a=True)
        self._VN = self._ZN @ (-A @ JB)
        self._IB = YB @ (A.T @ self._VN) + JB
        self._is_solved = True
        
        # Adjust result for branches with ideal current sources
        for branch, element in self._net.items():
            if isinstance(element, CurrentSourceIdeal):
                self._IB[branch - 1] = element.J

    def thevenin(self, node_1, node_2):
        """ Returns the Thevenin equivalent circuit between two nodes
            (voltage source, and series impedance).
            Exception NetworkNotSolved is raised if the network is not solved.
            Exception NetworkMissingNode is raised if 'node_1' or 'node_2'
            do not exist.
        """
        if not self._is_solved: raise NetworkNotSolved
        if node_1 not in self._list_nodes: raise NetworkMissingNode({node_1})
        if node_2 not in self._list_nodes: raise NetworkMissingNode({node_2})
    
        match (node_1, node_2):
            case (n1, n2) if n1 == n2:
                return 0, 0
            case (0, n):
                return (-self._VN[n - 1],
                        self._ZN[n - 1, n - 1])
            case (n, 0):
                return (self._VN[n - 1],
                        self._ZN[n - 1, n - 1])
            case (n1, n2):
                return (self._VN[n1 - 1] - self._VN[n2 - 1],
                        self._ZN[n1 - 1, n1 - 1] +
                        self._ZN[n2 - 1, n2 - 1] -
                        self._ZN[n1 - 1, n2 - 1] -
                        self._ZN[n2 - 1, n1 - 1])
    
    def norton(self, node_1, node_2):
        """ Returns the Norton equivalent circuit between two nodes
            (current source, and parallel admitance).
            May raise the same exceptions as the function thevenin.
        """
        eth, zth = self.thevenin(node_1, node_2)
        match (eth, zth):
            case (0, 0):
                return 0, float('inf')
            case (eth, 0):
                return float('inf'), float('inf')
            case _:
                return eth/zth, 1/zth
           
    def node_voltage(self, node):
        """ Returns the voltage of a node (with respect to node 0)
            May raise the same exceptions as the function thevenin.
        """
        vn, _ = self.thevenin(node, 0)
        return vn
    
    def branch_voltage(self, branch):
        """ Returns the voltage across a branch.
            Exception NetworkNotSolved is raised if the network is not solved.
            Exception NetworkMissingBranch is raised if 'branch' does not exist.
        """
        if not self._is_solved: raise NetworkNotSolved
        if branch not in self._list_branches: raise NetworkMissingBranch({branch})
    
        element = self._net[branch]
        vn, _ = self.thevenin(element.from_node, element.to_node)
        return vn
    
    def branch_current(self, branch):
        """ Returns the current through a branch.
            Exception NetworkNotSolved is raised if the network is not solved.
            Exception NetworkMissingBranch is raised if 'branch' does not exist.
        """
        if not self._is_solved: raise NetworkNotSolved
        if diff := {branch} - set(self._list_branches): raise NetworkMissingBranch(diff)
    
        return self._IB[branch - 1]
    
    def results(self, polar=False, u_fmt='', u_unit='', u_scale=1,
                i_fmt='', i_unit='', i_scale=1):
        """ Returns a multi-line string with the voltages and currents of all 
            the branches.
            A different scale factor can be applied to voltages and currents.
            Polar display can be chosen, and different formatting for voltages
            and currents are available. Also, different  units can be added at
            the end of every voltage and at the end of every current.
            Exception NetworkNotSolved is raised if the network is not solved.
        """
        if not self._is_solved: raise NetworkNotSolved
        
        if u_unit: u_unit = ' ' + u_unit
        if i_unit: i_unit = ' ' + i_unit
        s = []
        for branch in range(1, self.num_branches + 1):
            u = self.branch_voltage(branch)*u_scale
            i = self.branch_current(branch)*i_scale
            if polar:
                u = ComplexD(u)
                i = ComplexD(i)
            s.append(f'U{branch} = ' + format(u, u_fmt) + u_unit + 
                     f',  I{branch} = ' + format(i, i_fmt) + i_unit)
        return '\n'.join(s)
   
    def save(self, filename):
        """ Saves the contents of the network to a binary file """
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
                if self._is_solved:
                    pickle.dump(self._ZN, fh, pickle.HIGHEST_PROTOCOL)
                    pickle.dump(self._VN, fh, pickle.HIGHEST_PROTOCOL)
                    pickle.dump(self._IB, fh, pickle.HIGHEST_PROTOCOL)
                return True
        except (OSError, pickle.PicklingError) as err:
            print(f'{filename}: save error: {err}')
            return False
    
    def load(self, filename):
        """ Loads the network with the contents of a binary file """
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
                if self._is_solved:
                    self._ZN = pickle.load(fh)
                    self._VN = pickle.load(fh)
                    self._IB = pickle.load(fh)
                return True
        except (OSError, pickle.UnpicklingError) as err:
            print(f'{filename}: load error: {err}')
            return False

        
if __name__ == '__main__':
    import doctest
    doctest.testmod()