__author__ = 'Josep Mollera Barriga'
__version__ = '15.1'
__name__ = 'qed.utils'
__all__ = ['Complex', 'DEGREE_SYMBOL', 'ANGLE_SYMBOL']

import math
import cmath

DEGREE_SYMBOL = '°'
ANGLE_SYMBOL = '∠'


class Complex(complex):
    """
    Handles polar complex numbers in degrees or radians.

    >>> z = Complex(4, 45)  # same as Complex(4, 45, in_rad=False)
    >>> z.mod
    4.0
    >>> z.arg
    45.0
    >>> z.in_rad
    False
    >>> z  # __repr__
    Complex(4.0, 45.0, in_rad=False)
    >>> str(z)  # __str__
    '4.0∠45.0°'
    >>> print(z)  # __str__
    4.0∠45.0°
    >>> print(f'z = {z:.3f}')  # __format__
    z = 4.000∠45.000°
    >>> print(f'z = {z:.1f/.4f}')  # __format__
    z = 4.0∠45.0000°
    >>> z = Complex(4, 0.785, in_rad=True)
    >>> z.mod
    4.0
    >>> z.arg
    0.785
    >>> z.in_rad
    True
    >>> z  # __repr__
    Complex(4.0, 0.785, in_rad=True)
    >>> str(z)  # __str__
    '4.0∠0.785'
    >>> print(z)  # __str__
    4.0∠0.785
    >>> print(f'z = {z:.3f}')  # __format__
    z = 4.000∠0.785
    >>> print(f'z = {z:.1f/.4f}')  # __format__
    z = 4.0∠0.7850
    """

    def __new__(cls, *args: int | float | complex | str,
                in_rad: bool = False) -> complex | None:
        """
        Create a polar complex number in degrees or radians.

        Parameters
        ----------
        args:   One value compatible with complex (int|float|complex|str),
                or two values: module (positive) and argument compatible
                with float (int|float|str).
        in_rad: True -> argument in radians, False -> argument in degrees.
                The default value is False (i.e., argument in degrees).
        """
        try:
            match args:
                case [z]:
                    return super().__new__(cls, z)
                case [mod, arg]:
                    if (fmod := float(mod)) < 0:
                        raise ValueError
                    farg = float(arg)
                    farg_rad = farg if in_rad else math.radians(farg)
                    return super().__new__(cls, cmath.rect(fmod, farg_rad))
                case _:
                    raise TypeError
        except (TypeError, ValueError):
            print(f'{cls.__name__} called with bad parameters.')
            return None

    def __init__(self, *args: int | float | complex | str,
                 in_rad: bool = False) -> None:
        """ __init__ is automatically called by __new__."""
        self.__in_rad = in_rad

    def __format_polar(self, fmt: str = '') -> str:
        """Format polar representation with a given format string.

        Parameters
        ----------
        fmt: A string using Python's Format Specification Mini-Language,
             with the addition of "/" to separate the magnitude format and
             the argument format. If "/" is not present, the same format is
             used for both the magnitude and the argument. The default value is ''.
             """
        if '/' in fmt:
            fmt_mod, fmt_arg = fmt.split(sep='/', maxsplit=1)
        else:
            fmt_mod = fmt_arg = fmt
        unit = '' if self.__in_rad else DEGREE_SYMBOL
        return f'{self.mod:{fmt_mod}}{ANGLE_SYMBOL}{self.arg:{fmt_arg}}{unit}'

    def __repr__(self) -> str:
        """Class representation function"""
        return f'{type(self).__name__}({self.mod!r}, {self.arg!r}, in_rad={self.__in_rad!r})'

    def __str__(self) -> str:
        """Class printing function"""
        return self.__format_polar()

    def __format__(self, fmt: str = '') -> str:
        """ Class format function."""
        return self.__format_polar(fmt)

    @property
    def mod(self) -> float:
        """Modulus of the complex number."""
        return abs(self)

    @property
    def arg(self) -> float:
        """Argument of the complex number."""
        phase = cmath.phase(self)
        return phase if self.__in_rad else math.degrees(phase)

    @property
    def in_rad(self) -> bool:
        """Angular mode of the complex number."""
        return self.__in_rad


if __name__ == '__main__':
    import doctest
    doctest.testmod()
