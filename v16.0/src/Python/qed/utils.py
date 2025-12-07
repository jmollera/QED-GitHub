__author__ = 'Josep Mollera Barriga'
__version__ = '16.0'
__name__ = 'qed.utils'
__all__ = ['Complex', 'ComplexR', 'ComplexD']

import math
import cmath


class Complex(complex):
    """
    Handle polar complex numbers in degrees or radians.

    >>> z = Complex(4, 45, in_rad=False)
    >>> z.mod
    4.0
    >>> z.arg
    45.0
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

    angle_symbol:  str = '∠'  # Symbol to be used to separate module and argument.
    degree_symbol: str = '°'  # Symbol to be used for degrees.
    rad_symbol:    str = ''   # Symbol to be used for radians.

    def __new__(cls, *args: int | float | complex | str,
                in_rad: bool) -> complex | None:
        """
        Create a polar complex number in degrees or radians.

        Parameters
        ----------
        args:   One value compatible with complex (int|float|complex|str),
                or two values: module (positive) and argument compatible
                with float (int|float|str).
        in_rad: True -> argument in radians, False -> argument in degrees.
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

    def __repr__(self) -> str:
        """Class representation function"""
        return f'{type(self).__name__}({self.mod!r}, {self.arg!r}, in_rad={self.__in_rad!r})'

    def __str__(self) -> str:
        """Class printing function"""
        return self.__format_polar(fmt='')

    def __format__(self, fmt: str = '') -> str:
        """ Class format function."""
        return self.__format_polar(fmt=fmt)

    def __format_polar(self, fmt: str) -> str:
        """Format polar representation with a given format string.

        Parameters
        ----------
        fmt: A string using Python's Format Specification Mini-Language,
             with the addition of "/" to separate the magnitude format and
             the argument format. If "/" is not present, the same format is
             used for both the magnitude and the argument.
             """
        if '/' in fmt:
            fmt_mod, fmt_arg = fmt.split(sep='/', maxsplit=1)
        else:
            fmt_mod = fmt_arg = fmt
        unit = self.rad_symbol if self.__in_rad else self.degree_symbol
        return f'{self.mod:{fmt_mod}}{self.angle_symbol}{self.arg:{fmt_arg}}{unit}'

    @property
    def mod(self) -> float:
        """Modulus of the complex number."""
        return abs(self)

    @property
    def arg(self) -> float:
        """Argument of the complex number using the defined angular mode."""
        phase = cmath.phase(self)
        return phase if self.__in_rad else math.degrees(phase)


class ComplexR(Complex):
    """
    Handle polar complex numbers in radians.
    Subclass of Complex where in_rad is always True.

    >>> z = ComplexR(4, 1.5)
    >>> z.mod
    4.0
    >>> z.arg
    1.5
    >>> z  # __repr__
    ComplexR(4.0, 1.5)
    >>> str(z)  # __str__
    '4.0∠1.5'
    >>> print(z)  # __str__
    4.0∠1.5
    >>> print(f'z = {z:.3f}')  # __format__
    z = 4.000∠1.500
    >>> print(f'z = {z:.1f/.4f}')  # __format__
    z = 4.0∠1.5000
    """

    def __new__(cls, *args: int | float | complex | str) -> complex | None:
        """Call the parent function with in_rad=True."""
        return super().__new__(cls, *args, in_rad=True)

    def __init__(self, *args: int | float | complex | str) -> None:
        """Call the parent function with in_rad=True."""
        super().__init__(*args, in_rad=True)

    def __repr__(self) -> str:
        """Overwrite the parent function so that not to include in_rad=True."""
        return f'{type(self).__name__}({self.mod!r}, {self.arg!r})'


class ComplexD(Complex):
    """
    Handle polar complex numbers in degrees.
    Subclass of Complex where in_rad is always False.

    >>> z = ComplexD(4, 45)
    >>> z.mod
    4.0
    >>> z.arg
    45.0
    >>> z  # __repr__
    ComplexD(4.0, 45.0)
    >>> str(z)  # __str__
    '4.0∠45.0°'
    >>> print(z)  # __str__
    4.0∠45.0°
    >>> print(f'z = {z:.3f}')  # __format__
    z = 4.000∠45.000°
    >>> print(f'z = {z:.1f/.4f}')  # __format__
    z = 4.0∠45.0000°
    """

    def __new__(cls, *args: int | float | complex | str) -> complex | None:
        """Call the parent function with in_rad=False."""
        return super().__new__(cls, *args, in_rad=False)

    def __init__(self, *args: int | float | complex | str) -> None:
        """Call the parent function with in_rad=False."""
        super().__init__(*args, in_rad=False)

    def __repr__(self) -> str:
        """Overwrite the parent function so that not to include in_rad=False."""
        return f'{type(self).__name__}({self.mod!r}, {self.arg!r})'


if __name__ == '__main__':
    import doctest
    doctest.testmod()
