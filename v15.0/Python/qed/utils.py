__author__ = 'Josep Mollera Barriga'
__version__ = '15.0'
__name__ = 'qed.utils'
__all__ = ['Complex']

import math
import cmath


class Complex(complex):
    """
    Class to handle polar complex numbers in degrees or radians.

    >>> z = Complex(4, 45)  # or Complex(4, 45, in_rad=False)
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
    >>> w = Complex(4, 0.785, in_rad=True)
    >>> w.mod
    4.0
    >>> w.arg
    0.785
    >>> w.in_rad
    True
    >>> w  # Calls __repr__
    Complex(4.0, 0.785, in_rad=True)
    >>> str(w)  # __str__
    '4.0∠0.785'
    >>> print(w)  # __str__
    4.0∠0.785
    >>> print(f'w = {w:.3f}')  # __format__
    w = 4.000∠0.785
    >>> print(f'w = {w:.1f/.4f}')  # __format__
    w = 4.0∠0.7850
    """

    def __new__(cls, *args: int | float | complex | str,
                in_rad: bool = False) -> complex | None:
        """
        Create a polar complex number in degrees or radians.

        Parameters
        ----------
        args:   One value compatible with complex (int|float|complex|str),
                or two values: module (positive) and argument, compatibles
                with float (int|float|str).
        in_rad: True -> argument in radians, False -> argument in degrees.
                Default value is False (i.e. argument in degrees).
        """
        try:
            match args:
                case [z]:
                    return super().__new__(cls, z)
                case [mod, arg]:
                    if (fmod := float(mod)) < 0:
                        raise ValueError
                    farg = float(arg)
                    if in_rad:
                        return super().__new__(cls, cmath.rect(fmod, farg))
                    else:
                        return super().__new__(cls, cmath.rect(fmod, math.radians(farg)))
                case _:
                    raise TypeError
        except (TypeError, ValueError):
            print(f'{cls.__name__} called with bad parameters.')
            return None

    def __init__(self, *args: int | float | complex | str,
                 in_rad: bool = False) -> None:
        """ __init__ is automatically called by __new__. """
        self.__arg_in_rad = True if in_rad else False

    def __repr__(self) -> str:
        """Class representation function"""
        return f'{type(self).__name__}({self.mod!r}, {self.arg!r}, in_rad={self.in_rad!r})'

    def __str__(self) -> str:
        """Class printing function"""
        unit = '' if self.in_rad else '°'
        return f'{self.mod}∠{self.arg}{unit}'

    def __format__(self, fmt: str = '') -> str:
        """
        Format a polar complex number.

        Parameters
        ----------
        fmt: A string using the Python's Format Specification Mini-Language,
             with the addition of "/" to separate the magnitude format and
             the argument format. If "/" is not present, the same format is
             used for both the magnitude and the argument. Default value is ''.
        """
        if (slash_position := fmt.find('/')) == -1:
            fmt_mod = fmt_arg = fmt
        else:
            fmt_mod = fmt[:slash_position]
            fmt_arg = fmt[slash_position + 1:]
        unit = '' if self.in_rad else '°'
        return f'{self.mod:{fmt_mod}}∠{self.arg:{fmt_arg}}{unit}'

    @property
    def mod(self) -> float:
        """Modulus of the complex number."""
        return abs(self)

    @property
    def arg(self) -> float:
        """Argument of the complex number."""
        if self.in_rad:
            return cmath.phase(self)
        else:
            return math.degrees(cmath.phase(self))

    @property
    def in_rad(self) -> bool:
        """Angular mode of the complex number."""
        return self.__arg_in_rad


if __name__ == '__main__':
    import doctest
    doctest.testmod()
