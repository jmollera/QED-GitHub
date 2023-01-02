__author__ = 'Josep Mollera Barriga'
__version__ = '14.2'
__all__ = ['ComplexR', 'ComplexD']

import math
import cmath


class ComplexR(complex):
    """
    Class to handle polar complex numbers in radians.

    Reimplemented complex built-in functions
    ----------------------------------------
    __repr__, __str__, __iter__, __format__

    Properties
    ----------
    mod : float
        Returns the module.
    arg : float
        Returns the argument in radians.

    Examples
    --------
    >>> z = ComplexR(4, 1.5)
    >>> z.mod
    4.0
    >>> z.arg
    1.5
    >>> m, a = z  # __iter__
    >>> print(f'mod={m}, arg={a} rad')
    mod=4.0, arg=1.5 rad
    >>> z  # __repr__
    ComplexR(mod: 4.0, arg: 1.5 rad)
    >>> str(z)  # __str__
    '4.0∠1.5'
    >>> print(z)  # __str__
    4.0∠1.5
    >>> print(f'z = {z:.3f}')  # __format__
    z = 4.000∠1.500
    >>> print(f'z = {z:.1f/.4f}')  # __format__
    z = 4.0∠1.5000
    """

    def __new__(cls, *args):
        """Create a polar complex number in radians."""
        try:
            match args:
                case [z]:
                    # Assuming int, float, complex, str
                    return super().__new__(cls, z)
                case [mod, arg]:
                    # Assuming module and argument in radians
                    return super().__new__(cls, cmath.rect(float(mod), float(arg)))
                case _:
                    raise TypeError
        except TypeError:
            print('Only 1 or 2 values allowed.',
                  ' 1 value:  int, float, complex, str',
                  ' 2 values: module and argument in radians',
                  sep='\n')
            return None

    def __repr__(self):
        return f'{type(self).__name__}(mod: {self.mod!r}, arg: {self.arg!r} rad)'

    def __str__(self):
        return f'{self.mod}∠{self.arg}'

    def __iter__(self):
        return (x for x in (self.mod, self.arg))

    def __format__(self, fmt=''):
        """Format a polar complex number in radians.

        The Python's Format Specification Mini-Language is used, with the
        addition of "/" to separate the magnitude and the argument formats.
        If "/" is not present, the same format in used for both magnitude
        and argument.
        """
        if (pos := fmt.find('/')) == -1:
            fmt_mod = fmt_arg = fmt
        else:
            fmt_mod = fmt[:pos]
            fmt_arg = fmt[pos+1:]
        return ''.join([format(self.mod, fmt_mod), '∠',
                        format(self.arg, fmt_arg)])

    @property
    def mod(self):
        """Modulus of the complex number."""
        return abs(self)

    @property
    def arg(self):
        """Argument of the complex number in radians."""
        return cmath.phase(self)


class ComplexD(complex):
    """
    Class to handlepolar complex numberst in degrees.

    Reimplemented complex built-in functions
    ----------------------------------------
    __repr__, __str__, __iter__, __format__

    Properties
    ----------
    mod : float
        Returns the module.
    arg : float
        Returns the argument in degrees.

    Examples
    --------
    >>> z = ComplexD(4, 45)
    >>> z.mod
    4.0
    >>> z.arg
    45.0
    >>> m, a = z  # __iter__
    >>> print(f'mod={m}, arg={a}°')
    mod=4.0, arg=45.0°
    >>> z  # __repr__
    ComplexD(mod: 4.0, arg: 45.0 deg)
    >>> str(z)  # __str__
    '4.0∠45.0°'
    >>> print(z)  # __str__
    4.0∠45.0°
    >>> print(f'z = {z:.3f}')  # __format__
    z = 4.000∠45.000°
    >>> print(f'z = {z:.1f/.4f}')  # __format__
    z = 4.0∠45.0000°
    """

    def __new__(cls, *args):
        """Create a polar complex number in degrees."""
        try:
            match args:
                case [z]:
                    # Assuming int, float, complex, str
                    return super().__new__(cls, z)
                case [mod, arg]:
                    # Assuming module and argument in degrees
                    return super().__new__(cls, cmath.rect(float(mod), math.radians(float(arg))))
                case _:
                    raise TypeError
        except TypeError:
            print('Only 1 or 2 values allowed',
                  '  1 value:  int, float, complex, str',
                  '  2 values: module and argument in degrees',
                  sep='\n')
            return None

    def __repr__(self):
        return f'{type(self).__name__}(mod: {self.mod!r}, arg: {self.arg!r} deg)'

    def __str__(self):
        return f'{self.mod}∠{self.arg}°'

    def __iter__(self):
        return (x for x in (self.mod, self.arg))

    def __format__(self, fmt=''):
        """Format a polar complex number in degrees.

        The Python's Format Specification Mini-Language is used, with the
        addition of "/" to separate the magnitude and the argument formats.
        If "/" is not present, the same format in used for both magnitude
        and argument.
        """
        if (pos := fmt.find('/')) == -1:
            fmt_mod = fmt_arg = fmt
        else:
            fmt_mod = fmt[:pos]
            fmt_arg = fmt[pos+1:]
        return ''.join([format(self.mod, fmt_mod), '∠',
                        format(self.arg, fmt_arg), '°'])

    @property
    def mod(self):
        """Modulus of the complex number."""
        return abs(self)

    @property
    def arg(self):
        """Argument of the complex number in degrees."""
        return math.degrees(cmath.phase(self))


if __name__ == '__main__':
    import doctest
    doctest.testmod()
