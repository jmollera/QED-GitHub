import unittest
import random
import cmath
import math
from qed.utils import Complex


class TestComplexClass(unittest.TestCase):
    def test_new(self):
        self.assertIsNone(Complex('1x'))
        self.assertIsNone(Complex('1x', 45))
        self.assertIsNone(Complex(1, '45x'))
        self.assertIsNone(Complex(1, 2, 3, 4))
        self.assertIsNone(Complex(-1, 45))
        self.assertIsNone(Complex('1x', in_rad=True))
        self.assertIsNone(Complex('1x', 45, in_rad=True))
        self.assertIsNone(Complex(1, '45x', in_rad=True))
        self.assertIsNone(Complex(1, 2, 3, in_rad=True))
        self.assertIsNone(Complex(-1, 45, in_rad=True))
        self.assertIsNone(Complex('1x', in_rad=False))
        self.assertIsNone(Complex('1x', 45, in_rad=False))
        self.assertIsNone(Complex(1, '45x', in_rad=False))
        self.assertIsNone(Complex(1, 2, 3, in_rad=False))
        self.assertIsNone(Complex(-1, 45, in_rad=False))
        for _ in range(0, 1000):
            z = complex(random.randint(-100, 100), random.randint(-100, 100))
            m = abs(z)
            a_rad = cmath.phase(z)
            a_deg = math.degrees(a_rad)
            self.assertAlmostEqual(z, Complex(z))
            self.assertAlmostEqual(z, Complex(z, in_rad=False))
            self.assertAlmostEqual(z, Complex(z, in_rad=True))
            self.assertAlmostEqual(z, Complex(m, a_deg))
            self.assertAlmostEqual(z, Complex(m, a_deg, in_rad=False))
            self.assertAlmostEqual(z, Complex(m, a_rad, in_rad=True))


if __name__ == '__main__':
    unittest.main()
