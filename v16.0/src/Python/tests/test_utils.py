import unittest
import random
import cmath
import math
from qed.utils import Complex, ComplexR, ComplexD


class TestComplexClass(unittest.TestCase):
    def setUp(self):
        Complex.angle_symbol = '∠'
        Complex.degree_symbol = '°'
        Complex.rad_symbol = ''

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
            self.assertAlmostEqual(z, Complex(m, a_deg))
            self.assertAlmostEqual(z, Complex(m, a_deg, in_rad=False))
            self.assertAlmostEqual(z, Complex(m, a_rad, in_rad=True))

    def test_str(self):
        self.assertEqual(str(Complex(1.0, 45)), '1.0∠45.0°')
        self.assertEqual(str(Complex(1.0, 45, in_rad=False)), '1.0∠45.0°')
        self.assertEqual(str(Complex(1.0, 0.5, in_rad=True)), '1.0∠0.5')
        Complex.rad_symbol = ' rad'
        Complex.degree_symbol = ' deg'
        Complex.angle_symbol = ' @ '
        self.assertEqual(str(Complex(1.0, 45)), '1.0 @ 45.0 deg')
        self.assertEqual(str(Complex(1.0, 45, in_rad=False)), '1.0 @ 45.0 deg')
        self.assertEqual(str(Complex(1.0, 0.5, in_rad=True)), '1.0 @ 0.5 rad')


class TestComplexRClass(unittest.TestCase):
    def setUp(self):
        ComplexR.angle_symbol = '∠'
        ComplexR.rad_symbol = ''

    def test_new(self):
        self.assertIsNone(ComplexR('1x'))
        self.assertIsNone(ComplexR('1x', 1.5))
        self.assertIsNone(ComplexR(1, '1.5x'))
        self.assertIsNone(ComplexR(1, 1.5, 3))
        self.assertIsNone(ComplexR(-1, 1.5))
        for _ in range(0, 1000):
            z = complex(random.randint(-100, 100), random.randint(-100, 100))
            m = abs(z)
            a_rad = cmath.phase(z)
            self.assertAlmostEqual(z, ComplexR(z))
            self.assertAlmostEqual(z, ComplexR(m, a_rad))

    def test_str(self):
        self.assertEqual(str(ComplexR(1.0, 0.5)), '1.0∠0.5')
        ComplexR.rad_symbol = ' rad'
        ComplexR.angle_symbol = ' @ '
        self.assertEqual(str(ComplexR(1.0, 0.5)), '1.0 @ 0.5 rad')


class TestComplexDClass(unittest.TestCase):
    def setUp(self):
        ComplexD.angle_symbol = '∠'
        ComplexD.degree_symbol = '°'

    def test_new(self):
        self.assertIsNone(ComplexD('1x'))
        self.assertIsNone(ComplexD('1x', 45))
        self.assertIsNone(ComplexD(1, '45x'))
        self.assertIsNone(ComplexD(1, 2, 3))
        self.assertIsNone(ComplexD(-1, 45))
        for _ in range(0, 1000):
            z = complex(random.randint(-100, 100), random.randint(-100, 100))
            m = abs(z)
            a_deg = math.degrees(cmath.phase(z))
            self.assertAlmostEqual(z, ComplexD(z))
            self.assertAlmostEqual(z, ComplexD(m, a_deg))

    def test_str(self):
        self.assertEqual(str(ComplexD(1.0, 45)), '1.0∠45.0°')
        ComplexD.degree_symbol = ' deg'
        ComplexD.angle_symbol = ' @ '
        self.assertEqual(str(ComplexD(1.0, 45)), '1.0 @ 45.0 deg')


if __name__ == '__main__':
    unittest.main()
