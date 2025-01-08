import unittest
import qed.eng_elec as ee


class TestFunctions(unittest.TestCase):
    def test_z_cable(self):
        with self.assertRaises(ee.EEInvalidArguments):
            # Invalid cable
            ee.z_cable(25, 100, 1)
        with self.assertRaises(ee.EEInvalidArguments):
            # Invalid section
            ee.z_cable(2500, 100, ee.CableType.A)

    def test_r_cable(self):
        with self.assertRaises(ee.EEInvalidArguments):
            # Invalid material
            ee.r_cable(25, 100, 15, 'Fe')

    def test_ezs_u(self):
        with self.assertRaises(ee.EENetNotSolvable):
            # No solution
            ee.ezs_u(0.4+0.3j, 10, 0.6+0.45j)

    def test_curve_51(self):
        with self.assertRaises(ee.EEInvalidArguments):
            # Invalid curve
            ee.curve_51(200, 1, 100)


class TestNetworkClass(unittest.TestCase):
    def test_add(self):
        net = ee.Network('Test network')
        with self.assertRaises(ee.EEInvalidArguments):
            # Wrong element
            net.add(25)
        with self.assertRaises(ee.EEInvalidArguments):
            # Wrong element
            net.add(25, from_to=(1, 0), branch=1)
        with self.assertRaises(ee.EEInvalidArguments):
            # Wrong element
            net.add(25, coupled_branches=(1, 2))
        with self.assertRaises(ee.EEInvalidArguments):
            # Wrong parameters
            net.add(ee.Impedance(1))
        with self.assertRaises(ee.EEInvalidArguments):
            # Wrong parameters
            net.add(ee.Impedance(1), from_to=(1, 2))
        with self.assertRaises(ee.EEInvalidArguments):
            # Wrong parameters
            net.add(ee.Impedance(1), branch=1)
        with self.assertRaises(ee.EEInvalidArguments):
            # Wrong parameters
            net.add(ee.Impedance(1), coupled_branches=(1, 2))
        with self.assertRaises(ee.EEInvalidArguments):
            # Wrong parameters
            net.add(ee.MutualCoupling(1j))
        with self.assertRaises(ee.EEInvalidArguments):
            # Wrong parameters
            net.add(ee.MutualCoupling(1j), from_to=(1, 2))
        with self.assertRaises(ee.EEInvalidArguments):
            # Wrong parameters
            net.add(ee.MutualCoupling(1j), branch=1)
        with self.assertRaises(ee.EEInvalidArguments):
            # Wrong parameters
            net.add(ee.MutualCoupling(1j), from_to=(1, 2), branch=1)
        with self.assertRaises(ee.EEInvalidArguments):
            # Node < 0
            net.add(ee.Impedance(1), from_to=(-1, 1), branch=1)
        with self.assertRaises(ee.EEInvalidArguments):
            # Node < 0
            net.add(ee.Impedance(1), from_to=(1, -1), branch=1)
        with self.assertRaises(ee.EEInvalidArguments):
            # Two equal nodes
            net.add(ee.Impedance(1), from_to=(1, 1), branch=1)
        with self.assertRaises(ee.EEInvalidArguments):
            # Branch < 1
            net.add(ee.Impedance(1), from_to=(1, 0), branch=0)
        with self.assertRaises(ee.EEInvalidArguments):
            # Wrong element
            net.add(ee.Impedance(1), from_to=(1,), branch=1)
        with self.assertRaises(ee.EEInvalidArguments):
            # Wrong element
            net.add(ee.Impedance(1), from_to=(1, 2, 3), branch=1)
        with self.assertRaises(ee.EEInvalidArguments):
            # Impedance is not pure imaginary
            net.add(ee.MutualCoupling(1 + 1j), coupled_branches=(1, 2))
        with self.assertRaises(ee.EEInvalidArguments):
            # Impedance is not pure imaginary
            net.add(ee.MutualCoupling(1), coupled_branches=(1, 2))
        with self.assertRaises(ee.EEInvalidArguments):
            # Branch < 1
            net.add(ee.MutualCoupling(1j), coupled_branches=(1, 0))
        with self.assertRaises(ee.EEInvalidArguments):
            # Branch < 1
            net.add(ee.MutualCoupling(1j), coupled_branches=(0, 1))
        with self.assertRaises(ee.EEInvalidArguments):
            # Two equal branches
            net.add(ee.MutualCoupling(1j), coupled_branches=(1, 1))
        with self.assertRaises(ee.EEInvalidArguments):
            # Two equal branches
            net.add(ee.MutualCoupling(1j), coupled_branches=(1,))
        with self.assertRaises(ee.EEInvalidArguments):
            # Two equal branches
            net.add(ee.MutualCoupling(1j), coupled_branches=(1, 2, 3))

    def test_remove(self):
        net = ee.Network('Test network')
        net.add(ee.Impedance(Z=10), from_to=(1, 2), branch=1)
        net.add(ee.MutualCoupling(X=3j), coupled_branches=(1, 2))
        with self.assertRaises(ee.EEInvalidArguments):
            # Two None arguments not allowed
            net.remove()
        with self.assertRaises(ee.EEInvalidArguments):
            # Two arguments different from None not allowed
            net.remove(branch=1, coupled_branches=(1, 2))
        with self.assertRaises(TypeError):
            # Arguments are keyword only
            net.remove(1)
        with self.assertRaises(TypeError):
            # Arguments are keyword only
            net.remove((1, 2))
        with self.assertRaises(TypeError):
            # Arguments are keyword only
            net.remove(1, (1, 2))
        with self.assertRaises(TypeError):
            # Arguments are keyword only
            net.remove((1, 2), branch=1)
        with self.assertRaises(TypeError):
            # Arguments are keyword only
            net.remove(1, coupled_branches=(1, 2))
        with self.assertRaises(ee.EENetMissingBranch):
            # Branch does not exist
            net.remove(branch=2)
        with self.assertRaises(ee.EENetMissingMutualCoupling):
            # Branch does not exist
            net.remove(coupled_branches=(1, 3))
        with self.assertRaises(ee.EENetMissingMutualCoupling):
            # Branch does not exist
            net.remove(coupled_branches=(3, 2))

    def test_solve(self):
        net = ee.Network('Test network')
        net.add(ee.Impedance(Z=10), from_to=(1, 0), branch=2)
        net.add(ee.Impedance(Z=10), from_to=(1, 2), branch=5)
        with self.assertRaises(ee.EENetMissingBranch):
            # Missing branches
            net.solve()
        net = ee.Network('Test network')
        net.add(ee.Impedance(Z=10), from_to=(1, 0), branch=1)
        net.add(ee.Impedance(Z=10), from_to=(1, 2), branch=2)
        net.add(ee.MutualCoupling(X=10j), coupled_branches=(1, 3))
        with self.assertRaises(ee.EENetMissingBranch):
            # Missing coupled branches
            net.solve()
        net = ee.Network('Test network')
        net.add(ee.MutualCoupling(X=10j), coupled_branches=(1, 2))
        with self.assertRaises(ee.EENetMissingBranch):
            # No branches
            net.solve()
        net = ee.Network('Test network')
        net.add(ee.Impedance(Z=10), from_to=(1, 0), branch=1)
        net.add(ee.Impedance(Z=10), from_to=(1, 3), branch=2)
        with self.assertRaises(ee.EENetMissingNode):
            # Missing nodes
            net.solve()
        net = ee.Network('Test network')
        with self.assertRaises(ee.EENetNotSolvable):
            # Empty network
            net.solve()
        net = ee.Network('Test network')
        net.add(ee.VoltageSource(E=10), from_to=(1, 0), branch=1)
        net.add(ee.ShortCircuit(), from_to=(1, 0), branch=2)
        with self.assertRaises(ee.EENetNotSolvable):
            # Ideal voltage source with a shortcircuit
            net.solve()
        net = ee.Network('Test network')
        net.add(ee.VoltageSource(E=10), from_to=(1, 0), branch=1)
        net.add(ee.VoltageSource(E=10), from_to=(1, 0), branch=2)
        with self.assertRaises(ee.EENetNotSolvable):
            # Two ideal voltage sources in parallel
            net.solve()
        net = ee.Network('Test network')
        net.add(ee.CurrentSource(J=10), from_to=(0, 1), branch=1)
        with self.assertRaises(ee.EENetNotSolvable):
            # Just one branch with an ideal current source
            net.solve()
        net = ee.Network('Test network')
        net.add(ee.CurrentSource(J=10), from_to=(0, 1), branch=1)
        net.add(ee.CurrentSource(J=10), from_to=(1, 2), branch=2)
        net.add(ee.Impedance(2), from_to=(2, 0), branch=3)
        with self.assertRaises(ee.EENetNotSolvable):
            # Two ideal current sources in series
            net.solve()

    def test_thevenin(self):
        net = ee.Network('Test network')
        net.add(ee.Impedance(Z=10), from_to=(1, 0), branch=1)
        with self.assertRaises(ee.EENetNotSolved):
            # Net not solved yet
            net.thevenin(1, 0)
        net.solve()
        with self.assertRaises(ee.EENetMissingNode):
            # Missing node
            net.thevenin(1, 2)
        self.assertEqual((0, 0), net.thevenin(1, 1))

    def test_norton(self):
        net = ee.Network('Test network')
        net.add(ee.Impedance(Z=10), from_to=(1, 0), branch=1)
        with self.assertRaises(ee.EENetNotSolved):
            # Net not solved yet
            net.norton(1, 0)
        net.solve()
        with self.assertRaises(ee.EENetMissingNode):
            # Missing node
            net.norton(1, 2)
        self.assertEqual((0, float('inf')), net.norton(1, 1))

    def test_voltage(self):
        net = ee.Network('Test network')
        net.add(ee.Impedance(Z=10), from_to=(1, 0), branch=1)
        with self.assertRaises(ee.EENetNotSolved):
            # Net not solved yet
            net.voltage(node=1)
        with self.assertRaises(ee.EENetNotSolved):
            # Net not solved yet
            net.voltage(branch=1)
        net.solve()
        with self.assertRaises(ee.EENetMissingBranch):
            # Missing branch
            net.voltage(branch=2)
        with self.assertRaises(ee.EENetMissingNode):
            # Missing node
            net.voltage(node=2)
        with self.assertRaises(ee.EEInvalidArguments):
            # Two None arguments not allowed
            net.voltage()
        with self.assertRaises(ee.EEInvalidArguments):
            # Two arguments different from None not allowed
            net.voltage(node=1, branch=1)
        with self.assertRaises(TypeError):
            # Arguments are keyword only
            net.voltage(1)
        with self.assertRaises(TypeError):
            # Arguments are keyword only
            net.voltage(1, 1)
        with self.assertRaises(TypeError):
            # Arguments are keyword only
            net.voltage(1, node=1)
        with self.assertRaises(TypeError):
            # Arguments are keyword only
            net.voltage(1, branch=1)

    def test_current(self):
        net = ee.Network('Test network')
        net.add(ee.Impedance(Z=10), from_to=(1, 0), branch=1)
        with self.assertRaises(ee.EENetNotSolved):
            # Net not solved yet
            net.voltage(branch=1)
        net.solve()
        with self.assertRaises(ee.EENetMissingBranch):
            # Missing branch
            net.current(branch=2)


class TestMotor3phClass(unittest.TestCase):
    def setUp(self):
        R_1 = 0.095  # Ω/phase
        X_1 = 0.680  # Ω/phase
        R_2 = 0.3  # Ω/phase
        X_2 = 0.672  # Ω/phase
        R_Fe = 620  # Ω/phase
        X_m = 18.7  # Ω/phase
        p = 4  # number of poles
        f = 50  # Hz
        self.motor = ee.Motor3ph(p, f, R_1, X_1, R_2, X_2, R_Fe, X_m)

    def test_speed(self):
        with self.assertRaises(ee.EEInvalidArguments):
            # Invalid from_to
            self.motor.convert(1450, from_to=1)
        # Starting
        self.assertAlmostEqual(1.0, self.motor.convert(0, ee.Speed.rad_per_s_to_slip))
        self.assertAlmostEqual(1.0, self.motor.convert(0, ee.Speed.rpm_to_slip))
        self.assertAlmostEqual(0.0, self.motor.convert(1, ee.Speed.slip_to_rad_per_s))
        self.assertAlmostEqual(0.0, self.motor.convert(0, ee.Speed.rpm_to_rad_per_s))
        self.assertAlmostEqual(0.0, self.motor.convert(1, ee.Speed.slip_to_rpm))
        self.assertAlmostEqual(0.0, self.motor.convert(0, ee.Speed.rad_per_s_to_rpm))
        # Synchronous speed
        self.assertAlmostEqual(0.0, self.motor.convert(self.motor.ω_m_sync, ee.Speed.rad_per_s_to_slip))
        self.assertAlmostEqual(0.0, self.motor.convert(self.motor.n_m_sync, ee.Speed.rpm_to_slip))
        self.assertAlmostEqual(self.motor.ω_m_sync, self.motor.convert(0, ee.Speed.slip_to_rad_per_s))
        self.assertAlmostEqual(self.motor.ω_m_sync, self.motor.convert(self.motor.n_m_sync, ee.Speed.rpm_to_rad_per_s))
        self.assertAlmostEqual(self.motor.n_m_sync, self.motor.convert(0, ee.Speed.slip_to_rpm))
        self.assertAlmostEqual(self.motor.n_m_sync, self.motor.convert(self.motor.ω_m_sync, ee.Speed.rad_per_s_to_rpm))


if __name__ == '__main__':
    unittest.main()
