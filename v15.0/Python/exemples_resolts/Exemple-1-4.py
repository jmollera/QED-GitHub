from qed.utils import Complex
from qed.eng_elec import Network, VoltageSource, CurrentSource, Impedance

circuit = Network()
circuit.add(VoltageSource(E=5, Z=4), from_to=(0, 1), branch=1)
circuit.add(CurrentSource(J=2), from_to=(0, 2), branch=2)
circuit.add(Impedance(Z=5j), from_to=(1, 2), branch=3)
circuit.add(Impedance(Z=-2j), from_to=(1, 0), branch=4)
circuit.solve()

I_c = Complex(circuit.current(branch=4))
print(f'Corrent pel condensador: {I_c:.3f/.2f} A')
