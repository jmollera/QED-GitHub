from qed.utils import ComplexD
from qed.eng_elec import Network, VoltageSource, CurrentSourceIdeal, Impedance

circuit = Network()
circuit.add(VoltageSource(from_to=(0, 1), E=5, Z=4), branch=1)
circuit.add(CurrentSourceIdeal(from_to=(0, 2), J=2), branch=2)
circuit.add(Impedance(from_to=(1, 2), Z=5j), branch=3)
circuit.add(Impedance(from_to=(1, 0), Z=-2j), branch=4)
circuit.solve()

I_c = ComplexD(circuit.current(branch=4))
print(f'Corrent pel condensador: {I_c:.3f/.2f} A')
