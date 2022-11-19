from qed.utils import ComplexD
from qed.eng_elec import Network, VoltageSource, CurrentSourceIdeal, Impedance

circuit = Network()
circuit.add(VoltageSource(from_node=0, to_node=1, E=5, Z=4), branch=1)
circuit.add(CurrentSourceIdeal(from_node=0, to_node=2, J=2), branch=2)
circuit.add(Impedance(from_node=1, to_node=2, Z=5j), branch=3)
circuit.add(Impedance(from_node=1, to_node=0, Z=-2j), branch=4)
circuit.solve()

I_c = ComplexD(circuit.branch_current(branch=4))
print(f'Corrent pel condensador: {I_c:.3f/.2f} A')