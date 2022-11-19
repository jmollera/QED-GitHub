from jmb.ee import Network
from jmb.utils import ComplexD

net = Network()

net.add_voltage_source(1, 0, 1, 1.1, 0.25j)
net.add_voltage_source(2, 0, 2, 1.05+0.1j, 0.2j)
net.add_voltage_source(3, 0, 3, 1.08+0.12j, 0.25j)
net.add_impedance(4, 3, 4, 0.1j)
net.add_impedance(5, 2, 4, 0.405j)
net.add_impedance(6, 2, 3, 0.5j)
net.add_impedance(7, 1, 2, 0.16j)
net.add_current_source(8, 4, 0, 2-0.9j, 1/(-25j))
net.add_mutual_coupling(5, 6, 0.05j)

net.solve()

print(net, "\n")

branches = range(1, net.num_branches+1)
for branch in branches:
    print(f"Branch {branch} : Voltage = {ComplexD(net.branch_voltage(branch)):.4f} V")
print()
for branch in branches:
    print(f"Branch {branch} : Current = {ComplexD(net.branch_current(branch)):.4f} A")
