from jmb.ee import Network

net = Network()

net.add_voltage_source(1, 0, 1, 5, 1)
net.add_impedance(2, 2, 1, 1)
net.add_impedance(3, 1, 0, 1)
net.add_impedance(4, 3, 0, 1)
net.add_impedance(5, 2, 3, 1)
net.add_voltage_source(6, 3, 2, 2, 1)
print(net, "\n")

net.solve()

branches = range(1, net.num_branches+1)
for branch in branches:
    print(f"Branch {branch} : Voltage = {net.branch_voltage(branch):.3f} V")
print()
for branch in branches:
    print(f"Branch {branch} : Current = {net.branch_current(branch):.3f} A")
