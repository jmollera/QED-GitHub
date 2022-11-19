from jmb.ee import Network
from jmb.utils import ComplexD

net = Network()

net.add_voltage_source(1, 0, 1, 200, 10)
net.add_voltage_source(2, 1, 2, -50, 20j)
net.add_impedance(3, 2, 0, 5j)
net.add_impedance(4, 1, 2, 20)
net.add_current_source(5, 0, 2, 4, 1/10)
print("\n", net, "\n")

net.solve()

branches = range(1, net.num_branches+1)
for branch in (2, 5):
    print(f"Branch {branch} : Current = {ComplexD(net.branch_current(branch)):.4f} A")



net.add_mutual_coupling(2, 3, 5j)
print("\n", net, "\n")
net.solve()

eth, zth = net.thevenin(1,2)
jno, yno = net.norton(1,2)
print(f"E_Th(1,2) = {ComplexD(eth):.4f} V")
print(f"Z_Th(1,2) = {ComplexD(zth):.4f} Ω")
print(f"J_No(1,2) = {ComplexD(jno):.4f} A")
print(f"Y_No(1,2) = {ComplexD(yno):.4f} S")
