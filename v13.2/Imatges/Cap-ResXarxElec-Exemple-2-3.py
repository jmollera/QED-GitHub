from jmb.ee import Network


net = Network()

net.add_voltage_source(1, 0, 1, 200, 10)
net.add_voltage_source(2, 1, 2, -50, 20j)
net.add_impedance(3, 2, 0, 5j)
net.add_impedance(4, 1, 2, 20)
net.add_current_source(5, 0, 2, 4, 1/10)

net.solve()

print(net, "\n")

branches = range(1, net.num_branches+1)
for branch in (2, 5):
    print(f"Branch {branch} : Current = {net.branch_current(branch):.4f} A")


net.add_mutual_coupling(2, 3, 5j)

net.solve()

print("\n")
print(net, "\n")

eth, zth = net.thevenin(1,2)
jno, yno = net.norton(1,2)
print(f"E_Th(1,2) = {eth:.4f} V")
print(f"Z_Th(1,2) = {zth:.4f} Ω")
print(f"J_No(1,2) = {jno:.4f} A")
print(f"Y_No(1,2) = {yno:.4f} S")
