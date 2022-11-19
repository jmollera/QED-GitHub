from jmb.ee import Network

net = Network()

net.add_voltage_source(1, 0, 1, 125.1, 0.034)
net.add_voltage_source(2, 0, 1, 124.8, 0.041)
net.add_voltage_source(3, 0, 1, 125.2, 0.029)
print(net, "\n")

net.solve()

eth, zth = net.thevenin(1,0)
jno, yno = net.norton(1,0)
print(f"E_Th = {eth} V")
print(f"Z_Th = {zth*1000} mΩ")
print(f"J_No = {jno/1000} kA")
print(f"Y_No = {yno} S")

rq = 50
iq = eth/(zth+rq)
uq = rq*iq
print(f"I_Q = {iq} A")
print(f"U_Q = {uq} V")