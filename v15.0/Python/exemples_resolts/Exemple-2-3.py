from qed.eng_elec import ezs_u

E = 0.4+0.3j
Z = 0.1j
S = 0.6+0.45j

U = ezs_u(E, Z, S)
I = S.conjugate()/U.conjugate()
Z_S = U/I
print(f'U   = {U:.4f}')
print(f'I   = {I:.4f}')
print(f'Z_S = {Z_S:.4f}')
