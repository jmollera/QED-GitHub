from qed.eng_elec import awg_to_mm2

AWG = 14

S = awg_to_mm2(AWG)
print(f'AWG {AWG} equival a {S:.1f} mm²')
