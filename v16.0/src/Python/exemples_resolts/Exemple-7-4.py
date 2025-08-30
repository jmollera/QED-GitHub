from qed.eng_elec import mm2_to_awg

S = 4  # mm²

AWG = mm2_to_awg(S)
print(f'{S} mm² equival a AWG {AWG}')
