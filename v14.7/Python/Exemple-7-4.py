import qed.eng_elec as ee

S = 4  # mm²

AWG = ee.mm2_to_awg(S)
print(f'{S} mm² equival a AWG {AWG}')
