import qed.eng_elec as ee

AWG = 14

S = ee.awg_to_mm2(AWG)
print(f'AWG {AWG} equival a {S:.1f} mm²')
