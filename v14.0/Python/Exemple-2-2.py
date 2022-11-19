import qed.eng_elec as ee

Z_A, Z_B, Z_C = ee.D_to_Y(10, -10j, -10j)
print(f'Z_A = {Z_A} Ω')
print(f'Z_B = {Z_B} Ω')
print(f'Z_C = {Z_C} Ω')