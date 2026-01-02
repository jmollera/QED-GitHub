from qed.eng_elec import D_to_Y

Z_A, Z_B, Z_C = D_to_Y(10, -10j, -10j)
print(f'Z_A = {Z_A} Ω')
print(f'Z_B = {Z_B} Ω')
print(f'Z_C = {Z_C} Ω')
