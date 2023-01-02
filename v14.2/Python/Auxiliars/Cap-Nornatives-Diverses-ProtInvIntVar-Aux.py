import numpy as np

def Ia(t):
    """t en segons. Ia en ampere."""
    Id = 1.98*np.exp(-t/0.023) + 3.89*np.exp(-t/0.475) + 1.46 - 1.459*(1 - np.exp(-t/0.475))
    Iq = 3.159*np.exp(-t/0.023) + 0.781*np.exp(-t/0.106)
    return 4368*np.hypot(Id, Iq)

def t51(I):
    """I en ampere. t51 en segons. Cal definir el paràmetre T."""
    return 0.6*(0.014/((I/6000)**0.022 - 1) + 0.0226)


t = np.linspace(0, 3, 500)  # segons   
Iat = Ia(t)  # ampere
try:
    with open("Cap-Nornatives-Diverses-ProtInvIntVar-Aux-1.txt", "w", encoding="utf8") as fh:
        for i in range(len(t)):
            fh.write(f"{Iat[i]:11.5f}  {t[i]:6.4f} \n")
except OSError as err:
    print(err)

I = np.linspace(6100, 10_0000, 500)  # ampere
t51I = t51(I)  # segons
try:
    with open("Cap-Nornatives-Diverses-ProtInvIntVar-Aux-2.txt", "w", encoding="utf8") as fh:
        for i in range(len(I)):
            fh.write(f"{I[i]:11.5f}  {t51I[i]:7.4f} \n")
except OSError as err:
    print(err)
    