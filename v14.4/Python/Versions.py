import tabulate
import sys
import numba
import mkl
import numpy
import scipy
import sympy
import matplotlib
import pandas
import pandapower
import qed.utils
import qed.eng_elec

s = [['Llibreria', 'Versió'], [numba.__name__, numba.__version__],
     [mkl.__name__, mkl.__version__], [numpy.__name__,  numpy.__version__],
     [scipy.__name__, scipy.__version__], [sympy.__name__, sympy.__version__],
     [matplotlib.__name__, matplotlib.__version__],
     [pandas.__name__, pandas.__version__],
     [pandapower.__name__, pandapower.__version__],
     [qed.utils.__name__, qed.utils.__version__],
     [qed.eng_elec.__name__, qed.eng_elec.__version__]]

print('Python:', sys.version, '\n')
print(tabulate.tabulate(s, headers='firstrow'))
