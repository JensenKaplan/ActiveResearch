import sys
sys.path.append('..')
from JensenTools import *


L = 3
S = .5

LS = 100.5
x =  0.03629536921151444
bpf = -0.6570713391739674
# Assigning the coefficients from grid search
# Enforcing cubic constraints as a start
# and including the B20 term which is needed for tetragonal symmetry
B40 = bpf
B60 = x*bpf
B44 = 5*B40
B64 = -21*B60
B20 = 0

stev = {'B40': B40, 'B60': B60, 'B44' : B44, 'B64' : B64, 'B20' : B20 }


Pr = cef.LS_CFLevels.Bdict(Bdict = stev, L = L, S = S, SpinOrbitCoupling=LS)
Pr.diagonalize()

print(Pr.eigenvalues)

Emeas = [0,168, 335, 385] # The measured INS magnetic modes
numlevels = 5

e = kmeansSort(Pr.eigenvalues,numlevels)[:numlevels-1] #Excluding the highest mode which we did not detect in our INS runs
print(e)