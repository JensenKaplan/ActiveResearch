import sys
sys.path.append('..')
from JensenTools import *

# PCOLig, Pr = cef.importCIF('Sr2PrO4.cif','Pr1',LS_Coupling = 64.5)
PCOLig, Pr = cef.importCIF('Li8PrO6.cif','Pr1')
Pr.diagonalize()
Pr.printEigenvectors()
# print(dir(Pr.Bdict))
print(Pr.BnmLabels)
energy = np.linspace(0,700,1000)
CalculatedSpectrum = Pr.neutronSpectrum(energy, Temp=4.76, Ei=700, ResFunc = lambda x: 3)
plt.plot(energy, CalculatedSpectrum, label='point charge model')
plt.legend()
plt.title("Sr2PrO4 With Ce central")
# plt.show()

# print(Pr.B)
