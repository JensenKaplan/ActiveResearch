import sys
sys.path.append('..')
from JensenTools import *
from scipy.integrate import simps
from matplotlib import rcParams
from matplotlib import patches

# Plot Formatting
#####################################################################################################################################################################
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
rcParams.update({'font.size': 15})
rcParams['font.weight'] = 'bold'
rcParams['axes.linewidth'] = 4
rcParams['xtick.direction'] = 'out'
rcParams['ytick.direction'] = 'out'
rcParams['xtick.top'] = False
rcParams['ytick.right'] = False
rcParams['xtick.major.size'] = 12.5
rcParams['ytick.major.size'] = 12.5
rcParams['xtick.minor.size'] = 7.5
rcParams['ytick.minor.size'] = 7.5
rcParams['xtick.major.width'] = 3
rcParams['ytick.major.width'] = 3
rcParams['xtick.minor.width'] = 3
rcParams['ytick.minor.width'] = 3
rcParams['xtick.minor.visible'] = True
rcParams['ytick.minor.visible'] = True
rcParams['legend.frameon'] = True
rcParams['legend.fontsize'] = 12
#####################################################################################################################################################################

# Define important things
#####################################################################################################################################################################
comp = 'Sr2PrO4'
ion = 'Ce3+'
who = 'MPMS'
LS_on = True
per = 'spin'
saveplots = False

molweight = molweight[comp]

# The S,L,J values are as follows for the Pr4+ ion
S = 0.5
L = 3
J = 5./2
#####################################################################################################################################################################


# SrPr PAPER BEST FIT
#####################################################################################################################################################################
pf  =  0.24697226666068609
B20  =  10.079050365865657
B40  =  -0.07424576811339292
B60  =  0.040280418062079895
B44  =  -0.9997183336100587
B64  =  0.07450220927671151
LS  =  55.63834774730741
#####################################################################################################################################################################

# Create Stevens Bdict, Create CFLevels object, diagonalize.
#####################################################################################################################################################################
stev = {'B20' : B20, 'B40' : B40, 'B60' : B60, 'B44' : B44, 'B64' : B64}

# Create the CFLevels object and diagonalize it
if LS_on:
	Pr = cef.LS_CFLevels.Bdict(Bdict = stev, L = L, S = S, SpinOrbitCoupling= LS)
	Pr.diagonalize()
else:
	Pr = cef.CFLevels.Bdict(Bdict = stev, ion = ion)
	Pr.diagonalize()
#####################################################################################################################################################################

# Tracking eigenvalues vs field
#####################################################################################################################################################################
eByField = {}
for i in range(14):
	eByField['E{}'.format(i)] = []

print(eByField.keys())

fields = np.linspace(0,10000,1000)
powder = True

for i in fields:
	
	if powder:
		testField = i
	else:
		testField = [0,0,i]

	k = 0
	for j in eByField.keys():
		eByField[j].append(Pr.zeeman(Field = testField)[k])
		k+=1

plt.figure()
for i in eByField.keys():
	# if i == 'E0' or i == 'E1':
	plt.plot(fields,eByField[i], label = i, linewidth = 2)

plt.xlabel('Applied Field (T)')
plt.ylabel('Energy (meV)')
plt.legend(bbox_to_anchor=(1.1, 1))
if powder:
	plt.title(comp + ' Powder Averaging')
else:
	plt.title(comp + ' Applied Field in Z Direction')
plt.show()


# print("\n\nEigenvalues before applying magnetic field:\n{}.".format(Pr.eigenvalues))

# print("Eigenvalues after applying a vector {} T magnetic field:\n{}.".format(testField,Pr.zeeman(Field = testField)))