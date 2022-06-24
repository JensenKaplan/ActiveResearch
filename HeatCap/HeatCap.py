import sys
sys.path.append('..')
from JensenTools import *
from matplotlib import rcParams
from matplotlib import patches

rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
rcParams.update({'font.size': 13})
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
rcParams['legend.frameon'] = False
rcParams['legend.fontsize'] = 12


comp = 'Ba2YbNbO6'
who = 'PPMS'
# comp = 'Li8PrO6'
# who = 'MPMS'
dataType = 'HC'
saveDir = getSaveDir('m', comp = comp, dataType = dataType)
DataViewer = False
avo = 6.0221409e+23 #spin/mol

runs = []
for i in os.listdir(saveDir):
    if i.endswith('.DAT') or i.endswith('.dat'):
        runs.append(i)

temp = []
for i in runs:
	# print()
	temp.append((i.split('_')[-1]).split('.')[0][:-1]) # this creates a list of just temperatures as read by the filename   
temp = np.argsort([float(i) for i in temp]) # Sort by temperature
runs = [runs[i] for i in temp] # newly sorted listed
# print(runs)

data = {}
for i in runs:
	T,C,CErr,field, mass = getData(i,saveDir, dataType = dataType, who = who)
	
	if mass == -1:
		C = C
	else:
		pass
		C = C*molweight[comp]/mass/(10e6)
		CErr = CErr*molweight[comp]/mass/(10e6)

	# PPMS when doing HC can either go from lowT -> highT, or vice verse
	# it also performs multiple iterations of the same temp step.
	# This leads to possible messiness when plotting. Points connected that shouldn't be.
	# The below sorts the dsata appropriately and hence fixes plots.
	zippedHC = zip(T,C)
	sortedHC = sorted(zippedHC)
	HCtuples = zip(*sortedHC)
	print(HCtuples)
	T, C = [list(tuple) for tuple in  HCtuples]

	data[field] = np.array(T),np.array(C),np.array(CErr)
	# data[field] = T,h,hErr

CPlt = plt.figure()
CAx = CPlt.add_subplot(1,1,1)
CTPlt = plt.figure()
CTAx = CTPlt.add_subplot(1,1,1)
for i in data.keys():
	if i == 'HC1T' or i == 'HC0T' or i == 'DRHC0T':
		pass

	else:
		T = data[i][0]
		C = data[i][1]# To Tesla
		CErr = data[i][2] # To Tesla
		CTErr = CErr/T
		CAx.errorbar(T,C, yerr = CErr, linestyle = '--', marker = 'o',markersize = 5, linewidth = 3, label = i)
		CTAx.errorbar(T, C/T, yerr = CTErr, linestyle = '--', marker = 'o',markersize = 5, linewidth = 3, label = i)

CAx.legend(fontsize = '13')
CAx.set_title(comp)
CAx.set_xlabel('Temperature (K)', fontweight = 'bold')
CAx.set_ylabel('Heat Capacity (J/K/mol)', fontweight = 'bold')
CTAx.legend(fontsize = '13')
CTAx.set_title(comp)
CTAx.set_xlabel('Temperature (K)', fontweight = 'bold')
CTAx.set_ylabel('C/T (J/K^2 mol^-1)', fontweight = 'bold')
plt.show()

# Sample Temp (Kelvin),Samp HC (µJ/K),Samp HC Err (µJ/K),

# data = {}
# for i in runs:
# 	pd.read_csv(saveDir+i)

