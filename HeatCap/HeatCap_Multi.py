import sys
sys.path.append('..')
from JensenTools import *
from matplotlib import rcParams
from matplotlib import patches

fig = plt.figure(figsize = (10,10))
gs = fig.add_gridspec(1, 1)
ax1 = plt.subplot(gs[0])
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
rcParams.update({'font.size': 28})
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
rcParams['legend.fontsize'] = 18


comp = ['Sr2CeO4','Sr2PrO4']
who = 'PPMS'
# comp = 'Li8PrO6'
# who = 'MPMS'
dataType = 'HC'
saveDirDict = {}
for j in comp:
	saveDirDict[j] = getSaveDir('m', comp = j, dataType = dataType)

compounds = {}

for j in comp:
	runs = []
	for i in os.listdir(saveDirDict[j]):
	    if i.endswith('.DAT') or i.endswith('.dat'):
	        runs.append(i)

	temp = []
	for i in runs:
		# print()
		temp.append(i.split('_')[-1].split('.')[0][:-1]) # this creates a list of just temperatures as read by the filename   
	temp = np.argsort([float(i) for i in temp]) # Sort by temperature
	runs = [runs[i] for i in temp] # newly sorted listed

	data = {}
	for i in runs:
		T,h,hErr,field, mass = getData(i,saveDirDict[j], dataType = dataType, who = who)
		h = h * molweight[j]/mass
		data[field] = T,h,hErr,field

	compounds[j] = data


plt.figure()
for j in compounds.keys():
	for i in compounds[j].keys():
		T = compounds[j][i][0][:]
		h = compounds[j][i][1][:] # To Tesla
		if j == 'Sr2CeO4':
			h = h*35
		plt.plot(T,h, label = j + ' ' + i,linestyle = '--', marker = 'o',markersize = 5, linewidth = 3)
plt.legend(fontsize = '30')
plt.xlabel('Temperature (K)', fontweight = 'bold')
plt.ylabel('Heat Capacity (J/K) mol^-1', fontweight = 'bold')
plt.show()


# plt.figure()
# for i in data.keys():
# 	T = data[i][0][:]
# 	h = data[i][1][:]*10E-6 # To Tesla
# 	hErr = data[i][2][:]*10E-6 # To Tesla
# 	# if i == '0T':
# 		# print(type(T[0]))
# 	plt.errorbar(T,h, yerr = hErr, linestyle = '--', marker = 'o',markersize = 5, linewidth = 3, label = i)
# plt.legend(fontsize = '30')
# plt.xlabel('Temperature (K)', fontweight = 'bold')
# plt.ylabel('Heat Capacity (J/K)', fontweight = 'bold')
# plt.show()

# Sample Temp (Kelvin),Samp HC (µJ/K),Samp HC Err (µJ/K),

# data = {}
# for i in runs:
# 	pd.read_csv(saveDir+i)

