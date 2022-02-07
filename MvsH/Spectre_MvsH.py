import sys
sys.path.append('..')
from JensenTools import *

comp = 'Ba2YbNbO6'
who = 'PPMS'
dataType = 'MH'
saveDir = getSaveDir('m', comp = comp)
MHDir = getSaveDir('m', comp = comp, dataType = dataType)
molweight = molweight[comp]
per = 'spin'
mass = .038

#####################################################################################################################################################################
file = saveDir + 'Spectre/' + 'Ba2YbNbO6_MvsH_1P8K.dat'
headerList = ['Field', 'Mx', 'My', 'Mz']
df = pd.read_csv(file, sep = '\s+', names = headerList)
# print(df.columns)
df = df.replace('D','E', regex=True)
# df = pd.to_numeric(df)
HSpec = pd.to_numeric(df['Field'])
Mz = pd.to_numeric(df['Mz'])

MSpec = []
for i in Mz:
	# temp = normalize(i,mass,molweight,per)
	# temp = bohrToEmu2(i)
	MSpec.append(i)
#####################################################################################################################################################################


# Get all runs from data directory
# Sort the runs by temperature (this makes for prettier plotting).
#####################################################################################################################################################################
runs = [] #A list of all the data file names
for i in os.listdir(MHDir):
    if i.endswith('.DAT') or i.endswith('.dat'): #This was a safeguard against a situation arising at an earlier implementation of my code.
        runs.append(i)

temp = [] # temporary list for sorting
for i in runs:
    temp.append(getTemp(i, who = who)) # this creates a list of just temperatures as read by the filename   
temp = np.argsort([int(i) for i in temp]) # Sort by temperature
runs = [runs[i] for i in temp] # newly sorted listed

MHData = {}
# Normalizes and stores data as well as plotting in Emu/Oe for all temperatures.
plt.figure()
for i in runs:
    M, H, Err, mass, T = getData(i,MHDir,who = who, dataType = 'MH')
    M = normalize(M,mass,molweight,per)
    Err = normalize(Err,mass,molweight, per)
    MHData[T] = [M,H,Err,mass,i]
    # plt.plot(H, M, label = T)
# plt.title('{} Magnetization'.format(comp))
# plt.ylabel('Moment (emu {}^-1)'.format(per))
# plt.xlabel('Field (Oe)')
# plt.legend()

#####################################################################################################################################################################

print(MHData.keys())
temp = 1.8
H = oeToTesla(MHData[temp][1])
M = emuToBohr2(MHData[temp][0])

plt.plot(HSpec,MSpec, label = 'Spectre')
plt.plot(H,M, label = 'Measured')

plt.title('{} Spectre Magnetization at 20K'.format(comp))
plt.xlabel('Field (T)')
plt.ylabel('Moment (uB spin^-1)')
plt.legend()
plt.show()