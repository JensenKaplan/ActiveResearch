import sys
sys.path.append('..')
from JensenTools import *

comp = 'Sr2PrO4'
who = 'Arun'
dataType = 'MT'
saveDir = getSaveDir('m', comp = comp)
MTDir = getSaveDir('m', comp = comp, dataType = dataType)
molweight = molweight[comp]
per = 'mol'
mass = .02701

#####################################################################################################################################################################
file = saveDir + 'Spectre/' + 'Sr2PrO4_MvsT_0p1T_zdir.dat'
headerList = ['Temperature', 'Xx', 'XyXz', 'Powder']
df = pd.read_csv(file, sep = '\s+', names = headerList)
# print(df.columns)
df = df.replace('D','E', regex=True)
# df = pd.to_numeric(df)
TSpec = pd.to_numeric(df['Temperature'])
Xx = pd.to_numeric(df['Xx'])
XyXz = pd.to_numeric(df['XyXz'])
XPowder = pd.to_numeric(df['Powder'])

XSpec = []
for i in range(len(Xx)):
    XSpec.append(XPowder[i])
	# temp = normalize(i,mass,molweight,per)
	# temp = bohrToEmu2(i)
	# MSpec.append(i)
XiSpec = 1/np.array(XSpec)
#####################################################################################################################################################################


# Get all runs from data directory
# Sort the runs by temperature (this makes for prettier plotting).
#####################################################################################################################################################################
runs = []
for i in os.listdir(MTDir):
    if i.endswith('.DAT') or i.endswith('.dat'):
        runs.append(i)       
MTData = {}
for i in runs:
    M,H,T,MErr,mass,measType = getData(i,MTDir, who = who, dataType = dataType)
    M = normalize(M,mass,molweight, per)
    Merr = normalize(M,mass,molweight, per)
    MTData[measType] = [M,H,T,MErr,mass]
#####################################################################################################################################################################

runtype = 'FC'
H = oeToTesla(MTData[runtype][1])
M = emuToBohr2(MTData[runtype][0])
H = MTData[runtype][1]
M = MTData[runtype][0]
X = M/H
Xi = 1/X


plt.plot(TSpec,XiSpec, label = 'Spectre')
plt.plot(T,Xi, label = 'Measured')

plt.title('{} X^-1'.format(comp))
plt.xlabel('Field (T)')
plt.ylabel('Moment (emu^-1 T mol)')
plt.legend()
plt.show()