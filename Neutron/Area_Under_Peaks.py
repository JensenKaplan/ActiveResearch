import sys
sys.path.append('..')
from JensenTools import *
from scipy.integrate import simps

# Basics
#####################################################################################################################################################################
comp = 'Sr2PrO4'
dataType = 'IE'
saveDir = getSaveDir('m', comp = comp, dataType = dataType)
#####################################################################################################################################################################

# List all Ei runs and select one to analyze
#####################################################################################################################################################################
dirList = os.listdir(saveDir)
runList = []
for i in dirList:
    if '.ascii' in i:
        runList.append(i)
print("All of the 1D cuts in the folder will be listed with their index. So run this once and fill in magrun.")
for i in range(len(runList)):
    print('{} {}'.format(runList[i], i))
#####################################################################################################################################################################

# Choose a run
magrun = runList[4]
nonmagrun = runList[6]
# Loading data from chosen file
#####################################################################################################################################################################
f = open(saveDir + magrun, 'r')
datamag = np.genfromtxt(f, comments = '#')
f.close()

E = datamag[:,0]
I = datamag[:,1]
IErr = datamag[:,2]

f = open(saveDir + nonmagrun, 'r')
datanonmag = np.genfromtxt(f, comments = '#')
f.close()

Enonmag = datanonmag[:,0]
Inonmag = datanonmag[:,1]
IErrnonmag = datanonmag[:,2]

sub = np.array(I) - .6*np.array(Inonmag)
#####################################################################################################################################################################

peak1Range = [306,361]
peak2Range = [363,410]
peak1New = []
peak2New = []

for i in range(len(E)):
    if (E[i] >= peak1Range[0] and E[i]<= peak1Range[1]):
        peak1New.append(E[i])
    if (E[i] >= peak2Range[0] and E[i]<= peak2Range[1]):
        peak2New.append(E[i])

area1 = simps(peak1New, dx=1)
print("Area under peak 1 =", area1)
area2 = simps(peak2New, dx=1)
print("Area under peak 1=", area2)


plt.figure()
plt.plot(E,I, label = 'Sr2PrO4', linestyle = '--')
plt.plot(Enonmag,Inonmag, label = 'Sr2CeO4', linestyle = '--')
plt.plot(E,sub, label = 'Subtracted')
plt.ylabel('Intensity (a.u.)')
plt.xlabel('Energy (meV)')
plt.legend()
plt.show()