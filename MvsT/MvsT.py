import sys
sys.path.append('..')
from JensenTools import *

#####################################################################################################################################################################
# comp = 'Sr2PrO4'
# who = 'Arun'
comp = 'Li8PrO6'
who = 'MPMS'
dataType = 'MT'
saveDir = getSaveDir('m', comp = comp, dataType = dataType)
molweight = molweight[comp]
per = 'spin'

# The S,L,J values are as follows for the Pr4+ ion
S = 0.5
L = 3
J = 5./2

massErr = .00005
#####################################################################################################################################################################

# Inverse Curie Law for LMFIT
# Constant calculation function
#####################################################################################################################################################################
# This is the inverse of the Curie-Weiss Law.
# We will be creating an LMFIT model with this and fitting for:
## Curie's constant 'c', Weiss' constant 'wc'
def Curiei(t,c,wc):
    return(t-wc)/c

# Takes in Curie's constant and the system's total angular momentum J value
# Returns effective moment (bohr magnetons) and effective g-factor
def calcConstants(c,J):
    ueff = np.sqrt(8*c)
    gj = ueff/np.sqrt(J*(J+1))
    return ueff, gj
#####################################################################################################################################################################

# Load all M vs H runs
# Normalize 
#####################################################################################################################################################################
runs = []
for i in os.listdir(saveDir):
    if i.endswith('.DAT') or i.endswith('.dat'):
        runs.append(i)       
data = {}

for i in runs:
    M,H,T,MErr,mass,measType = getData(i,saveDir, who = who, dataType = dataType)
    M = normalize(M,mass,molweight, per)
    Merr = normalize(M,mass,molweight, per)
    data[measType] = [M,H,T,MErr,mass]

#Choose here
name = '0.1T_ZFC'
M,H,T,MErr,samplemass = data[name]
MBohr = emuToBohr2(M)
HTes = oeToTesla(H)
#####################################################################################################################################################################

# Calculate susceptibility (Emu/Oe).
# Choose a temp range to fit over (K)
#####################################################################################################################################################################
# X,XErr,Xi,XiErr = normSusc(M,H,MErr,molweight,samplemass,massErr)
X = M/H
Xi = 1 / X
XBohr = MBohr/HTes
XiBohr = 1/XBohr


tr = [0,300] #temprange = [low,high]
newT = []
newXi = []
newErr = []
for i in range(len(T)):
    if (T[i] >= tr[0] and T[i]<= tr[1]):
        newT.append(T[i])
        newXi.append(Xi[i])
        # newErr.append(XiErr[i])
#####################################################################################################################################################################

cmodeli =  Model(Curiei, independent_vars = ['t'])
params = cmodeli.make_params()
params['wc'].set(value = -10)
params['c'].set(value = 10)   

resulti = cmodeli.fit(newXi, params, t = newT) #fit
# resulti = cmodeli.fit(newXi, params, t = newT, weights = newErr) #fit

fullLine = []
for i in T:
    fullLine.append(Curiei(i,resulti.params['c'].value,resulti.params['wc'].value))

#####################################################################################################################################################################
plt.figure()
plt.plot(T,Xi,label = 'Measured 1/X')
# plt.errorbar(T,Xi,yerr = XiErr,label = 'Measured 1/X')
plt.plot(T,fullLine,'orange', linestyle = '--', label = 'Fitted X^-1')
plt.title("{} {} fitted over T = [{},{}]".format(comp,name,tr[0],tr[1]), fontsize = 20)
plt.xlabel('Temperature (K)', fontsize = 13)
plt.ylabel('X^-1 (emu ^-1 Oe)', fontsize = 13)
plt.legend()

plt.figure()
plt.plot(T,XiBohr,label = 'Measured X^-1')
# plt.plot(T,fullLine,'orange', linestyle = '--', label = 'Fitted 1/X')
plt.title("{} {}".format(comp,name), fontsize = 20)
plt.xlabel('Temperature (K)', fontsize = 13)
plt.ylabel('X^-1 (uB^-1 T)', fontsize = 13)
plt.legend()
print('The Weiss constant = {:.2f} K\nThe Curie constant = {:.3f}'.format(resulti.params['wc'].value,resulti.params['c'].value))

ueff, gj = calcConstants(resulti.params['c'].value,J)

print('Effective moment for {:} is {:.3f} bohr magnetons, with J={} -> gj factor = {:.3f}'.format(comp,ueff,J,gj))
plt.show()

#####################################################################################################################################################################