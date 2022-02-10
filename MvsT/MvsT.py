import sys
sys.path.append('..')
from JensenTools import *

<<<<<<< Updated upstream
#####################################################################################################################################################################
comp = 'Sr2PrO4'
who = 'MPMS'
# comp = 'Li8PrO6'
# who = 'MPMS'
dataType = 'MT'
saveDir = getSaveDir('m', comp = comp, dataType = dataType)
molweight = molweight[comp]
per = 'mol'

fit = False

# The S,L,J values are as follows for the Pr4+ ion
S = 0.5
L = 3
J = 5./2

massErr = .00005
#####################################################################################################################################################################
=======
comp = 'Ba2YbNbO6'
who = 'PPMS'
dataType = 'MT'
saveDir = getSaveDir('m', comp = comp, dataType = dataType)
molweight = molweight[comp]
J = 1/2
>>>>>>> Stashed changes

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

<<<<<<< Updated upstream
# Load all M vs H runs
# Normalize 
=======

massErr = .00005

>>>>>>> Stashed changes
#####################################################################################################################################################################
runs = []
for i in os.listdir(saveDir):
    if i.endswith('.DAT') or i.endswith('.dat'):
<<<<<<< Updated upstream
        runs.append(i)       
data = {}
=======
        runs.append(i)     
data = {}

for i in runs:
    M,H,T,MErr,samplemass,measType = getData(i,saveDir, who = who, dataType = dataType)
    data[measType] = [M,H,T,MErr,samplemass]

#####################################################################################################################################################################

>>>>>>> Stashed changes

for i in runs:
    M,H,T,MErr,mass,measType = getData(i,saveDir, who = who, dataType = dataType)
    M = normalize(M,mass,molweight, per)
    Merr = normalize(M,mass,molweight, per)
    data[measType] = [M,H,T,MErr,mass]


byField = {}
#Choose here
for i in data.keys():
    name = i
    fieldStr = i.split('_')[0]
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
    XT = X*T
    if fieldStr in byField.keys():
        byField[fieldStr].append([name,T,X,Xi,XT])
    else:
        byField[fieldStr] = [name,T,X,Xi,XT]
    XBohr = MBohr/HTes
    XiBohr = 1/XBohr
#####################################################################################################################################################################

if fit:
    tr = [0,300] #temprange = [low,high]
    newT = []
    newXi = []
    newErr = []
    for i in range(len(T)):
        if (T[i] >= tr[0] and T[i]<= tr[1]):
            newT.append(T[i])
            newXi.append(Xi[i])
            # newErr.append(XiErr[i])

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
for i in byField.keys():
    plt.plot(byField[i][1],byField[i][3],label = byField[i][0])
    # print(byField[i,2])
    plt.title("{} {}".format(comp, 'X^-1'), fontsize = 20)
    plt.xlabel('Temperature (K)', fontsize = 13)
    plt.ylabel('X^-1 (emu ^-1 Oe)', fontsize = 13)
    plt.legend()    
plt.show()

plt.figure()
for i in byField.keys():
    plt.plot(byField[i][1],byField[i][2],label = byField[i][0])
    plt.title("{} {}".format(comp, 'X'), fontsize = 20)
    plt.xlabel('Temperature (K)', fontsize = 13)
    plt.ylabel('X (emu Oe^-1)', fontsize = 13)
    plt.legend()   
plt.show()

plt.figure()
for i in byField.keys():
    plt.plot(byField[i][1],byField[i][4],label = byField[i][0])
    plt.title("{} {}".format(comp, 'X*T'), fontsize = 20)
    plt.xlabel('Temperature (K)', fontsize = 13)
    plt.ylabel('X*T (emu K Oe^-1)', fontsize = 13)
    plt.legend()   
plt.show()

# fig,ax = plt.subplots()

# ax.plot(T,Xi, color ='b',label = 'Measured X^-1')
# ax.set_xlabel('Temperature (K)', fontsize = 13)
# ax.set_ylabel('X^-1 (emu ^-1 Oe)', fontsize = 13, color = 'b')
# ax2 = ax.twinx()
# ax2.plot(T,X, color = 'orange', label = 'Measured X')
# ax2.set_ylabel('X (emu Oe^-1)', fontsize = 13, color = 'orange')

# ax2.plot(T, T*X, color = 'red', label = 'X*T')

# # plt.errorbar(T,Xi,yerr = XiErr,label = 'Measured 1/X')
# if fit:
#     plt.plot(T,fullLine,'orange', linestyle = '--', label = 'Fitted X^-1')
# ax.set_title("{} {}".format(comp,name), fontsize = 20)
# fig.legend()
# plt.show()

# plt.figure()
# plt.plot(T,Xi,label = 'Measured X^-1')
# plt.plot(T,X, label = 'Measured X')
# plt.plot(T, T*X, label = 'X*T')
# # plt.errorbar(T,Xi,yerr = XiErr,label = 'Measured 1/X')
# if fit:
#     plt.plot(T,fullLine,'orange', linestyle = '--', label = 'Fitted X^-1')
# plt.title("{} {}".format(comp,name), fontsize = 20)
# plt.xlabel('Temperature (K)', fontsize = 13)
# plt.ylabel('X^-1 (emu ^-1 Oe)', fontsize = 13)
# plt.legend()

# plt.figure()
# plt.plot(T,XiBohr,label = 'Measured X^-1')
# if fit:
#     plt.plot(T,fullLine,'orange', linestyle = '--', label = 'Fitted 1/X')
# plt.title("{} {}".format(comp,name), fontsize = 20)
# plt.xlabel('Temperature (K)', fontsize = 13)
# plt.ylabel('X^-1 (uB^-1 T)', fontsize = 13)
# plt.legend()
# if fit:
#     print('The Weiss constant = {:.2f} K\nThe Curie constant = {:.3f}'.format(resulti.params['wc'].value,resulti.params['c'].value))
#     ueff, gj = calcConstants(resulti.params['c'].value,J)

#     print('Effective moment for {:} is {:.3f} bohr magnetons, with J={} -> gj factor = {:.3f}'.format(comp,ueff,J,gj))
# plt.show()

#####################################################################################################################################################################