import sys
sys.path.append('..')
from JensenTools import *

# Important stuff
#####################################################################################################################################################################
comp = 'Ba2DyNbO6'
who = 'PPMS'
dataType = 'MT'
saveDir = getSaveDir('m', comp = comp, dataType = dataType)
MTDir = getSaveDir('m', comp = comp, dataType = dataType)
molweight = molweight[comp]
J = 7/2
massErr = .00005
fit = False
savepic = True
per = 'mol'
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
        # print(i)    
           
data = {}
for i in runs:
    M,H,T,MErr,mass,measType = getData(i,saveDir, who = who, dataType = dataType)
    M = normalize(M,mass,molweight, per)
    Merr = normalize(M,mass,molweight, per)
    data[measType] = [M,H,T,MErr,mass]

XiPlt = plt.figure()
XiAx = XiPlt.add_subplot(1,1,1)
XPlt = plt.figure()
XAx = XPlt.add_subplot(1,1,1)
XTPlt = plt.figure()
XTAx = XTPlt.add_subplot(1,1,1)
# Plot all
for i in data.keys():
    M,H,T,MErr,samplemass = data[i]
    MBohr = emuToBohr2(M)
    HTes = oeToTesla(H)
    MBohrErr = emuToBohr2(MErr)


    # X,XErr,Xi,XiErr = normSusc(M,H,MErr,molweight,samplemass,massErr)
    X = M/H
    Xi = 1/X
    XBohr = MBohr/HTes
    XiBohr = 1/XBohr
    XErr = MErr/H
    XBohrErr = MBohrErr/HTes
    XiBohrErr = POEXi(MBohr,MBohrErr,HTes,samplemass,massErr,comp,per)
    XiErr = POEXi(M, MErr, H, samplemass, massErr, comp, per)


    XiAx.errorbar(T, Xi, yerr = XiErr, label = 'Measured 1/X {}'.format(i), marker = '.', linestyle = 'none')
    XiAx.set_title("{}".format(comp), fontsize = 15)
    XiAx.set_xlabel('Temperature (K)', fontsize = 13)
    XiAx.set_ylabel('1/X (emu^-1 Oe {})'.format(per), fontsize = 13)
    XiAx.legend()

    XAx.errorbar(T, X, yerr = XiErr, label = 'Measured 1/X {}'.format(i), marker = '.', linestyle = 'none')
    XAx.set_title("{}".format(comp), fontsize = 15)
    XAx.set_xlabel('Temperature (K)', fontsize = 13)
    XAx.set_ylabel('X (emu Oe^-1 {})'.format(per), fontsize = 13)
    XAx.legend()

    XTAx.errorbar(T, X*T, yerr = XiErr, label = 'Measured 1/X {}'.format(i), marker = '.', linestyle = 'none')
    XTAx.set_title("{}".format(comp), fontsize = 15)
    XTAx.set_xlabel('Temperature (K)', fontsize = 13)
    XTAx.set_ylabel('X*T (emu K Oe^-1 {})'.format(per), fontsize = 13)
    XTAx.legend()


if savepic:
    XiPlt.savefig(MTDir+'{}_XivsT_emu_Oe.pdf'.format(comp))
    XPlt.savefig(MTDir+'{}_XvsT_emu_Oe.pdf'.format(comp))
    XTPlt.savefig(MTDir+'{}_XTvsT_emu_Oe.pdf'.format(comp))
plt.show()

XiPlt = plt.figure()
XiAx = XiPlt.add_subplot(1,1,1)
XPlt = plt.figure()
XAx = XPlt.add_subplot(1,1,1)
XTPlt = plt.figure()
XTAx = XTPlt.add_subplot(1,1,1)

for i in data.keys():
    M,H,T,MErr,samplemass = data[i]
    MBohr = emuToBohr2(M)
    HTes = oeToTesla(H)
    MBohrErr = emuToBohr2(MErr)
    XErr = MErr/H

    X = M/H
    Xi = 1/X
    XBohr = MBohr/HTes
    XiBohr = 1/XBohr
    XiBohrErr = POEXi(MBohr,MBohrErr,HTes,samplemass,massErr,comp,per)
    XiErr = POEXi(M, MErr, H, samplemass, massErr, comp, per)

    XiAx.errorbar(T, Xi, yerr = XiErr, label = 'Measured 1/X {}'.format(i), marker = '.', linestyle = 'none')
    XiAx.set_title("{}".format(comp), fontsize = 15)
    XiAx.set_xlabel('Temperature (K)', fontsize = 13)
    XiAx.set_ylabel('1/X (uB^-1 T {})'.format(per), fontsize = 13)
    XiAx.legend()

    XAx.errorbar(T, X, yerr = XiErr, label = 'Measured 1/X {}'.format(i), marker = '.', linestyle = 'none')
    XAx.set_title("{}".format(comp), fontsize = 15)
    XAx.set_xlabel('Temperature (K)', fontsize = 13)
    XAx.set_ylabel('X (uB T^-1 {})'.format(per), fontsize = 13)
    XAx.legend()

    XTAx.errorbar(T, X*T, yerr = XiErr, label = 'Measured 1/X {}'.format(i), marker = '.', linestyle = 'none')
    XTAx.set_title("{}".format(comp), fontsize = 15)
    XTAx.set_xlabel('Temperature (K)', fontsize = 13)
    XTAx.set_ylabel('X*T (uB K T^-1 {})'.format(per), fontsize = 13)
    XTAx.legend()

if savepic:
    XiPlt.savefig(MTDir+'{}_XivsT_uB_T.pdf'.format(comp))
    XPlt.savefig(MTDir+'{}_XvsT_uB_T.pdf'.format(comp))
    XTPlt.savefig(MTDir+'{}_XTvsT_uB_T.pdf'.format(comp))
plt.show()
#####################################################################################################################################################################

# FITTING
#####################################################################################################################################################################
if fit:
    tr = [200,300] #temprange = [low,high]
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
#####################################################################################################################################################################
    fullLine = []
    for i in T:
        fullLine.append(Curiei(i,resulti.params['c'].value,resulti.params['wc'].value))

#####################################################################################################################################################################
# plt.figure()
# plt.errorbar(T, Xi, yerr = XiErr, label = 'Measured 1/X', marker = '.', linestyle = 'none')
# # plt.errorbar(T,Xi,yerr = XiErr,label = 'Measured 1/X')

# if fit:
#     plt.plot(T,fullLine,'orange', linestyle = '--', label = 'Fitted 1/X')
#     plt.title("{} {} fitted over T = [{},{}]".format(comp,measType,tr[0],tr[1]), fontsize = 15)
# else:
#     plt.title("{} {}".format(comp,measType), fontsize = 15) 

# plt.xlabel('Temperature (K)', fontsize = 13)
# plt.ylabel('1/X (emu ^-1 Oe {})'.format(per), fontsize = 13)
# plt.legend()

# plt.figure()
# plt.errorbar(T,XiBohr, yerr = XiBohrErr, label = 'Measured 1/X', marker = '.', linestyle = 'none')
# # plt.plot(T,fullLine,'orange', linestyle = '--', label = 'Fitted 1/X')
# plt.title("{} {}".format(comp,measType), fontsize = 15)
# plt.xlabel('Temperature (K)', fontsize = 13)
# plt.ylabel('1/X (uB^-1 T {})'.format(per), fontsize = 13)
# plt.legend()

# if fit:
#     print('The Weiss constant = {:.2f} K\nThe Curie constant = {:.3f}'.format(resulti.params['wc'].value,resulti.params['c'].value))
#     ueff, gj = calcConstants(resulti.params['c'].value,J)
#     print('Effective moment for {:} is {:.3f} bohr magnetons, with J={} -> gj factor = {:.3f}'.format(comp,ueff,J,gj))

plt.show()

#####################################################################################################################################################################