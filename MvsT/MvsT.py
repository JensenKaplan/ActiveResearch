import numpy as np
import matplotlib.pyplot as plt
import os
from lmfit import Model
import pandas as pd
import sys
sys.path.append('..')
from JensenTools import *

comp = 'Sr2PrO4'
who = 'Arun'
dataType = 'MT'
saveDir = getSaveDir('m', comp = comp, dataType = dataType)
J = 1/2


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

molweight = 380.15
massE = .00005

#####################################################################################################################################################################
runs = []
for i in os.listdir(saveDir):
    if i.endswith('.DAT'):
        runs.append(i)       
data = {}
for i in runs:
    M,H,T,Err,samplemass,measType = getData(i,saveDir, who = who, dataType = dataType)
    data[measType] = [M,H,T,Err,samplemass]
#####################################################################################################################################################################


M,H,T,Err,samplemass = data['ZFC']

#####################################################################################################################################################################
SuscNorm = []
SuscErrorNorm = []
for i in range(len(M)):
    SuscNorm.append(M[i]/H[i]*molweight/samplemass/1)
    SuscErrorNorm.append(Err[i]/H[i]*molweight/samplemass + M[i]/H[i]*molweight*massE/(samplemass**2))

Xi = 1/np.array(SuscNorm)
XiErr = np.array(SuscErrorNorm)/np.array(SuscNorm)**2

tr = [1,300] #temprange = [low,high]
newT = []
newXi = []
newErr = []
for i in range(len(T)):
    if (T[i] >= tr[0] and T[i]<= tr[1]):
        newT.append(T[i])
        newXi.append(Xi[i])
        newE.append(XiErr[i])
#####################################################################################################################################################################

cmodeli =  Model(Curiei, independent_vars = ['t'])
params = cmodeli.make_params()
params['wc'].set(value = 10)
params['c'].set(value = 10)   
resulti = cmodeli.fit(newXi, params, t = newT, weights = newErr) #fit

fullLine = []
for i in T:
    fullLine.append(Curiei(i,resulti.params['c'].value,resulti.params['wc'].value))

#####################################################################################################################################################################
plt.errorbar(T,Xi,yerr = XiErr,label = 'Measured 1/X')
plt.plot(T,fullLine,'orange', linestyle = '--', label = 'Fitted 1/X')
plt.title("{} {} fitted over T = [{},{}]".format(comp,measType,tr[0],tr[1]), fontsize = 20)
plt.xlabel('Temperature (K)', fontsize = 13)
plt.ylabel('1/X (spin/emu)', fontsize = 13)
plt.legend()
plt.show()

print('The Weiss constant = {:.2f} K\nThe Curie constant = {:.3f}'.format(resulti.params['wc'].value,resulti.params['c'].value))

ueff, gj = calcConstants(resulti.params['c'].value,J)

print('Effective moment for {:} is {:.3f} bohr magnetons, with J={} -> gj factor = {:.3f}'.format(comp,ueff,J,gj))
#####################################################################################################################################################################