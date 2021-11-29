import numpy as np
import matplotlib.pyplot as plt
import os
from lmfit.models import LinearModel
import pandas as pd
import sys
sys.path.append('..')
from JensenTools import *

who = 'Arun'
comp = 'Sr2PrO4'
MHDir = getSaveDir('m', comp = comp, dataType = 'MH')
molweight = 380.15

runs = [] #A list of all the data file names
for i in os.listdir(MHDir):
    if i.endswith('.DAT'): #This was a safeguard against a situation arising at an earlier implementation of my code.
        runs.append(i)


mass = getMass(runs[0], who = who)
T = getTemp(runs[0], who = who)


temp = [] #The newly sorted list
for i in runs:
    temp.append(getTemp(i, who = who)) #this creates a list of just temperatures as read by the filename
   
temp = np.argsort([int(i) for i in temp]) #Sort by temperature
runs = [runs[i] for i in temp]


MHdata = {}
plt.figure()
for i in runs:
    M, H, Err, mass, T = getData(i,MHDir,who = who, dataType = 'MH')
    # M = emuToBohr(M,mass,molweight)
    # H = oeToTesla(H)
    # Err = emuToBohr(Err,mass,molweight)
    M = normalize(M,mass,molweight,'spin')
    Err = normalize(Err,mass,molweight,'spin')
    MHdata[T] = [M,H,Err,mass,i]
    plt.errorbar(H, M, yerr = Err, label = T)
plt.title(comp)
plt.ylabel('Moment (emu) spin^-1')
plt.xlabel('Field (Oe)')
plt.legend()
# plt.show()

plt.figure()
for i in MHdata.keys():
    M,H,Err,Mass,T = MHdata[i]
    # M= normalize(M,mass,molweight,'spin')
    M = emuToBohr2(M)
    H = oeToTesla(H)
    Err = emuToBohr2(Err)
    plt.errorbar(H,M, yerr = Err,label = i)

plt.title(comp)
plt.ylabel('Moment (\N{GREEK SMALL LETTER MU}B) spin^-1')
plt.xlabel('Field (T)')
plt.legend()
plt.show()



temp = '20.0K'
curRun = MHdata[temp] #loading the data from my current run

M = curRun[0]
H = curRun[1]
Err = curRun[2]



fieldRange = [0,14]
newH = []
newM = []
newErr = []
for i in range(len(curRun[0])):
    if (H[i] >= fieldRange[0] and H[i] <= fieldRange[1]):
        newH.append(H[i])
        newM.append(M[i])
        newErr.append(Err[i])



linModel = LinearModel()
params = linModel.guess(newM, x = newH)
fitted = linModel.fit(newM, x = newH, weights = newErr)

MLine = []
for i in H:
    MLine.append(fitted.params['slope'].value*i + fitted.params['intercept'].value)



plt.plot(H,M, label = temp)
plt.plot(H,MLine, linestyle = '--', label = 'Fitted')
plt.xlabel('Field (T)')
plt.ylabel('Magnetization (bohr magneton)')
plt.legend()
plt.title(comp)
plt.show()

print('Saturation magnetization =  {:.3f} bohr magneton'.format(fitted.params['intercept'].value))

