#!/usr/bin/env python
# coding: utf-8

# In[72]:


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


# # Load the M vs H data files
# ## Manually name the sample in the variable "comp"
# ## Manually enter the molecular weight (g/mol)
# ## Load the sample mass (grams) as read in by one of the data files' name
# 
# Molecular weight calculated by https://www.lenntech.com/calculators/molecular/molecular-weight-calculator.htm

# In[74]:


runs = [] #A list of all the data file names
for i in os.listdir(MHDir):
    if i.endswith('.DAT'): #This was a safeguard against a situation arising at an earlier implementation of my code.
        runs.append(i)

molweight = 380.15
mass = getMass(runs[0], who = who)
T = getTemp(runs[0], who = who)


# # Order the files by temperature
# This guarantees nicer plotting

# In[75]:


temp = [] #The newly sorted list
for i in runs:
    temp.append(getTemp(i, who = who)) #this creates a list of just temperatures as read by the filename
   
temp = np.argsort([int(i) for i in temp]) #Sort by temperature
runs = [runs[i] for i in temp]



# This dictionary easily allows me to access the data from a specific run
## {'Temperature': [H,M,Err]} 
MHdata = {}
plt.figure()
for i in runs:
    M, H, Err, mass, T = getData(i,MHDir,who = who, dataType = 'MH')
    M = emuToBohr(M,mass,molweight)
    H = oeToTesla(H)
    Err = emuToBohr(Err,mass,molweight)
    MHdata[T] = [M,H,Err,mass,i]
    plt.errorbar(H,M, yerr = Err, label = T)

# MHData = {}
# for i in runs:
#     M, H, Err, mass, T = getData(i,MHDir,who = who, dataType = 'MH')
#     MHdata[T] = [M,H,Err,mass,i]

plt.title(comp)
plt.ylabel('Moment (\N{GREEK SMALL LETTER MU}B)')
plt.xlabel('Field (T)')
plt.legend()
plt.show()


temp = '20K'
curRun = MHdata[temp] #loading the data from my current run
H = curRun[0]
M = curRun[1]
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

