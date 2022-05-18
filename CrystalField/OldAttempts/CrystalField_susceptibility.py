#!/usr/bin/env python
# coding: utf-8

# # Crystal Electric Field (CF) Analysis of the Pr4+ ion
# Here the goal is to use collected inelastic neutron scattering (INS) data to determine the ground state crystal field Hamiltonian. The Stevens' coefficients first guess comes from the grid search. I use PyCrystalField in conjunction with LMFIT to find the true values of the Stevens' Coefficients by fitting predicted energy levels to measured one. I then calculate the CF Hamiltonian and predict the compounds thermodynamic properties. Magnetization (M vs H) and Susceptibility (M vs T) data. 
# 
# It's important to note that I don't actually use Pr4+ as my central ion, I use Ce3+. The correct orbital values for the Ln4+ oxidized states have not been calculated and tabulated. Since the operators depend on electronic orbital states, I can use an electrically equivalent ion, Ce3+, instead of the actual Pr4+ ion.

# In[15]:

import sys
sys.path.append('..')
from JensenTools import *
from lmfit import Model
import PyCrystalField as cef
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn
seaborn.set()

#Simple function for pritting my parameters after fitting so that I can copy paste the values from output for further iterations
def paramPrint(fittedparams):
	print('\nEasy copy pasting:')
	for i in fittedparams:
		# print(i, ' = ', i.value)
		print(i, ' = ',fittedparams[i].value )


# ##  PCF Energy Calculation Function To Be Used As An LMFIT Model

# In[16]:


def energyFit(B20, B40, B44, B60,B64, LS, numlevels, **kwargs ):
	Stev = {} #Creating the Stevens' Coefficients dictionary and assigning values
	Stev['B20'] = B20
	Stev['B40'] = B40
	Stev['B60'] = B60
	Stev['B44'] = B44
	Stev['B64'] = B64

	if kwargs['LS_on']: #If LS Basis
		Pr = cef.LS_CFLevels.Bdict(Bdict=Stev, L=3, S=0.5, SpinOrbitCoupling = LS) #Create CF_Levels obejct wtih the given coefficients.
		Pr.diagonalize()
		if kwargs['Kmeans']: #If using Kmeans sorting
			e = kmeansSort(Pr.eigenvalues,numlevels)[0:3] #Excluding the highest mode which we did not detect in our INS runs
			# e.append(e[2]/e[1]) #A ratio of energy levels, to allow for more features to vary.
		else: #If not Kmeans then storing all eigenvalues
			e = Pr.eigenvalues
            
	else: #If J Basis
		Pr = cef.CFLevels.Bdict(Bdict = Stev, ion = kwargs['ion'])
		Pr.diagonalize()
		if kwargs['Kmeans']: #If Kmeans sorting
			e = kmeansSort(Pr.eigenvalues,numlevels) #Excluding the highest mode which we did not detect in our INS runs
		else: #If not Kmeans then storing all eigenvalues
			e =  Pr.eigenvalues
	return e


# ## Define kwargs and other important things

# In[17]:


#kwargs
#####################################################################################################################################################################
comp = 'Sr2PrO4'
ion = 'Ce3+'
who = 'Arun'
LS_on = True
Kmeans = True
LSValue = 100
numlevels = 4
#####################################################################################################################################################################

# The L,S values are as follows for the Pr4+ ion
L = 3
S = 0.5
#The measured INS magnetic modes. The last entry being a ratio.
Emeas = [168, 335,385] 


# ## Load the project's data directories from my DropBox

# In[18]:


saveDir = getSaveDir('m',comp = comp) #General Directory for the project
MHDir = getSaveDir('m',comp = comp, dataType = 'MH') #MvsH data
MTDir = getSaveDir('m',comp = comp, dataType = 'MT') #MvsT data


# ## LMFIT Model Creation and Fitting

# In[19]:


eModel = Model(energyFit, independent_vars = ['numlevels'])
params = eModel.make_params()

#From GridSearch For LS
#####################################################################################################################################################################
# LS = LSValue
# x =  0.03629536921151444
# bpf = -0.6570713391739674
# # Assigning the coefficients from grid search
# # Enforcing cubic constraints as a start
# # and including the B20 term which is needed for tetragonal symmetry
# B40 = bpf
# B60 = x*bpf
# B44 = 5*B40
# B64 = -21*B60
# B20 = 0
#####################################################################################################################################################################

#Best fit for LS after a few iterations.
#####################################################################################################################################################################
B40  =  -0.6568663783690575
B60  =  -0.02328250024945387
LS  =  100.00007580463522
B44  =  -3.1415463304732714
B64  =  0.504906552605772
B20  =  0.4858075931009187
#####################################################################################################################################################################

# Since we only have 3 training points, only 3 parameters can vary at once.
params['B20'].set(value = B20, vary = False)
params['B40'].set(value=B40, vary=True)
params['B60'].set(value=B60, vary=True)
params['B44'].set(value = B44, vary = True )
params['B64'].set(value = B64, vary = False )
params['LS'].set(value=LS, vary=False)


#Fit model to data
fitted = eModel.fit(Emeas,params, numlevels = 4, LS_on = LS_on, Kmeans = Kmeans)

#Create a dictionary of the fitted parameters (stevens coefficients)
stev = { 'B20' :fitted.params['B20'].value, 'B40': fitted.params['B40'].value, 'B44' : fitted.params['B44'].value, 'B60': fitted.params['B60'].value, 'B64' : fitted.params['B64'].value }
#####################################################################################################################################################################


# ## Printing fitted parameters and reduced chi-squared.

# In[20]:


#Print the parameters and reduced chi sqr value
print('\n\nFitted parameters:')
fitted.params.pretty_print()
print('\nReduced Chi Sqr = {}\n'.format(fitted.result.redchi))
#Uncomment to print out in easy copy paste format
paramPrint(fitted.params)


# ## Use fitted B parameters to obtain the CF Hamiltonian

# In[21]:


#CF Analysis
#####################################################################################################################################################################
#Create the CFLevels object and diagonalize it
if LS_on:
	Pr = cef.LS_CFLevels.Bdict(Bdict = stev, L = L, S = S, SpinOrbitCoupling=fitted.params['LS'].value)
	Pr.diagonalize()
else:
	Pr = cef.CFLevels.Bdict(Bdict = stev, ion = 'Ce3+')
	Pr.diagonalize()
	

#Print final matrix
print('\n\nEnergy values as measured by INS (meV): {}'.format(Emeas[:-1]))
Pr.printEigenvectors()
#####################################################################################################################################################################


# ## PCF G Tensor

# In[22]:


#Calculate and neatly print G-Tensor using Pandas
gt = Pr.gtensor()
rows = ['gx','gy','gz']
df = pd.DataFrame(gt, columns = rows, index = rows)
print(df)


# ## Magnetization (M vs H)

# In[23]:


# ## PCF Magnetization
#####################################################################################################################################################################
runs = []
for i in os.listdir(MHDir):
    if i.endswith('.DAT'): #This was a safeguard against a situation arising at an earlier implementation of my code.
        runs.append(i)

mass = getMass(runs[0],who = who)
molweight = 380.15
MHdata = {}
# plt.figure()
for i in runs:
    H, M, Err, T = getData(i,MHDir,who = who, dataType = 'MH')
    # M = emuToBohr(M,mass,molweight)
    # H = oeToTesla(H)
    # Err = emuToBohr(Err,mass,molweight)
    MHdata[T] = [H,M,Err,i]
    # plt.errorbar(H,M, yerr = Err, label = name)

# print(MHdata.keys())
T = '20K'

# samplemass = getMass(MHdata[T][3])
Temp = getTemp(MHdata[T][3], who = who)
#####################################################################################################################################################################
H, M, Err, filename = MHdata[T]
MNorm = normalizeSpin(M,mass,molweight)
HTes = oeToTesla(H)
MBohr = emuToBohr(M,mass,molweight)
ErrBohr = emuToBohr(Err,mass,molweight)

#Generate a magnetization curve for comparing results to experiment
magCalcBohr = []
# fieldT = np.linspace(0.01,14,1000)

for i in HTes:
	if LS_on:
		magCalcBohr.append(Pr.magnetization( Temp = Temp, Field = [0, 0, i])[2])		
	else:
		magCalcBohr.append(Pr.magnetization( Temp = Temp, Field = [0, 0, i], ion = ion)[2])


print(MBohr[100])
tester = bohrToEmu2(MBohr[100])
print(tester)
tester = normalizeSpin(tester,mass,molweight)
print(tester)
print(M[100])
# M = normalizeSpin(M,mass,molweight)
magCalcEmu = bohrToEmu2(magCalcBohr)
magCalcEmu = normalizeSpin(magCalcEmu,mass,molweight)
#####################################################################################################################################################################
# plt.figure()
# plt.plot(H,-1*np.array(magCalcBohr), label = 'PCF')
# plt.xlabel('Field (Oe)')
# plt.ylabel('Magnetization (\N{GREEK SMALL LETTER MU}B mol^-1)')
# plt.title('PCF Magnetization at {} K'.format(Temp))
# plt.legend()
# # plt.show()

# plt.figure()
# plt.plot(H,-1*np.array(magCalcEmu), label = 'PCF')
# plt.xlabel('Field (Oe)')
# plt.ylabel('Magnetization (emu mol^-1)')
# plt.title('PCF Magnetization at {} K'.format(Temp))
# plt.legend()
# plt.show()

plt.figure()
plt.plot(H,-1.*np.array(magCalcEmu), label = 'PCF Calculated')
# plt.errorbar(H,M, yerr = Err, label = 'Measured')
plt.plot(H,MNorm,label = 'Measured')
plt.xlabel('Field (oe)')
plt.ylabel('Magnetization (emu mol^-1)')
plt.title('{} Magnetization at {} K'.format(comp,Temp))
plt.legend()
plt.show()

plt.figure()
plt.plot(HTes,-1.*np.array(magCalcBohr), label = 'PCF Calculated')
plt.errorbar(HTes,MBohr, yerr = ErrBohr, label = 'Measured')
plt.xlabel('Field (T)')
plt.ylabel('Magnetization \N{GREEK SMALL LETTER MU}B/spin')
plt.title('{} Magnetization at {} K'.format(comp,Temp))
plt.legend()
plt.show()
#####################################################################################################################################################################


# ## Susceptibility (M vs T)

# In[24]:


#Loading data for M vs T 
#####################################################################################################################################################################
runs = []
for i in os.listdir(MTDir):
    if i.endswith('.DAT'): #This was a safeguard against a situation arising at an earlier implementation of my code.
        runs.append(i)

MTdata = {}
# plt.figure()
for i in runs:
    M,H,T,MErr,samplemass,measType = getData(i,MTDir, who = who, dataType = 'MT')
    MTdata[measType] = [M,H,T,MErr,samplemass]

#Either 'ZFC' or 'FC'
M,H,T,MErr,samplemass = MTdata['FC']
#####################################################################################################################################################################

massErr = .00005

fieldT = 0.1
deltaField = .00001
if LS_on:
	XCalcBohr = Pr.susceptibility(Temps = T, Field = fieldT, deltaField = deltaField)
else:
	XCalcBohr = Pr.susceptibility(Temps = T, Field = fieldT, deltaField = deltaField, ion = ion)
# for i in Temp:
	# suscCalc.append(Pr.susceptibility(i,fieldT,deltaField))
XiCalcBohr = 1/XCalcBohr

XEmu,XErrEmu,XiEmu,XiErrEmu = normSusc(M,H,MErr,molweight,samplemass,massErr)

# XBohr = emuToBohr(XEmu, samplemass,molweight)
# XErrBohr = emuToBohr(XErrEmu, samplemass,molweight)

# XiBohr = 1/np.array(XBohr)
# XiErrBohr = 1/np.array(XErrBohr)

XCalcEmu = bohrToEmu2(XCalcBohr)
XCalcEmu = normalizeMol(XCalcEmu,samplemass,molweight)
XiCalcEmu = 1/XCalcEmu

# plt.figure()
# plt.plot(T,XCalcEmu)
# plt.xlabel('Temperature (K)')
# plt.ylabel('Susceptibility (unsure unit)')
# plt.title('PCF Susceptibility with Scalar {} T field'.format(fieldT))

plt.figure()
plt.plot(T,-1*XiCalcEmu)
plt.xlabel('Temperature (K)')
plt.ylabel('1/Susceptibility (unsure unit)')
plt.title('PCF Inverse Susceptibility with Scalar {} T field'.format(fieldT))
# plt.show()

plt.figure()
plt.errorbar(T,XiEmu, yerr = XiErrEmu, label = 'Measured')
plt.plot(T,-1*XiCalcEmu, label = 'PCF')
plt.xlabel('Temperature (K)')
plt.ylabel('1/Susceptibility (Oe*mol/emu)')
plt.title('{} Inverse Susceptibility Experiment with Scalar {} T field'.format(comp,fieldT))
plt.legend()
plt.show()

# plt.figure()
# plt.errorbar(T,XBohr, yerr = XiErrBohr, label = 'Measured')
# plt.plot(T,-1*XiCalcBohr, label = 'PCF')
# plt.xlabel('Temperature (K)')
# plt.ylabel('1/Susceptibility (Tesla*mol/BohrMagneton)')
# plt.title('{} Inverse Susceptibility Experiment and Model with Scalar {} T field'.format(comp,fieldT))
# plt.legend()
# plt.show()
#####################################################################################################################################################################


# ## PCF Neutron Spectrum
# Define: an energy array that will be used to calculate lineshape; neutron incident energy, Ei (mev); temperature (K); instrument resolution (meV).
# 
# Instrument resolution is typically ~3% - 5% of the Ei.
# 
# Lineshape fitting would be ideal for determining the CF Hamiltonian; however, nuances arise. If not all magnetic modes are seen at a single Ei then certain scaling has to be taken into account to adjust the intensities of the signal. For this reason we decided to go down the route of fitting to energy levels.

# In[25]:


# Define these values for neutron lineshape calculation
#####################################################################################################################################################################
energy = np.linspace(.01,700,1000)
Ei = 700
Temp = 4.82
res = 9
#####################################################################################################################################################################

# Calculate neutron spectrum and plot it.
#####################################################################################################################################################################
CalculatedSpectrum = Pr.neutronSpectrum(energy, Temp=Temp, Ei=Ei, ResFunc = lambda x: res )
# ResFunc = lambda x: 9 if (energy < 200) else 21

plt.figure()
plt.plot(energy,CalculatedSpectrum)
plt.ylabel('Intensity (arb. units)')
plt.xlabel('Energy (meV)')
plt.title('PCF Spectrum: Ei = {}meV, Temp = {}K, Res = {}'.format(Ei,Temp,res))
plt.show()
#####################################################################################################################################################################


# ## Convert our fitted Stevens' to Wybourne

# In[26]:


# wyb = cef.StevensToWybourne('Ce3+',stev, LS=LS_on)
# print(wyb)


# In[27]:


Pr.printLaTexEigenvectors()


# In[ ]:




