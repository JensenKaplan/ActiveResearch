#!/usr/bin/env python
# coding: utf-8

# # Crystal Electric Field (CF) Analysis of the Pr4+ ion
# ### Here the goal is to use collected inelastic neutron scattering (INS) data to determine the ground state crystal field Hamiltonian. I use PyCrystalField in conjunction with LMFIT to find the values of the Stevens' Coefficients by fitting predicted energy levels to measured one. I then calculate the CF Hamiltonian and predict the compounds thermodynamic properties. Magnetization (M vs H) and Susceptibility (M vs T) data. 
# 
# ### It's important to note that I don't actually use Pr4+ as my central ion, I use Ce3+. The correct orbital values for the Ln4+ oxidized states have not been calculated and tabulated. Since the operators depend on electronic orbital states, I can use an electrically equivalent ion, Ce3+, instead of the actual Pr4+ ion.

# In[2]:


# get_ipython().run_line_magic('reload_ext', 'autoreload')

import sys
sys.path.append('..')
from JensenTools import *

from lmfit import Model
import PyCrystalField as cef
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd



#Simple function for pritting my parameters after fitting so that I can copy paste the values from output for further iterations
def paramPrint(fittedparams):
	print()
	for i in fittedparams:
		# print(i, ' = ', i.value)
		print(i, ' = ',fittedparams[i].value )


# ##  My Energy Calculation Function To Be Used As An LMFIT Model

# In[3]:


#energy level calculation used with LMFIT.
#calculates 4 energies, the eigenvalues as given by PCF
#Returns 3 energies (chosen because we only have 3 observed modes from INS) as well as a ratio of the 3rd energy to the 2nd energy.
#this ratio comes from the fact that those two are visible in one Ei from INS. It was thought that it may improve fitting, but ultimately doesn't change result
def energyCalcKFit(B40,B60,LS, numlevels, B44, B64, B20):
	numlevels = numlevels
	Stev = {} #Creating the Stevens' Coefficients dictionary and assigning values
	Stev['B20'] = B20
	Stev['B40'] = B40
	Stev['B60'] = B60
	Stev['B44'] = B44
	Stev['B64'] = B64

	Pr = cef.LS_CFLevels.Bdict(Bdict=Stev, L=3, S=0.5, SpinOrbitCoupling=LS) #Create CF_Levels obejct wtih the given coefficients.
	Pr.diagonalize()
	e = kmeansSort(Pr.eigenvalues,numlevels)[:3] #Excluding the highest mode which we did not detect in our INS runs
	e.append(e[2]/e[1]) #The aforementioned ratio
    
	return e

def energyFit(B40,B60, numlevels, B44, B64, B20, **kwargs ):
	numlevels = numlevels
	Stev = {} #Creating the Stevens' Coefficients dictionary and assigning values
	Stev['B20'] = B20
	Stev['B40'] = B40
	Stev['B60'] = B60
	Stev['B44'] = B44
	Stev['B64'] = B64

	if kwargs['LS_on']:
		Pr = cef.LS_CFLevels.Bdict(Bdict=Stev, L=3, S=0.5, SpinOrbitCoupling = kwargs['LSValue']) #Create CF_Levels obejct wtih the given coefficients.
		Pr.diagonalize()
		if kwargs['Kmeans']:
			e = JT.kmeansSort(Pr.eigenvalues,numlevels)[:3] #Excluding the highest mode which we did not detect in our INS runs
			e.append(e[2]/e[1]) #The aforementioned ratio
		else: 
			e = Pr.eigenvalues
	else:
		Pr = cef.CFLevels.Bdict(Bdict = Stev, ion = kwargs['ion'])
		Pr.diagonalize()
		if kwargs['Kmeans']:   	
			e = JT.kmeansSort(Pr.eigenvalues,numlevels) #Excluding the highest mode which we did not detect in our INS runs
		else:
			e =  Pr.eigenvalues
	return e


# ## L and S values for our ion, energy levels as measured by INS, and Stevens' Starting Point From Grid Search

# In[4]:

comp = 'Sr2PrO4'
who = 'Arun'
LS_on = True
Kmeans = True
LSValue = 100

saveDir = getSaveDir('m',comp = comp)
MHDir = getSaveDir('m',comp = comp, dataType = 'MH')
MTDir = getSaveDir('m',comp = comp, dataType = 'MT')


# The L,S values are as follows for the Pr4+ ion
L = 3
S = 0.5
Emeas = [168, 335,385,385/335] #The measured INS magnetic modes. The last entry being a ratio

# The starting values from the grid search
LS = LSValue
x =  0.03629536921151444
bpf = -0.6570713391739674

# Assigning the coefficients from grid search
# Enforcing cubic constraints as a start
# and including the B20 term which is needed for tetragonal symmetry
B40 = bpf
B60 = x*bpf
B44 = 5*B40
B64 = -21*B60
B20 = 0


# ## LMFIT Model Creation

# In[5]:


# Create the LMFIT model using my energy calc function
# Let all values be dependent params except for number of levels
# We can aleways lock or unlock the varying of each param (making them effectively independent/dependent)
eModel = Model(energyCalcKFit, independent_vars = ['numlevels'])
params = eModel.make_params()


## Uncomment below to make sure the proper variables are set
# print('parameter names: {}'.format(eModel.param_names))
# print('independent variables: {}'.format(eModel.independent_vars))


# ## Best fit results.
# ### Since there are more features than training points, I let a few coefficients vary freely, lock in those values, then let the others vary. The values below are after a few iterations,  produced some of the best results.

# In[6]:


B40  =  -0.6568663783690575
B60  =  -0.02328250024945387
LS  =  100.00007580463522
B44  =  -3.1415463304732714
B64  =  0.504906552605772
B20  =  0.4858075931009187


# ## Set the LMFIT parameters and use it to fit to the previously given energies.

# In[7]:


# Since we only have 4 training points, only 4 parameters can vary at once.
params['B20'].set(value = B20, vary = True)
params['B40'].set(value=B40, vary=True)
params['B60'].set(value=B60, vary=True)
params['B44'].set(value = B44, vary = True )
params['B64'].set(value = B64, vary = False )
params['LS'].set(value=LS, vary=False)

#Fit model to data
fitted = eModel.fit(Emeas,params, numlevels = 4 , LS_on = LS_on, Kmeans = Kmeans)

#Print the parameters and reduced chi sqr value
print('\n\nFitted parameters:')
fitted.params.pretty_print()
print('\nReduced Chi Sqr = {}'.format(fitted.result.redchi))

#Uncomment to print out in easy copy paste format
paramPrint(fitted.params)


# ## Use our fitted B parameters to obtain the CF Hamiltonian

# In[8]:


#Create a dictionary of the fitted parameters (stevens coefficients)
stev = {'B40': fitted.params['B40'].value, 'B60': fitted.params['B60'].value, 'B44' : fitted.params['B44'].value, 'B64' : fitted.params['B64'].value, 'B20' :fitted.params['B20'].value }

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


# Pr.printLaTexEigenvectors()


# ## PCF G Tensor

# In[9]:


#Calculate and neatly print G-Tensor using Pandas
gt = Pr.gtensor()
rows = ['gx','gy','gz']
df = pd.DataFrame(gt, columns = rows, index = rows)
print(df)


# ## PCF Magnetization

# In[10]:
runs = []
for i in os.listdir(MHDir):
    if i.endswith('.DAT'): #This was a safeguard against a situation arising at an earlier implementation of my code.
        runs.append(i)

mass = getMass(runs[0])
molweight = 380.15
MHdata = {}
# plt.figure()
for i in runs:
    H, M, Err, T = getData(i,MHDir,who = 'Arun')
    M = emuToBohr(M,mass,molweight)
    H = oeToTesla(H)
    Err = emuToBohr(Err,mass,molweight)
    MHdata[T] = [H,M,Err,i]
    # plt.errorbar(H,M, yerr = Err, label = name)

T = '20K'

# samplemass = getMass(MHdata[T][3])
Temp = getTemp(MHdata[T][3])
H, M, Err, filename = MHdata[T]


#Generate a magnetization curve for comparing results to experiment
magCalc = []
fieldT = np.linspace(0.01,14,1000)

for i in H:
	magCalc.append(Pr.magnetization( Temp = Temp, Field = [0, 0, i])[2])


plt.plot(H,magCalc)
plt.xlabel('Field (T)')
plt.ylabel('Magnetization \N{GREEK SMALL LETTER MU}B')
plt.title('PCF Magnetization at {} K'.format(Temp))

plt.figure()
plt.plot(H,-1.*np.array(magCalc), label = 'PCF Calculated')
plt.errorbar(H,M, yerr = Err, label = 'Measured')
plt.xlabel('Field (T)')
plt.ylabel('Magnetization \N{GREEK SMALL LETTER MU}B')
plt.title('Magnetization at {} K'.format(Temp))
plt.legend()
# plt.show()
#####################################################################################################################################################################

# ## PCF Susceptibility

# In[11]:
runs = []
for i in os.listdir(MTDir):
    if i.endswith('.DAT'): #This was a safeguard against a situation arising at an earlier implementation of my code.
        runs.append(i)


MTdata = {}
# plt.figure()
for i in runs:
    M,H,T,Err,samplemass,measType = getData(i,MTDir, who = who, dataType = 'MT')
    MTdata[measType] = [M,H,T,Err,samplemass]
    # plt.errorbar(H,M, yerr = Err, label = name)

M,H,T,Err,samplemass = MTdata['ZFC']
# samplemass = JT.getMass(MHdata['ZFC'][4])

suscCalc = []
Temp = np.linspace(.1,400,1000)
fieldT = 0.01
deltaField = 0.0001
suscCalc = Pr.susceptibility(Temps = Temp, Field = fieldT, deltaField = deltaField)
# for i in Temp:
	# suscCalc.append(Pr.susceptibility(i,fieldT,deltaField))

suscCalcI = []
for i in suscCalc:
	suscCalcI.append(1/i)
plt.figure()
plt.plot(Temp,suscCalc)
plt.xlabel('Temperature (K)')
plt.ylabel('Susceptibility (unsure unit)')
plt.title('PCF Susceptibility with Scalar {} T field'.format(fieldT))

plt.figure()
plt.plot(Temp,suscCalcI)
plt.xlabel('Temperature (K)')
plt.ylabel('1/Susceptibility (unsure unit)')
plt.title('PCF Inverse Susceptibility with Scalar {} T field'.format(fieldT))
plt.show()




# ## PCF Neutron Spectrum
# Define: an energy array that will be used to calculate lineshape; neutron incident energy, Ei (mev); temperature (K); instrument resolution (meV).
# 
# Instrument resolution is typically ~3% - 5% of the Ei.
# 
# Lineshape fitting would be ideal for determining the CF Hamiltonian; however, nuances arise. If not all magnetic modes are seen at a single Ei then certain scaling has to be taken into account to adjust the intensities of the signal. For this reason we decided to go down the route of fitting to energy levels.

# In[12]:


energy = np.linspace(.01,700,1000)
Ei = 700
Temp = 4.82
res = 9

CalculatedSpectrum = Pr.neutronSpectrum(energy, Temp=Temp, Ei=Ei, ResFunc = lambda x: res )
# ResFunc = lambda x: 9 if (energy < 200) else 21
plt.figure()
plt.plot(energy,CalculatedSpectrum)
plt.ylabel('Intensity (arb. units)')
plt.xlabel('Energy (meV)')
plt.title('PCF Spectrum: Ei = {}meV, Temp = {}K, Res = {}'.format(Ei,Temp,res))
plt.show()


# ## Convert our fitted Stevens' to Wybourne

# In[13]:


# wyb = cef.StevensToWybourne('Ce3+',stev, LS=True)
# print(wyb)

