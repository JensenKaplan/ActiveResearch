#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import sys
sys.path.append('..')
from JensenTools import *


# In[ ]:


# Define important things
#####################################################################################################################################################################
comp = 'Sr2PrO4'
ion = 'Ce3+'
who = 'Arun'
LS_on = True
Kmeans = True
per = 'spin'
molweight = molweight[comp]
LSValue = 100
# The L,S values are as follows for the Pr4+ ion
L = 3
S = 0.5
#####################################################################################################################################################################


# In[ ]:


# Function to be made into an LMFIT model.
def magFit(B40,B60, B44, B64, B20, LS, HTes, **kwargs ):
    Stev = {} #Creating the Stevens' Coefficients dictionary and assigning values
    Stev['B20'] = B20
    Stev['B40'] = B40
    Stev['B60'] = B60
    Stev['B44'] = B44
    Stev['B64'] = B64
    mag = []
    if kwargs['LS_on']:
        Pr = cef.LS_CFLevels.Bdict(Bdict=Stev, L=3, S=0.5, SpinOrbitCoupling = LS) #Create CF_Levels obejct wtih the given coefficients.
        Pr.diagonalize()
        for i in HTes:
            mag.append((Pr.magnetization(Temp = Temp, Field = [i, 0, 0])[0] + Pr.magnetization(Temp = Temp, Field = [0, i, 0])[1] + Pr.magnetization(Temp = Temp, Field = [0, 0, i])[2])/3)
    else:
        Pr = cef.CFLevels.Bdict(Bdict = Stev, ion = kwargs['ion'])
        Pr.diagonalize()
        for i in HTes:
            mag.append((Pr.magnetization(Temp = Temp, Field = [i, 0, 0], ion = ion)[0] + Pr.magnetization(Temp = Temp, Field = [0, i, 0], ion = ion)[1] + Pr.magnetization(Temp = Temp, Field = [0, 0, i], ion = ion)[2])/3)

    return -1*np.array(mag)

# Function to be made into an LMFIT model.
def magFit2(Pr,a,b,c, HTes, **kwargs ):
    # Stev = {} #Creating the Stevens' Coefficients dictionary and assigning values
    # Stev['B20'] = B20
    # Stev['B40'] = B40
    # Stev['B60'] = B60
    # Stev['B44'] = B44
    # Stev['B64'] = B64
    mag = []
    if kwargs['LS_on']:
        # Pr = cef.LS_CFLevels.Bdict(Bdict=Stev, L=3, S=0.5, SpinOrbitCoupling = LS) #Create CF_Levels obejct wtih the given coefficients.
        # Pr.diagonalize()
        for i in HTes:
            mag.append(a*Pr.magnetization(Temp = Temp, Field = [i, 0, 0])[0] + b*Pr.magnetization(Temp = Temp, Field = [0, i, 0])[1] + c*Pr.magnetization(Temp = Temp, Field = [0, 0, i])[2])
    else:
        Pr = cef.CFLevels.Bdict(Bdict = Stev, ion = kwargs['ion'])
        Pr.diagonalize()
        for i in HTes:
            mag.append((Pr.magnetization(Temp = Temp, Field = [i, 0, 0], ion = ion)[0] + Pr.magnetization(Temp = Temp, Field = [0, i, 0], ion = ion)[1] + Pr.magnetization(Temp = Temp, Field = [0, 0, i], ion = ion)[2])/3)

    return -1*np.array(mag)

# In[ ]:


#Best Fit LS
#####################################################################################################################################################################
if LS_on:	
	B40  =  -0.6568663783690575
	B60  =  -0.02328250024945387
	LS  =  100.00007580463522
	B44  =  -3.1415463304732714
	B64  =  0.504906552605772
	B20  =  0.4858075931009187
#####################################################################################################################################################################

#Fix B20 to different values, check g tensor 

# Best Fit J
#####################################################################################################################################################################
if not LS_on:
	# Red Chi = ~5
	# B40  =  -0.5572886105373519
	# B60  =  0.4673
	# B44  =  -3.0342208316734602
	# B64  =  -9.8133
	# B20  =  12.606195910392971

	# # Red Chi = ~.01
	B40  =  -0.5572886105373519
	B60  =  0.4673
	B44  =  -3.0946858584804335
	B64  =  -9.8133
	B20  =  12.606195720794622
#####################################################################################################################################################################


# In[ ]:


saveDir = getSaveDir('m',comp = comp) #General Directory for the project
MHDir = getSaveDir('m',comp = comp, dataType = 'MH') #MvsH data

stev = { 'B20' :B20, 'B40': B40, 'B44' : B44, 'B60': B60, 'B64' : B64 }

#Create the CFLevels object and diagonalize it
if LS_on:
	Pr = cef.LS_CFLevels.Bdict(Bdict = stev, L = L, S = S, SpinOrbitCoupling=LS)
	Pr.diagonalize()
else:
	Pr = cef.CFLevels.Bdict(Bdict = stev, ion = ion)
	Pr.diagonalize()


# In[ ]:


# Loading data for M vs H
#####################################################################################################################################################################
runs = []
for i in os.listdir(MHDir):
    if i.endswith('.DAT'): #This was a safeguard against a situation arising at an earlier implementation of my code.
        runs.append(i)
MHdata = {}
for i in runs:
    M, H, Err, mass, T = getData(i,MHDir,who = who, dataType = 'MH')
    M = normalize(M,mass,molweight, per)
    Err = normalize(Err,mass,molweight, per)
    MHdata[T] = [M,H,Err,mass,i]
    
# Choosing 20K run
T = '20K'
Temp = getTemp(MHdata[T][-1], who = who)
M, H, Err, mass, filename = MHdata[T]

MBohr = emuToBohr2(M)
HTes = oeToTesla(H)
#####################################################################################################################################################################


# In[ ]:


# Make LMFIT model and fit
# Create stevens coefficients dictionary from fitted parameters
#####################################################################################################################################################################
magModel = Model(magFit2, independent_vars = ['HTes', 'Pr'])
params = magModel.make_params()

# Since we only have 4 training points, only 4 parameters can vary at once.
# params['B20'].set(value = B20, vary = True)
# params['B40'].set(value=B40, vary=True)
# params['B60'].set(value=B60, vary=True)
# params['B44'].set(value = B44, vary = True )
# params['B64'].set(value = B64, vary = True )

# if LS_on:
# 	params['LS'].set(value=LS, vary=False)
    

params['a'].set(value = 1/3, min = .0001, max = 1 )
params['b'].set(value = 1/3, min = .0001, max = 1  )
params['c'].set(value = 1/3, min = .0001, max = 1  )
# params['a'].set(value = 1/3)
# params['b'].set(value = 1/3)
# params['c'].set(value = 1/3)
# Fit model to data
fitted = magModel.fit(MBohr,params,Pr = Pr, HTes = HTes, LS_on = LS_on, ion = ion)
# print(len(fitted))
# Create a dictionary of the fitted parameters (stevens coefficients)
# stev = {'B40': fitted.params['B40'].value, 'B60': fitted.params['B60'].value, 'B44' : fitted.params['B44'].value, 'B64' : fitted.params['B64'].value, 'B20' :fitted.params['B20'].value }
#####################################################################################################################################################################

fitted.params.pretty_print()
# In[ ]:


magCalcBohr = []
magCalcBohrPowder = []
for i in HTes:
    if LS_on:
        magCalcBohr.append(Pr.magnetization( Temp = Temp, Field = [0, 0, i])[2])
        magCalcBohrPowder.append((fitted.params['a'].value*Pr.magnetization(Temp = Temp, Field = [i, 0, 0])[0] + fitted.params['b'].value*Pr.magnetization(Temp = Temp, Field = [0, i, 0])[1] + fitted.params['c'].value*Pr.magnetization(Temp = Temp, Field = [0, 0, i])[2]))

    else:
        magCalcBohr.append(Pr.magnetization( Temp = Temp, Field = [0, 0, i], ion = ion)[2])
        magCalcBohrPowder.append((Pr.magnetization(Temp = Temp, Field = [i, 0, 0], ion = ion)[0] + Pr.magnetization(Temp = Temp, Field = [0, i, 0], ion = ion)[1] + Pr.magnetization(Temp = Temp, Field = [0, 0, i], ion = ion)[2])/3)
        
magCalcBohr = np.array(magCalcBohr)
magCalcBohrPowder = np.array(magCalcBohrPowder)


# In[ ]:


plt.figure()
# plt.plot(HTes,-1*magCalcBohrPowder, label = 'PCF Calculated')
plt.plot(HTes,MBohr, label = 'Measured')
plt.plot(HTes, fitted.best_fit, label = 'Powder Averaging Fitted')
plt.xlabel('Field (T)')
plt.ylabel('Magnetization (uB / spin)')
plt.title('Powder {} Magnetization at {} K'.format(comp,Temp))
plt.legend()

plt.show()

