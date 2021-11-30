import sys
sys.path.append('..')
from JensenTools import *
from lmfit import Model
import PyCrystalField as cef
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


#####################################################################################################################################################################
#Simple function for pritting my parameters after fitting so that I can copy paste the values from output for further iterations
def paramPrint(fittedparams):
	print()
	for i in fittedparams:
		# print(i, ' = ', i.value)
		print(i, ' = ',fittedparams[i].value )

def energyFit(B40,B60, B44, B64, B20, numlevels, LS, **kwargs ):
	numlevels = numlevels
	Stev = {} #Creating the Stevens' Coefficients dictionary and assigning values
	Stev['B20'] = B20
	Stev['B40'] = B40
	Stev['B60'] = B60
	Stev['B44'] = B44
	Stev['B64'] = B64

	if kwargs['LS_on']:
		Pr = cef.LS_CFLevels.Bdict(Bdict=Stev, L=3, S=0.5, SpinOrbitCoupling = LS) #Create CF_Levels obejct wtih the given coefficients.
		Pr.diagonalize()
		if kwargs['Kmeans']:
			e = kmeansSort(Pr.eigenvalues,numlevels)[:numlevels-1] #Excluding the highest mode which we did not detect in our INS runs
			# e.append(e[2]/e[1]) #The aforementioned ratio
		else: 
			e = Pr.eigenvalues
	else:
		Pr = cef.CFLevels.Bdict(Bdict = Stev, ion = kwargs['ion'])
		Pr.diagonalize()
		if kwargs['Kmeans']:   	
			e = kmeansSort(Pr.eigenvalues,numlevels)[:numlevels-1] #Excluding the highest mode which we did not detect in our INS runs
			e.append(e[1]/e[0]) #The aforementioned ratio
		else:
			e =  Pr.eigenvalues
	return e
#####################################################################################################################################################################


#Define important things
#####################################################################################################################################################################
comp = 'Sr2PrO4'
ion = 'Ce3+'
who = 'Arun'
LS_on = False
Kmeans = True
molweight = molweight[comp]
LSValue = 100
if LS_on:
	numlevels = 4
	Emeas = [168, 335,385] #The measured INS magnetic modes
else:
	numlevels = 3
	Emeas = [168, 335, 335/168] #The measured INS magnetic modes, only first 2 for J basis
# The L,S values are as follows for the Pr4+ ion
L = 3
S = 0.5
#####################################################################################################################################################################

saveDir = getSaveDir('m',comp = comp)
MHDir = getSaveDir('m',comp = comp, dataType = 'MH')
MTDir = getSaveDir('m',comp = comp, dataType = 'MT')

# From GridSearch For LS
#####################################################################################################################################################################
if LS_on:
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
#####################################################################################################################################################################

# From GridSearch For J
#####################################################################################################################################################################
if not LS_on:
	x = -1.0000
	bpf = -0.4673
	# Assigning the coefficients from grid search
	# Enforcing cubic constraints as a start
	# and including the B20 term which is needed for tetragonal symmetry	
	B40 = bpf
	B60 = x*bpf
	B44 = 5*B40
	B64 = -21*B60
	B20 = 0
#####################################################################################################################################################################

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

# Make LMFIT model and fit
# Create stevens coefficients dictionary from fitted parameters
#####################################################################################################################################################################
eModel = Model(energyFit, independent_vars = ['numlevels'])
params = eModel.make_params()

# Since we only have 4 training points, only 4 parameters can vary at once.
params['B20'].set(value = B20, vary = False)
params['B40'].set(value=B40, vary=False)
params['B60'].set(value=B60, vary=False)
params['B44'].set(value = B44, vary = False )
params['B64'].set(value = B64, vary = False )
if LS_on:
	params['LS'].set(value=LS, vary=False)
# Fit model to data
fitted = eModel.fit(Emeas,params, numlevels = numlevels, LS_on = LS_on, Kmeans = Kmeans, ion = ion)
# Create a dictionary of the fitted parameters (stevens coefficients)
stev = {'B40': fitted.params['B40'].value, 'B60': fitted.params['B60'].value, 'B44' : fitted.params['B44'].value, 'B64' : fitted.params['B64'].value, 'B20' :fitted.params['B20'].value }
#####################################################################################################################################################################

# Print the parameters and reduced chi sqr value
#####################################################################################################################################################################
print('\n\nFitted parameters:')
fitted.params.pretty_print()
print('\nReduced Chi Sqr = {}'.format(fitted.result.redchi))
#Uncomment to print out in easy copy paste format
paramPrint(fitted.params)
#####################################################################################################################################################################

# CF Analysis
#####################################################################################################################################################################
# Create the CFLevels object and diagonalize it
if LS_on:
	Pr = cef.LS_CFLevels.Bdict(Bdict = stev, L = L, S = S, SpinOrbitCoupling=fitted.params['LS'].value)
	Pr.diagonalize()
else:
	Pr = cef.CFLevels.Bdict(Bdict = stev, ion = 'Ce3+')
	Pr.diagonalize()
	

# Print final matrix
print('\n\nEnergy values as measured by INS (meV): {}'.format(Emeas[:-1]))
Pr.printEigenvectors()

# Calculate and neatly print G-Tensor using Pandas
gt = Pr.gtensor()
rows = ['gx','gy','gz']
df = pd.DataFrame(gt, columns = rows, index = rows)
print(df)
#####################################################################################################################################################################

# ## PCF Magnetization
#####################################################################################################################################################################
runs = []
for i in os.listdir(MHDir):
    if i.endswith('.DAT'): #This was a safeguard against a situation arising at an earlier implementation of my code.
        runs.append(i)

# mass = getMass(runs[0],who = who)

MHdata = {}
# plt.figure()
for i in runs:
    M, H, Err,  mass, T = getData(i,MHDir,who = who, dataType = 'MH')
    M = normalize(M,mass,molweight, 'spin')
    Err = normalize(Err,mass,molweight, 'spin')
    MHdata[T] = [M,H,Err,mass,i]
    # plt.errorbar(H,M, yerr = Err, label = name)

T = '20K'
Temp = getTemp(MHdata[T][-1], who = who)
M, H, Err, mass, filename = MHdata[T]
M = emuToBohr2(M)
Err = emuToBohr2(Err)
H = oeToTesla(H)

#Generate a magnetization curve for comparing results to experiment
magCalc = []
for i in H:
	if LS_on:
		magCalc.append(Pr.magnetization( Temp = Temp, Field = [0, 0, i])[2])		
	else:
		magCalc.append(Pr.magnetization( Temp = Temp, Field = [0, 0, i], ion = ion)[2])


plt.figure()
plt.plot(H,-1.*np.array(magCalc), label = 'PCF Calculated')
plt.errorbar(H,M, yerr = Err, label = 'Measured')
plt.xlabel('Field (T)')
plt.ylabel('Magnetization \N{GREEK SMALL LETTER MU}B')
plt.title('Magnetization at {} K'.format(Temp))
plt.legend()
plt.show()
#####################################################################################################################################################################

# ## PCF Susceptibility
#####################################################################################################################################################################
runs = []
for i in os.listdir(MTDir):
    if i.endswith('.DAT') or i.endswith('.dat'): #This was a safeguard against a situation arising at an earlier implementation of my code.
        runs.append(i)

MTdata = {}
# plt.figure()
for i in runs:
    M,H,T,Err,samplemass,measType = getData(i,MTDir, who = who, dataType = 'MT')
    M = normalize(M,mass,molweight, 'spin')
    Err = normalize(Err,mass,molweight, 'spin')    
    MTdata[measType] = [M,H,T,Err,samplemass]
    # plt.errorbar(H,M, yerr = Err, label = name)

M,H,T,Err,samplemass = MTdata['ZFC']
M = emuToBohr2(M)
Err = emuToBohr2(Err)
H = oeToTesla(H)
# samplemass = JT.getMass(MHdata['ZFC'][4])

suscCalc = []
Temp = np.linspace(.1,400,1000)
fieldT = 0.1
deltaField = 0.0001
if LS_on:
	suscCalc = Pr.susceptibility(Temps = Temp, Field = fieldT, deltaField = deltaField)
else:
	suscCalc = Pr.susceptibility(Temps = Temp, Field = fieldT, deltaField = deltaField, ion = ion)
# for i in Temp:
	# suscCalc.append(Pr.susceptibility(i,fieldT,deltaField))

suscCalcI = 1/suscCalc

plt.figure()
plt.plot(Temp,-1*np.array(suscCalc))
plt.xlabel('Temperature (K)')
plt.ylabel('X (uB T^-1 Spin^-1)')
plt.title('PCF Susceptibility with Scalar {} T field'.format(fieldT))

plt.figure()
plt.plot(Temp,-1*np.array(suscCalcI))
plt.xlabel('Temperature (K)')
plt.ylabel('1/X (uB ^-1 T Spin)')
plt.title('PCF Inverse Susceptibility with Scalar {} T field'.format(fieldT))
# plt.show()
#####################################################################################################################################################################



# ## PCF Neutron Spectrum
#####################################################################################################################################################################
Ei = 700
Temp = 4.82
res = 9
energy = np.linspace(.01,Ei,1000)

CalculatedSpectrum = Pr.neutronSpectrum(energy, Temp=Temp, Ei=Ei, ResFunc = lambda x: res )
# ResFunc = lambda x: 9 if (energy < 200) else 21
plt.figure()
plt.plot(energy,CalculatedSpectrum)
plt.ylabel('Intensity (arb. units)')
plt.xlabel('Energy (meV)')
plt.title('PCF Spectrum: Ei = {}meV, Temp = {}K, Res = {}'.format(Ei,Temp,res))
plt.show()
#####################################################################################################################################################################

print()
Pr.printLaTexEigenvectors()

# wyb = cef.StevensToWybourne('Ce3+',stev, LS=True)
# print(wyb)

