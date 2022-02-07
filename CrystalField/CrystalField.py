import sys
sys.path.append('..')
from JensenTools import *

#####################################################################################################################################################################
# Function to be made into an LMFIT model.
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
			# e.append(e[2]/e[1]) #The aforementioned ratio
		else:
			e =  Pr.eigenvalues
	return e

def energyFit2(B40,B60, B44, B64, B20, numlevels, LS, **kwargs ):
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
			e = kmeansSort2(Pr.eigenvalues,numlevels)[:3] #Excluding the highest mode which we did not detect in our INS runs
			e.append(e[2]/e[1]) #The aforementioned ratio
		else: 
			e = Pr.eigenvalues
	else:
		Pr = cef.CFLevels.Bdict(Bdict = Stev, ion = kwargs['ion'])
		Pr.diagonalize()
		if kwargs['Kmeans']:   	
			e = kmeansSort2(Pr.eigenvalues,numlevels) #Excluding the highest mode which we did not detect in our INS runs
		else:
			e =  Pr.eigenvalues
	return e

# Function to be made into an LMFIT model.
def energyFitAll(B20, B21, B22, B40, B41, B42, B43, B44, B60, B61, B62, B63, B64, B65, B66, numlevels, LS, **kwargs ):
	numlevels = numlevels
	Stev = {} #Creating the Stevens' Coefficients dictionary and assigning values
	Stev['B20'] = B20
	Stev['B21'] = B21
	Stev['B22'] = B22
	Stev['B40'] = B40
	Stev['B40'] = B41
	Stev['B40'] = B42
	Stev['B40'] = B43
	Stev['B40'] = B44
	Stev['B60'] = B60
	Stev['B61'] = B61
	Stev['B62'] = B62
	Stev['B63'] = B63
	Stev['B64'] = B64
	Stev['B65'] = B65
	Stev['B66'] = B66

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
			# e.append(e[2]/e[1]) #The aforementioned ratio
		else:
			e =  Pr.eigenvalues
	return e

# Simple function for pritting my parameters after fitting so that I can copy paste the values from output for further iterations.
def paramPrint(fittedparams):
	print()
	for i in fittedparams:
		# print(i, ' = ', i.value)
		print(i, ' = ',fittedparams[i].value )
#####################################################################################################################################################################

# Define important things
#####################################################################################################################################################################
comp = 'Ba2YbNbO6'
ion = 'Yb3+'
who = 'Arun'
LS_on = False
Kmeans = True
molweight = molweight[comp]
LSValue = 100
per = 'mol'

# The S,L,J values are as follows for the Pr4+ ion
S = 0.5
L = 3
J = 5./2

if LS_on:
	numlevels = 4
	Emeas = [168, 335,385, 385/335] # The measured INS magnetic modes
else:
	numlevels = 4
	Emeas = [168, 335, 335/168] # The measured INS magnetic modes, only first 2 for J basis
#####################################################################################################################################################################

#Load necessary directories
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

# Loading data for Susceptibility (M vs T) and Magnetization (M vs H) data
#####################################################################################################################################################################
runs = []
for i in os.listdir(MTDir):
    if i.endswith('.DAT') or i.endswith('.dat'): #This was a safeguard against a situation arising at an earlier implementation of my code.
        runs.append(i)
# Normalizes and stores data in MTdata dict
MTdata = {}
for i in runs:
    M,H,T,MErr,mass,measType = getData(i,MTDir, who = who, dataType = 'MT')
    M = normalize(M,mass,molweight,per)
    Merr = normalize(MErr,mass,molweight,per)
    MTdata[measType] = [M,H,T,MErr,mass]

runs = []
for i in os.listdir(MHDir):
    if i.endswith('.DAT') or i.endswith('.dat'): #This was a safeguard against a situation arising at an earlier implementation of my code.
        runs.append(i)       
# Normalizes and stores data in MHdata dict
MHdata = {}
for i in runs: 
    M, H, Err, mass, T = getData(i,MHDir,who = who, dataType = 'MH')
    M = normalize(M,mass,molweight,per)
    Err = normalize(Err,mass,molweight,per)
    MHdata[T] = [M,H,Err,mass,i]
#####################################################################################################################################################################


# Make LMFIT model and fit
# Create stevens coefficients dictionary from fitted parameters
#####################################################################################################################################################################
eModel = Model(energyFit2, independent_vars = ['numlevels'])
params = eModel.make_params()

# Since we only have 4 training points, only 4 parameters can vary at once.
params['B20'].set(value = B20, vary = False)
# params['B21'].set(value = 0, vary = False)
# params['B22'].set(value = 0, vary = False)
params['B40'].set(value=B40, vary=False)
# params['B41'].set(value=0, vary=False)
# params['B42'].set(value=0, vary=False)
# params['B43'].set(value=0, vary=False)
params['B44'].set(value=B44, vary=False)
params['B60'].set(value=B60, vary=False)
# params['B61'].set(value=0, vary=False)
# params['B62'].set(value=0, vary=False)
# params['B63'].set(value=0, vary=False)
params['B64'].set(value=B64, vary=False)
# params['B65'].set(value=0, vary=False)
# params['B66'].set(value=0, vary=False)

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
	Pr = cef.CFLevels.Bdict(Bdict = stev, ion = ion)
	Pr.diagonalize()

# Print final matrix
print('\n\nEnergy values as measured by INS (meV): {}'.format(Emeas))
Pr.printEigenvectors()

# Calculate and neatly print G-Tensor using Pandas
gt = Pr.gtensor()
rows = ['gx','gy','gz']
df = pd.DataFrame(gt, columns = rows, index = rows)
print(df)
#####################################################################################################################################################################

# ## PCF Magnetization
#####################################################################################################################################################################
T = '20K'
Temp = getTemp(MHdata[T][-1], who = who)
M, H, Err, mass, filename = MHdata[T]
MBohr = emuToBohr2(M)
ErrBohr = emuToBohr2(Err)
HTes = oeToTesla(H)

#Generate a magnetization curve for comparing results to experiment
magCalc = []
for i in HTes:
	if LS_on:
		magCalc.append((Pr.magnetization( Temp = Temp, Field = [i, 0, 0])[0] + Pr.magnetization( Temp = Temp, Field = [0, i, 0])[1] + Pr.magnetization( Temp = Temp, Field = [0, 0, i])[2])/3)		
	else:
		magCalc.append((Pr.magnetization( Temp = Temp, Field = [i, 0, 0], ion = ion)[0] + Pr.magnetization( Temp = Temp, Field = [0, i, 0], ion = ion)[1] + Pr.magnetization( Temp = Temp, Field = [0, 0, i], ion = ion)[2])/3)		

magCalcEmu = bohrToEmu2(magCalc)*6.02214076*10**23

plt.figure()
plt.plot(H,-1.*np.array(magCalcEmu), label = 'PCF Powder Average')
plt.errorbar(H,M, yerr = Err, label = 'Measured')
plt.xlabel('Field (Oe)')
plt.ylabel('Magnetization (emu {}^-1)'.format(per))
plt.title('{} Magnetization at {} K'.format(comp, Temp))
plt.legend()


# plt.show()
#####################################################################################################################################################################

# ## PCF Susceptibility
#####################################################################################################################################################################
M,H,T,Err,samplemass = MTdata['ZFC']
X = M/H
Xi  = 1/X

MBohr = emuToBohr2(M)
ErrBohr = emuToBohr2(Err)
HTes = oeToTesla(H)

XBohr = MBohr/HTes
XBohrI = 1/XBohr

XCalc = []
# Temp = np.linspace(.001,T,1000)
fieldT = 0.1
deltaField = 0.0001
if LS_on:
	XCalc = Pr.susceptibility(Temps = T, Field = fieldT, deltaField = deltaField)
else:
	XCalc = Pr.susceptibility(Temps = T, Field = fieldT, deltaField = deltaField, ion = ion)
XCalc = XCalc*6.02214076*10**23

XCalcI = 1/XCalc
    
XCalcEmu = bohrToEmu2(XCalc)/10000
XCalcEmuI = 1/XCalcEmu

plt.figure()
plt.plot(T, -1*np.array(XCalcEmuI), label = 'PCF Powder Average')
plt.plot(T, Xi, label = 'Measured')
plt.xlabel('Temperature (K)')
plt.ylabel('1/X (emu^-1 T {})'.format(per))
plt.title('{} Inverse Susceptibility with Scalar {} T field'.format(comp, fieldT))
plt.legend()

plt.show()

# plt.figure()
# plt.plot(T,-1*np.array(XCalc), label = 'PCF')
# plt.plot(T,XBohr, label = 'Masured')
# plt.xlabel('Temperature (K)')
# plt.ylabel('X (uB T^-1 Spin^-1)')
# plt.title('Susceptibility with Scalar {} T field'.format(fieldT))
# plt.legend()
# plt.figure()

# plt.plot(T, -1*np.array(XCalcEmu), label = 'PCF')
# plt.plot(T,X, label = 'Measured')
# plt.xlabel('Temperature (K)')
# plt.ylabel('X (emu oe^-1 Spin^-1)')
# plt.title('Susceptibility with Scalar {} T field'.format(fieldT))
# plt.legend()
# plt.figure()

# plt.plot(T, -1*np.array(XCalcEmuI),label = 'PCF')
# plt.plot(T, Xi, label = 'Measured')
# plt.xlabel('Temperature (K)')
# plt.ylabel('1/X (emu ^-1 Oe Spin)')
# plt.title('Inverse Susceptibility with Scalar {} T field'.format(fieldT))
# plt.legend()

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
# plt.show()
#####################################################################################################################################################################



print()
# Pr.printLaTexEigenvectors()
print()

wyb = cef.StevensToWybourne('Ce3+',stev, LS=True)
print("Fitted coefficients in Wybourne's")
print(wyb)