import sys
sys.path.append('..')
from JensenTools import *

# Define important things
#####################################################################################################################################################################
comp = 'Sr2PrO4'
ion = 'Ce3+'
who = 'Arun'
LS_on = True
per = 'spin'
molweight = molweight[comp]
LSValue = 100.5

# The S,L,J values are as follows for the Pr4+ ion
S = 0.5
L = 3
J = 5./2

Emeas = [0, 168, 335, 385, 385/335] # The measured INS magnetic modes
numlevels = 6
#####################################################################################################################################################################

#LMFIT Models
#####################################################################################################################################################################
# Fits the concatenated X^-1 and magnetization
def fullFit(B20, B40,B60, B44, B64, LS, TempX, FieldX, TempM, FieldM, **kwargs ):
    deltaField = .0001
    numlevels = 4

    Stev = {} #Creating the Stevens' Coefficients dictionary and assigning values
    Stev['B20'] = B20
    Stev['B40'] = B40
    Stev['B60'] = B60
    Stev['B44'] = B44
    Stev['B64'] = B64
    
    M = []
    if kwargs['LS_on']:
        Pr = cef.LS_CFLevels.Bdict(Bdict=Stev, L=3, S=0.5, SpinOrbitCoupling = LS) #Create CF_Levels obejct wtih the given coefficients.
        Pr.diagonalize()
        X = Pr.susceptibility(Temps = TempX, Field = FieldX, deltaField = deltaField)
        for i in FieldM:
            M.append(1/3*Pr.magnetization(Temp = TempM, Field = [i, 0, 0])[0] + 1/3*Pr.magnetization(Temp = TempM, Field = [0, i, 0])[1] + 1/3*Pr.magnetization(Temp = TempM, Field = [0, 0, i])[2])
    
    else:
        Pr = cef.CFLevels.Bdict(Bdict = Stev, ion = kwargs['ion'])
        Pr.diagonalize()
        X = Pr.susceptibility(Temps = TempX, Field = FieldX, deltaField = deltaField, ion = ion)
        for i in FieldM:
            M.append(1/3*Pr.magnetization(Temp = TempM, Field = [i, 0, 0], ion=ion)[0] + 1/3*Pr.magnetization(Temp = TempM, Field = [0, i, 0], ion = ion)[1] + 1/3*Pr.magnetization(Temp = TempM, Field = [0, 0, i], ion = ion)[2])
    
    e = kmeansSort(Pr.eigenvalues,numlevels)[:numlevels-1] #Excluding the highest mode which we did not detect in our INS runs
    M = -1*np.array(M)
    Xi = -1/X
    total = np.concatenate((e,Xi,M), axis = None)
    return total

def fullFit2(B20, B21, B22, B40, B41, B42, B43, B44, B60, B61, B62, B63, B64, B65, B66, LS, TempX, FieldX, TempM, FieldM, **kwargs ):
    numlevels = 4
    deltaField = .0001

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

    M = []
    if kwargs['LS_on']:
        Pr = cef.LS_CFLevels.Bdict(Bdict=Stev, L=3, S=0.5, SpinOrbitCoupling = LS) #Create CF_Levels obejct wtih the given coefficients.
        Pr.diagonalize()
        X = Pr.susceptibility(Temps = TempX, Field = FieldX, deltaField = deltaField)
        for i in FieldM:
            M.append(1/3*Pr.magnetization(Temp = TempM, Field = [i, 0, 0])[0] + 1/3*Pr.magnetization(Temp = TempM, Field = [0, i, 0])[1] + 1/3*Pr.magnetization(Temp = TempM, Field = [0, 0, i])[2])
    else:
        Pr = cef.CFLevels.Bdict(Bdict = Stev, ion = kwargs['ion'])
        Pr.diagonalize()
        X = Pr.susceptibility(Temps = TempX, Field = FieldX, deltaField = deltaField, ion = ion)
        for i in FieldM:
            M.append(1/3*Pr.magnetization(Temp = TempM, Field = [i, 0, 0], ion=ion)[0] + 1/3*Pr.magnetization(Temp = TempM, Field = [0, i, 0], ion = ion)[1] + 1/3*Pr.magnetization(Temp = TempM, Field = [0, 0, i], ion = ion)[2])

    e = kmeansSort(Pr.eigenvalues,numlevels)[:numlevels-1] #Excluding the highest mode which we did not detect in our INS runs
    
    M = -1*np.array(M)
    Xi = -1/X
    total = np.concatenate((e,Xi,M), axis = None)
    return total

# Fits the concatenated X^-1 and magnetization
def thermoFit(B20, B40,B60, B44, B64, LS, TempX, FieldX, TempM, FieldM, **kwargs ):
    deltaField = .0001
    
    Stev = {} #Creating the Stevens' Coefficients dictionary and assigning values
    Stev['B20'] = B20
    Stev['B40'] = B40
    Stev['B60'] = B60
    Stev['B44'] = B44
    Stev['B64'] = B64
    
    M = []
    if kwargs['LS_on']:
        Pr = cef.LS_CFLevels.Bdict(Bdict=Stev, L=3, S=0.5, SpinOrbitCoupling = LS) #Create CF_Levels obejct wtih the given coefficients.
        Pr.diagonalize()
        X = Pr.susceptibility(Temps = TempX, Field = FieldX, deltaField = deltaField)
        for i in FieldM:
            M.append(1/3*Pr.magnetization(Temp = TempM, Field = [i, 0, 0])[0] + 1/3*Pr.magnetization(Temp = TempM, Field = [0, i, 0])[1] + 1/3*Pr.magnetization(Temp = TempM, Field = [0, 0, i])[2])
    else:
        Pr = cef.CFLevels.Bdict(Bdict = Stev, ion = kwargs['ion'])
        Pr.diagonalize()
        X = Pr.susceptibility(Temps = TempX, Field = FieldX, deltaField = deltaField, ion = ion)
        for i in FieldM:
            M.append(1/3*Pr.magnetization(Temp = TempM, Field = [i, 0, 0], ion=ion)[0] + 1/3*Pr.magnetization(Temp = TempM, Field = [0, i, 0], ion = ion)[1] + 1/3*Pr.magnetization(Temp = TempM, Field = [0, 0, i], ion = ion)[2])
    M = -1*np.array(M)
    Xi = -1/X
    total = np.concatenate((Xi,M), axis = None)
    return total

def thermoFit2(B20, B21, B22, B40, B41, B42, B43, B44, B60, B61, B62, B63, B64, B65, B66, LS, TempX, FieldX, TempM, FieldM, **kwargs ):
    deltaField = .0001

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

    M = []
    if kwargs['LS_on']:
        Pr = cef.LS_CFLevels.Bdict(Bdict=Stev, L=3, S=0.5, SpinOrbitCoupling = LS) #Create CF_Levels obejct wtih the given coefficients.
        Pr.diagonalize()
        X = Pr.susceptibility(Temps = TempX, Field = FieldX, deltaField = deltaField)
        for i in FieldM:
            M.append(1/3*Pr.magnetization(Temp = TempM, Field = [i, 0, 0])[0] + 1/3*Pr.magnetization(Temp = TempM, Field = [0, i, 0])[1] + 1/3*Pr.magnetization(Temp = TempM, Field = [0, 0, i])[2])
    else:
        Pr = cef.CFLevels.Bdict(Bdict = Stev, ion = kwargs['ion'])
        Pr.diagonalize()
        X = Pr.susceptibility(Temps = TempX, Field = FieldX, deltaField = deltaField, ion = ion)
        for i in FieldM:
            M.append(1/3*Pr.magnetization(Temp = TempM, Field = [i, 0, 0], ion=ion)[0] + 1/3*Pr.magnetization(Temp = TempM, Field = [0, i, 0], ion = ion)[1] + 1/3*Pr.magnetization(Temp = TempM, Field = [0, 0, i], ion = ion)[2])
    M = -1*np.array(M)
    Xi = -1/X
    total = np.concatenate((Xi,M), axis = None)
    return total

# Fits X^-1
def susFit(B20, B40,B60, B44, B64, LS, TempX, FieldX, TempM, FieldM, **kwargs ):
    deltaField = .0001
    
    Stev = {} #Creating the Stevens' Coefficients dictionary and assigning values
    Stev['B20'] = B20
    Stev['B40'] = B40
    Stev['B60'] = B60
    Stev['B44'] = B44
    Stev['B64'] = B64
    
    if kwargs['LS_on']:
        Pr = cef.LS_CFLevels.Bdict(Bdict=Stev, L=3, S=0.5, SpinOrbitCoupling = LS) #Create CF_Levels obejct wtih the given coefficients.
        Pr.diagonalize()
        X = Pr.susceptibility(Temps = TempX, Field = FieldX, deltaField = deltaField)
    else:
        Pr = cef.CFLevels.Bdict(Bdict = Stev, ion = kwargs['ion'])
        Pr.diagonalize()
        X = Pr.susceptibility(Temps = TempX, Field = FieldX, deltaField = deltaField, ion = ion)
    Xi = -1/X
    return Xi

# Fits magnetization
def magFit(B20, B40,B60, B44, B64, LS, TempX, FieldX, TempM, FieldM, **kwargs ):
    deltaField = .0001
    
    Stev = {} #Creating the Stevens' Coefficients dictionary and assigning values
    Stev['B20'] = B20
    Stev['B40'] = B40
    Stev['B60'] = B60
    Stev['B44'] = B44
    Stev['B64'] = B64
    
    M = []
    if kwargs['LS_on']:
        Pr = cef.LS_CFLevels.Bdict(Bdict=Stev, L=3, S=0.5, SpinOrbitCoupling = LS) #Create CF_Levels obejct wtih the given coefficients.
        Pr.diagonalize()
        for i in FieldM:
            M.append(1/3*Pr.magnetization(Temp = TempM, Field = [i, 0, 0])[0] + 1/3*Pr.magnetization(Temp = TempM, Field = [0, i, 0])[1] + 1/3*Pr.magnetization(Temp = TempM, Field = [0, 0, i])[2])
        # M = 1/3*Pr.magnetization(Temp = TempM, Field = [FieldM, 0, 0])[0] + 1/3*Pr.magnetization(Temp = TempM, Field = [0, FieldM, 0])[1] + 1/3*Pr.magnetization(Temp = TempM, Field = [0, 0, FieldM])[2]
 
    else:
        Pr = cef.CFLevels.Bdict(Bdict = Stev, ion = kwargs['ion'])
        Pr.diagonalize()
        for i in FieldM:
            M.append(1/3*Pr.magnetization(Temp = TempM, Field = [i, 0, 0], ion=ion)[0] + 1/3*Pr.magnetization(Temp = TempM, Field = [0, i, 0], ion = ion)[1] + 1/3*Pr.magnetization(Temp = TempM, Field = [0, 0, i], ion = ion)[2])
    M = -1*np.array(M)
    return M

# Fitting to eigenvalues
def energyFit(B40, B60, B44, B64, B20, LS, TempX, FieldX, TempM, FieldM, **kwargs ):
    numlevels = kwargs['numlevels']
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
            e.append(e[:-1]/e[:-2]) #The aforementioned ratio
        else:
            e =  Pr.eigenvalues
    return e
#####################################################################################################################################################################


#Best Fit LS
#####################################################################################################################################################################
if LS_on:	
	B40  =  -0.6568663783690575
	B60  =  -0.02328250024945387
	LS   =  LSValue
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


# # From GridSearch For LS
# #####################################################################################################################################################################
# if LS_on:
# 	LS = LSValue
# 	x =  0.03629536921151444
# 	bpf = -0.6570713391739674
# 	# Assigning the coefficients from grid search
# 	# Enforcing cubic constraints as a start
# 	# and including the B20 term which is needed for tetragonal symmetry
# 	B40 = bpf
# 	B60 = x*bpf
# 	B44 = 5*B40
# 	B64 = -21*B60
# 	B20 = 0
# #####################################################################################################################################################################

# # From GridSearch For J
# #####################################################################################################################################################################
# if not LS_on:
# 	x = -1.0000
# 	bpf = -0.4673
# 	# Assigning the coefficients from grid search
# 	# Enforcing cubic constraints as a start
# 	# and including the B20 term which is needed for tetragonal symmetry	
# 	B40 = bpf
# 	B60 = x*bpf
# 	B44 = 5*B40
# 	B64 = -21*B60
# 	B20 = 0
# #####################################################################################################################################################################



saveDir = getSaveDir('m',comp = comp) #General Directory for the project
MTDir = getSaveDir('m',comp = comp, dataType = 'MT') #MvsT data
MHDir = getSaveDir('m',comp = comp, dataType = 'MH') #MvsT data


# Loading data for M vs T 
#####################################################################################################################################################################
runs = []
for i in os.listdir(MTDir):
    if i.endswith('.DAT') or i.endswith('.dat'): #This was a safeguard against a situation arising at an earlier implementation of my code.
        runs.append(i)
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
# Normalizes and stores data as well as plotting in Emu/Oe for all temperatures.
MHdata = {}
for i in runs: 
    M, H, Err, mass, T = getData(i,MHDir,who = who, dataType = 'MH')
    M = normalize(M,mass,molweight,per)
    Err = normalize(Err,mass,molweight,per)
    MHdata[T] = [M,H,Err,mass,i]
#####################################################################################################################################################################

#Either 'ZFC' or 'FC'
M,H,TempX,MErr,mass = MTdata['FC']
MBohr = emuToBohr2(M)
HTes = oeToTesla(H)
XBohr = MBohr/HTes
XiBohr = 1/XBohr
XEmu = M/H
XiEmu = 1/XEmu

# Choosing 20K run
Tmh = '20K'
TempM = getTemp(MHdata[Tmh][-1], who = who)
M, H, Err, mass, filename = MHdata[Tmh]
MBohr = emuToBohr2(M)
HTes = oeToTesla(H)
# X = M/H
# Xi = 1/X

# total = np.concatenate((Emeas,Xi,M), axis = None)
total = np.concatenate((XiBohr,MBohr), axis = None)

# ENorm = 1/89/len(Emeas)*np.ones(len(Emeas))
# XiNorm = 44/89/len(Xi)*np.ones(len(Xi))
# MNorm = 44/89/len(M)*np.ones(len(M))
# error = np.concatenate((ENorm,XiNorm,MNorm),axis = None)


# Make LMFIT model and fit
# Create stevens coefficients dictionary from fitted parameters
#####################################################################################################################################################################
myModel = Model(susFit, independent_vars = ['TempX', 'FieldX', 'TempM', 'FieldM'])
params = myModel.make_params()

# Since we only have 4 training points, only 4 parameters can vary at once.
params['B20'].set(value=B20, vary=True)
# params['B21'].set(value = 0, vary = True)
# params['B22'].set(value = 0, vary = True)
params['B40'].set(value=B40, vary=True)
# params['B41'].set(value=0, vary=True)
# params['B42'].set(value=0, vary=True)
# params['B43'].set(value=0, vary=True)
params['B44'].set(value=B44, vary=True)
params['B60'].set(value=B60, vary=True)
# params['B61'].set(value=0, vary=True)
# params['B62'].set(value=0, vary=True)
# params['B63'].set(value=0, vary=True)
params['B64'].set(value=B64, vary=True)
# params['B65'].set(value=0, vary=True)
# params['B66'].set(value=0, vary=True)

if LS_on:
	params['LS'].set(value=LS, vary=False)
    
# Fit model to data
fitted = myModel.fit(XiBohr,params, TempX = TempX, FieldX = .1, TempM = TempM, FieldM = HTes, LS_on = LS_on, ion = ion, numlevels = numlevels, Kmeans = True)
# fitted = myModel.fit(Xi,params, TempX = TempX, FieldX = .1, TempM = TempM, FieldM = FieldM, LS_on = LS_on, ion = ion, weights = error)

# Create a dictionary of the fitted parameters (stevens coefficients)
stev = {'B20' :fitted.params['B20'].value, 'B40': fitted.params['B40'].value, 'B44': fitted.params['B44'].value, 'B60': fitted.params['B60'].value,'B64': fitted.params['B64'].value}
# stev = {'B20' :fitted.params['B20'].value, 'B21' :fitted.params['B21'].value, 'B22' :fitted.params['B22'].value, 'B40': fitted.params['B40'].value,  'B41': fitted.params['B41'].value,  'B42': fitted.params['B42'].value,  'B43': fitted.params['B43'].value,  'B44': fitted.params['B44'].value, 'B60': fitted.params['B60'].value, 'B61': fitted.params['B61'].value, 'B62': fitted.params['B62'].value, 'B63': fitted.params['B63'].value, 'B64': fitted.params['B64'].value, 'B65': fitted.params['B65'].value, 'B66': fitted.params['B66'].value}

# Create the CFLevels object and diagonalize it
if LS_on:
	Pr = cef.LS_CFLevels.Bdict(Bdict = stev, L = L, S = S, SpinOrbitCoupling=fitted.params['LS'].value)
	Pr.diagonalize()
else:
	Pr = cef.CFLevels.Bdict(Bdict = stev, ion = ion)
	Pr.diagonalize()
#####################################################################################################################################################################

params.pretty_print()
magCalcPowderBohr = []

if LS_on:
    XCalcPowderBohr = Pr.susceptibility(Temps = TempX, Field = .1, deltaField = .0001)
    for i in HTes:
        magCalcPowderBohr.append((Pr.magnetization(Temp = TempM, Field = [i, 0, 0])[0] + Pr.magnetization(Temp = TempM, Field = [0, i, 0])[1] + Pr.magnetization(Temp = TempM, Field = [0, 0, i])[2])/3)

else:
    XCalcPowderBohr = Pr.susceptibility(Temps = TempX, Field = .1, deltaField = .0001, ion = ion)
    for i in HTes:
        magCalcPowderBohr.append((Pr.magnetization(Temp = TempM, Field = [i, 0, 0], ion = ion)[0] + Pr.magnetization(Temp = TempM, Field = [0, i, 0], ion = ion)[1] + Pr.magnetization(Temp = TempM, Field = [0, 0, i], ion = ion)[2])/3)


XCalcPowderEmu = bohrToEmu2(XCalcPowderBohr)*6.02214076*10**23/10000
XiCalcPowderEmu = -1/XCalcPowderEmu

XEmu = XEmu*6.02214076*10**23
XiEmu = 1/XEmu

# XiCalcPowderBohr = -1/XCalcPowderBohr
magCalcPowderBohr = -1*np.array(magCalcPowderBohr)
magCalcPowderEmu = bohrToEmu2(magCalcPowderBohr)*6.02214076*10**23

M = M*6.02214076*10**23

plt.figure()
plt.plot(TempX,XiEmu, label = 'Measured')
plt.plot(TempX,XiCalcPowderEmu, label = 'PCF Powder Average')
plt.xlabel('Temperature (K)')
plt.ylabel('X^-1 (emu^-1 T mol)')
plt.legend()
plt.title('{} Inverse Susceptbility Applied {} Tesla '.format(comp,.1))

plt.figure()
plt.plot(H,M, label = 'Measured')
plt.plot(H,magCalcPowderEmu, label = 'PCF Powder Average')
plt.xlabel('Field (oe)')
plt.ylabel('M (emu mol^-1)')
plt.legend()
plt.title('{} Magnetization at {} K '.format(comp,TempM))

Pr.printEigenvectors()

Pr.printLaTexEigenvectors()
plt.show()
