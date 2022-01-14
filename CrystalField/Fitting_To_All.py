import sys
sys.path.append('..')
from JensenTools import *

# Define important things
#####################################################################################################################################################################
comp = 'Ba2YbNbO6'
ion = 'Yb3+'
who = 'PPMS'
LS_on = False
per = 'spin'
molweight = molweight[comp]
LSValue = 100.5

# The L,S values are as follows for the Pr4+ ion
L = 3
S = 0.5

# Emeas = [168, 335, 385] # The measured INS magnetic modes
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

# Fits the concatenated X^-1 and magnetization
def thermoFit(B40,B60, B44, B64, LS, TempX, FieldX, TempM, FieldM, **kwargs ):
    deltaField = .0001
    
    Stev = {} #Creating the Stevens' Coefficients dictionary and assigning values
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

# Fits X^-1
def susFit( B40,B60, B44, B64, LS, TempX, FieldX, TempM, FieldM, **kwargs ):
    deltaField = .0001
    
    Stev = {} #Creating the Stevens' Coefficients dictionary and assigning values
    Stev['B40'] = B40
    Stev['B60'] = B60
    Stev['B44'] = B44
    Stev['B64'] = B64
    
#     M = []
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
#####################################################################################################################################################################


saveDir = getSaveDir('m',comp = comp) #General Directory for the project
PCOLig, Yb = cef.importCIF(saveDir + 'Ba2YbNbO6.cif','Yb1')
MTDir = getSaveDir('m',comp = comp, dataType = 'MT') #MvsT data
MHDir = getSaveDir('m',comp = comp, dataType = 'MH') #MvsT data


B40 = Yb.B[0]
B44 = Yb.B[1]
B60 = Yb.B[2]
B64 = Yb.B[3]

# print(MHdata.keys())
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
X = MBohr/HTes
Xi = 1/X


# print(MHdata.keys())
# Choosing 20K run
Tmh = 1.8
TempM = getTemp(MHdata[Tmh][-1], who = who)
M, FieldM, Err, mass, filename = MHdata[Tmh]
M = emuToBohr2(M)
FieldM = oeToTesla(FieldM)

total = np.concatenate((-Xi,-M), axis = None)

# ENorm = 1/7/len(Emeas)*np.ones(len(Emeas))
# XiNorm = 3/7/len(Xi)*np.ones(len(Xi))
# MNorm = 3/7/len(M)*np.ones(len(M))
# error = np.concatenate((ENorm,XiNorm,MNorm),axis = None)

# Make LMFIT model and fit
# Create stevens coefficients dictionary from fitted parameters
#####################################################################################################################################################################
susModel = Model(thermoFit, independent_vars = ['TempX', 'FieldX', 'TempM', 'FieldM'])
params = susModel.make_params()

# Since we only have 4 training points, only 4 parameters can vary at once.
# params['B20'].set(value = B20, vary = False)
params['B40'].set(value = B40, vary=True)
params['B60'].set(value = B60, vary=True)
params['B44'].set(value = B44, vary = True)
params['B64'].set(value = B64, vary = True)

if LS_on:
	params['LS'].set(value=LS, vary=True)


# Fit model to data
fitted = susModel.fit(total,params, TempX = TempX, FieldX = .1, TempM = TempM, FieldM = FieldM, LS_on = LS_on, ion = ion)

# Create a dictionary of the fitted parameters (stevens coefficients)
stev = {'B40': fitted.params['B40'].value, 'B60': fitted.params['B60'].value, 'B44' : fitted.params['B44'].value, 'B64' : fitted.params['B64'].value}

# Create the CFLevels object and diagonalize it
if LS_on:
	Pr = cef.LS_CFLevels.Bdict(Bdict = stev, L = L, S = S, SpinOrbitCoupling=fitted.params['LS'].value)
	Pr.diagonalize()
else:
	Pr = cef.CFLevels.Bdict(Bdict = stev, ion = ion)
	Pr.diagonalize()
#####################################################################################################################################################################

params.pretty_print()
magCalcBohrPowder = []

if LS_on:
    XCalcPowder = Pr.susceptibility(Temps = TempX, Field = .1, deltaField = .0001)
    for i in FieldM:
        magCalcBohrPowder.append((Pr.magnetization(Temp = TempM, Field = [i, 0, 0])[0] + Pr.magnetization(Temp = TempM, Field = [0, i, 0])[1] + Pr.magnetization(Temp = TempM, Field = [0, 0, i])[2])/3)

else:
    XCalcPowder = Pr.susceptibility(Temps = TempX, Field = .1, deltaField = .0001, ion = ion)
    for i in FieldM:
        magCalcBohrPowder.append((Pr.magnetization(Temp = TempM, Field = [i, 0, 0], ion = ion)[0] + Pr.magnetization(Temp = TempM, Field = [0, i, 0], ion = ion)[1] + Pr.magnetization(Temp = TempM, Field = [0, 0, i], ion = ion)[2])/3)

XiCalcPowder = -1/XCalcPowder
magCalcBohrPowder = -1*np.array(magCalcBohrPowder)


plt.figure()
plt.plot(TempX,Xi, label = 'Measured')
plt.plot(TempX,XiCalcPowder, label = 'PCF Powder Average')
plt.xlabel('Temperature (K)')
plt.ylabel('X^-1 uB^-1 Tesla spin')
plt.legend()
plt.title('Inverse Susceptbility Applied {} Tesla '.format(.1))

plt.figure()
plt.plot(FieldM,M, label = 'Measured')
plt.plot(FieldM,magCalcBohrPowder, label = 'PCF Powder Average')
plt.xlabel('Field (T)')
plt.ylabel('X^-1 uB^-1 Tesla spin')
plt.legend()
plt.title('Magnetization at {} K '.format(TempM))

Pr.printEigenvectors()
plt.show()