import sys
sys.path.append('..')
from JensenTools import *

# Define important things
#####################################################################################################################################################################
comp = 'Sr2PrO4'
ion = 'Ce3+'
who = 'MPMS'
LS_on = True
per = 'spin'
molweight = molweight[comp]
LSValue = 100.5

# The S,L,J values are as follows for the Pr4+ ion
S = 0.5
L = 3
J = 5./2

Emeas = [0, 168, 336, 384.9, 384.9/336,336/168] # The measured INS magnetic modes
EMeasErr = [0, 0.1231, 0.3966, 0.993, 0.993/.3966, 0.3966/0.1231]
# Emeas = [0,168,336,384.9]
# EMeasErr = [0, 0.1231, 0.3966, 0.993]
numlevels = 7
# Emeas = [168, 335, 385] # The measured INS magnetic modes
# numlevels = 4
#####################################################################################################################################################################

#LMFIT Models
#####################################################################################################################################################################
# Fits the concatenated X^-1 and magnetization
def fullFit(B20, B40,B60, B44, B64, LS, TempX, FieldX, TempM, FieldM, **kwargs ):
    deltaField = .0001
    numlevels = kwargs['numlevels']

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
    
    e = kmeansSort(Pr.eigenvalues,numlevels)[:numlevels-3]
    # print(e) #Excluding the highest mode which we did not detect in our INS runs
    e.append(e[3]/e[2]) #The aforementioned ratio
    e.append(e[2]/e[1]) #The aforementioned ratio
    # print(e)
    M = -1*np.array(M)
    Xi = -1/X
    total = np.concatenate((e,X), axis = None)
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
    total = np.concatenate((-X,M), axis = None)
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
            e = kmeansSort(Pr.eigenvalues,numlevels)[:numlevels-3] #Excluding the highest mode which we did not detect in our INS runs
            e.append(e[3]/e[2]) #The aforementioned ratio
            e.append(e[2]/e[1]) #The aforementioned ratio
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

# Simple function for pritting my parameters after fitting so that I can copy paste the values from output for further iterations.
def paramPrint(fittedparams):
    print()
    for i in fittedparams:
        # print(i, ' = ', i.value)
        print(i, ' = ',fittedparams[i].value )
#####################################################################################################################################################################


# Best Fit Params from [E0,E1,E2,E3,E3/E2,E2/E1]

#####################################################################################################################################################################
B40  =  -0.5012606800352563
B60  =  0.035606510173406256
B44  =  -2.3362611638838318
B64  =  0.16051908499496623
B20  =  2.786665787837128
LS  =  78.52277597643308
#####################################################################################################################################################################

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
    M, H, MErr, mass, T = getData(i,MHDir,who = who, dataType = 'MH')
    M = normalize(M,mass,molweight,per)
    MErr = normalize(MErr,mass,molweight,per)
    MHdata[T] = [M,H,MErr,mass,i]
#####################################################################################################################################################################

# Susceptibility MvsT
# Either 'ZFC' or 'FC' (usually ZFC)
Mx,Hx,TempX,MErrxEmu,mass = MTdata['3T_ZFC']
FieldX = 3.
MBohr = emuToBohr2(Mx)
HTes = oeToTesla(Hx)
MErrBohr = emuToBohr2(MErrxEmu)

XEmu = Mx/Hx
XErrEmu = MErrxEmu/Hx
XiEmu = 1/XEmu

XBohr = MBohr/HTes
XErrBohr = MErrBohr/HTes
XiBohr = 1/XBohr

# Magnetization MvsH
# Choosing 50K run
Tmh = '50K'
TempM = getTemp(MHdata[Tmh][-1], who = who)
Mm, Hm, MErrmEmu, mass, filename = MHdata[Tmh]
MBohr = emuToBohr2(Mm)
HTes = oeToTesla(Hm)
MErrmBohr = emuToBohr2(MErrmEmu)

# total = np.concatenate((Emeas,XBohr,M), axis = None)
total = np.concatenate((XBohr,MBohr), axis = None)
# total = np.concatenate((Emeas,XBohr), axis = None)


XNorm = 2/3/len(XErrBohr)*XErrBohr
MNorm = 1/3/len(MErrmBohr)*MErrmBohr
# EMeasErrNorm = 4/5/len(EMeasErr)*np.array(EMeasErr)


print("X data length {}. M data length {}. Total data length {}".format(len(XBohr),len(MBohr), len(total)))
# print("X error length {}. M error length {}. Total error length {}".format(len(XNorm),len(MNorm), len(XNorm) + len(MNorm)))


# error = np.concatenate((EMeasErrNorm,XNorm),axis = None)
error = np.concatenate((XNorm,MBohr),axis = None)


# Make LMFIT model and fit
# Create stevens coefficients dictionary from fitted parameters
#####################################################################################################################################################################
myModel = Model(energyFit, independent_vars = ['TempX', 'FieldX', 'TempM', 'FieldM'])
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
	params['LS'].set(value=LS, vary=True)
# Fit model to data
fitted = myModel.fit(Emeas,params, TempX = TempX, FieldX = FieldX, TempM = TempM, FieldM = HTes, LS_on = LS_on, ion = ion, numlevels = numlevels, Kmeans = True, weights = EMeasErr)
# fitted = myModel.fit(total,params, TempX = TempX, FieldX = FieldX, TempM = TempM, FieldM = HTes, LS_on = LS_on, ion = ion, numlevels = numlevels, weights = error)

# Create a dictionary of the fitted parameters (stevens coefficients)
stev = {'B20' :fitted.params['B20'].value, 'B40': fitted.params['B40'].value, 'B44': fitted.params['B44'].value, 'B60': fitted.params['B60'].value,'B64': fitted.params['B64'].value}
# stev = {'B20' :fitted.params['B20'].value, 'B21' :fitted.params['B21'].value, 'B22' :fitted.params['B22'].value, 'B40': fitted.params['B40'].value,  'B41': fitted.params['B41'].value,  'B42': fitted.params['B42'].value,  'B43': fitted.params['B43'].value,  'B44': fitted.params['B44'].value, 'B60': fitted.params['B60'].value, 'B61': fitted.params['B61'].value, 'B62': fitted.params['B62'].value, 'B63': fitted.params['B63'].value, 'B64': fitted.params['B64'].value, 'B65': fitted.params['B65'].value, 'B66': fitted.params['B66'].value}


paramPrint(fitted.params)

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
    XCalcPowderBohr = Pr.susceptibility(Temps = TempX, Field = FieldX, deltaField = .0001)
    for i in HTes:
        magCalcPowderBohr.append((Pr.magnetization(Temp = TempM, Field = [i, 0, 0])[0] + Pr.magnetization(Temp = TempM, Field = [0, i, 0])[1] + Pr.magnetization(Temp = TempM, Field = [0, 0, i])[2])/3)

else:
    XCalcPowderBohr = Pr.susceptibility(Temps = TempX, Field = FieldX, deltaField = .0001, ion = ion)
    for i in HTes:
        magCalcPowderBohr.append((Pr.magnetization(Temp = TempM, Field = [i, 0, 0], ion = ion)[0] + Pr.magnetization(Temp = TempM, Field = [0, i, 0], ion = ion)[1] + Pr.magnetization(Temp = TempM, Field = [0, 0, i], ion = ion)[2])/3)


XCalcPowderEmu = bohrToEmu2(XCalcPowderBohr)*6.02214076*10**23/10000
XiCalcPowderEmu = -1/XCalcPowderEmu

XEmu = XEmu*6.02214076*10**23
XiEmu = 1/XEmu

# XiCalcPowderBohr = -1/XCalcPowderBohr
magCalcPowderBohr = -1*np.array(magCalcPowderBohr)
magCalcPowderEmu = bohrToEmu2(magCalcPowderBohr)*6.02214076*10**23

Mm = Mm*6.02214076*10**23

plt.figure()
plt.plot(TempX, XiEmu, label = 'Measured')
plt.plot(TempX,XiCalcPowderEmu, label = 'PCF Powder Average')
plt.xlabel('Temperature (K)')
plt.ylabel('X^-1 (emu^-1 T mol)')
plt.legend()
plt.title('{} Inverse Susceptbility Applied {} Tesla '.format(comp,FieldX))

plt.figure()
plt.errorbar(Hm, Mm, yerr =MErrmEmu, label = 'Measured')
plt.plot(Hm,magCalcPowderEmu, label = 'PCF Powder Average')
plt.xlabel('Field (oe)')
plt.ylabel('M (emu mol^-1)')
plt.legend()
plt.title('{} Magnetization at {} K '.format(comp,TempM))

Pr.printEigenvectors()

Pr.printLaTexEigenvectors()
plt.show()


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
# # From GridSearch For LS
# #####################################################################################################################################################################
# if LS_on:
#   LS = LSValue
#   x =  0.03629536921151444
#   bpf = -0.6570713391739674
#   # Assigning the coefficients from grid search
#   # Enforcing cubic constraints as a start
#   # and including the B20 term which is needed for tetragonal symmetry
#   B40 = bpf
#   B60 = x*bpf
#   B44 = 5*B40
#   B64 = -21*B60
#   B20 = 0
# #####################################################################################################################################################################

# # From GridSearch For J
# #####################################################################################################################################################################
# if not LS_on:
#   x = -1.0000
#   bpf = -0.4673
#   # Assigning the coefficients from grid search
#   # Enforcing cubic constraints as a start
#   # and including the B20 term which is needed for tetragonal symmetry    
#   B40 = bpf
#   B60 = x*bpf
#   B44 = 5*B40
#   B64 = -21*B60
#   B20 = 0
# #####################################################################################################################################################################


# #Best Fit LS
# #####################################################################################################################################################################
# if LS_on: 
#   B40  =  -0.6568663783690575
#   B60  =  -0.02328250024945387
#   LS   =  LSValue
#   B44  =  -3.1415463304732714
#   B64  =  0.504906552605772
#   B20  =  0.4858075931009187
# #####################################################################################################################################################################

#Fix B20 to different values, check g tensor 

# Best Fit J
# #####################################################################################################################################################################
# if not LS_on:
#   # Red Chi = ~5
#   # B40  =  -0.5572886105373519
#   # B60  =  0.4673
#   # B44  =  -3.0342208316734602
#   # B64  =  -9.8133
#   # B20  =  12.606195910392971

#   # # Red Chi = ~.01
#   B40  =  -0.5572886105373519
#   B60  =  0.4673
#   B44  =  -3.0946858584804335
#   B64  =  -9.8133
#   B20  =  12.606195720794622
# #####################################################################################################################################################################

