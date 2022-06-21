import sys
sys.path.append('..')
from JensenTools import *
from scipy.integrate import simps
from matplotlib import rcParams
from matplotlib import patches


avo = 6.02214076*10**23

#####################################################################################################################################################################
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
rcParams.update({'font.size': 28})
rcParams['font.weight'] = 'bold'
rcParams['axes.linewidth'] = 4
rcParams['xtick.direction'] = 'out'
rcParams['ytick.direction'] = 'out'
rcParams['xtick.top'] = False
rcParams['ytick.right'] = False
rcParams['xtick.major.size'] = 12.5
rcParams['ytick.major.size'] = 12.5
rcParams['xtick.minor.size'] = 7.5
rcParams['ytick.minor.size'] = 7.5
rcParams['xtick.major.width'] = 3
rcParams['ytick.major.width'] = 3
rcParams['xtick.minor.width'] = 3
rcParams['ytick.minor.width'] = 3
rcParams['xtick.minor.visible'] = True
rcParams['ytick.minor.visible'] = True
rcParams['legend.frameon'] = False
rcParams['legend.fontsize'] = 18
#####################################################################################################################################################################

Ei = 700
Temp = 4.82
res = 6
massErr = .00005


# Define important things
#####################################################################################################################################################################
comp = 'Sr2PrO4'
ion = 'Ce3+'

who = 'Arun'
LS_on = True
per = 'spin'
molweight = molweight[comp]
LSValue = 100.5

# The L,S values are as follows for the Pr4+ ion
L = 3
S = 0.5


Emeas = [168, 335, 385] # The measured INS magnetic modes
#####################################################################################################################################################################

#LMFIT Models
#####################################################################################################################################################################
# Fits the concatenated X^-1 and magnetization
def fullFit(B20, B40,B60, B44, B64, LS, TempX, FieldX, TempM, FieldM, **kwargs ):
    deltaField = .0001
    numlevels = 4
=======
who = 'MPMS'
LS_on = True
per = 'spin'

molweight = molweight[comp]
LSValue = 100.5

# The S,L,J values are as follows for the Pr4+ ion
S = 0.5
L = 3
J = 5./2
#####################################################################################################################################################################

#March Meeting Best Fit
#####################################################################################################################################################################
pf  =  0.31466005532972374
B20  =  10.76366288180152
B40  =  -0.06594735005189663
B60  =  0.03827417917473945
B44  =  -1.1393836672476119
B64  =  0.0017325741601411152
LS  =  64.52581063639214
#####################################################################################################################################################################


#LMFIT Models
#####################################################################################################################################################################
# Fits the concatenated X^-1 and magnetization
def lineshapeFit(pf, B20, B40,B60, B44, B64, LS, TempX, FieldX, TempM, FieldM, energy, **kwargs ):

    Stev = {} #Creating the Stevens' Coefficients dictionary and assigning values
    Stev['B20'] = B20
    Stev['B40'] = B40
    Stev['B60'] = B60
    Stev['B44'] = B44
    Stev['B64'] = B64
    
    if kwargs['LS_on']:
        Pr = cef.LS_CFLevels.Bdict(Bdict=Stev, L=3, S=0.5, SpinOrbitCoupling = LS) #Create CF_Levels obejct wtih the given coefficients.
        Pr.diagonalize()
        CalculatedSpectrum = Pr.neutronSpectrum(energy, Temp=Temp, Ei=Ei, ResFunc = lambda x: res )
        
    else:
        Pr = cef.CFLevels.Bdict(Bdict = Stev, ion = kwargs['ion'])
        Pr.diagonalize()
        CalculatedSpectrum = Pr.neutronSpectrum(energy, Temp=Temp, Ei=Ei, ResFunc = lambda x: res )

    return pf*CalculatedSpectrum

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
    
    e = kmeansSort(Pr.eigenvalues,numlevels)[:numlevels-1] #Excluding the highest mode which we did not detect in our INS runs
    M = -1*np.array(M)
    Xi = -1/X
    total = np.concatenate((e,Xi,M), axis = None)
    return total

    e = kmeansSort(Pr.eigenvalues,numlevels)[:numlevels-3]
    # print(e) #Excluding the highest mode which we did not detect in our INS runs
    e.append(e[3]/e[2]) #The aforementioned ratio
    e.append(e[2]/e[1]) #The aforementioned ratio
    # print(e)
    M = -1*np.array(M)
    Xi = -1/X
    total = np.concatenate((e,X), axis = None)
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

    total = np.concatenate((-X,M), axis = None)
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

# Fitting to eigenvalues
def energyIntensityFit(B40, B60, B44, B64, B20, LS, TempX, FieldX, TempM, FieldM, **kwargs ):

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
            e = kmeansSort(Pr.eigenvalues,numlevels)[:numlevels-2] #Excluding the highest mode which we did not detect in our INS runs

            E = Pr.neutronSpectrum(energy, Temp=Temp, Ei=Ei, ResFunc = lambda x: res )
            peak1Range = [308,358]
            peak2Range = [360,405]
            peak1New = []
            peak2New = []

            for i in range(len(energy)):
                if (energy[i] >= peak1Range[0] and energy[i]<= peak1Range[1]):
                    peak1New.append(E[i])
                if (energy[i] >= peak2Range[0] and energy[i]<= peak2Range[1]):
                    peak2New.append(E[i])

            area1 = simps(peak1New, dx=.01)
            area2 = simps(peak2New, dx=.01)

            e.append(area2/area1)
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
< HEAD
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

# Choosing 20K run
Tmh = '20K'
TempM = getTemp(MHdata[Tmh][-1], who = who)
M, FieldM, Err, mass, filename = MHdata[Tmh]
M = emuToBohr2(M)
FieldM = oeToTesla(FieldM)

total = np.concatenate((Emeas,Xi,M), axis = None)

ENorm = 1/7/len(Emeas)*np.ones(len(Emeas))
XiNorm = 3/7/len(Xi)*np.ones(len(Xi))
MNorm = 3/7/len(M)*np.ones(len(M))

error = np.concatenate((ENorm,XiNorm,MNorm),axis = None)

# Normalizes and stores data as well as plotting in Emu/Oe for all temperatures.
MHdata = {}
for i in runs: 
    M, H, MErr, mass, T = getData(i,MHDir,who = who, dataType = 'MH')
    M = normalize(M,mass,molweight,'mol')
    MErr = normalize(MErr,mass,molweight,'mol')
    MHdata[T] = [M,H,MErr,mass,i]
#####################################################################################################################################################################


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

XiErr = POEXi(Mx,MErrxEmu,Hx,mass,massErr,comp,per)

# Magnetization MvsH
# Choosing 50K run
# Tmh = '10K'
# TempM = getTemp(MHdata[Tmh][-1], who = who)
# Mm, Hm, MErrmEmu, mass, filename = MHdata[Tmh]
# MBohr = emuToBohr2(Mm)
# HTes = oeToTesla(Hm)
# MErrmBohr = emuToBohr2(MErrmEmu)

# Tmh = '15K'
# TempM = getTemp(MHdata[Tmh][-1], who = who)
# Mm, Hm, MErrmEmu, mass, filename = MHdata[Tmh]
# MBohr = emuToBohr2(Mm)
# HTes = oeToTesla(Hm)
# MErrmBohr = emuToBohr2(MErrmEmu)

# total = np.concatenate((Emeas,XBohr,M), axis = None)
# total = np.concatenate((XBohr,MBohr), axis = None)
# total = np.concatenate((Emeas,XBohr), axis = None)


# XNorm = 2/3/len(XErrBohr)*XErrBohr
# MNorm = 1/3/len(MErrmBohr)*MErrmBohr
# EMeasErrNorm = 4/5/len(EMeasErr)*np.array(EMeasErr)

# error = np.concatenate((EMeasErrNorm,XNorm),axis = None)
# error = np.concatenate((XNorm,MBohr),axis = None)
#####################################################################################################################################################################


# Make LMFIT model and fit
# Create stevens coefficients dictionary from fitted parameters
#####################################################################################################################################################################
< HEAD
susModel = Model(fullFit, independent_vars = ['TempX', 'FieldX', 'TempM', 'FieldM'])
params = susModel.make_params()

# Since we only have 4 training points, only 4 parameters can vary at once.
params['B20'].set(value = B20, vary = True)
params['B40'].set(value = B40, vary=True)
params['B60'].set(value = B60, vary=True)
params['B44'].set(value = B44, vary = True)
params['B64'].set(value = B64, vary = True)

if LS_on:
	params['LS'].set(value=LS, vary=True)
    
# Fit model to data
fitted = susModel.fit(total,params, TempX = TempX, FieldX = .1, TempM = TempM, FieldM = FieldM, LS_on = LS_on, ion = ion, weights = error)

# Create a dictionary of the fitted parameters (stevens coefficients)
stev = {'B40': fitted.params['B40'].value, 'B60': fitted.params['B60'].value, 'B44' : fitted.params['B44'].value, 'B64' : fitted.params['B64'].value, 'B20' :fitted.params['B20'].value }

myModel = Model(lineshapeFit, independent_vars = ['TempX', 'FieldX', 'TempM', 'FieldM', 'energy'])
params = myModel.make_params()


# Since we only have 4 training points, only 4 parameters can vary at once.
# params['B20'].set(value=B20, min = .25*B20, max = 1.75*B20, vary=True)
# params['B20'].set(value=0, vary=True)
# params['B40'].set(value=B40, vary = False)
# params['B44'].set(value=5*B40, min = .5*5*params['B40'], max = 1.5*5*params['B40'], vary = False)
# # params['B44'].set(value=5*B40, vary = True)
# params['B60'].set(value=B60, vary = False)
# params['B64'].set(value=-21*B64, min = .5*-21*params['B60'], max = 1.5*-21*params['B60'], vary = False)
# # params['B64'].set(value=-21*B60, vary = True)
# params['pf'].set(value = pf, vary = True)

params['B20'].set(value = B20, vary = False)
params['B40'].set(value = B40, vary = False)
params['B44'].set(value = B44, vary = False)
params['B60'].set(value = B60, vary = False)
params['B64'].set(value = B64, vary = False)
params['pf'].set(value = pf,   vary = False)

if LS_on:
	params['LS'].set(value = LS, vary = False)
# Fit model to data

f ='/Users/jensenkaplan/Dropbox (GaTech)/Jensen/Sr2PrO4/Sr2PrO4_gaussians.csv' # We need to re-open the file
df = pd.read_csv(f)
energy = df['Energy']
intensity = df['Intensity']
lineErr = df['Error']

TempM = '50K'
fitted = myModel.fit(intensity, params, TempX = TempX, FieldX = FieldX, TempM = TempM, FieldM = HTes, energy = energy, LS_on = LS_on, ion = ion, Kmeans = True, weights = lineErr)
# fitted = myModel.fit(EmeasIntensity,params, TempX = TempX, FieldX = FieldX, TempM = TempM, FieldM = HTes, LS_on = LS_on, ion = ion, numlevels = numlevels, Kmeans = True, weights = EmeasIntensityErr)
# fitted = myModel.fit(total,params, TempX = TempX, FieldX = FieldX, TempM = TempM, FieldM = HTes, LS_on = LS_on, ion = ion, numlevels = numlevels, weights = error)

# Create a dictionary of the fitted parameters (stevens coefficients)
stev = {'B20' : fitted.params['B20'].value, 'B40' : fitted.params['B40'].value, 'B44': fitted.params['B44'].value, 'B60': fitted.params['B60'].value,'B64': fitted.params['B64'].value}
# stev = {'B20' :fitted.params['B20'].value, 'B21' :fitted.params['B21'].value, 'B22' :fitted.params['B22'].value, 'B40': fitted.params['B40'].value,  'B41': fitted.params['B41'].value,  'B42': fitted.params['B42'].value,  'B43': fitted.params['B43'].value,  'B44': fitted.params['B44'].value, 'B60': fitted.params['B60'].value, 'B61': fitted.params['B61'].value, 'B62': fitted.params['B62'].value, 'B63': fitted.params['B63'].value, 'B64': fitted.params['B64'].value, 'B65': fitted.params['B65'].value, 'B66': fitted.params['B66'].value}

paramPrint(fitted.params)
fitted.params.pretty_print()

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


# Use best fit to make thermo and neutron predictions
#####################################################################################################################################################################
magCalcPowderBohr10 = []
Tmh = '10K'
TempM = getTemp(MHdata[Tmh][-1], who = who)
Mm, Hm, MErrmEmu, mass, filename = MHdata[Tmh]
MBohr10 = emuToBohr2(Mm)/avo
HTes10 = oeToTesla(Hm)
MErrmBohr10 = emuToBohr2(MErrmEmu)/avo

if LS_on:
    XCalcPowderBohr = Pr.susceptibility(Temps = TempX, Field = FieldX, deltaField = .0001)
    for i in HTes10:
        magCalcPowderBohr10.append((Pr.magnetization(Temp = TempM, Field = [i, 0, 0])[0] + Pr.magnetization(Temp = TempM, Field = [0, i, 0])[1] + Pr.magnetization(Temp = TempM, Field = [0, 0, i])[2])/3)
else:
    XCalcPowderBohr = Pr.susceptibility(Temps = TempX, Field = FieldX, deltaField = .0001, ion = ion)
    for i in HTes:
        magCalcPowderBohr.append((Pr.magnetization(Temp = TempM, Field = [i, 0, 0], ion = ion)[0] + Pr.magnetization(Temp = TempM, Field = [0, i, 0], ion = ion)[1] + Pr.magnetization(Temp = TempM, Field = [0, 0, i], ion = ion)[2])/3)

magCalcPowderBohr15 = []
Tmh = '15K'
TempM = getTemp(MHdata[Tmh][-1], who = who)
Mm, Hm, MErrmEmu, mass, filename = MHdata[Tmh]
MBohr15 = emuToBohr2(Mm)/avo
HTes15 = oeToTesla(Hm)
MErrmBohr15 = emuToBohr2(MErrmEmu)/avo

if LS_on:
    XCalcPowderBohr = Pr.susceptibility(Temps = TempX, Field = FieldX, deltaField = .0001)
    for i in HTes15:
        magCalcPowderBohr15.append((Pr.magnetization(Temp = TempM, Field = [i, 0, 0])[0] + Pr.magnetization(Temp = TempM, Field = [0, i, 0])[1] + Pr.magnetization(Temp = TempM, Field = [0, 0, i])[2])/3)

magCalcPowderBohr50 = []
Tmh = '50K'
TempM = getTemp(MHdata[Tmh][-1], who = who)
Mm, Hm, MErrmEmu, mass, filename = MHdata[Tmh]
MBohr50 = emuToBohr2(Mm)/avo
HTes50 = oeToTesla(Hm)
MErrmBohr50 = emuToBohr2(MErrmEmu)/avo

if LS_on:
    XCalcPowderBohr = Pr.susceptibility(Temps = TempX, Field = FieldX, deltaField = .0001)
    for i in HTes50:
        magCalcPowderBohr50.append((Pr.magnetization(Temp = TempM, Field = [i, 0, 0])[0] + Pr.magnetization(Temp = TempM, Field = [0, i, 0])[1] + Pr.magnetization(Temp = TempM, Field = [0, 0, i])[2])/3)


XCalcPowderEmu = bohrToEmu2(XCalcPowderBohr)/10000
XiCalcPowderEmu = -1/XCalcPowderEmu

XEmu = XEmu
XiEmu = 1/XEmu

# XiCalcPowderBohr = -1/XCalcPowderBohr
magCalcPowderBohr10 = -1*np.array(magCalcPowderBohr10)
magCalcPowderBohr15 = -1*np.array(magCalcPowderBohr15)
magCalcPowderBohr50 = -1*np.array(magCalcPowderBohr50)
# magCalcPowderEmu = bohrToEmu2(magCalcPowderBohr)*6.02214076*10**23

# Mm = Mm*avo


df = pd.DataFrame()
df['TempX'] = TempX
df['X'] = -XCalcPowderEmu
# df.to_csv('/Users/jensenkaplan/Dropbox (GaTech)/Jensen/Sr2PrO4/Sr2PrO4_fitted_X.csv')

df = pd.DataFrame()
df['Hm'] = HTes10
df['Mm'] = magCalcPowderBohr10

# df.to_csv('/Users/jensenkaplan/Dropbox (GaTech)/Jensen/Sr2PrO4/Sr2PrO4_fitted_M.csv')



#####################################################################################################################################################################
plt.figure()
plt.errorbar(TempX, XiEmu, yerr = XiErr,linestyle = 'none', marker = 'o',color='magenta',label='Data',markersize = 5)
plt.plot(TempX,XiCalcPowderEmu, linestyle = '-',color='black',label='PCF Prediction',linewidth = 4)
plt.xlabel('Temperature (K)', fontweight = 'bold')
plt.ylabel('X^-1 (emu^-1 T mol)', fontweight = 'bold')
plt.legend(fontsize = 30)


plt.figure()
plt.errorbar(TempX, XBohr*TempX, yerr = XErrEmu,linestyle = 'none', marker = 'o',color='magenta',label='Data',markersize = 5)
plt.plot(TempX,-XCalcPowderBohr*TempX,  linestyle = '-',color='black',label='PCF Prediction',linewidth = 4)
plt.xlabel('Temperature (K)', fontweight = 'bold')
plt.ylabel('X(T)*T (uB T^-1 K spin^-1)', fontweight = 'bold')
plt.legend(fontsize = 30)


plt.figure()
plt.errorbar(HTes10, MBohr10, yerr = MErrmBohr10, color = 'blue',linestyle = 'none', marker = 'o',label='10K Data',markersize = 5)
plt.plot(HTes10,magCalcPowderBohr10,linestyle = '-',color='blue',label='10K PCF Prediction',linewidth = 4)
plt.xlabel('Field (T)', fontweight = 'bold')
plt.ylabel('M (uB per Pr4+ Ion)', fontweight = 'bold')
plt.errorbar(HTes15, MBohr15, yerr = MErrmBohr15, color = 'orange',linestyle = 'none', marker = 'o',label='15K Data',markersize = 5)
plt.plot(HTes15,magCalcPowderBohr15, color = 'orange',linestyle = '-',label='15K PCF Prediction',linewidth = 4)
plt.errorbar(HTes50, MBohr50, yerr = MErrmBohr50, color = 'red', linestyle = 'none', marker = 'o',label='50K Data',markersize = 5)
plt.plot(HTes50,magCalcPowderBohr50, color = 'red',linestyle = '-',label='50K PCF Prediction',linewidth = 4)
plt.legend(fontsize = 25)

Pr.printEigenvectors()

Pr.printLaTexEigenvectors()
# plt.show()


# energy = np.linspace(.01,Ei,1000)

CalculatedSpectrum = fitted.params['pf'].value*Pr.neutronSpectrum(energy, Temp=Temp, Ei=Ei, ResFunc = lambda x: res )
# ResFunc = lambda x: 9 if (energy < 200) else 21
plt.figure()
plt.plot(energy, CalculatedSpectrum, linestyle = '--', label = 'Fitted')
plt.plot(energy, intensity,  label = 'Generated Gaussian')
plt.legend()
plt.ylabel('Intensity (arb. units)')
plt.xlabel('Energy (meV)')
plt.title('PCF Spectrum: Ei = {}meV, Temp = {}K'.format(Ei,Temp))


print(fitted.params['B44']/fitted.params['B40'])
print(fitted.params['B64']/fitted.params['B60'])


wyb = cef.StevensToWybourne('Ce3+',stev, LS=True)
print(wyb)
print(Pr.gtensor())

plt.show()

