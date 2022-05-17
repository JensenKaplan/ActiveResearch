import sys
sys.path.append('..')
from JensenTools import *
from scipy.integrate import simps
from matplotlib import rcParams
from matplotlib import patches

# Plot Formatting
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

# Define important things
#####################################################################################################################################################################
comp = 'Li8PrO6'
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

avo = 6.02214076*10**23
Ei = 730
Temp = 4.88
res = 9
massErr = .00005
#####################################################################################################################################################################

# Save directories
#####################################################################################################################################################################
saveDir = getSaveDir('m',comp = comp) #General Directory for the project
MTDir = getSaveDir('m',comp = comp, dataType = 'MT') #MvsT data
MHDir = getSaveDir('m',comp = comp, dataType = 'MH') #MvsT data
#####################################################################################################################################################################

#Working on LiPr Non Cubic
#####################################################################################################################################################################
# BEST FIT ALL POS G TENSOR
B20  =  -0.018083295642218118
B40  =  -0.06404195129061385
B60  =  0.05112605185045689
B44  =  -0.012768518174137693
B64  =  0.7817977687376129
LS  =  64.86031437157135
#####################################################################################################################################################################

#LiPr Cubic
#####################################################################################################################################################################
B40  =  -0.08609482451631531
B60  =  0.04517067265624785
B44  =  -0.4304741225815766
B64  =  -0.9485841257812049
LS  =  74.7127360214523
#####################################################################################################################################################################

# LiPr new coeffs
#####################################################################################################################################################################
B20 = -0.018083295642218118
B40 = -0.08609482451631531
B60 = 0.04517067265624785
B66 = 0
LS = 74.7127360214523

B20  =  -0.018740828951674882
B40  =  -0.1704440647933064
B60  =  0.054276299350620526
B66  =  0.09633732483573106
LS  =  70.33588463748987

B20  =  8.241047680844808
B40  =  -0.06387402046718269
B60  =  0.10318589428804337
B66  =  0.8609904803371305
LS  =  137.84187718821397

B20  =  9.207738801541156
B40  =  -0.07725596766051734
B60  =  0.09244137940696769
B66  =  0.8410725396428753
LS  =  137.27527299160724

B20  =  11.541484761784544
B40  =  -0.0834225822292725
B60  =  0.08128305194622058
B66  =  0.9522533873874438
LS  =  134.11855121466655

B20  =  11.541484761784544
B40  =  -0.0834225822292725
B60  =  0.08128305194622058
B66  =  0.9522533873874438
LS  =  134.11855121466655

B20  =  11.51368298045048
B40  =  -0.08669828608853813
B60  =  0.08074019367214985
B66  =  0.9499231173304348
LS  =  133.64434694201643

B20  =  11.516850634214837
B40  =  -0.08561415828325324
B60  =  0.08053823288925545
B66  =  0.9461708916174897
LS  =  133.85691322020156

B20  =  11.380975494350661
B40  =  -0.09048303760312179
B43 = 0
B4_3 = 0
B60  =  0.0804037998942952
B63 = 0
B6_3 = 0
B66  =  0.9371634330523263
B6_6 = 0
LS  =  133.43733674784053
pf = 1

#####################################################################################################################################################################


#####################################################################################################################################################################
# Energy sus fit

# B43  =  -0.003218934440639337
# B4_3  =  -0.03404382028798658
# B63  =  0.06824800910345545
# B6_3  =  0.05965159806334571
# B6_6  =  0.010291956528105042
# B20  =  11.677782967441276
# B40  =  -0.09042648876618348
# B60  =  0.07788772043563691
# B66  =  0.9249125685061678
# LS  =  136.24027589368032

# Sus fit, this used to create ACII for Arun 05/01
B43  =  0.8789232900581587
B4_3  =  -1.0661345501645305
B63  =  0.467032745043109
B6_3  =  -0.3910696658845045
B6_6  =  0.2964874388843928
B20  =  20.486467198021053
B40  =  -0.0508870118070716
B60  =  0.039534782166764655
B66  =  1.332718654111215
LS  =  61.6389049134499

# best all fit, includes energy mode from ab initio theory
B43  =  1.0789394702894177
B4_3  =  -2.0050105349422513
B63  =  0.7979019717442609
B6_3  =  -0.6966748477923841
B6_6  =  0.2583814708645062
B20  =  31.109329674116104
B40  =  -0.14507762644744152
B60  =  0.022378111026331526
B66  =  1.8590224185193926
LS  =  52.33471773829955
#####################################################################################################################################################################

#LMFIT Models
#####################################################################################################################################################################
# Fits the concatenated X^-1 and magnetization
def lineshapeFit(pf, B43, B4_3, B63, B6_3, B6_6,B20, B40,B60, B66, LS, TempX, FieldX, TempM, FieldM, energy, **kwargs ):

    Stev = {} #Creating the Stevens' Coefficients dictionary and assigning values
    Stev['B20'] = B20
    Stev['B40'] = B40
    Stev['B60'] = B60
    Stev['B66'] = B66

    Stev['B43'] = B43
    Stev['B4-3'] = B4_3
    Stev['B63'] = B63
    Stev['B6-3'] = B6_3
    Stev['B6-6'] = B6_6

    if kwargs['LS_on']:
        Pr = cef.LS_CFLevels.Bdict(Bdict=Stev, L=3, S=0.5, SpinOrbitCoupling = LS) #Create CF_Levels obejct wtih the given coefficients.
        Pr.diagonalize()
        CalculatedSpectrum = Pr.neutronSpectrum(energy, Temp=Temp, Ei=Ei, ResFunc = lambda x: res )
        
    else:
        Pr = cef.CFLevels.Bdict(Bdict = Stev, ion = kwargs['ion'])
        Pr.diagonalize()
        CalculatedSpectrum = Pr.neutronSpectrum(energy, Temp=Temp, Ei=Ei, ResFunc = lambda x: res )

    return pf*CalculatedSpectrum


# Fits X^-1
def susFit( B43, B4_3, B63, B6_3, B6_6,B20, B40,B60, B66, LS, TempX, FieldX, TempM, FieldM, **kwargs ):
    deltaField = .0001
    
    Stev = {} #Creating the Stevens' Coefficients dictionary and assigning values
    Stev['B20'] = B20
    Stev['B40'] = B40
    Stev['B60'] = B60
    Stev['B66'] = B66

    Stev['B43'] = B43
    Stev['B4-3'] = B4_3
    Stev['B63'] = B63
    Stev['B6-3'] = B6_3
    Stev['B6-6'] = B6_6
    
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

    # Fits X^-1
def energysusFit(B43, B4_3, B63, B6_3, B6_6,B20, B40,B60, B66, LS, TempX, FieldX, TempM, FieldM, **kwargs ):
    deltaField = .0001
    numlevels = kwargs['numlevels']
    
    Stev = {} #Creating the Stevens' Coefficients dictionary and assigning values
    Stev['B20'] = B20
    Stev['B40'] = B40
    Stev['B60'] = B60
    Stev['B66'] = B66

    Stev['B43'] = B43
    Stev['B4-3'] = B4_3
    Stev['B63'] = B63
    Stev['B6-3'] = B6_3
    Stev['B6-6'] = B6_6

    if kwargs['LS_on']:
        Pr = cef.LS_CFLevels.Bdict(Bdict=Stev, L=3, S=0.5, SpinOrbitCoupling = LS) #Create CF_Levels obejct wtih the given coefficients.
        Pr.diagonalize()
        X = Pr.susceptibility(Temps = TempX, Field = FieldX, deltaField = deltaField)
        if kwargs['Kmeans']:
            e = kmeansSort(Pr.eigenvalues,numlevels) #Excluding the highest mode which we did not detect in our INS runs
            print(e)
            # e.append(e[3]/e[1]) #The aforementioned ratio
            # e.append(e[2]/e[1]) #The aforementioned ratio
            total = np.concatenate((e,-1/X))
        else: 
            e = Pr.eigenvalues[0:10]
            total = np.concatenate((e,-1/X))
    else:
        Pr = cef.CFLevels.Bdict(Bdict = Stev, ion = kwargs['ion'])
        Pr.diagonalize()
        X = Pr.susceptibility(Temps = TempX, Field = FieldX, deltaField = deltaField, ion = ion)

        if kwargs['Kmeans']:
            e = kmeansSort(Pr.eigenvalues,numlevels) #Excluding the highest mode which we did not detect in our INS runs
            print(e)
            # e.append(e[3]/e[1]) #The aforementioned ratio
            # e.append(e[2]/e[1]) #The aforementioned ratio
            total = np.concatenate((e,-1/X))
        else: 
            e = Pr.eigenvalues
    return total

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
def energyFit(B43, B4_3, B63, B6_3, B6_6,B20, B40,B60, B66, LS, TempX, FieldX, TempM, FieldM, **kwargs ):
    numlevels = kwargs['numlevels']

    Stev = {} #Creating the Stevens' Coefficients dictionary and assigning values
    Stev['B20'] = B20
    Stev['B40'] = B40
    Stev['B60'] = B60
    Stev['B66'] = B66

    Stev['B43'] = B43
    Stev['B4-3'] = B4_3
    Stev['B63'] = B63
    Stev['B6-3'] = B6_3
    Stev['B6-6'] = B6_6

    if kwargs['LS_on']:
        Pr = cef.LS_CFLevels.Bdict(Bdict=Stev, L=3, S=0.5, SpinOrbitCoupling = LS) #Create CF_Levels obejct wtih the given coefficients.
        Pr.diagonalize()
        if kwargs['Kmeans']:
            e = kmeansSort(Pr.eigenvalues,numlevels) #Excluding the highest mode which we did not detect in our INS runs
        else: 
            e = Pr.eigenvalues[0:10]
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



# Loading data for MT and MH
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
    M = normalize(M,mass,molweight,'mol')
    MErr = normalize(MErr,mass,molweight,'mol')
    MHdata[T] = [M,H,MErr,mass,i]
#####################################################################################################################################################################

# Choose susceptibility and mag runs
#####################################################################################################################################################################
# Susceptibility MvsT
# Either 'ZFC' or 'FC' (usually ZFC)
# The measured ones are in emu oe
Mx,Hx,TempX,MErrxEmu,mass = MTdata['3T_ZFC']
FieldX = 3.

# Convert measurements to tesla and uB
MBohr = emuToBohr2(Mx)
HTes = oeToTesla(Hx)
MErrBohr = emuToBohr2(MErrxEmu)

# Susc in emu
XEmu = Mx/Hx
XErrEmu = MErrxEmu/Hx
XiEmu = 1/XEmu

# Susc in uB
XBohr = MBohr/HTes
XErrBohr = MErrBohr/HTes
XiBohr = 1/XBohr

#Inverse susceptibility prop of error
XiErr = POEXi(Mx,MErrxEmu,Hx,mass,massErr,comp,per)

# Magnetization MvsH
# Choosing 50K run
Tmh = '50K'
TempM = 50
# TempM = getTemp(MHdata[Tmh][-1], who = who)
Mm, Hm, MErrmEmu, mass, filename = MHdata[Tmh]

# Convert to uB Tes
MBohr = emuToBohr2(Mm)
HTes = oeToTesla(Hm)
MErrmBohr = emuToBohr2(MErrmEmu)

# total = np.concatenate((XiBohr,MBohr), axis = None)

# XNorm = 3/4/len(XErrBohr)*XErrBohr
# MNorm = 1/4/len(MErrmBohr)*MErrmBohr
# errorNorm = np.concatenate((XNorm,MNorm),axis = None)
# error = np.concatenate((XiErr,MErrmBohr),axis = None)


# total = np.concatenate((Emeas,XBohr), axis = None)
# XNorm = 2/3/len(XErrBohr)*XErrBohr
# MNorm = 1/3/len(MErrmBohr)*MErrmBohr
# EMeasErrNorm = 4/5/len(EMeasErr)*np.array(EMeasErr)
# error = np.concatenate((EMeasErrNorm,XNorm),axis = None)
#####################################################################################################################################################################


# Make LMFIT model and fit
# Create stevens coefficients dictionary from fitted parameters
#####################################################################################################################################################################
# myModel = Model(susFit, independent_vars = ['TempX', 'FieldX', 'TempM', 'FieldM'])
# myModel = Model(energyFit, independent_vars = ['TempX', 'FieldX', 'TempM', 'FieldM'])
myModel = Model(energysusFit, independent_vars = ['TempX', 'FieldX', 'TempM', 'FieldM'])
# myModel = Model(lineshapeFit, independent_vars = ['TempX', 'FieldX', 'TempM', 'FieldM', 'energy'])
params = myModel.make_params()

params['B20'].set(value = B20, vary = True)
params['B40'].set(value = B40, vary = True)
params['B60'].set(value = B60, vary = True)
params['B66'].set(value = B66, vary = True)
params['B43'].set(value = B43, vary = True)

params['B4_3'].set(value = B4_3, vary = True)
params['B63'].set(value = B63, vary = True)
params['B6_3'].set(value = B6_3, vary = True)
params['B6_6'].set(value = B6_6, vary = True)
# params['pf'].set(value = pf,  vary = True)

# params.add('B44', expr = '5.0*B40')
# params.add('B64', expr = '-21.0*B60')


if LS_on:
	params['LS'].set(value = LS, vary = True)

# Open created gaussian lineshape
f ='/Users/jensenkaplan/Dropbox (GaTech)/Jensen/Li8PrO6/Li8PrO6_gaussians.csv' # We need to re-open the file
df = pd.read_csv(f)
# Store lineshape in variables
energy = df['Energy']
intensity = df['Intensity']
lineErr = df['Error']

# plt.figure()
# plt.plot(energy,intensity)
# plt.show()

#Measured eigenvalues from INS
energies = [0,0,269,269,269,269,396,396,676,676]
# energies = [0,269,396,676]tea
x = 2/3
XiNorm = (1-x)/len(XiBohr)*np.ones(len(XiBohr))
energiesNorm = x/len(energies)*np.ones(len(energies))
errorNorm = np.concatenate((energiesNorm,XiNorm),axis = None)
eXi = np.concatenate((energies,XiBohr), axis = None)
eE = np.concatenate((energies,intensity), axis = None)

# fitted = myModel.fit(energies, params, fit_kws={'maxfev': 10000}, TempX = TempX, FieldX = FieldX, TempM = TempM, FieldM = HTes, LS_on = LS_on, energy = energy, ion = ion, Kmeans = True, numlevels = len(energies))
# fitted = myModel.fit(XiBohr, params,fit_kws={'maxfev': 10000}, TempX = TempX, FieldX = FieldX, TempM = TempM, FieldM = HTes, LS_on = LS_on, ion = ion, Kmeans = True, weights = XiErr)
fitted = myModel.fit(eXi, params,fit_kws={'maxfev': 10000}, TempX = TempX, FieldX = FieldX, TempM = TempM, FieldM = HTes, LS_on = LS_on, ion = ion, Kmeans = False, numlevels = len(energies), weights = errorNorm)
# fitted = myModel.fit(intensity, params, fit_kws={'maxfev': 10000}, TempX = TempX, FieldX = FieldX, TempM = TempM, FieldM = HTes, energy = energy, LS_on = LS_on, ion = ion, Kmeans = True, weights = lineErr)



# Create a dictionary of the fitted parameters (stevens coefficients)
# stev = {'B20' : fitted.params['B20'].value, 'B40': fitted.params['B40'].value, 'B60': fitted.params['B60'].value,'B66': fitted.params['B66'].value}
stev = {'B20' : fitted.params['B20'].value, 'B40': fitted.params['B40'].value, 'B60': fitted.params['B60'].value,'B66': fitted.params['B66'].value, 'B43': fitted.params['B43'].value, 'B4-3' :  fitted.params['B4_3'].value,'B63' :  fitted.params['B63'].value,'B6-3' :  fitted.params['B6_3'].value, 'B6-6' :  fitted.params['B6_6'].value,}

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

# Plot neutron spectrum and best fit
#####################################################################################################################################################################
# CalculatedSpectrum = fitted.params['pf'].value*Pr.neutronSpectrum(energy, Temp=Temp, Ei=Ei, ResFunc = lambda x: res )
CalculatedSpectrum = Pr.neutronSpectrum(energy, Temp=Temp, Ei=Ei, ResFunc = lambda x: res )
# ResFunc = lambda x: 9 if (energy < 200) else 21
plt.figure()
plt.plot(energy, CalculatedSpectrum, linestyle = '--', label = 'Fitted')
plt.plot(energy, intensity,  label = 'Generated Gaussian')
plt.legend()
plt.ylabel('Intensity (arb. units)')
plt.xlabel('Energy (meV)')
plt.title('PCF Spectrum: Ei = {}meV, Temp = {}K'.format(Ei,Temp))
#####################################################################################################################################################################

# Susceptibility prediction from fit
#####################################################################################################################################################################
if LS_on:
    XCalcPowderBohr = Pr.susceptibility(Temps = TempX, Field = FieldX, deltaField = .0001)
else:
    XCalcPowderBohr = Pr.susceptibility(Temps = TempX, Field = FieldX, deltaField = .0001, ion = ion)

XCalcPowderEmu = -1*bohrToEmu2(XCalcPowderBohr)/10000
XiCalcPowderEmu = 1/XCalcPowderEmu

df = pd.DataFrame()
df['Temperature (K)'] = TempX
df['Susceptibility (uB/mol/T)'] = XCalcPowderBohr
df['Susceptibility (uB/mol/T)'] = XCalcPowderEmu
# df.to_csv('/Users/jensenkaplan/Dropbox (GaTech)/Jensen/Li8PrO6/Final Stretch/XvsT/Li8PrO6_XvsT_3T.csv')
#####################################################################################################################################################################

# Use best fit to make thermo and neutron predictions
# For 10K
#####################################################################################################################################################################
magCalcPowderBohr10 = []
Tmh = '10K'
TempM = getTemp(MHdata[Tmh][-1], who = who)
Mm, Hm, MErrmEmu, mass, filename = MHdata[Tmh]
MBohr10 = emuToBohr2(Mm)/avo
HTes10 = oeToTesla(Hm)
MErrmBohr10 = emuToBohr2(MErrmEmu)/avo

if LS_on:
    for i in HTes10:
        magCalcPowderBohr10.append((Pr.magnetization(Temp = TempM, Field = [i, 0, 0])[0] + Pr.magnetization(Temp = TempM, Field = [0, i, 0])[1] + Pr.magnetization(Temp = TempM, Field = [0, 0, i])[2])/3)
else:
    for i in HTes:
        magCalcPowderBohr10.append((Pr.magnetization(Temp = TempM, Field = [i, 0, 0], ion = ion)[0] + Pr.magnetization(Temp = TempM, Field = [0, i, 0], ion = ion)[1] + Pr.magnetization(Temp = TempM, Field = [0, 0, i], ion = ion)[2])/3)

magCalcPowderBohr10 = -1*np.array(magCalcPowderBohr10)

df = pd.DataFrame()
df['Field (T)'] = HTes10
df['Magnetization (uB/mol)'] = magCalcPowderBohr10
# df.to_csv('/Users/jensenkaplan/Dropbox (GaTech)/Jensen/Li8PrO6/Final Stretch/MvsH/Li8PrO6_MvsH_10K.csv')
#####################################################################################################################################################################

#For 15K
#####################################################################################################################################################################
magCalcPowderBohr15 = []
Tmh = '15K'
TempM = getTemp(MHdata[Tmh][-1], who = who)
Mm, Hm, MErrmEmu, mass, filename = MHdata[Tmh]
MBohr15 = emuToBohr2(Mm)/avo
HTes15 = oeToTesla(Hm)
MErrmBohr15 = emuToBohr2(MErrmEmu)/avo

if LS_on:
    for i in HTes15:
        magCalcPowderBohr15.append((Pr.magnetization(Temp = TempM, Field = [i, 0, 0])[0] + Pr.magnetization(Temp = TempM, Field = [0, i, 0])[1] + Pr.magnetization(Temp = TempM, Field = [0, 0, i])[2])/3)
else:
    for i in HTes:
        magCalcPowderBohr15.append((Pr.magnetization(Temp = TempM, Field = [i, 0, 0], ion = ion)[0] + Pr.magnetization(Temp = TempM, Field = [0, i, 0], ion = ion)[1] + Pr.magnetization(Temp = TempM, Field = [0, 0, i], ion = ion)[2])/3)

magCalcPowderBohr15 = -1*np.array(magCalcPowderBohr15)

df = pd.DataFrame()
df['Field (T)'] = HTes15
df['Magnetization (uB/mol)'] = magCalcPowderBohr15
# df.to_csv('/Users/jensenkaplan/Dropbox (GaTech)/Jensen/Li8PrO6/Final Stretch/MvsH/Li8PrO6_MvsH_15K.csv')
#####################################################################################################################################################################

# For 50K
#####################################################################################################################################################################
magCalcPowderBohr50 = []
Tmh = '50K'
TempM = getTemp(MHdata[Tmh][-1], who = who)
Mm, Hm, MErrmEmu, mass, filename = MHdata[Tmh]
MBohr50 = emuToBohr2(Mm)/avo
HTes50 = oeToTesla(Hm)
MErrmBohr50 = emuToBohr2(MErrmEmu)/avo

if LS_on:
    for i in HTes50:
        magCalcPowderBohr50.append((Pr.magnetization(Temp = TempM, Field = [i, 0, 0])[0] + Pr.magnetization(Temp = TempM, Field = [0, i, 0])[1] + Pr.magnetization(Temp = TempM, Field = [0, 0, i])[2])/3)
else:
    for i in HTes:
        magCalcPowderBohr50.append((Pr.magnetization(Temp = TempM, Field = [i, 0, 0], ion = ion)[0] + Pr.magnetization(Temp = TempM, Field = [0, i, 0], ion = ion)[1] + Pr.magnetization(Temp = TempM, Field = [0, 0, i], ion = ion)[2])/3)

magCalcPowderBohr50 = -1*np.array(magCalcPowderBohr50)

df = pd.DataFrame()
df['Field (T)'] = HTes50
df['Magnetization (uB/mol)'] = magCalcPowderBohr50
# df.to_csv('/Users/jensenkaplan/Dropbox (GaTech)/Jensen/Li8PrO6/Final Stretch/MvsH/Li8PrO6_MvsH_50K.csv')
#####################################################################################################################################################################




# df = pd.DataFrame()
# df['TempX'] = TempX
# df['X'] = -XCalcPowderEmu
# # df.to_csv('/Users/jensenkaplan/Dropbox (GaTech)/Jensen/Sr2PrO4/Sr2PrO4_fitted_X.csv')
# df = pd.DataFrame()
# df['Hm'] = HTes10
# df['Mm'] = magCalcPowderBohr10
# # df.to_csv('/Users/jensenkaplan/Dropbox (GaTech)/Jensen/Sr2PrO4/Sr2PrO4_fitted_M.csv')


# Plotting thermo results
#####################################################################################################################################################################
# Xi Emu
plt.figure()
plt.errorbar(TempX, XiEmu, yerr = XiErr,linestyle = 'none', marker = 'o',color='magenta',label='Data',markersize = 5)
plt.plot(TempX,XiCalcPowderEmu, linestyle = '-',color='black',label='PCF Prediction',linewidth = 4)
plt.xlabel('Temperature (K)', fontweight = 'bold')
plt.ylabel('X^-1 (emu^-1 T mol)', fontweight = 'bold')
plt.title('LS = {:.1f}'.format(LS))
plt.legend(fontsize = 30)

# X*T Bohr
plt.figure()
plt.errorbar(TempX, XBohr*TempX, yerr = XErrEmu,linestyle = 'none', marker = 'o',color='magenta',label='Data',markersize = 5)
plt.plot(TempX,-XCalcPowderBohr*TempX,  linestyle = '-',color='black',label='PCF Prediction',linewidth = 4)
plt.xlabel('Temperature (K)', fontweight = 'bold')
plt.ylabel('X(T)*T (uB T^-1 K spin^-1)', fontweight = 'bold')
plt.title('LS = {:.1f}'.format(LS))
plt.legend(fontsize = 30)

# M vs H @ 10K,15K, 50k
plt.figure()
plt.errorbar(HTes10, MBohr10, yerr = MErrmBohr10, color = 'blue',linestyle = 'none', marker = 'o',label='10K Data',markersize = 5)
plt.plot(HTes10,magCalcPowderBohr10,linestyle = '-',color='blue',label='10K PCF Prediction',linewidth = 4)
plt.xlabel('Field (T)', fontweight = 'bold')
plt.ylabel('M (uB per Pr4+ Ion)', fontweight = 'bold')
plt.errorbar(HTes15, MBohr15, yerr = MErrmBohr15, color = 'orange',linestyle = 'none', marker = 'o',label='15K Data',markersize = 5)
plt.plot(HTes15,magCalcPowderBohr15, color = 'orange',linestyle = '-',label='15K PCF Prediction',linewidth = 4)
plt.errorbar(HTes50, MBohr50, yerr = MErrmBohr50, color = 'red', linestyle = 'none', marker = 'o',label='50K Data',markersize = 5)
plt.plot(HTes50,magCalcPowderBohr50, color = 'red',linestyle = '-',label='50K PCF Prediction',linewidth = 4)
plt.title('LS = {:.1f}'.format(LS))
plt.legend(fontsize = 25)
#####################################################################################################################################################################

Pr.printEigenvectors()
# Pr.printLaTexEigenvectors()


# wyb = cef.StevensToWybourne('Ce3+',stev, LS=True)
# print(wyb)
print(Pr.gtensor())
plt.show()