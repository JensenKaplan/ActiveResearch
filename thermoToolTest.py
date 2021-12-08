from JensenTools import *


# Define important things
#####################################################################################################################################################################
comp = 'Sr2PrO4'
ion = 'Ce3+'
who = 'Arun'
LS_on = False
per = 'spin'
molweight = molweight[comp]
# The L,S values are as follows for the Pr4+ ion
L = 3
S = 0.5
#####################################################################################################################################################################

# Best Fit LS
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
	# # Red Chi = ~.01
	B40  =  -0.5572886105373519
	B60  =  0.4673
	B44  =  -3.0946858584804335
	B64  =  -9.8133
	B20  =  12.606195720794622
#####################################################################################################################################################################

# Create Pr CFLevels object using the best fit coefficients.
#####################################################################################################################################################################
stev = { 'B20' :B20, 'B40': B40, 'B44' : B44, 'B60': B60, 'B64' : B64 }
#Create the CFLevels object and diagonalize it
if LS_on:
	Pr = cef.LS_CFLevels.Bdict(Bdict = stev, L = L, S = S, SpinOrbitCoupling=LS)
	Pr.diagonalize()
else:
	Pr = cef.CFLevels.Bdict(Bdict = stev, ion = ion)
	Pr.diagonalize()
#####################################################################################################################################################################

# Loading data for M vs H
#####################################################################################################################################################################
saveDir = getSaveDir('m',comp = comp) #General Directory for the project
MHDir = getSaveDir('m',comp = comp, dataType = 'MH') #MvsH data

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




# TX = np.linspace(0,300,100)
# HX = .1
# TM = 20
# HM = np.linspace(0,14,100)
# M,X = thermoDiagnostic(Pr,TX,HX,TM,HM, LS_on = LS_on, ion = ion, comp = comp )



MPowder = []
for i in HTes:
	MPowder.append(powder_average(i,10,Pr,Temp))


#Generate a magnetization curve for comparing results to experiment
magCalcBohr = []
magCalcBohrPowder = []
for i in HTes:
    if LS_on:
        magCalcBohr.append(Pr.magnetization( Temp = Temp, Field = [0, 0, i])[2])
        magCalcBohrPowder.append((Pr.magnetization(Temp = Temp, Field = [i, 0, 0])[0] + Pr.magnetization(Temp = Temp, Field = [0, i, 0])[1] + Pr.magnetization(Temp = Temp, Field = [0, 0, i])[2])/3)

    else:
        magCalcBohr.append(Pr.magnetization( Temp = Temp, Field = [0, 0, i], ion = ion)[2])
        magCalcBohrPowder.append((Pr.magnetization(Temp = Temp, Field = [i, 0, 0], ion = ion)[0] + Pr.magnetization(Temp = Temp, Field = [0, i, 0], ion = ion)[1] + Pr.magnetization(Temp = Temp, Field = [0, 0, i], ion = ion)[2])/3)

magCalcBohrPowder = -1*np.array(magCalcBohrPowder)

# print(np.size(MPowder))
# print(MPowder)
MPowder = -1*np.array(MPowder)
plt.figure()
plt.plot(HTes,MPowder, label = 'Zhilling Powder Average')
plt.plot(HTes, MBohr, label = 'Measured')
plt.plot(HTes,magCalcBohrPowder, label = 'PCF 1/3 Powder Averaged')
plt.xlabel('Field (T)')
plt.ylabel('Moment (uB)')
plt.legend()
plt.show()

