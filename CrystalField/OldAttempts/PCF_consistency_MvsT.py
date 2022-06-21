import sys
sys.path.append('..')
from JensenTools import *
import PyCrystalField as cef
import numpy as np
import matplotlib.pyplot as plt

# kwargs
#####################################################################################################################################################################
comp = 'Sr2PrO4'
ion = 'Ce3+'
who = 'Arun'
LS_on = True
Kmeans = True
LSValue = 100
numlevels = 4
# The L,S values are as follows for the Pr4+ ion
L = 3
S = 0.5
massErr = .00005
molweight = 380.15
#####################################################################################################################################################################

# Make LS_CFLevels object with best fit parameters
#####################################################################################################################################################################
B20  =  0.4858075931009187
B40  =  -0.6568663783690575
B44  =  -3.1415463304732714
B60  =  -0.02328250024945387
B64  =  0.504906552605772
LS  =  100.00007580463522

saveDir = getSaveDir('m',comp = comp) #General Directory for the project
MTDir = getSaveDir('m',comp = comp, dataType = 'MT') #MvsT data

stev = { 'B20' :B20, 'B40': B40, 'B44' : B44, 'B60': B60, 'B64' : B64 }

#Create the CFLevels object and diagonalize it
if LS_on:
	Pr = cef.LS_CFLevels.Bdict(Bdict = stev, L = L, S = S, SpinOrbitCoupling=LS)
	Pr.diagonalize()
else:
	Pr = cef.CFLevels.Bdict(Bdict = stev, ion = 'Ce3+')
	Pr.diagonalize()
#####################################################################################################################################################################

# Loading data for M vs T 
#####################################################################################################################################################################
runs = []
for i in os.listdir(MTDir):
    if i.endswith('.DAT') or i.endswith('.dat'): #This was a safeguard against a situation arising at an earlier implementation of my code.
        runs.append(i)
MTdata = {}
for i in runs:
    M,H,T,MErr,mass,measType = getData(i,MTDir, who = who, dataType = 'MT')
    M = normalize(M,mass,molweight,'spin')
    MTdata[measType] = [M,H,T,MErr,mass]

#Either 'ZFC' or 'FC'
M,H,T,MErr,mass = MTdata['FC']
#####################################################################################################################################################################

# Use PCF function to calculate Magnetization (uB/spin) at fixed Field (.1T) over temop range
# Divide this magnetization by field to get Susceptibiltiy in (uB/T/spin)
#####################################################################################################################################################################
Ha = .1
magCalcBohr = []
for i in T:
	if LS_on:
		magCalcBohr.append(Pr.magnetization( Temp = i, Field = [0, 0, Ha])[2])		
	else:
		magCalcBohr.append(Pr.magnetization( Temp = i, Field = [0, 0, Ha], ion = ion)[2])

sus = np.array(magCalcBohr)/Ha
#####################################################################################################################################################################

# Use PCF function to calculate Susceptibility (uB/Tesla/Spin) over T
#####################################################################################################################################################################
fieldT = 0.1 #To match Arun's measurements
deltaField = .00001
if LS_on:
	XCalcBohr = Pr.susceptibility(Temps = T, Field = fieldT, deltaField = deltaField)
else:
	XCalcBohr = Pr.susceptibility(Temps = T, Field = fieldT, deltaField = deltaField, ion = ion)
#####################################################################################################################################################################


PCFRatio = XCalcBohr/sus

# Plotting
#####################################################################################################################################################################
# Plot ratio of calculated PCFs
plt.figure()
plt.plot(T,PCFRatio, label = 'Ratio')
plt.xlabel('Temperature (K)')
plt.ylabel('Ratio (PCF susc)/(PCF M/H)')
plt.legend()
plt.title('Ratio of Susceptibilities {} Tesla'.format(Ha))
plt.show()


# # Plot PCF M/H
# plt.figure()
# plt.plot(T,-1*np.array(sus), label = 'PCF')
# plt.xlabel('Temperature (K)')
# plt.ylabel('Susceptibility (uB) Tesla^-1 spin^-1')
# plt.legend()
# plt.title('PCF M/H {} Tesla'.format(Ha))
# # plt.show()

# #Plot PCF Susc
# plt.figure()
# plt.plot(T,-1*XCalcBohr, label = 'PCF')
# plt.xlabel('Temperature (K)')
# plt.ylabel('Susceptibility (uB) Tesla^-1 spin^-1')
# plt.legend()
# plt.title('PCF Susceptbility {} Tesla '.format(Ha))
# # plt.show()

M = emuToBohr2(M)
H = oeToTesla(H)
# Plot experimental susc (M/H)
# plt.figure()
# plt.plot(T,M/H, label = 'Measured')
# plt.xlabel('Temperature (K)')
# plt.ylabel('Susceptibility (uB) Tesla^-1 spin^-1')
# plt.legend()
# plt.title('Measured Susceptbility {} Tesla '.format(Ha))

# Overplot all 3
plt.figure()
plt.plot(T,M/H, label = 'Measured M/H')
plt.plot(T,-1*XCalcBohr, label = 'PCF Susceptibility')
plt.plot(T,-1*np.array(sus), label = 'PCF M/H')
plt.xlabel('Temperature (K)')
plt.ylabel('Susceptibility (uB) Tesla^-1 spin^-1')
plt.legend()
plt.title('Susceptbility All Methods Plotted {} Tesla '.format(Ha))

plt.show()
#####################################################################################################################################################################
