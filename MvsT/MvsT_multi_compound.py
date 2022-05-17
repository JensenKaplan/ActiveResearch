import sys
sys.path.append('..')
from JensenTools import *
from matplotlib import rcParams
from matplotlib import patches

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
comp = ['ErOI','ErOBr']
who = 'PPMS'
# comp = 'Li8PrO6'
# who = 'MPMS'
dataType = 'MT'
saveDirDict = {}
for j in comp:
    saveDirDict[j] = getSaveDir('m', comp = j, dataType = dataType)
molweightDict = {}
for j in comp:
    molweightDict[j] = molweight[j]
per = 'mol'

fit = False

# The S,L,J values are as follows for the Pr4+ ion
S = 0.5
L = 3
J = 5./2

massErr = .00005

# Inverse Curie Law for LMFIT
# Constant calculation function
#####################################################################################################################################################################
# This is the inverse of the Curie-Weiss Law.
# We will be creating an LMFIT model with this and fitting for:
## Curie's constant 'c', Weiss' constant 'wc'
def Curiei(t,c,wc):
    return(t-wc)/c

# Takes in Curie's constant and the system's total angular momentum J value
# Returns effective moment (bohr magnetons) and effective g-factor
def calcConstants(c,J):
    ueff = np.sqrt(8*c)
    gj = ueff/np.sqrt(J*(J+1))
    return ueff, gj
#####################################################################################################################################################################

compounds = {}

for j in comp:
#####################################################################################################################################################################
    runs = []
    for i in os.listdir(saveDirDict[j]):
        if i.endswith('.DAT') or i.endswith('.dat'):
            runs.append(i)       


    # print(data.keys())
    temp = [] # temporary list for sorting
    for i in runs:
        temporary = i.split('_')[-2][:-1].replace('p','.')
        temporary = temporary.replace('P','.')
        temp.append(temporary) # this creates a list of just temperatures as read by the filename   


    temp = np.argsort([float(i) for i in temp]) # Sort by temperature
    runs = [runs[i] for i in temp] # newly sorted listed

    data = {}
    for i in runs:
        fieldStr = i.split('_')[-2][:-1].replace('p','.')
        fieldStr = fieldStr.replace('P','.')
        M,H,T,MErr,samplemass,measType = getData(i,saveDirDict[j], who = who, dataType = dataType)
        data[j + ' ' + measType + ' ' + fieldStr + 'T'] = [M,H,T,MErr,samplemass]
    print(data.keys())
    #####################################################################################################################################################################
    for i in runs:
        fieldStr =i.split('_')[-2][:-1].replace('p','.')
        fieldStr =fieldStr.replace('P','.')
        M,H,T,MErr,mass,measType = getData(i,saveDirDict[j], who = who, dataType = dataType)
        M = normalize(M,mass,molweightDict[j], per)
        Merr = normalize(M,mass,molweightDict[j], per)
        data[j + ' ' + measType + ' ' + fieldStr + 'T'] = [M,H,T,MErr,mass]


    byField = {}
    #Choose here
    for i in data.keys():
        name = i
        # print(name)
        fieldStr = i.split('_')[0]
        # print(fieldStr)
        M,H,T,MErr,samplemass = data[name]
        MBohr = emuToBohr2(M)
        HTes = oeToTesla(H)
        #####################################################################################################################################################################

        # Calculate susceptibility (Emu/Oe).
        # Choose a temp range to fit over (K)
        #####################################################################################################################################################################
        # X,XErr,Xi,XiErr = normSusc(M,H,MErr,molweight,samplemass,massErr)
        X = M/H
        Xi = 1 / X
        XT = X*T
        XiErr = POEXi(M,MErr,H,samplemass,.0005,j,per)
        if fieldStr in byField.keys():
            byField[fieldStr].append([name,T,X,Xi,XT,XiErr])
        else:
            byField[fieldStr] = [name,T,X,Xi,XT,XiErr]
        XBohr = MBohr/HTes
        XiBohr = 1/XBohr
    compounds[j] = byField
    #####################################################################################################################################################################

    # print(byField.keys())

    if fit:
        name,T, X, Xi, XT,XiErr = byField['3T'][0],byField['3T'][1],byField['3T'][2],byField['3T'][3],byField['3T'][4],byField['3T'][5],
        tr = [1,20] #temprange = [low,high]
        newT = []
        newXi = []
        newErr = []
        for i in range(len(T)):
            if (T[i] >= tr[0] and T[i]<= tr[1]):
                newT.append(T[i])
                newXi.append(Xi[i])
                # newErr.append(XiErr[i])

        cmodeli =  Model(Curiei, independent_vars = ['t'])
        params = cmodeli.make_params()
        params['wc'].set(value = -10)
        params['c'].set(value = 10)   

        resulti = cmodeli.fit(newXi, params, t = newT) #fit
        # resulti = cmodeli.fit(newXi, params, t = newT, weights = newErr) #fit

        fullLine = []
        for i in T:
            fullLine.append(Curiei(i,resulti.params['c'].value,resulti.params['wc'].value))


#####################################################################################################################################################################

# plt.figure()
# plt.errorbar(T,Xi,yerr = XiErr,label = 'Measured Data', marker = 'o', color = 'magenta', linestyle = 'none')
# if fit:
#     plt.plot(T,fullLine,'black', linestyle = '-', label = 'Curie Weiss Fit', linewidth = 4)
# # ax.set_title("{} {}".format(comp,name), fontsize = 20)
# plt.xlabel('Temperature (K)', fontweight = 'bold')
# plt.ylabel('X^-1 (emu ^-1 Oe)', fontweight = 'bold')
# plt.legend(fontsize = 30)
# plt.title(comp)
# # plt.show()
# if fit:
#     print('Curie: ', resulti.params['c'].value )
#     ueff, gj = calcConstants(resulti.params['c'].value, J =5./2)
#     print('ueff = {}, gj = {}'.format(ueff,gj))
#     resulti.params.pretty_print()

print()
print(compounds.keys())
print()

print()
print(compounds['ErOI'].keys())
print()

print()
print(compounds['ErOI']['ErOI FC 0.1T'])
print()

plt.figure()
for j in compounds.keys():
    for i in compounds[j].keys():
        print(i)
        if (i == 'ErOI ZFC 0.1T'):
            continue

        plt.errorbar(compounds[j][i][1],compounds[j][i][3], compounds[j][i][5], label = compounds[j][i][0], marker = 'o', linestyle = 'none')
        # print(byField[i,2])

plt.title("{} {}".format(comp, 'X^-1'), fontsize = 20)
plt.xlabel('Temperature (K)', fontsize = 13)
plt.ylabel('X^-1 (emu ^-1 Oe)', fontsize = 13)
plt.legend(fontsize = 30) 
plt.show()


# plt.figure()
# for i in byField.keys():
#     plt.plot(byField[i][1],byField[i][3],label = byField[i][0])
#     # print(byField[i,2])
#     plt.title("{} {}".format(comp, 'X^-1'), fontsize = 20)
#     plt.xlabel('Temperature (K)', fontsize = 13)
#     plt.ylabel('X^-1 (emu ^-1 Oe)', fontsize = 13)
#     plt.legend(fontsize = 30)    
#     plt.title(comp)

# plt.show()

# plt.figure()
# for i in byField.keys():
#     plt.plot(byField[i][1],byField[i][2],label = byField[i][0])
#     plt.title("{} {}".format(comp, 'X'), fontsize = 20)
#     plt.xlabel('Temperature (K)', fontsize = 13)
#     plt.ylabel('X (emu Oe^-1)', fontsize = 13)
#     plt.legend(fontsize = 30)   
#     plt.title(comp)

# plt.show()

# plt.figure()
# for i in byField.keys():
#     plt.plot(byField[i][1],byField[i][4],label = byField[i][0], linewidth = 2)
#     plt.title("{} {}".format(comp, 'X*T(T)'))
#     plt.xlabel('Temperature (K)')
#     plt.ylabel('X*T (emu K Oe^-1)')
#     plt.legend(fontsize = 30)   
#     plt.title(comp)

# plt.show()

# fig,ax = plt.subplots()

# ax.plot(T,Xi, color ='b',label = 'Measured X^-1')
# ax.set_xlabel('Temperature (K)', fontsize = 13)
# ax.set_ylabel('X^-1 (emu ^-1 Oe)', fontsize = 13, color = 'b')
# ax2 = ax.twinx()
# ax2.plot(T,X, color = 'orange', label = 'Measured X')
# ax2.set_ylabel('X (emu Oe^-1)', fontsize = 13, color = 'orange')

# ax2.plot(T, T*X, color = 'red', label = 'X*T')



# plt.figure()
# plt.plot(T,Xi,label = 'Measured X^-1')
# plt.plot(T,X, label = 'Measured X')
# plt.plot(T, T*X, label = 'X*T')
# # plt.errorbar(T,Xi,yerr = XiErr,label = 'Measured 1/X')
# if fit:
#     plt.plot(T,fullLine,'orange', linestyle = '--', label = 'Fitted X^-1')
# plt.title("{} {}".format(comp,name), fontsize = 20)
# plt.xlabel('Temperature (K)', fontsize = 13)
# plt.ylabel('X^-1 (emu ^-1 Oe)', fontsize = 13)
# plt.legend()

# plt.figure()
# plt.plot(T,XiBohr,label = 'Measured X^-1')
# if fit:
#     plt.plot(T,fullLine,'orange', linestyle = '--', label = 'Fitted 1/X')
# plt.title("{} {}".format(comp,name), fontsize = 20)
# plt.xlabel('Temperature (K)', fontsize = 13)
# plt.ylabel('X^-1 (uB^-1 T)', fontsize = 13)
# plt.legend()
# if fit:
#     print('The Weiss constant = {:.2f} K\nThe Curie constant = {:.3f}'.format(resulti.params['wc'].value,resulti.params['c'].value))
#     ueff, gj = calcConstants(resulti.params['c'].value,J)

#     print('Effective moment for {:} is {:.3f} bohr magnetons, with J={} -> gj factor = {:.3f}'.format(comp,ueff,J,gj))
# plt.show()

#####################################################################################################################################################################