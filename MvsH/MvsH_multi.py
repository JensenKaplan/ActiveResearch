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
# who = 'Arun'
# comp = 'Li8PrO6'
who = 'PPMS'
dataType =  'MH'

saveDirDict = {}
for j in comp:
    saveDirDict[j] = getSaveDir('m', comp = j, dataType = dataType)
molweightDict = {}
for j in comp:
    molweightDict[j] = molweight[j]
per = 'spin'

# MHDir = getSaveDir('m', comp = comp, dataType = 'MH')

fit = False

# The S,L,J values are as follows for the Pr4+ ion
S = 0.5
L = 3
J = 5./2

massErr = .00005

#####################################################################################################################################################################

compounds = {}

# Get all runs from data directory
# Sort the runs by temperature (this makes for prettier plotting).
#####################################################################################################################################################################
for j in comp:
    runs = [] #A list of all the data file names
    for i in os.listdir(saveDirDict[j]):
        if i.endswith('.DAT') or i.endswith('.dat'): #This was a safeguard against a situation arising at an earlier implementation of my code.
            runs.append(i)

    temp = [] # temporary list for sorting
    for i in runs:
        temp.append(getTemp(i, who = who)) # this creates a list of just temperatures as read by the filename   
    temp = np.argsort([int(i) for i in temp]) # Sort by temperature
    runs = [runs[i] for i in temp] # newly sorted listed
    #####################################################################################################################################################################

    ## Get mass and temp from one of the filenames. Doesn't matter since they should all have the same name/mass.
    mass = getMass(runs[0], who = who)
    T = getTemp(runs[0], who = who)

    # Normalize measurement to spin^-1
    # Plot the data in both Emu/Oe and uB/T
    #####################################################################################################################################################################
    MHdata = {} #Data will be normalized and stored in Emu/Oe. Conversions to uB/T are done as needed.

    # Normalizes and stores data as well as plotting in Emu/Oe for all temperatures.
    # plt.figure()
    for i in runs:
        M, H, Err, mass, T = getData(i,saveDirDict[j],who = who, dataType = 'MH')
        M = normalize(M,mass,molweight[j],per)
        Err = normalize(Err,mass,molweight[j], per)
        MHdata[T] = [M,H,Err,mass,i]
    #     plt.plot(H, M, label = T)
    # plt.title('{} Magnetization'.format(comp))
    # plt.ylabel('Moment (emu {}^-1)'.format(per))
    # plt.xlabel('Field (Oe)')
    # plt.legend()
    compounds[j] = MHdata
#####################################################################################################################################################################


#Plot in uB/T for all temperatures for all compounds
#####################################################################################################################################################################
plt.figure()
for j in compounds.keys():
    for i in compounds[j].keys():
        # if j == 'ErOI':
        # if i == 3.0:
        M,H,Err,Mass,T = compounds[j][i]
        M = emuToBohr2(M)
        H = oeToTesla(H)
        Err = emuToBohr2(Err)
        plt.errorbar(H,M,yerr = Err, label = j + ' ' + str(i) + 'K')
plt.title('Magnetization')
plt.ylabel('Moment (\N{GREEK SMALL LETTER MU}B {}^-1)'.format(per))
plt.xlabel('Field (T)')
plt.legend()
plt.show()
#####################################################################################################################################################################



if fit:
# Choose a run to find the saturation magnetization
# Perform a linear fit over a chosen field range (Tesla)
#####################################################################################################################################################################
    temp = '20K' #Select a temperature to analyze
    curRun = MHdata[temp] #loading the data from my current run

    # Converto to uB/T
    M = curRun[0]
    H = curRun[1]
    MBohr = emuToBohr2(curRun[0])
    HTes = oeToTesla(curRun[1])
    ErrBohr = emuToBohr2(curRun[2])


    # Choose the field range to fit over (Tesla)
    fieldRange = [13.5,14]
    newH = []
    newM = []
    newErr = []
    for i in range(len(M)):
        if (H[i] >= fieldRange[0] and H[i] <= fieldRange[1]):
            newH.append(HTes[i])
            newM.append(MBohr[i])
            newErr.append(ErrBohr[i])

    # Create LMFIT Linear Model
    linModel = LinearModel()
    params = linModel.guess(newM, x = newH)
    fitted = linModel.fit(newM, x = newH, weights = newErr)

    #Drawing a full line from the fitted results
    MLine = []
    for i in H:
        MLine.append(fitted.params['slope'].value*i + fitted.params['intercept'].value)

    # Plot the data and the fit
    plt.figure()
    plt.plot(H,M, label = temp)
    plt.xlabel('Field (Oe)')
    plt.ylabel('Magnetization (emu mol^-1)')
    plt.legend(fontsize = 30)
    plt.title(comp)

    plt.figure()
    # Plot the data and the fit
    plt.plot(HTes,MBohr, label = temp)
    if fit:
        plt.plot(HTes,MLine, linestyle = '--', label = 'Fitted')
    plt.xlabel('Field (T)')
    plt.ylabel('Magnetization (uB)')
    plt.legend(fontsize = 30)
    plt.title(comp)
    # plt.show()

if fit:
    print('Saturation magnetization =  {:.3f} uB'.format(fitted.params['intercept'].value))
#####################################################################################################################################################################
