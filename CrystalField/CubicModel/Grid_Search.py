import sys
sys.path.append('../../')
from JensenTools import *

#kwargs
#####################################################################################################################################################################
comp = 'Sr2PrO4' #Compound name
LS_on = False
ion = 'Ce3+'
saveDir = getSaveDir('m', comp = comp, dataType = 'grid')
tol = .01 #tolerance allowed between measured and calculated energy.
Emeas = [168, 335,385] # Our measured levels of Sr2PrO4
#####################################################################################################################################################################

# Getting correct directory
#####################################################################################################################################################################
LSDir = '/cubic_matrix_LS/' 
JDir = '/cubic_matrix_J/'
if LS_on: 
    saveDir = saveDir + LSDir
else:
    saveDir = saveDir + JDir
grid = 'J_Grid_2_levels.mat' # grid name, needed for JBasis
#####################################################################################################################################################################

# Param Searching
#####################################################################################################################################################################
print('Energies as measured by paper (meV):  ', Emeas)
# LS Basis
#####################################################################################################################################################################
if(LS_on):
    LSNames, EList, data = loadMatrix(saveDir, LS_on = LS_on) #Load in all created 800x800 grids

    for c in range(len(LSNames)):

        #Loading the x,bpf, and LS of each file.
        x = data[c]['X'][0]
        bpf = data[c]['B'][0]
        LS = data[c]['LS'][0][0]

        # plotContours(data[c],EList[c], LS_on = LS_on, LSName =LSNames[c]) #Contour plotting for 4 E levels


        #Choose which bands to look for compatibilities.
        #For since we only measure 3 magnetic modes, only search for compatibilities with energies [E1,E2,E3].
        index = [1,2]
        Eindex = []
        EListindex = []
        for i in index:
            Eindex.append(Emeas[i-1])
            EListindex.append(EList[c][i-1])
            
        #Function call that searches for compatible (x,bpf) coordinates.
        coords = paramFinder(data[c],EListindex,Eindex,tol,comp,LSName = LSNames[c],LS_on = LS_on)

        #Printing results
        if len(coords) !=0:
            for j in [coords[0]]:
                print('!!! Compatibilities Found !!!')
                print('With x = ', x[j[0]], ' and bpf = ', bpf[j[1]])
                count = 1
                for i in EList[c]:
                    print('E%i = '%count, data[c][i][j[0]][j[1]], 'meV')
                    count += 1
                print()
        else:
            print('No compatibilities found')

        
        #If there is a compatibility then print an example of the matrix generated by PCF with cubic constraints.
        if(len(coords) != 0):
            print('\nFor ', LSNames[c])
            xind,bind = coords[0][0], coords[0][1]

            print('\nFor ', comp, ' at x[%i] = %.4f and bpf[%i] = %.4f'%(xind,x[xind],bind,bpf[bind]))
            print('Using these values lets construct the CF Hamiltonian\n')
            printPCFEigens(x[xind],bpf[bind],LS = LS, LS_on = LS_on)

# J Basis
#####################################################################################################################################################################
else:
    EList, data = loadMatrix(saveDir, LS_on = LS_on, grid = grid)
    print(EList)
    #Loading the x,bpf, and LS of each file.
    x = data['X'][0]
    bpf = data['B'][0]
    print("Length of x = {}, length of bpf = {}".format(len(x),len(bpf)))
    # plotContours(data,EList, LS_on = LS_on) #Contour plotting for 4 E levels


    #Choose which bands to look for compatibilities.
    #For since we only measure 3 magnetic modes, only search for compatibilities with energies [E1,E2,E3].
    index = [1]
    Eindex = []
    EListindex = []
    for i in index:
        Eindex.append(Emeas[i-1])
        EListindex.append(EList[i-1])

    coords = paramFinder(data,EListindex,Eindex,tol,comp,LS_on = LS_on)

    #Printing results
    if len(coords) !=0:
        for j in [coords[0]]:
            print('!!! Compatibilities Found !!!')
            print('With x = ', x[j[0]], ' and bpf = ', bpf[j[1]])
            count = 1
            for i in EList:
                print('E%i = '%count, data[i][j[0]][j[1]], 'meV')
                count += 1
            print()
    else:
        print('No compatibilities found')

    
    #If there is a compatibility then print an example of the matrix generated by PCF with cubic constraints.
    if(len(coords) != 0):
        print('\nFor J Basis')
        xind,bind = coords[0][0], coords[0][1]

        print('\nFor ', comp, ' at x[%i] = %.4f and bpf[%i] = %.4f'%(xind,x[xind],bind,bpf[bind]))
        print('Using these values lets construct the CF Hamiltonian\n')
        printPCFEigens(x[xind],bpf[bind], LS_on = LS_on, ion = ion)
#####################################################################################################################################################################

plt.show()
