import numpy as np
import matplotlib.pyplot as plt
import os
from lmfit import Model
from lmfit.models import LinearModel, LorentzianModel, GaussianModel, VoigtModel, PolynomialModel
from matplotlib import rcParams
from matplotlib import patches

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


def labelMakerFinal(run):
    x = run.split('.')[0].split('_')
    return (x[0] + ' Ei = ' + x[1].replace("p",".") + " meV, T = " + x[2].replace("p",".") + " K" )
def labelMaker(run, energyIntegration, qIntegration):
    x = run.split('.')[0].split('_')
    return (x[0] + ' Ei = ' + x[1].replace("p",".") + " meV, T = " + x[2].replace("p",".") + " K" + " -- " + energyIntegration + " -- " + qIntegration)

runDir = '/Users/jensenkaplan/Dropbox (GaTech)/Jensen/IE_runs/'
# runDir = '/Users/jensenkaplan/Dropbox (GaTech)/Jensen/March Meeting/slice/'
dirList = os.listdir(runDir)
runList = []

for i in dirList:
    if '.ascii' in i:
        runList.append(i)
print("All of the 1D cuts in the folder will be listed with their index. So run this once and fill in magrun and nonmagrun.")
for i in range(len(runList)):
    print(runList[i], i)

# Ei = 700meV SrPr
# nonmagrun = runList[3]
# magrun = runList[2]
# Ei = 500 meV SrPr
nonmagrun = runList[0]
magrun = runList[1]
# Ei = 500 meV LiPr
nonmagrun = runList[11]
magrun = runList[18]

f = open(runDir + nonmagrun, 'r')
datanonmag = np.genfromtxt(f, comments = '#')
f.close()
f = open(runDir + magrun, 'r')
datamag = np.genfromtxt(f, comments = '#')
f.close()

enonmag = datanonmag[:,0]
signalnonmag = datanonmag[:,1]
errornonmag = datanonmag[:,2]

emag = datamag[:,0]
signalmag = datamag[:,1]
errormag = datamag[:,2]

tr = [120,700]

newE = []
newSignalMag = []
newSignalNonMag = []
newErrorMag = []
newErrorNonMag = []
for i in range(len(emag)):
    if (emag[i] >= tr[0] and emag[i]<= tr[1]):
        newE.append(emag[i])
        newSignalMag.append(signalmag[i])
        newErrorMag.append(errormag[i])
        newSignalNonMag.append(signalnonmag[i])
        newErrorNonMag.append(errornonmag[i])
        
        
enonmag = newE
emag = newE
signalnonmag = newSignalNonMag
signalmag = newSignalMag
errornonmag = newErrorNonMag
errormag = newErrorMag

# SrPr
#####################################################################################################################################################################
bkg = PolynomialModel(prefix = 'bg_', degree = 3)
params = bkg.guess(signalmag,x = emag)
bkgparams = bkg.make_params()
bkgparams['bg_c0'].set(value = 3.5e-6, vary = False )
# bkgparams['bg_c1'].set(value =  -2.526e-07  )
# bkgparams['bg_c2'].set(value = 2.926e-10 )

g1 = GaussianModel(prefix = 'g1_')
g2 = GaussianModel(prefix = 'g2_')
params.update(g1.make_params())
params.update(g2.make_params())

myModel = bkg + g1 + g2

# Ei = 700meV
params['g1_center'].set(value = 330)
params['g1_sigma'].set(value = 5)
params['g2_center'].set(value = 385)
params['g2_sigma'].set(value = 5)
params['g2_amplitude'].set(value = .15, min = 0)


# # ## Ei = 500 meV
g3 = GaussianModel(prefix = 'g3_')
params.update(g3.make_params())
params['g1_center'].set(value = 168)
params['g1_sigma'].set(value = 7)
params['g2_center'].set(value = 220)
params['g2_sigma'].set(value = 2)
params['g3_center'].set(value = 330)
params['g3_sigma'].set(value = 5)
myModel = bkg + g1 + g2 + g3
####################################################################################################################################################################


#####################################################################################################################################################################
# bkg = PolynomialModel(prefix = 'bg_', degree = 3)
params = bkg.guess(signalmag,x = emag)
bkgparams = bkg.make_params()
bkgparams['bg_c0'].set(value = 3.5e-6, vary = False )
# bkgparams['bg_c1'].set(value =  -2.526e-07  )
# bkgparams['bg_c2'].set(value = 2.926e-10 )

g1 = GaussianModel(prefix = 'g1_')
# g2 = GaussianModel(prefix = 'g2_')
params.update(g1.make_params())
# params.update(g2.make_params())

myModel = bkg + g1

# # Ei = 700meV
# # params['g1_center'].set(value = 330)
# # params['g1_sigma'].set(value = 5)
# # params['g2_center'].set(value = 385)
# # params['g2_sigma'].set(value = 5)
# # params['g2_amplitude'].set(value = .15, min = 0)


# ## Ei = 500 meV
# g3 = GaussianModel(prefix = 'g3_')
# params.update(g3.make_params())
params['g1_center'].set(value = 270)
params['g1_sigma'].set(value = 7)
# params['g2_center'].set(value = 220)
# params['g2_sigma'].set(value = 2)
# params['g3_center'].set(value = 330)
# params['g3_sigma'].set(value = 5)
myModel = bkg + g1
#####################################################################################################################################################################

# bkgresult = bkg.fit(signalmag,bkgparams, x = emag, weights = errormag)
result = myModel.fit(signalmag,params, x = emag)
result.params.pretty_print()
# bkgresult.params.pretty_print()

comps = result.eval_components()

print(comps.keys())
plt.figure()
# plt.errorbar(emag, signalmag, yerr = errormag, fmt = 'o', linestyle = 'none',alpha=.5,label = labelMakerFinal(magrun))
# plt.errorbar(emag,signalnonmag, yerr = errornonmag, fmt = 'o', linestyle = 'none',alpha=.5,label = labelMakerFinal(nonmagmagrun))
plt.errorbar(emag, signalmag, linestyle = '--',color='magenta',marker = 'o',markersize = 5,linewidth = 4, label = 'Sr2PrO4')
plt.errorbar(enonmag, signalnonmag,linestyle = '--',color='green',marker = 'o',markersize = 5,linewidth = 4, label = 'Sr2CeO4')
plt.plot(emag, result.best_fit, linestyle = '-',color='black',label='Fitted',linewidth = 4)
# plt.title('Ei = 500 meV, T = 4.81K')
plt.xlabel('Energy (meV)', fontweight = 'bold')
plt.ylabel('Intensity (Arb. Units)', fontweight = 'bold')
# plt.plot(emag, result.init_fit)
# plt.plot(emag,comps['bg_'],emag,comps['g1_'],emag,comps['g2_'])
plt.legend(fontsize = 30)
plt.show()
print(labelMakerFinal(nonmagrun))
print(labelMakerFinal(magrun))
# print(result.params['g1_amplitude']/result.params['g3_amplitude'])
# print(result.params['g2_amplitude']/result.params['g1_amplitude'])
