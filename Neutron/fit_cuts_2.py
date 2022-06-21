import numpy as np
import matplotlib.pyplot as plt
import os
from lmfit import Model
from lmfit.models import LinearModel, LorentzianModel, GaussianModel, VoigtModel, PolynomialModel
from matplotlib import rcParams
from matplotlib import patches


def paramPrint(fittedparams):
    print()
    for i in fittedparams:
        # print(i, ' = ', i.value)
        print(i, ' = ',fittedparams[i].value )

#####################################################################################################################################################################
fig = plt.figure(figsize = (10,10))
gs = fig.add_gridspec(1, 1)
ax1 = plt.subplot(gs[0])
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
rcParams.update({'font.size': 15})
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
def plotSaver(run):
    x = run.split('.')[0].split('_')
    return (x[0] + '_Ei_' + x[1] + "_meV.png")
def labelMaker(run, energyIntegration, qIntegration):
    x = run.split('.')[0].split('_')
    return (x[0] + ' Ei = ' + x[1].replace("p",".") + " meV, T = " + x[2].replace("p",".") + " K" + " -- " + energyIntegration + " -- " + qIntegration)

runDir = '/Users/jensenkaplan/Dropbox (GaTech)/Jensen/Sr2PrO4/IE_runs/'
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
# nonmagrun = runList[14]
# magrun = runList[12]
# Ei = 500 meV SrPr
# nonmagrun = runList[2]
# magrun = runList[9]
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

#SrPr 700mev
tr = [225,500]
#SrPr 500mev
tr = [0,1000]
#LiPr 500meV
tr = [0,1000]


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

emag = np.array(newE)
signalmag = np.array(newSignalMag)
errormag = np.array(newErrorMag)

newE = []
for i in range(len(enonmag)):
    if (enonmag[i] >= tr[0] and enonmag[i]<= tr[1]):
        newE.append(enonmag[i])
        newSignalNonMag.append(signalnonmag[i])
        newErrorNonMag.append(errornonmag[i])
                
enonmag = np.array(newE)
signalnonmag = np.array(newSignalNonMag)
errornonmag = np.array(newErrorNonMag)

# plt.plot(enonmag,signalnonmag,emag,signalmag)
# plt.show()
# SrPr
#####################################################################################################################################################################
# bkg = PolynomialModel(prefix = 'bg_', degree = 3)
# params = bkg.guess(signalmag,x = emag)
# bkgparams = bkg.make_params()
# bkgparams['bg_c0'].set(value = 3.5e-6, vary = True )
# # bkgparams['bg_c1'].set(value =  -2.526e-07  )
# # bkgparams['bg_c2'].set(value = 2.926e-10 )

# g1 = GaussianModel(prefix = 'g1_')
# g2 = GaussianModel(prefix = 'g2_')
# params.update(g1.make_params())
# params.update(g2.make_params())

# myModel = bkg + g1 + g2

# # Ei = 700meV
# params['g1_center'].set(value = 334.5)
# params['g1_sigma'].set(value = 5)
# params['g2_center'].set(value = 388.8, vary = True)
# params['g2_sigma'].set(value = 3)
# params['g2_amplitude'].set(value = .15, min = 0)


# # # # # ## Ei = 500 meV
# # g3 = GaussianModel(prefix = 'g3_')
# # params.update(g3.make_params())
# # params['g1_center'].set(value = 168)
# # params['g1_sigma'].set(value = 7)
# # params['g2_center'].set(value = 220)
# # params['g2_sigma'].set(value = 2)
# # params['g3_center'].set(value = 330)
# # params['g3_sigma'].set(value = 5)
# # myModel = bkg + g1 + g2 + g3
# ####################################################################################################################################################################

#LiPr
# #####################################################################################################################################################################
bkg = PolynomialModel(prefix = 'bg_', degree = 3)
params = bkg.guess(signalmag,x = emag)
bkgparams = bkg.make_params()
bkgparams['bg_c0'].set(value = 3.5e-6, vary = True )
# bkgparams['bg_c1'].set(value =  -2.526e-07  )
# bkgparams['bg_c2'].set(value = 2.926e-10 )

g1 = GaussianModel(prefix = 'g1_')
# # g2 = GaussianModel(prefix = 'g2_')
params.update(g1.make_params())
# # params.update(g2.make_params())


# ## Ei = 500 meV
params['g1_center'].set(value = 270)
params['g1_sigma'].set(value = 7)
myModel = bkg + g1
#####################################################################################################################################################################

print('\n\n\n')
# bkgresult = bkg.fit(signalmag,bkgparams, x = emag, weights = errormag)
result = myModel.fit(signalmag,params, x = emag)
result.params.pretty_print()
paramPrint(result.params)
# bkgresult.params.pretty_print()

# comps = result.eval_components()

# print(comps.keys())
# plt.figure()
# plt.errorbar(emag, signalmag, yerr = errormag, fmt = 'o', linestyle = 'none',alpha=.5,label = labelMakerFinal(magrun))
# plt.errorbar(emag,signalnonmag, yerr = errornonmag, fmt = 'o', linestyle = 'none',alpha=.5,label = labelMakerFinal(nonmagmagrun))
plt.errorbar(emag, signalmag, linestyle = '--',color='magenta',marker = 'o',markersize = 5,linewidth = 4, label = labelMakerFinal(magrun))
plt.errorbar(enonmag, signalnonmag,linestyle = '--',color='green',marker = 'o',markersize = 5,linewidth = 4, label = labelMakerFinal(nonmagrun))
plt.plot(emag, result.best_fit, linestyle = '-',color='black',label='Fitted',linewidth = 4)
# plt.title('Ei = 500 meV, T = 4.81K')
plt.xlabel('Energy (meV)', fontweight = 'bold')
plt.ylabel('Intensity (Arb. Units)', fontweight = 'bold')
# plt.plot(emag, result.init_fit)
# plt.plot(emag,comps['bg_'],emag,comps['g1_'],emag,comps['g2_'])
plt.legend()
plt.savefig(runDir + 'pics/' +plotSaver(magrun))
plt.show()

# print(labelMakerFinal(nonmagrun))
# print(labelMakerFinal(magrun))
# print(result.params['g1_amplitude']/result.params['g3_amplitude'])
# print(result.params['g2_amplitude']/result.params['g1_amplitude'])

#LiPr For Paper
#####################################################################################################################################################################
bg_c0  =  0.00018505104444559132
bg_c1  =  -1.5019691570118256e-06
bg_c2  =  4.323248077115866e-09
bg_c3  =  -4.257487047838987e-12
g1_amplitude  =  0.0002769855396335671
g1_center  =  274.2196764501529
g1_sigma  =  9.101311129042994
g1_fwhm  =  21.431949472893024
g1_height  =  1.214124500101291e-05
#####################################################################################################################################################################

#SrPr For Paper
#####################################################################################################################################################################
#500 meV
bg_c0  =  0.00018455782658675292
bg_c1  =  -1.4718141769515983e-06
bg_c2  =  4.1180944701962e-09
bg_c3  =  -3.91981609197623e-12
g1_amplitude  =  0.00010025238254331347
g1_center  =  167.86803392064175
g1_sigma  =  8.28491499033998
g1_fwhm  =  19.509483517552393
g1_height  =  4.827438316378923e-06
g2_amplitude  =  0.00021057360331221768
g2_center  =  227.71534226725424
g2_sigma  =  14.189659743957016
g2_fwhm  =  33.414094558264864
g2_height  =  5.920277098993856e-06
g3_amplitude  =  0.00011474894377347192
g3_center  =  334.8342674591908
g3_sigma  =  6.138000561491445
g3_fwhm  =  14.453886482211287
g3_height  =  7.4581628158789436e-06

#700 meV
bg_c0  =  0.00016963603713744712
bg_c1  =  -1.2593727621092305e-06
bg_c2  =  3.101749296944818e-09
bg_c3  =  -2.4733272591912303e-12
g1_amplitude  =  0.00022882035627084753
g1_center  =  335.3525716811439
g1_sigma  =  13.90613734716213
g1_fwhm  =  32.74645034784433
g1_height  =  6.564448267594624e-06
g2_amplitude  =  5.277865776909252e-05
g2_center  =  382.2838548009725
g2_sigma  =  8.998765891237555
g2_fwhm  =  21.190473896004022
g2_height  =  2.3398363037555323e-06
#####################################################################################################################################################################



# # # Ei = 700meV
# # params['g1_center'].set(value = 330)
# # params['g1_sigma'].set(value = 5)
# # params['g2_center'].set(value = 385)
# # params['g2_sigma'].set(value = 5)
# # params['g2_amplitude'].set(value = .15, min = 0)