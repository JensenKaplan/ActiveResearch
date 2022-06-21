import sys
sys.path.append('..')
from JensenTools import *
from matplotlib import rcParams
from matplotlib import patches
import pandas as pd

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

def gaussian(mu,sigma,amp,x):
    g = []
    for i in x:
        g.append((amp/sigma/np.sqrt(2*np.pi))*np.exp(-(i-mu)**2/2/sigma**2))
    return np.array(g)




# #Create the 3 gaussians with parameters determine by LMFIT
# g1 = gaussian(mu = 168.22215251994697,sigma = 4.095279662587834,amp = .0002701750011412199, x = x) #Ei = 300meV, T = 4.88K
# g2 = gaussian(mu = 336.00897026613654,sigma = 15.304438676288882, amp = 0.00032946653900736716, x = x) #Ei = 700meV, T=4.76K
# g3 = gaussian(mu = 384.9032448977165,sigma = 10.061567665493945, amp = 7.359747321640221e-05, x = x) #Ei = 700meV, T=4.76K
# #Forcing widths to be the same
# g1 = gaussian(mu = 168.22215251994697,sigma = 3,amp = .0002701750011412199, x = x) #Ei = 300meV, T = 4.88K
# g2 = gaussian(mu = 336.00897026613654,sigma = 3, amp = 0.00032946653900736716, x = x) #Ei = 700meV, T=4.76K
# g3 = gaussian(mu = 384.9032448977165,sigma = 3, amp = 7.359747321640221e-05, x = x) #Ei = 700meV, T=4.76K


# g1 = gaussian(mu = 168.1, sigma = .2276, amp = 2.444e-7, x = x) #Ei = 500meV, T = 4.88K
# g2 = gaussian(mu = 334.4, sigma =  5.779, amp = .0001288, x = x) #Ei = 500meV, T = 4.88K
# g3 = gaussian(mu = 332.8, sigma =  23.93, amp = .0005951, x = x) #Ei = 700meV, T = 4.76K
# g4 = gaussian(mu = 388, sigma =  24.48, amp = .000292, x = x) #Ei = 700meV, T = 4.76K

# SrPr
# #####################################################################################################################################################################
# g1amp = .000324
# g2amp = 8.938e-5
# g22amp = .0001005 
# g3amp = 7.9943e-6
# g1g2 = 0.7057918532469233
# g3g2 = .25

# print('g1amp/g2amp = ', g1g2 )
# print('g3amp/g2amp = ', g3g2 )

# sigma = 3
# width = sigma/2.35
# x = np.linspace(30,500,1000) # Generate energy X-axis. 0-450 should do since highest peak is ~380


# g1 = gaussian(mu = 166.2, sigma = width, amp = g1g2, x = x) #Ei = 500meV, T = 4.88K
# g2 = gaussian(mu = 334.5, sigma =  width, amp = 1, x = x) #Ei = 500meV, T = 4.88K
# # g22 = gaussian(mu = 332.8, sigma =  3, amp = .0005951, x = x) #Ei = 700meV, T = 4.76K
# g3 = gaussian(mu = 388.8, sigma =  width, amp = g3g2, x = x) #Ei = 700meV, T = 4.76K



# curve = np.array(g1) + np.array(g2)  + np.array(g3)
#####################################################################################################################################################################

#LiPr 274 mode
#####################################################################################################################################################################
# g1_amplitude = 0.0003207 
# g1_center = 274.3 
# g1_fwhm  = 23.41 
# g1_height = 1.287e-05 
# g1_sigma  = 9.941    
#####################################################################################################################################################################

#SRPR FOR PAPER
#####################################################################################################################################################################
# 500 meV
g1_amplitude  =  0.00010025238254331347
g1_center  =  167.86803392064175
g1_sigma  =  8.28491499033998
g1_fwhm  =  19.509483517552393
g1_height  =  4.827438316378923e-06
g2_amplitude  =  0.00011474894377347192
g2_center  =  334.8342674591908
g2_sigma  =  6.138000561491445
g2_fwhm  =  14.453886482211287
g2_height  =  7.4581628158789436e-06


#700 meV
g2_amplitude  =  0.00022882035627084753
g2_center  =  335.3525716811439
g2_sigma  =  13.90613734716213
g2_fwhm  =  32.74645034784433
g2_height  =  6.564448267594624e-06
g3_amplitude  =  5.277865776909252e-05
g3_center  =  382.2838548009725
g3_sigma  =  8.998765891237555
g3_fwhm  =  21.190473896004022
g3_height  =  2.3398363037555323e-06

sigma = 6
width = sigma/2.35
x = np.linspace(30,700,1000) # Generate energy X-axis. 0-450 should do since highest peak is ~380

# g1_cen = 166.2
# g2_cen = 334.5
# g3_cen = 388.8

g1 = gaussian(mu = g1_center, sigma = width, amp = g1_amplitude/g2_amplitude, x = x) # Ei = 500meV, T = 4.88K
g2 = gaussian(mu = g2_center, sigma = width, amp = 1, x = x) # From Optics
g3 = gaussian(mu = g3_center, sigma = width, amp = g3_amplitude/g2_amplitude, x = x) # From Optics

curve = np.array(g1) + np.array(g2) + np.array(g3)
#####################################################################################################################################################################




lineErr = .00005*np.ones(len(curve))
# lineErr = 1*(curve)


# plt.figure()
plt.errorbar(x,curve,yerr = lineErr)
plt.xlabel('Energy (meV)')
plt.ylabel('Intensity (a.u.)')
plt.savefig('/Users/jensenkaplan/Dropbox (GaTech)/Jensen/Sr2PrO4/Sr2PrO4_Produced_Gaussians_FINALPAPER.png')
plt.show()

# np.savetxt('/Users/jensenkaplan/Dropbox (GaTech)/Jensen/Sr2PrO4/Sr2PrO4_scaled.ascii',np.transpose([x,curve]), header = 'Energy (meV), Intensity (arb. units)')
df = pd.DataFrame()
df['Energy'] = x
df['Intensity'] = curve
df['Error'] = lineErr

df.to_csv('/Users/jensenkaplan/Dropbox (GaTech)/Jensen/Sr2PrO4/Sr2PrO4_Produced_Gaussians_FINALPAPER.csv')


# f ='/Users/jensenkaplan/Dropbox (GaTech)/Jensen/Sr2PrO4/Sr2PrO4_scaled.csv' # We need to re-open the file
# df = pd.read_csv(f)

# energy = df['Energy']
# intensity = df['Intensity']


# plt.figure()
# plt.plot(x,curve)
# plt.show()