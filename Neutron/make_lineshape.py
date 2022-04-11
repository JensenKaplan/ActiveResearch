import sys
sys.path.append('..')
from JensenTools import *

def gaussian(mu,sigma,amp,x):
    g = []
    for i in x:
        g.append((amp/sigma/np.sqrt(2*np.pi))*np.exp(-(i-mu)**2/2/sigma**2))
    return np.array(g)


x = np.linspace(30,500,1000) # Generate energy X-axis. 0-450 should do since highest peak is ~380

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



g1amp = .000324
g2amp = 8.938e-5
g22amp = .0001005
g3amp = 7.9943e-6
g1g2 = 0.7057918532469233
g3g2 = .25
print('g1amp/g2amp = ', g1g2 )
print('g3amp/g2amp = ', g3g2 )

sigma = 6
width = sigma/2.35

g1 = gaussian(mu = 166.2, sigma = width, amp = g1g2, x = x) #Ei = 500meV, T = 4.88K
g2 = gaussian(mu = 334.5, sigma =  width, amp = 1, x = x) #Ei = 500meV, T = 4.88K
# g22 = gaussian(mu = 332.8, sigma =  3, amp = .0005951, x = x) #Ei = 700meV, T = 4.76K
g3 = gaussian(mu = 388.8, sigma =  width, amp = g3g2, x = x) #Ei = 700meV, T = 4.76K


curve = np.array(g1) + np.array(g2)  + np.array(g3)


lineErr = .003*np.ones(len(curve))
# lineErr = .01*curve


plt.figure()
plt.errorbar(x,curve,yerr = lineErr)
plt.show()

# np.savetxt('/Users/jensenkaplan/Dropbox (GaTech)/Jensen/Sr2PrO4/Sr2PrO4_scaled.ascii',np.transpose([x,curve]), header = 'Energy (meV), Intensity (arb. units)')
df = pd.DataFrame()
df['Energy'] = x
df['Intensity'] = curve
df['Error'] = lineErr


df.to_csv('/Users/jensenkaplan/Dropbox (GaTech)/Jensen/Sr2PrO4/Sr2PrO4_gaussians.csv')


# f ='/Users/jensenkaplan/Dropbox (GaTech)/Jensen/Sr2PrO4/Sr2PrO4_scaled.csv' # We need to re-open the file
# df = pd.read_csv(f)

# energy = df['Energy']
# intensity = df['Intensity']


# plt.figure()
# plt.plot(x,curve)
# plt.show()