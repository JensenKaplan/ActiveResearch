from JensenTools import *

comp = 'Sr2PrO4'
saveDir = getSaveDir('m', comp = comp)

wyb = {'B20': -311.4421417091645, 'B40': -2664.8979343417955, 'B44': -1318.0240415089413, 'B60': 709.0624635353416, 'B64': -1209.5640959298416}


x = '/Users/jensenkaplan/Dropbox (GaTech)/Jensen/Sr2PrO4/Sr2PrO4_fitted_X.csv' # We need to re-open the file
m = '/Users/jensenkaplan/Dropbox (GaTech)/Jensen/Sr2PrO4/Sr2PrO4_fitted_M.csv'
dfX = pd.read_csv(x)
dfM = pd.read_csv(m)
TempX = dfX['TempX']
X = dfX['X']
Hm = dfM['Hm']
Mm = dfM['Mm']


Xfile = saveDir + 'Spectre/' + 'Sr2PrO4_XvsT_3T.dat'
headerList = ['Temperature', 'Xx', 'XyXz', 'Powder']
df = pd.read_csv(Xfile, sep = '\s+', names = headerList)
df = df.replace('D','E', regex=True)
TSpec = pd.to_numeric(df['Temperature'])
Xx = pd.to_numeric(df['Xx'])
XyXz = pd.to_numeric(df['XyXz'])
XPowder = pd.to_numeric(df['Powder'])

Mfile = saveDir + 'Spectre/' + 'Sr2PrO4_MvsH_50K.dat'
headerList = ['Field', 'Mx', 'My', 'Mz']
df = pd.read_csv(Mfile, sep = '\s+', names = headerList)
df = df.replace('D','E', regex=True)
HSpec = pd.to_numeric(df['Field'])
Mz = pd.to_numeric(df['Mz'])

TempM = 50
FieldX = 3.

plt.figure()
plt.plot(Hm,Mm, label = 'PCF Powder Average')
plt.plot(HSpec, Mz, label = 'SPECTRE Mz')
plt.xlabel('Field (oe)')
plt.ylabel('M (emu mol^-1)')
plt.legend()
plt.title('{} Magnetization at {} K '.format(comp,TempM))


plt.figure()
plt.plot(TempX,X*TempX, label = 'PCF Powder Average')
plt.plot(TSpec,XPowder*TSpec, label = 'SPECTRE Powder Average')
plt.xlabel('Temperature (K)')
plt.ylabel('X*T (emu T^-1 mol^-1 K)')
plt.legend()
plt.title('{} X(T)*T Applied {} Tesla '.format(comp,FieldX))

plt.show()