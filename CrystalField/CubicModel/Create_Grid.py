import sys
sys.path.append('../../')
from JensenTools import *


comp = 'Sr2PrO4'
LS_on = False
numlevels = 1
saveDir = getSaveDir('m',dataType = 'grid', comp = comp)
ion = 'Ce3+'
Kmeans = True

LSDir = '/cubic_matrix_LS_test/'
JDir = '/cubic_matrix_J/'
if LS_on:
	saveDir = saveDir + LSDir
else:
	saveDir = saveDir + JDir

grid = 'J_Grid_2_levels.mat'

# grid = 'LS_100.mat'
LS = [100,110]


xmin, xmax = -1,1
bpfmin, bpfmax =  -1,1
numx, numbpf = 200, 200


if __name__ == '__main__':
		saveMatrixPar(xmin,xmax,numx,bpfmin,bpfmax,numbpf,saveDir,LS_on = LS_on, LSList = LS, numlvls = numlevels,  L = 3, S = .5, ion = ion, grid = grid, Kmeans = Kmeans)

