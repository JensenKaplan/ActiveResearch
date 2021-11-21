#!/usr/bin/env python
# coding: utf-8

# # The Creation of the Cubic System Grid Search
# You can specify the bounds of the coefficient ratio (x) and the pre-factor (bpf) as well as the grid density. The code has been parallelized for optimal results however this is still an intensive calculation. In the example below I only allow a 20x20 grid for 3 different LS values. Even this takes a few minutes to compute and create the .mat files. In this directory you'll find 800x800 grids that I created. These took over a day to compute. The idea was to create super-fine mesh grids once and be able to search through those .mat files as opposed to grid creation each time.

# In[1]:


# get_ipython().run_line_magic('reload_ext', 'autoreload')
import sys
sys.path.append('../../')
from JensenTools import *

saveDir = getSaveDir('m')


# JDir = 'cubic_matrix_J/'
# saveDir = saveDir + JDir
# grid = 'J_Grid_test_6lvls.mat'

LSDir = 'cubic_matrix_LS_test/'
saveDir = saveDir + LSDir
# grid = 'LS_100.mat'
LS = [100,110]

# In[2]:


numlevels = 4 #PCF produces a number of eigenvalues by diagonalizing the CF Hamiltonian. I specify the number of unique levels to cluster.
xmin, xmax = -1,1
bpfmin, bpfmax =  -1,1
numx, numbpf = 20, 20


# In[3]:


if __name__ == '__main__':
		saveMatrixPar(xmin,xmax,numx,bpfmin,bpfmax,numbpf,saveDir,LS_on = True, LSList = LS, numlvls = numlevels,  L = 3, S = .5)

