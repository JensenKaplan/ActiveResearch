#!/usr/bin/env python
# coding: utf-8

# # The Creation of the Cubic System Grid Search
# You can specify the bounds of the coefficient ratio (x) and the pre-factor (bpf) as well as the grid density. The code has been parallelized for optimal results however this is still an intensive calculation. In the example below I only allow a 20x20 grid for 3 different LS values. Even this takes a few minutes to compute and create the .mat files. In this directory you'll find 800x800 grids that I created. These took over a day to compute. The idea was to create super-fine mesh grids once and be able to search through those .mat files as opposed to grid creation each time.

# In[1]:


# get_ipython().run_line_magic('reload_ext', 'autoreload')
import sys
sys.path.append('../../')
from JensenTools import *
saveDir = getSaveDir()
dataMat ='cubic_matrix_J/'
grid = 'J_Grid.mat'
saveDir = saveDir + dataMat + grid

# In[2]:


xmin, xmax = -1,1
bpfmin, bpfmax =  -1,1
numx, numbpf = 200,200


# In[3]:



saveMatrixParJ(xmin,xmax,numx,bpfmin,bpfmax,numbpf,saveDir,numlevels, LS = False)

