#!/usr/bin/env python
# coding: utf-8

# # Here is a real use of the cubic grid search. I use the grid to find starting parameters for the compound Sr2PrO4 which has tetragonal symmetry. The starting parameters will then be used in the file "Crystal Field Analysis".

# In[1]:


# get_ipython().run_line_magic('reload_ext', 'autoreload')
import sys
sys.path.append('../../')
from JensenTools import *
import numpy as np
import matplotlib.pyplot as plt
import PyCrystalField as cef
import os
import scipy.io as sio
from functools import reduce
# import time
saveDir = getSaveDir('w')
saveDir = saveDir + 'cubic_matrix_J/'
grid = 'J_Grid_test_6lvls.mat'
saveDir = saveDir + grid


# ### Define the measured energy levels (from INS data) and define an allowable tolerance between calculated and measured energy.

# In[2]:


tol = .1 #tolerance allowed between measured and calculated energy.
Emeas = [168, 335] # Our measured levels of Sr2PrO4
comp = 'Sr2PrO4' #Compound name

# ### In the following section we scan through all LS grids and find the (x,bpf) points that create matching energy levels.

# In[3]:


print('Energies as measured by paper (meV):  ', Emeas)

EList, data = loadMatrixJ(saveDir) #Load in all created 800x800 grids
 


#Loading the x,bpf, and LS of each file.
x = data['X'][0]
bpf = data['B'][0]
print('Size is: {} x {}'.format(len(x),len(bpf)))
print('Print the bands {}'.format(EList))
plotContoursJ(data,EList)

# plotContoursJ(data,EList) #Contour plotting for 4 E levels

plt.show()
