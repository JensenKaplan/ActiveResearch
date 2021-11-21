#!/usr/bin/env python
# coding: utf-8

# In[5]:


get_ipython().run_line_magic('reload_ext', 'autoreload')
import sys
sys.path.append('..')
import JensenTools as JT
from lmfit import Model
import PyCrystalField as cef
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn
seaborn.set()


# In[6]:


B40  =  -0.6568663783690575
B60  =  -0.02328250024945387
LS  =  100.00007580463522
B44  =  -3.1415463304732714
B64  =  0.504906552605772
B20  =  0.4858075931009187


# In[7]:


stev = {'B20': B20, 'B40': B40, 'B44': B44, 'B60' : B60 , 'B64' : B64}
Pr = cef.CFLevels.Bdict(Bdict = stev, ion = 'Ce3+')
Pr.diagonalize()
Pr.printEigenvectors()


# In[ ]:




