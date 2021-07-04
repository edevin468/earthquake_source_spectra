#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  2 11:11:43 2021

@author: emmadevin
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

path = '/Users/emmadevin/Work/USGS 2021/Data/Prelim/record_spectra/38443183/CI_AVM_HHNE__38443183.out'

data = np.genfromtxt(path, dtype = float, comments = '#', delimiter = None, usecols = (0,1,2)) #only read in first two cols
df = pd.DataFrame(data)

freq = df[0]
spectra = df[1]

plt.plot(freq, spectra)
plt.xscale('log')
plt.yscale('log')
