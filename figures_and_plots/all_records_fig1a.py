#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  2 11:11:43 2021

@author: emmadevin
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob

fig = plt.figure(figsize = (3,5))
plt.style.use('classic')
fig.patch.set_facecolor('white')



 
    
path = glob.glob('/Users/emmadevin/Work/USGS 2021/Data/Prelim_filtered/record_spectra/*/*.out')

for i in range(len(path)):
    data = np.genfromtxt(path[i], dtype = float, comments = '#', delimiter = None, usecols = (0,1,2)) #only read in first two cols

    freq = data.T[0]
    spectra = data.T[1]
    
    
    plt.plot(freq, spectra, c= 'k')
    plt.title('(a)', loc='left')
    plt.xlim(0.04, 50)
    # plt.ylim(10**-4,10**2)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('frequency (Hz)')
    plt.ylabel('velocity amplitude (m)')
    
    
# path = glob.glob('/Users/emmadevin/Work/USGS 2021/Data/Prelim_filtered/Andrews_inversion/Events/*.out')

# for i in range(len(path)):
#     data = np.genfromtxt(path[i], dtype = float, comments = '#', delimiter = None, usecols = (0,1,2)) #only read in first two cols

#     freq = data.T[0]
#     spectra = data.T[1]
    
    
#     plt.plot(freq, spectra, c= 'blue')
#     plt.title('(b)', loc='left')





# x = np.linspace(100,120,10)
# y = x

# plt.plot(x,y, c='grey', label='station spectra')
# plt.plot(x,y, c='blue', label='event spectra')
# plt.legend(loc = 'lower center', fontsize = 11)