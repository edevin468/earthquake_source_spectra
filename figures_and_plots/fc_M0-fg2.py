#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  1 11:01:18 2021

@author: emmadevin
"""


import numpy as np
import glob
import os.path as path
import pandas as pd
import matplotlib.pyplot as plt


working_dir =  '/Users/emmadevin/Work/USGS 2021/Data/Prelim_filtered'
outfile_path = working_dir + '/stress_drops'
df = pd.read_csv(outfile_path + '/stress_drops_prelim_filtered.csv')

y = df['best-fit M0']
x = df['best-fit fc']

y_err = df['std M0']
x_err = df['std fc']

sigma1 = 1e6
sigma2 = 10e6
beta = 3500
fc = np.logspace(-1,1, 10)
Mo1 = sigma1*(0.42*beta/fc)**3
Mo2 = sigma2*(0.42*beta/fc)**3

fig = plt.figure(figsize = (4,4))
plt.style.use('classic')
fig.patch.set_facecolor('white')
plt.scatter(x, y, edgecolors = 'none', facecolors = 'springgreen', s = 40)
plt.errorbar(x, y,xerr = x_err,yerr=y_err, fmt='.', c='k')
plt.plot(fc, Mo1, c='k')
plt.plot(fc, Mo2, c='k')
plt.xlim(0.3,10)
plt.ylim(10**13,10**17)
plt.text(1, 4*10**16, r'$\sigma =  10 \mathrm{MPa}$')
plt.text(0.36, 4*10**15, r'$\sigma =  1\mathrm{MPa}$')
plt.xscale('log')
plt.yscale('log')
plt.title('(c)', loc = 'left')
plt.ylabel(r'best fit $M_0$ (dyn-cm)')
plt.xlabel(r'best fit $f_c$ (Hz)')