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

x = df['Brune fc']
y = df['best-fit fc']

x1 = [0.03, 10]
y1 = x1

fig = plt.figure(figsize = (4,4))
plt.style.use('classic')
fig.patch.set_facecolor('white')
plt.scatter(x, y, edgecolors = 'none', facecolors = 'mediumblue', s = 40)
plt.plot(x1, y1, c='k')
plt.xlim(0.03,10)
plt.ylim(0.03,10)
plt.xscale('log')
plt.yscale('log')
plt.title('(e)', loc = 'left')
plt.ylabel(r'best fit $f_c$')
plt.xlabel(r'theoretical $f_c$')