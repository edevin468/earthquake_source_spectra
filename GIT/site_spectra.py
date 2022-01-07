#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 17 11:16:40 2021

compare our site spectra with linear site reponse model from Bayless and Abrahamson (2019) 

@author: emmadevin
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob
import os
import os.path as path

working_dir = '/Users/emmadevin/Work/USGS 2021/Data/Prelim_filtered/'

station_vs30 = pd.read_csv(working_dir + 'Station_info/station_locs.csv')
station_ids = station_vs30['network']+station_vs30['station']
vs30s = station_vs30['VS30']
vs30_dict = {station_ids[i]: vs30s[i] for i in range(len(station_ids))}

coeffs =pd.read_csv('/Users/emmadevin/Work/USGS 2021/Site Spectra/Bayless_ModelCoefs.csv')

site_spectra = glob.glob(working_dir + 'Andrews_inversion/Stations/*.out')

for i in range(len(site_spectra)):
    site = path.basename(site_spectra[i]).split('.')[0]
    vs30 =  vs30_dict[site]

    
    data = np.genfromtxt(site_spectra[i], dtype = float, comments = '#', delimiter = None, usecols = (0,1,2)) #only read in first two cols

    freq = data.T[0]
    spectra = data.T[1]
    
    c8 = coeffs['c8']
    f = coeffs['f (Hz)']
    fsl = np.exp(c8 * np.log(min(vs30, 1000)/1000))
    
    fig = plt.figure(figsize = (5,5))
    plt.style.use('classic')
    fig.patch.set_facecolor('white')
    
    plt.plot(freq, spectra, c= 'grey')
    plt.plot(f, fsl, c= 'r')
    # plt.xlim(0.04, 50)
    # plt.ylim(10**-4,10**2)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('frequency (Hz)')
    plt.ylabel('velocity amplitude (m)')
    
    
