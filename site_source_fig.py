#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 11:14:52 2021

@author: emmadevin
"""
import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from obspy import read


fig = plt.figure(figsize=(8,6))
fig.patch.set_facecolor('white')
plt.style.use('classic')
working_dir = '/Users/emmadevin/Work/USGS 2021/Data/Prelim/Andrews_inversion/'

event = working_dir+'Events/38443183.out'
station = working_dir+'Stations/CIAVM.out'
trace = '/Users/emmadevin/Work/USGS 2021/Data/Prelim/RC_beta/38443183/AVM.CI.HHZ..2019.185.173349.38443183.ms'


source = np.genfromtxt(event, dtype = float, comments = '#', delimiter = None)
source = pd.DataFrame(source)
f_source = source[0]
a_source = source[1]

site = np.genfromtxt(station, dtype = float, comments = '#', delimiter = None)
site = pd.DataFrame(site)
f_site = site[0]
a_site = site[1]
 
st = read(trace)
tr = st[0] 
t = np.arange(0, tr.stats.npts / tr.stats.sampling_rate, tr.stats.delta) 
 


a1 = plt.subplot(211)
a1.plot(t,tr, lw=0.5, c='k')
# a1.set_yticks([])
a1.set_yticklabels([])
a1.set_xticks([0,100,200,300])
a1.set_title('(a) ', loc='left')
a1.set_xlabel('time (s)')
a1.text(305,1.1*10**6,'event ID: 38443183 \nnetwork: CI \nstation: AVM', bbox={'facecolor':'w'}, fontsize = 11)
a1.text(10,1.9*10**6,'06 JULY 2019 03:19:53',bbox={'facecolor':'w'}, fontsize = 11)

a2 = plt.subplot(212)
a2.plot(f_source,a_source, c = 'b',label='source spectra')
a2.plot(f_site,a_site, c = 'r', label='site spectra')
plt.legend(loc='lower left', ncol = 2, fontsize = 11)
a2.set_xscale('log')
a2.set_yscale('log')
a2.set_yticklabels([])
a2.set_title('(b)', loc = 'left')
a2.set_xlabel('frequency (Hz)')
a2.set_xlim(10**-2,60)
a2.set_ylim(10**-2,10**4)
a2.text(7.4,1*10**2,'event ID: 38443183 \nnetwork: CI \nstation: AVM', bbox={'facecolor':'w'}, fontsize = 11)
a2.text(0.012,1.4*10**3,'06 JULY 2019 03:19:53', bbox={'facecolor':'w'}, fontsize = 11)
plt.tight_layout()

plt.savefig('/Users/emmadevin/Work/USGS 2021/Figures/SCEC/example_spectra.pdf')

     