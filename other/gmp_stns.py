#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  3 09:57:16 2022

@author: emmadevin
"""

import obspy as op
from obspy import read
import os
import os.path as path
import glob
import pandas as pd

# station ids from SDVS dataset
stn_ids = pd.read_csv( '/Users/emmadevin/Work/USGS 2021/Data/Prelim_filtered/Station_info/stations.csv')
og_stns = stn_ids['network'] + '.' + stn_ids['station']

# stations downloaded with gmprocess
working_dir = '/Users/emmadevin/Work/USGS 2021/Data/gmprocess_test2/data/ci38443095/raw'

files = glob.glob(working_dir + '/*.xml')

gm_stns = []
for file in files:
    filename = path.basename(file).split('/')[-1]
    net = filename.split('.')[0]
    stat = filename.split('.')[1]
    gm_stn = net +'.'+stat
    gm_stns.append(gm_stn)
    
both = set(og_stns).intersection(gm_stns)

print('orignal dataset: ' + str(len(og_stns)) + ' stations')
print('gm dataset: ' + str(len(gm_stns)) + ' stations')
print('interesection: ' + str(len(both)) + ' stations')
    