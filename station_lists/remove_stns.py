#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 19 12:13:57 2022

@author: emmadevin
"""
import pandas as pd
import os 
import os.path as path
import glob

ch = 'HN'

working_dir = '/Users/emmadevin/Work/USGS 2021/Data/gmprocess/station_downloads/'+ch+'/data/ci38548295/raw/'
stn_files = glob.glob(working_dir + '*.xml')

stns = []
for stn in stn_files: 
    name = path.basename(stn).split('/')[-1]
    name = name.split('.')[0] + '.' + name.split('.')[1]
    stns.append(name)

df = pd.read_csv('/Users/emmadevin/Work/USGS 2021/Data/stns_channels/stns_'+ch+'.csv')
ids = list(df['network']+'.'+df['station'])

for stn in stns:
    if stn not in ids: 
        os.remove(working_dir + stn + '.xml')
        


