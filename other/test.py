#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 11 09:36:08 2021

@author: emmadevin
"""

import glob
import os.path as path
import numpy as np
import obspy
from obspy import read    
import time
import random
import pandas as pd
import dread

#read in the cut and corrected spectra = records
#records are in m/s
#should be binned frequencies and amplitudes
#make list of records and the corresponding events and stations
working_dir = '/Users/emmadevin/Work/USGS 2021/Data/Prelim'
outfile_path = working_dir + '/Andrews_inversion'

# df with station locations
stations = pd.read_csv(working_dir + '/station_locs.csv')

#list of record files
ev = working_dir + '/record_spectra/38443183/CI_AVM_HHNE__38443183.out'


   

# get lists of station ids and station locations
stn_list = (stations['network'] + stations['station']).tolist()
# stn_list = stn_list.tolist()
stn_lat = stations['latitude']
stn_lon = stations['longitude']

    
    
eventidlist = []
event_lat = []
event_lon = []
event_depth = []
record_freq = []
record_spec = []
record_std = []

record = (ev.split('/')[-1])
base = path.basename(record)
event = base.split('_')[-1]
event = event.split('.')[0]
ntwk = base.split('_')[0]
stn = base.split('_')[1]
stn_id = ntwk + stn
print(event)

phase_file = working_dir + '/RC_phase_beta/' + event + '.phase'
phase = pd.read_csv(phase_file, sep = '\s+', index_col=0, nrows = 0).columns.tolist()

evlat = float(phase[4])
evlon = float(phase[5])
evdepth = float(phase[6])

index = stn_list.index(stn_id)

stlat = stn_lat[stn_list.index(stn_id)]


# evlat = phase_df[]

        
