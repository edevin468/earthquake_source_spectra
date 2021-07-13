#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 10:31:36 2021

determine if 1/r is appropriate geometrical spreading correction

@author: emmadevin
"""

import glob
import os.path as path
import numpy as np
import obspy
from obspy import read    
import dread
import time
import random
import pandas as pd
import matplotlib.pyplot as plt



#read in the cut and corrected spectra = records
#records are in m/s
#should be binned frequencies and amplitudes
#make list of records and the corresponding events and stations

# working directory and outfil path for inversion
working_dir = '/Users/emmadevin/Work/USGS 2021/Data/Prelim'
outfile_path = working_dir + '/Andrews_inversion'

# df with station locations
stations = pd.read_csv(working_dir + '/station_locs.csv')

# list of event folders
ev = glob.glob(working_dir + '/record_spectra/*')

# list of files
   
    
# get lists of station ids and station locations
stn_list = (stations['network']+stations['station']).tolist()
stn_lat = stations['latitude']
stn_lon = stations['longitude']


df = pd.DataFrame()

events = []

t1 = time.time()

for event in ev:
    
    distance = []
    amp1 = []
    amp2 = []
    amp3 = []
    amp4 = []
    amp5 = []
        
    evt = event.split('/')[-1]
    print(evt)
    events.append(evt)
    base = path.basename(evt)
    print(working_dir + '/' + evt + '/*')
    
    filelist = glob.glob(working_dir + '/record_spectra/' + evt + '/*')
    
    for file in filelist: 
        
        base = file.split('/')[-1]
        eventid = base.split('_')[-1]
        eventid = eventid.split('.')[0]
        ntwk = base.split('_')[0]
        stn = base.split('_')[1]
        stn_id = ntwk + stn
        
        # print some updates to track progress
        print('Event: ', eventid)
        print('Station: ', stn)
        
        # read phase file fo get event locations
        phase_file = working_dir + '/RC_phase_beta/' + eventid + '.phase'
        phase = pd.read_csv(phase_file, sep = '\s+', index_col=0, nrows = 0).columns.tolist()
        
        # assign event coordinates and depth
        evlat = float(phase[4])
        evlon = float(phase[5])
        evdepth = float(phase[6])
        
        # assign station corrdinates and depth
        stlat = stn_lat[stn_list.index(stn_id)]
        stlon = stn_lon[stn_list.index(stn_id)]
        stdepth = 0
        
        #find distance between event and station
        dist =  dread.compute_rrup(evlon, evlat, evdepth, stlon, stlat, stdepth) #in km
        
        #read in spectra file
        data = np.genfromtxt(file, dtype = float, comments = '#', delimiter = None, usecols = (0,1,2))  #only read in first two cols
       
        amp1.append(data[0,1])
        amp2.append(data[27,1])
        amp3.append(data[57,1])
        amp4.append(data[69,1])
        amp5.append(data[74,1])
        
        distance.append(dist)
    
        
    df[str(evt)] = [distance, amp1, amp2, amp3, amp4, amp5]
    
    
x = np.logspace(1,3,100)
y = 1/x

    
plt.plot(x,y, c = 'grey', label = '~1/r')   
plt.plot(x,y/10**1, c = 'grey') 
plt.plot(x,y/10**2, c = 'grey') 
plt.plot(x,y/10**3, c = 'grey') 
plt.plot(x,y/10**4, c = 'grey') 
plt.plot(x,y/10**5, c = 'grey') 
plt.plot(x,y/10**6, c = 'grey') 
plt.plot(x,y/10**7, c = 'grey') 

f = 5

plt.scatter(df.loc[0,events[0]],df.loc[f,events[0]], c = 'b', label = events[0])
plt.scatter(df.loc[0,events[1]],df.loc[f,events[1]], c = 'r', label = events[1])
plt.scatter(df.loc[0,events[2]],df.loc[f,events[2]], c = 'm', label = events[2])
plt.scatter(df.loc[0,events[3]],df.loc[f,events[3]], c = 'c', label = events[3])
plt.scatter(df.loc[0,events[4]],df.loc[f,events[4]], c = 'k', label = events[4])
plt.scatter(df.loc[0,events[5]],df.loc[f,events[5]], c = 'dimgrey', label = events[5])
plt.scatter(df.loc[0,events[6]],df.loc[f,events[6]], c = 'purple', label = events[6])
plt.scatter(df.loc[0,events[7]],df.loc[f,events[7]], c = 'orange', label = events[7])
# plt.scatter(distance, amp2, c = 'r', label = 'f = 0.057 Hz')
# plt.scatter(distance, amp3, c = 'm', label = 'f = 1.55 Hz')
# plt.scatter(distance, amp4, c = 'k', label = 'f = 5.78 Hz')
# plt.scatter(distance, amp5, c = 'c', label = 'f = 10 Hz')
plt.ylim(10**-10,10**-2)
plt.xlim(10**1,20**2)
plt.xlabel('distance (km)')
plt.ylabel('amplitude')
plt.title('f = 10Hz', c = 'dimgrey')
plt.xscale('log')
plt.yscale('log')
leg = plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), facecolor = 'w', ncol = 2, fontsize = 10)

for text in leg.get_texts():
    plt.setp(text, color = 'dimgrey')

    
