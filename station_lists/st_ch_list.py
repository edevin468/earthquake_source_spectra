#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 09:19:54 2021

Create *.csv file in each event folder with list of stations 

@author: emmadevin
"""

import obspy as op
from obspy import read, Stream
import os
import os.path as path
import glob
import pandas as pd


working_dir = '/Users/emmadevin/Work/USGS 2021/Data/All_M3+'
event_dirs = glob.glob(working_dir + '/RC_M3above/*')

    
data = []
i = 0    
for event_dir in event_dirs: 
    i +=1 
    event = path.basename(event_dir).split('/')[-1]
    print(event + '  (' + str(i) + ' of ' + str(len(event_dirs)) + ')') 
    
    file_list = glob.glob(event_dir + '/*.ms')
    
    for file in file_list: 
       
        st = read(file)
        filename = path.basename(file)
       
        ntwk = filename.split('.')[1]
        
        stn = filename.split('.')[0]
        
        ch = filename.split('.')[2]
        ch = ch[:2]
        
        t = (ntwk, stn, ch)
        
        if t not in data: 
            data.append(t)
        
df = pd.DataFrame(data, columns = ['network', 'station', 'channel'])  
df.to_csv('/Users/emmadevin/Work/USGS 2021/Data/stns_ch.csv')