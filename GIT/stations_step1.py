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

# working directory
working_dir = '/Users/emmadevin/Work/USGS 2021/Data/Prelim+'

# event directories and outpath
event_dirs = glob.glob(working_dir + '/RC_beta/*')


# create list of event directory names
events = []
for i in range(len(event_dirs)):
    events.append(path.basename(event_dirs[i]))
    

# loop through *.ms files and apply instrument corrections and save the resulting files         
for event in events: 
    
    # event directory
    event_dir = working_dir + '/RC_beta/' + event
    
    # create list of all files in event directory
    file_list = glob.glob(event_dir + '/*.ms')
    
    # create lists for stations, networks, and files we cannot do instrument corrections for
    stns = []
    ntwks = []
    stn_ids = []
    uncorrected = []
    for file in file_list: 
        # read in file and determine station 
        st = read(file)
        filename = path.basename(file)
        ntwk = filename.split('.')[1]
        stn = filename.split('.')[0]
        stn_id = ntwk + '|' + stn
        
        if stn_id not in stn_ids:
            stn_ids.append(stn_id)
    

    df = pd.DataFrame()
    df['station_id'] = stn_ids
                
    df.to_csv(working_dir + '/RC_beta/'+event+'/stations.csv')