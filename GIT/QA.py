#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 09:50:10 2021

@author: emmadevin
"""

from progress_bar import progress_bar
import obspy as op
from obspy import read
import os
import os.path as path
import glob
import pandas as pd
import warnings
import time
import QA_functions as qa
from obspy.signal.trigger import plot_trigger


# working directory
working_dir = '/Users/emmadevin/Work/USGS 2021/Data/QA_testing'

# event directories and outpath
event_dirs = glob.glob(working_dir + '/All/*')

outpath = working_dir + '/QA_auto'

# create list of event directory names
events = []
for i in range(len(event_dirs)):
    events.append(path.basename(event_dirs[i]))
    

# create out directories for QA filtered events   
for i in range(len(events)):
    if not path.exists(outpath + '/' + events[i]):
        os.makedirs(outpath + '/'  + events[i])
        

# examine files from each event and determine whether or not they pass QA criteria.
uncorrected = []      
for event in events: 
    
    print('\nEvent: '+event)
    # event directory
    event_dir = working_dir + '/All/' + event
    
    # create list of all files in event directory
    file_list = glob.glob(event_dir + '/*.ms')
    
    # loop through all files for event
    for file in file_list: 
        filename = path.basename(file)
        
        st = read(file)
        
        # perform QA checks on file
        checklist = []
        sta_lta = qa.check_sta_lta(st, sta_length=1.0, lta_length=20.0, threshold=3.0) 
        
        # detrend (linear), demean, check 0 crossings, sta/lta
        
        tr = st[0]
        if sta_lta == False:
            color = 'red'
        else: color = 'blue'
        
        tr.plot(color = color)
        
        checklist.append(sta_lta)
        
        start = qa.signal_split(st, origin, model=None,picker_config=None,config=None)
        # checklist.append(qa.check_zero_crossings(st, min_crossings=0.1))
        # cft, tr = qa.check_zero_crossings(st, min_crossings=0.1)
        # plot_trigger(tr, cft, 1.5, 0.5)
        
        # # if file passes QA check save to  
        # if False not in checklist:
        #     filename = filename.replace('ms', 'mseed')
        #     print(filename)
        #     # st.write(outpath +'/'+ event + '/' + filename)
        
        
        
        
        
        
    
    
    
    
    
    
    
