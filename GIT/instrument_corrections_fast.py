#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 28 14:15:31 2021

Applies instrument correction to all of the preliminary dataset

@author: emmadevin
"""

import obspy as op
from obspy import read
import os
import os.path as path
import glob
import pandas as pd
import warnings
import time
from progress_bar import progress_bar

warnings.filterwarnings("ignore", category=UserWarning)

# working directory
working_dir = '/Users/emmadevin/Work/USGS 2021/Data/Prelim_filtered'

# event directories and outpath
event_dirs = glob.glob(working_dir + '/RC_beta/*')
outpath = working_dir + '/corrected/'

# station ids
stn_ids = pd.read_csv( '/Users/emmadevin/Work/USGS 2021/Data/Prelim_filtered/Station_info/stations.csv')

# create list of event directory names
events = []
for i in range(len(event_dirs)):
    events.append(path.basename(event_dirs[i]))
    
# create out directories for corrected events   
# for i in range(len(events)):
#     if not path.exists(outpath + '/' + events[i]):
#         os.makedirs(outpath + '/'  + events[i])


        

# loop through *.ms files and apply instrument corrections and save the resulting files   
uncorrected = []   

for event in events: 
    
    print('\nEvent: '+event)
   
    # event directory
    event_dir = working_dir + '/RC_beta/' + event
    
    # create list of all files in event directory
    file_list = glob.glob(event_dir + '/*.ms')
    
    # read in *.csv that says which response file to use
    r = pd.read_csv(event_dir + '/station_inv.csv')
    
    # load inventory files
    inv1 = op.read_inventory(event_dir + '/' + event + '.scedc.xml')
    inv2 = op.read_inventory(event_dir + '/' + event + '.ncedc.xml')
    inv3 = op.read_inventory(event_dir + '/' + event + '.iris.xml')

    #initilize progress bar
    l = len(file_list)     
    bar = progress_bar(l) 
    i = 0   
    
    for file in file_list: 
        bar.get_progress(i)
        i+=1
        
        
        t1 = time.time()
        
        # read in file and determine station 
        st = read(file)
        filename = path.basename(file)
        ntwk = filename.split('.')[1]
        stn = filename.split('.')[0]
        stn_id = ntwk + '|' + stn
        
        # determine which *.xml file to use, if none exists, move on and add to list of uncorrected files
        row = r.loc[r['station_id']==stn_id]
        row = row.reset_index()
        
        if row.loc[0,'type']=='none':
            uncorrected.append(filename)
        elif row.loc[0,'type']=='scedc':
            inv = inv1
        elif row.loc[0,'type']=='ncedc':
            inv = inv2
        elif row.loc[0,'type']=='iris':
            inv = inv3
        
            
        t2 = time.time() 
        
        # print('Time to determine which inv file to use: ', (t2-t1))
        
         
        t1 = time.time()
        
        # extract trace and make copy of it
        tr = st[0] 
        tr_corrected=tr.copy()
        
        # determine sampling rate, define nyquist frequency, and define filter band
        sample_rate= tr.stats.sampling_rate
        nyquist = sample_rate*0.5
        pre_filt = [0, 0.001, 35, nyquist]
        
        # detrend
        tr_corrected.detrend('linear')
        
        t2 = time.time() 
        # print('Time extract trace and detrend: ', (t2-t1))
        
        
        t1 = time.time()
        
        # remove instrument response, if any exception occurs, move on and add to list of uncorrected
        try: 
            tr_corrected.remove_response(inventory=inv,output='VEL',pre_filt=pre_filt,water_level=60,taper=False,plot=False)
        except Exception:
            uncorrected.append(filename)
            continue
        
        t2 = time.time()
        
        # print('Time remove instrument response: ', (t2-t1))
        
        # write miniseed files  
        # filename = filename.replace('ms', 'mseed')
        # tr_corrected.write(outpath + event + '/' + filename)
        
    print('\n---------------------------------------------------------------------------------------')