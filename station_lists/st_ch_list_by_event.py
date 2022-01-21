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


working_dir = '/Users/emmadevin/Work/USGS 2021/Data/gmprocess/qa_processing/data/ci38443095/raw'

    
data = []
file_list = glob.glob(working_dir + '/*.ms')

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
df.to_csv('/Users/emmadevin/Work/USGS 2021/Data/gmprocess/qa_processing/data/ci38443095/stns_ch.csv')