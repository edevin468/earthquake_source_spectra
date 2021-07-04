#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  2 13:02:10 2021

@author: emmadevin
"""

import obspy as op
from obspy import read, Stream
import os
import os.path as path
import glob
import pandas as pd


# working directory
working_dir = '/Users/emmadevin/Work/USGS 2021/Data/Prelim'

# event directories and outpath
event_dirs = glob.glob(working_dir + '/RC_beta/*')

# create list of event directory names
events = []
for i in range(len(event_dirs)):
    events.append(path.basename(event_dirs[i]))
    
# loop through events and check which station appears in which response file      
for event in events: 
    
    # event directory
    event_dir = working_dir + '/RC_beta/' + event
    
    # create list of all files in event directory
    file_list = glob.glob(event_dir + '/*.txt')
    
    for file in file_list:
         