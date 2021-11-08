#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 20 11:35:28 2021

@author: emmadevin
"""
import obspy as op
from obspy import read
import os
import os.path as path
import glob
import pandas as pd
from random import randint

working_dir = '/Users/emmadevin/Work/USGS 2021/Data/Prelim_filtered/RC_beta'

eventlist = glob.glob(working_dir + '/*')


for event in eventlist:
    
    filelist = glob.glob(event +'/*.ms')
    
    random1 = randint(0, len(filelist))
    random2 = randint(0, len(filelist))
    random3 = randint(0, len(filelist))
    
    st1 = read(filelist[random1])
    st2 = read(filelist[random2])
    st3 = read(filelist[random3])
    
    
    st1.plot()
    st2.plot()
    st3.plot()
    
    print(event)
      
    n = input("Continue to next event [y/n]?")
        
    if n == 'y': continue
       