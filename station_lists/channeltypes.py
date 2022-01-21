#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 17 09:45:56 2022

@author: emmadevin
"""

import glob
import pandas as pd
import os
import os.path as path

df1 = pd.DataFrame()
df2 = pd.DataFrame()

working_dir = '/Users/emmadevin/Work/USGS 2021/Data/gmprocess/'
data1 = 'gmprocess_us/data/ci38443095/raw/*.ms'
data2 = 'download3/data/ci38443095/raw/*.mseed'

list1 = glob.glob(working_dir + data1)
list2 = glob.glob(working_dir + data2)

stns1 = []
ntwks1 = []
chs1 = []

for e in list1:
    name = path.basename(e).split('/')[-1]
    stn = name.split('.')[0]
    ntwk = name.split('.')[1]
    ch = name.split('.')[2]
    
    stns1.append(stn)
    ntwks1.append(ntwk)
    chs1.append(ch)

stns2 = []
ntwks2 = []
chs2 = []    

for e in list2:
    name = path.basename(e).split('/')[-1]
    stn = name.split('.')[1]
    ntwk = name.split('.')[0]
    ch_time = name.split('.')[3]
    ch = ch_time.split('__')[0]
    
    stns2.append(stn)
    ntwks2.append(ntwk)
    chs2.append(ch)
    
    
df1['stn'] = stns1
df1['ntwk'] = ntwks1
df1['ch1'] = chs1

df2['stn'] = stns2
df2['ntwk'] = ntwks2
df2['ch2'] = chs2

merge = pd.merge(df1,df2, how = 'inner', on = ['stn','ntwk'])

    
    