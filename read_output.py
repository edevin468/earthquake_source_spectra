#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 11:59:01 2022

@author: emmadevin
"""

import os
import pkg_resources
from gmprocess.io.asdf.stream_workspace import StreamWorkspace


data_path =  '/Users/emmadevin/Work/USGS 2021/Data/gmprocess/gmprocess_python/data/ci38443095/workspace.h5'
outpath = '/Users/emmadevin/Work/USGS 2021/Data/gmprocess/gmprocess_python/data/ci38443095/plots/'
workspace = StreamWorkspace.open(data_path)

ds = workspace.dataset
stn_list = ds.waveforms.list()

for stn in stn_list:
    
    sc = workspace.getStreams(
      'ci38443095', stations=[stn], labels=['default'])

    sc.describe()
    
    sta_st = sc[0]
    check = sta_st.passed
    print(check)
    
    if check == True: c = 'b'
    elif check == False: c = 'r'
    
    st = ds.waveforms[stn]['ci38443095_default']
    st.plot(color = c, outfile = outpath+stn+'.png' )
    