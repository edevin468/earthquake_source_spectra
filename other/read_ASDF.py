#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 29 11:40:25 2021

@author: emmadevin
"""

import pyasdf as pf
working_dir = '/Users/emmadevin/Work/USGS 2021/Data/gmprocess_test/data/ci38443095'

df = pf.ASDFDataSet(working_dir + '/workspace.h5')

sts = df.waveforms.list()

for st in sts:
    print(st)
    
    sta = df.waveforms[st]
    print(sta.get_waveform_tags())
    


for obj in df.auxiliary_data.TraceProcessingParameters: 
    print(obj)
    
for obj in df.auxiliary_data.Cache: 
    print(obj)
    # net = obj.split('.')[0]
    # stat = obj.split('.')[1]
    # obj.name = net +'_'+stat
    
    # print(df.auxiliary_data.TraceProcessingParameters.obj)