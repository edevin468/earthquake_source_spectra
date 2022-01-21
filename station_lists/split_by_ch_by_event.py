#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 19 10:17:56 2022

@author: emmadevin
"""

import pandas as pd


df = pd.read_csv('/Users/emmadevin/Work/USGS 2021/Data/gmprocess/qa_processing/data/ci38443095/stns_ch.csv')

data = [tuple(r) for r in df.to_numpy().tolist()]


d_hn = []
d_hh = []
d_eh = []

for e in data: 
    if 'HN' in e:
        d_hn.append(e)
    elif 'HH' in e:
        d_hh.append(e)
    elif 'EH' in e:
        d_eh.append(e)
    else: print('Something wrong: ' + e)
    
df_HN = pd.DataFrame(d_hn, columns = ['', 'network', 'station', 'channel'])
df_HH = pd.DataFrame(d_hh, columns = ['', 'network', 'station', 'channel'])
df_EH = pd.DataFrame(d_eh, columns = ['', 'network', 'station', 'channel'])

df_HN = df_HN.drop(columns='')
df_HH = df_HH.drop(columns='')
df_EH = df_EH.drop(columns='')

df_HN.to_csv('/Users/emmadevin/Work/USGS 2021/Data/gmprocess/qa_processing/data/ci38443095/stns_HN.csv')
df_HH.to_csv('/Users/emmadevin/Work/USGS 2021/Data/gmprocess/qa_processing/data/ci38443095/stns_HH.csv')
df_EH.to_csv('/Users/emmadevin/Work/USGS 2021/Data/gmprocess/qa_processing/data/ci38443095/stns_EH.csv')