#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 10 09:33:45 2021

@author: emmadevin
"""

import pandas as pd

small = pd.read_csv('/Users/emmadevin/Work/USGS 2021/Data/Prelim+/Station_info/station_counts.csv')
big = pd.read_csv('/Users/emmadevin/Work/USGS 2021/Data/Prelim+/Station_info/station_locs.csv')


stn_list = small['station'].tolist()
big_list = big['station'].tolist()
networks = big['network'].tolist()
latitudes = big['latitude']
longitudes = big['longitude']

small_lats = []
small_lons = []
ntwks = []

for stn in stn_list: 

    index = big_list.index(stn)
    lat = latitudes[index]
    lon = longitudes[index]
    ntwk = networks[index]
    
    small_lats.append(lat)
    small_lons.append(lon)
    ntwks.append(ntwk)
    
small['latitude'] = small_lats
small['longitude'] = small_lons
small['network'] = ntwks

small.to_csv('/Users/emmadevin/Work/USGS 2021/Data/Prelim+/Station_info/station_counts.csv')
    
    