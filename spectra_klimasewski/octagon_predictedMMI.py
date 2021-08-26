#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 25 12:21:13 2021

@author: emmadevin
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from math import sin, cos, sqrt, atan2, radians

def haversine_distance(lat1, lon1, lat2, lon2):
   r = 6371
   phi1 = np.radians(lat1)
   phi2 = np.radians(lat2)
   delta_phi = np.radians(lat2 - lat1)
   delta_lambda = np.radians(lon2 - lon1)
   a = np.sin(delta_phi / 2)**2 + np.cos(phi1) * np.cos(phi2) *   np.sin(delta_lambda / 2)**2
   res = r * (2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a)))
   return np.round(res, 2)

df = pd.read_csv('/Users/emmadevin/Work/USGS Summer 2020/USGS PLUM/EEW Example/ShakeAlert_Northridge_Predicted_MMI.csv')

eventlat = 34.213
eventlon = -118.53

mmi = df['MMI'].round(1)

keep = mmi == 3.5

lat = np.array(df['lat'])
lon = np.array(df['lon'])

lat = lat[keep]
lon = lon[keep]

dist = haversine_distance(eventlat,eventlon,lat, lon)

avgdist = np.mean(dist)

x_s, y_s = [], []
angles = [0,np.pi/4,np.pi/2,np.pi/4*3,np.pi,np.pi/4*5,np.pi/2*3,np.pi/4*7]

for angle in angles:
    x_dist = avgdist*np.cos(angle)   
    y_dist = avgdist*np.sin(angle)
    
    x_dist_deg = x_dist*360/(40075*np.cos(eventlat))
    y_dist_deg = y_dist/111.32
    
    x_s.append(eventlon + x_dist_deg)
    y_s.append(eventlat + y_dist_deg)


plt.scatter(x_s,y_s)

df2 = pd.DataFrame()
df2['lat']=y_s
df2['lon']=x_s

df2.to_csv('/Users/emmadevin/Work/USGS Summer 2020/USGS PLUM/EEW Example/oct_pred_35.csv')
