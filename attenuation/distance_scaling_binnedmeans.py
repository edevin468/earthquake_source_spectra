#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 16 12:11:25 2021

@author: emmadevin
"""

import glob
import os.path as path
import numpy as np
import dread
import time
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt

fig, (ax1,ax2) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [3.3,4]},figsize=(7,4))
plt.style.use('classic')
fig.patch.set_facecolor('white')

#read in the cut and corrected spectra = records
#records are in m/s
#should be binned frequencies and amplitudes
#make list of records and the corresponding events and stations

# working directory and outfil path for inversion
working_dir = '/Users/emmadevin/Work/USGS 2021/Data/Prelim'
outfile_path = working_dir + '/Andrews_inversion'

# df with station locations
stations = pd.read_csv(working_dir + '/station_locs.csv')

# list of record files
record_list = glob.glob(working_dir + '/record_spectra/*/*')

# list of files
   
    
# get lists of station ids and station locations
stn_list = (stations['network']+stations['station']).tolist()
stn_lat = stations['latitude']
stn_lon = stations['longitude']


df = pd.DataFrame()

amp1 = []
amp2 = []
amp3 = []
amp4 = []
amp5 = []
f1 = 0.01
f2 = 0.1
f3 = 1.0
f4 = 10
distlist = []
maglist = []


for record in record_list: 
    
    base = record.split('/')[-1]
    eventid = base.split('_')[-1]
    eventid = eventid.split('.')[0]
    ntwk = base.split('_')[0]
    stn = base.split('_')[1]
    stn_id = ntwk + stn
    
    # print some updates to track progress
    print('Event: ', eventid)
    print('Station: ', stn)
    
    # read phase file fo get event locations
    phase_file = working_dir + '/RC_phase_beta/' + eventid + '.phase'
    phase = pd.read_csv(phase_file, sep = '\s+', index_col=0, nrows = 0).columns.tolist()
    
    # assign event coordinates and depth
    evlat = float(phase[4])
    evlon = float(phase[5])
    evdepth = float(phase[6])
    mag = int(float(phase[7]))
    
    # assign station corrdinates and depth
    stlat = stn_lat[stn_list.index(stn_id)]
    stlon = stn_lon[stn_list.index(stn_id)]
    stdepth = 0
    
    #find distance between event and station
    dist =  dread.compute_rrup(evlon, evlat, evdepth, stlon, stlat, stdepth) #in km
    distlist.append(dist)
    
    #read in record
    data = np.genfromtxt(record, dtype = float, comments = '#', delimiter = None, usecols = (0,1,2))  #only read in first two cols
    
    amp1.append(data[7,1])
    amp2.append(data[30,1])
    amp3.append(data[52,1])
    amp4.append(data[74,1])
    
    maglist.append(mag)
    
sizes = np.array(maglist)*10
# create list of colors
colors = []
for mag in maglist:
    if mag == 1:   colors.append('r')  
    elif mag == 2:  colors.append('b')  
    elif mag == 3:  colors.append('g')  
    elif mag == 4:  colors.append('k')  
    elif mag == 5:  colors.append('orange')  
    elif mag == 6:  colors.append('cyan')  
    elif mag == 7:  colors.append('m')  
    elif mag == 8:  colors.append('purple')  
    elif mag == 9:  colors.append('grey')  
  
df['record'] = record_list
df['amp:f1'] = amp1
df['amp:f2'] = amp2
df['amp:f3'] = amp3
df['amp:f4'] = amp4
df['mag'] = maglist

df['dist'] = distlist

edges = np.logspace(1,2,10)

bin_means1, bin_edges1, binnumber1 = stats.binned_statistic(df['dist'], df['amp:f1'], statistic='mean', bins=edges)
bin_width1 = (bin_edges1[1] - bin_edges1[0])
bin_centers1 = bin_edges1[1:] - bin_width1/2

bin_means2, bin_edges2, binnumber2 = stats.binned_statistic(df['dist'], df['amp:f2'], statistic='mean', bins=edges)
bin_width2 = (bin_edges2[1] - bin_edges2[0])
bin_centers2 = bin_edges2[1:] - bin_width2/2

bin_means3, bin_edges3, binnumber3 = stats.binned_statistic(df['dist'], df['amp:f3'], statistic='mean', bins=edges)
bin_width3 = (bin_edges3[1] - bin_edges3[0])
bin_centers3 = bin_edges3[1:] - bin_width3/2

bin_means4, bin_edges4, binnumber4 = stats.binned_statistic(df['dist'], df['amp:f4'], statistic='mean', bins=edges)
bin_width4 = (bin_edges4[1] - bin_edges4[0])
bin_centers4 = bin_edges4[1:] - bin_width4/2



ax1.scatter(distlist,amp2, c = maglist,edgecolors='none', cmap = 'winter',alpha = 0.6, marker = 'o',s = sizes)
# plt.scatter(bin_centers1,bin_means1, edgecolors = 'k', facecolors = 'darkgrey',s=50,label = 'binned means')


x = np.logspace(0.1,3,100)
y1 = 1/10*1/x
y2 = 1/100*1/x
y3 = 1/1000*1/x
y4 = 1/10000*1/x
y5 = 1/100000*1/x
y6 = 1/1000000*1/x
y7 = 1/10000000*1/x


ax1.plot(x,y1,ls ='--', c = 'k',label = '~1/r')
ax1.plot(x,y2,ls ='--', c = 'k')
ax1.plot(x,y3,ls ='--', c = 'k')
ax1.plot(x,y4,ls ='--', c = 'k')
ax1.plot(x,y5, ls ='--',c = 'k')
ax1.plot(x,y6, ls ='--',c = 'k')
# ax1.plot(x,y7, c = 'k')
ax1.legend(fontsize = 10, loc='lower left')

ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlim(6,11**2)
ax1.set_ylim(5*10**-9, 5*10**-2)
ax1.set_xlabel('distance (km)')
ax1.set_ylabel('amplitude')
ax1.set_title('(a) f = '+str(f2)+' Hz', loc= 'left', fontsize= 12)

# leg = plt.legend(loc='upper center', bbox_to_anchor=(1.1, -0.2), facecolor = 'w', ncol = 4, fontsize = 10)



im = ax2.scatter(distlist,amp3, c = maglist,edgecolors='none', cmap = 'winter',alpha = 0.6, marker = 'o',s = sizes)
# plt.scatter(bin_centers3,bin_means3, edgecolors = 'k', facecolors = 'darkgrey',s=50,label = 'binned means')

x = np.logspace(0.1,3,100)
y1 = 1/10*1/x
y2 = 1/100*1/x
y3 = 1/1000*1/x
y4 = 1/10000*1/x
y5 = 1/100000*1/x
y6 = 1/1000000*1/x
y7 = 1/10000000*1/x


# ax2.plot(x,y1, c = 'k',label = '~1/r' )
ax2.plot(x,y2,ls ='--', c = 'k',label = '~1/r')
ax2.plot(x,y3, ls ='--',c = 'k')
ax2.plot(x,y4, ls ='--',c = 'k')
ax2.plot(x,y5,ls ='--', c = 'k')
ax2.plot(x,y6, ls ='--',c = 'k')
# ax2.plot(x,y7, c = 'k')
ax2.legend(fontsize = 10, loc='lower left')

ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlim(6,11**2)
ax2.set_ylim(10**-9, 10**-2)
ax2.set_xlabel('distance (km)')
# ax2.set_yticks([])
cbar = fig.colorbar(im)

cbar.set_label('magnitude', rotation=270, labelpad=20)

ax2.set_title('(b) f = '+str(f3)+' Hz', loc= 'left', fontsize= 12)



# leg = plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), facecolor = 'w', ncol = 2, fontsize = 10)

plt.savefig('/Users/emmadevin/Work/USGS 2021/Figures/SCEC/attenuation.pdf', bbox_inches='tight')

