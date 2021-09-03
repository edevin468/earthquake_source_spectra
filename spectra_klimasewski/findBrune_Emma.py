#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 09:53:15 2017

@author: Alexis Klimasewski

compute “B” for a M3 earthquake, and then for every ~M3 earthquake
in your dataset (call it “A”), find A/B for every frequency range you inverted for.  
Then sum up A/B over all frequency ranges and find which earthquake has the minimum value for that
because then that would suggest that earthquake is closest to a brune spectrum
"""

import glob
import numpy as np
import obspy
import matplotlib.pyplot as plt
import pandas as pd

#properties for unit
beta = 3500. #3500m/s
stressdrop = 5e3 #5e6 pascals
U = 0.38#0.63
rho = 2750. #2750kg/m^3

working_dir =  '/Users/emmadevin/Work/USGS 2021/Data/Prelim'
event_spectra_dir = working_dir + '/Andrews_inversion/'
event_list = glob.glob(event_spectra_dir + '3*.out')

writefile = 'yes'


    
#compute Brune spectra for all of the events in directory
cf_list = []
cf2_list = []
Brune_list = []
spec_list = []
ev_list = []
mag_list = []

spec_demean_list = []
Brune_demean_list = []

spec_demean_list_log = []
Brune_demean_list_log = []

event_spectra = []

for event in event_list:
    
    filename = (event.split('/')[-1])
    event_id = filename.split('.')[0]
    
    # obtain event magnitudes
    phase_file = working_dir + '/RC_phase_beta/' + event_id + '.phase'
    phase = pd.read_csv(phase_file, sep = '\s+', index_col=0, nrows = 0).columns.tolist()
    
    mag = float(phase[7])
    
    m_l = 0
    m_u = 10
    
    if mag <= m_u and mag >= m_l:
        keep = True
    else:
        keep = False
        
    if keep == True:
        event_spectra.append(event)
    else:
        continue
    

    
for event in event_spectra:
    
    filename = (event.split('/')[-1])
    event_id = filename.split('.')[0]
    print(event_id)
    
    # obtain event magnitudes
    phase_file = working_dir + '/RC_phase_beta/' + event_id + '.phase'
    phase = pd.read_csv(phase_file, sep = '\s+', index_col=0, nrows = 0).columns.tolist()
    
    mag = float(phase[7])
    print(mag)
    mag_list.append(mag)
    
    
    data = np.genfromtxt(event, dtype = float, comments = '#', delimiter = None, usecols = (0,1))#only read in first two cols
    freq = data[:,0]
    spec = (data[:,1])  # these record spectra are in m


    # if less than 3, convert local magnitude to moment magnitude
    if mag < 3.0:
        M = 0.884 + 0.754*mag  # 0.884 + 0.667*ml, 754
    else:
        M = mag
        
    # compute Brune in SI units
    # moment from the moment magnitude
    M0 = 10.**((3./2.)*M + 9.1)
    
    #corner frequency
    fc = beta*(stressdrop/(8.47*M0))**(1./3.)
    omega0 = (M0*U)/(4.*rho*np.pi*(beta**(3.0)))
    
    #brune spectra over all frequencies
    Brune = (2.*np.pi*(freq)*omega0)/(1.+((1./fc)*freq)**2.)
    
    #stay in meters
    shift1 = np.mean(Brune[27:74])
    shift2 = np.mean(spec[27:74])
    
    cf_list.append(np.log10(spec/shift2)-np.log10(Brune/shift1))

    Brune_list.append(Brune)
    spec_list.append(spec)
    

#for each event, find A/B for all other events
#sum up all A/B over freqencies we are fitting
cfarray = np.array(cf_list)

# ind = cfarray.index(0.5)
sum_list =list(map(sum,cfarray[:,np.arange(27,74)]**2.0)) # found the best fit from 1-32.7Hz

# find the minimum in log space
ind = sum_list.index(min(sum_list))
print(event_spectra[ind])
print(min(sum_list))        

print(mag_list[ind])
fig = plt.figure(figsize = (12,10))
plt.style.use('classic')
fig.patch.set_facecolor('white')


plt.ylabel('Velocity amplitude (m)', fontsize = 16)
plt.xlim(0.5,70)
plt.loglog(freq , spec_list[ind], color = 'green', label = 'event spectra')
plt.grid()
plt.loglog(freq, Brune_list[ind], color = 'blue', label = 'Brune spectra')
plt.legend(loc = 'lower right', fontsize = 16)
plt.xlabel('Frequency (Hz)', fontsize = 16)
plt.title(event_spectra[ind], fontsize = 16)
plt.tick_params(axis='both', which='major', labelsize=15)
plt.tick_params(axis='both', which='both', length = 5, width = 1)
# plt.text(0.7, .1, 'Median log(diff) 1-32.7 Hz (demeaned): ' + str(round(sum_list[ind],3)), fontsize = 16)
plt.show()



# #write the constraint file in linear space to agree with the event and station spectra
# if writefile == 'yes':
#     outfile = open(working_dir + 'test_codes/constraint_' + ev_list[ind] + '.out', 'w')
#     out = (np.array([freq, (10.**(cf_list[ind]))]).T)
#     outfile.write('#freq_bins \t cf_m \n')
#     np.savetxt(outfile, out, fmt=['%E', '%E'], delimiter='\t')
#     outfile.close()

