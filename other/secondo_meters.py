#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 11:36:09 2017

@author: Alexis Klimasewski
Edited by Emma Devin for Ridgecrest dataset

inputs: reads in the text files of NE average spectra from record_spectra

method: reads in record spectra and corrects for distances using dread function, 
uses the Andrews 86 method and svd to compute station and event spectra from record spectra, 
plots the event and station spectra, then if a constraint station or event file is specified, 
the output files are rewritten for the constraint and plots are generated again, 
if no constraint, comment out that part of the code

outputs: writes spectra text files for each event and station in event_site_spectra directory
"""
import glob
import os.path as path
import numpy as np
import obspy
from obspy import read    
import dread
import time
import random
import pandas as pd

#read in the cut and corrected spectra = records
#records are in m/s
#should be binned frequencies and amplitudes
#make list of records and the corresponding events and stations

# working directory and outfil path for inversion
working_dir = '/Users/emmadevin/Work/USGS 2021/Data/Prelim+'
outfile_path = working_dir + '/Andrews_inversion'

# df with station locations
stations = pd.read_csv(working_dir + 'Station_info/station_locs.csv')

#list of record files
ev = glob.glob(working_dir + '/record_spectra/*/*')


#===================================================================#

def secondo(record_path, out_file_path):
    print('Number of records: ', len(record_path))
   
    
    # get lists of station ids and station locations
    stn_list = (stations['network']+stations['station']).tolist()
    stn_lat = stations['latitude']
    stn_lon = stations['longitude']
    
    
    eventidlist = []
    event_lat = []
    event_lon = []
    event_depth = []
    record_freq = []
    record_spec = []
    record_std = []
    
    t1 = time.time()

    for i in range(len(record_path)):
        
        # read record filename and extract event and station identifiers
        record = (record_path[i].split('/')[-1])
        base = path.basename(record)
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
        
        # assign station corrdinates and depth
        stlat = stn_lat[stn_list.index(stn_id)]
        stlon = stn_lon[stn_list.index(stn_id)]
        stdepth = 0
        
        #find distance between event and station
        dist =  dread.compute_rrup(evlon, evlat, evdepth, stlon, stlat, stdepth) #in km
        
        #km to m
        dist = dist*1000.

        #read in spectra file
        data = np.genfromtxt(record_path[i], dtype = float, comments = '#', delimiter = None, usecols = (0,1,2))  #only read in first two cols
        record_freq.append(data[:,0])
        
        # data is NE spectra; square for power spectra; correct for geometrical spreading 
        record_spec.append((data[:,1]*dist)**2.)
        
        # propagation here for errors going lin to log power
        std_prop = np.abs(2*(data[:,2])/(data[:,1]*np.log(10)))
        record_std.append(std_prop) #in log power here
        
        #if event info not part of lists yet add
        if eventid not in eventidlist:
            eventidlist.append(eventid)
            event_lat.append(evlat)
            event_lon.append(evlon)
            event_depth.append(evdepth)
            
    t2 = time.time()
    
    print('Time to read and distance correct all records: ', (t2-t1)/60.)
    print('Number of records (frequency): ', len(record_freq))
    print('Number of records (amplitude): ',len(record_spec))

    freq_list = record_freq[0]
    print(freq_list)
    F_bins = len(freq_list)
    print(F_bins)
    
    rows = len(record_path) #testing first 10
    print(rows)
    
    index_matrix = [[0 for j in range(3)] for i in range(rows)]
    
    for i in range(len(record_path)):
    # for i in range(rows):
    # for eventid in eventidlist:
        record = record_path[i].split('/')[-1]
        base = path.basename(record)
        
        print(base)
        
        network = base.split('_')[0]
        station = base.split('_')[1]
        stnid = network+station
        eventid = base.split('_')[-1]
        eventid = eventid.split('.')[0]
        
        print(network, station)
        print(eventid)
        
        
        
        #make a tuple of record, event, station so indices can be assigned
        index_matrix[i] = [base, eventidlist.index(eventid), stn_list.index(stnid)]
    
    print(eventidlist[0])
    print(stn_list)
    
    I = len(eventidlist)#events
    J = len(stn_list)#stations
    K = len(record_path)#records
    K = rows
    
    print('Number of events: ', I, ' Number of stations: ', J)
    print('Number of rows (records): ', K, ' Number of cols (events+stations): ', I+J)
    
    #make the G matrix of 1s and 0s and R matrix of records
    G1 = np.zeros((K,I))
    G2 = np.zeros((K,J))

    
    for k in range(K):#for all records
        G1[k][index_matrix[k][1]] = 1 #record row, eventid col
        G2[k][index_matrix[k][2]] = 1 #record row, station col
    
    
    G = np.concatenate((G1,G2), axis = 1)
    
    print(G)
    
    R = np.zeros((K,F_bins))
    cov = np.zeros((K,F_bins))
    
    #populate R matrix with the log of the record spectra (power record spectra)
    for k in range(K):#for all records
        #each row is a record, and col is a frequency band
        #set row equal to the that spectral array
        #here we take the log
        R[k:,] = np.log10(record_spec[k])
        cov[k:,] = record_std[k]
        
    m1 = np.zeros((I+J, F_bins))
    m_cov = np.zeros((I+J, F_bins))
    
    #do the inversion for each freq
    for f in range(F_bins):
        t1 = time.time()
        d = R[:,f]#record for given frequency col
        dT = d.T
        print('inverting for frequency: ', f, freq_list[f])
        G_inv = np.linalg.pinv(G, rcond=1e-13)
        covd = np.diag(cov[:,f])

        covm = np.dot((np.dot(G_inv, covd)), G_inv.T)
        m1[:,f] = np.dot(G_inv,dT)
        m_cov[:,f]= covm.diagonal()
        t2 = time.time()
        print('time for inversion: (min) ', round((t2-t1)/60., 4))
    
    
    print(m1.shape)
    #now split m into an event matrix and a station matrix
    event = m1[0:I,:] #take first I rows
    station = m1[I:I+J,:]
    event_cov = m_cov[0:I,:]
    station_cov = m_cov[I:I+J,:]
    print(event.shape, station.shape)
    print(event_cov.shape, station_cov.shape)

    for i in range(I):#for each event
        #go from the log of the power spectra to the regular spectra in m
        amp = np.sqrt(np.power(10.0, event[i,:]))

        std = (np.sqrt(np.abs(event_cov[i,:])/2.)*((amp)*(np.log(10))))
        
        outfile = open(outfile_path + '/' + eventidlist[i] + '.out', 'w')
        out = (np.array([freq_list, amp, std])).T
        outfile.write('#freq_bins \t vel_spec_NE_m \t stdev_m \n')
        np.savetxt(outfile, out, fmt=['%E', '%E', '%E'], delimiter='\t')
        outfile.close()


    print(outfile_path)
    for i in range(J):#for each station
        amp = np.sqrt(np.power(10.0, station[i,:]))

        std1 = np.sqrt((station_cov[i,:]))
        std = np.abs((std1/2.)*(amp)*(np.log(10)))
        outfile = open(outfile_path + '/' + stn_list[i] + '.out', 'w')
        out = (np.array([freq_list, amp, std])).T
        outfile.write('#freq_bins \t vel_spec_NE_m \t stdev_m \n')
        np.savetxt(outfile, out, fmt=['%E', '%E', '%E'], delimiter='\t')
        outfile.close()

# ===================================================================#


secondo(record_path = ev, out_file_path = outfile_path)
