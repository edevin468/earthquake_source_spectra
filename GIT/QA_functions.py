#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 10:15:00 2021

funtions for QA processing
code mostly borrowed from gmprocess
modified somewhat

@author: emmadevin
"""

from obspy.signal.trigger import classic_sta_lta
import obspy as op
from obspy import read
import os
import os.path as path
import glob
import warnings
import time
from progress_bar import progress_bar
import logging

import numpy as np
import pandas as pd

from openquake.hazardlib.gsim.base import RuptureContext
from openquake.hazardlib import const
from openquake.hazardlib import imt

from obspy.geodetics.base import gps2dist_azimuth

from gmprocess.waveform_processing.phase import (
    pick_power, pick_ar, pick_baer, pick_kalkan, pick_travel)
from gmprocess.utils.config import get_config
from gmprocess.metrics.station_summary import StationSummary
from gmprocess.utils.models import load_model


def check_free_field(st, reject_non_free_field=False):
    """
    Checks free field status of stream.
    Args:
        st (obspy.core.stream.Stream):
            Stream of data.
        reject_non_free_field (bool):
            Should non free-field stations be failed?
    Returns:
        Stream that has been checked for free field status.
    """
    l = []
    for trace in st:
        if not trace.free_field and reject_non_free_field:
            l.append(False)
        else: 
            l.append(True)
            continue

    if False in l:
        result = False
    else: 
        result = True
        
        
    return result





def check_sta_lta(st, sta_length=1.0, lta_length=20.0, threshold=5.0):
    '''
    Checks that the maximum STA/LTA ratio for AT LEAST ONE of the stream's
    traces is above a certain threshold.
    Args:
        st (obspy.core.stream.Stream):
            Stream of data.
        sta_length (float):
            Length of time window for STA (seconds).
        lta_length (float):
            Length of time window for LTA (seconds).
        threshold (float):
            Required maximum STA/LTA ratio to pass the test.
    Returns:
        Stream that has been checked for sta/lta requirements.
    '''
    l = []
    for tr in st:
        sr = tr.stats.sampling_rate
        nlta = lta_length * sr + 1
        
        
        if len(tr) >= nlta:
            sta_lta = classic_sta_lta(tr.data, sta_length * sr + 1, nlta)
            if sta_lta.max() < threshold:
                result = False
                # print('Failed sta/lta check because threshold sta/lta '
                #         'is not exceeded.')
            else: result = True
        else:
            result = False
            # print('Failed sta/lta check because record length is shorter '
            #         'than lta length.')
            
        l.append(result)
    
    if False in l:
        return False
    else: 
        return True

def signal_split(
        st, origin, model=None,
        picker_config=None,
        config=None):
    """
    This method tries to identifies the boundary between the noise and signal
    for the waveform. The split time is placed inside the
    'processing_parameters' key of the trace stats.
    The P-wave arrival is used as the split between the noise and signal
    windows. Multiple picker methods are suppored and can be configured in the
    config file
    '~/.gmprocess/picker.yml
    Args:
        st (StationStream):
            Stream of data.
        origin (ScalarEvent):
            ScalarEvent object.
        model (TauPyModel):
            TauPyModel object for computing travel times.
        picker_config (dict):
            Dictionary containing picker configuration information.
        config (dict):
            Dictionary containing system configuration information.
    Returns:
        trace with stats dict updated to include a
        stats['processing_parameters']['signal_split'] dictionary.
    """
    if picker_config is None:
        picker_config = get_config(section='pickers')
    if config is None:
        config = get_config()

    loc, mean_snr = pick_travel(st, origin, model)
    if loc > 0:
        tsplit = st[0].stats.starttime + loc
        preferred_picker = 'travel_time'
    else:
        pick_methods = ['ar', 'baer', 'power', 'kalkan']
        columns = ['Stream', 'Method', 'Pick_Time', 'Mean_SNR']
        df = pd.DataFrame(columns=columns)
        for pick_method in pick_methods:
            try:
                if pick_method == 'ar':
                    loc, mean_snr = pick_ar(
                        st, picker_config=picker_config, config=config)
                elif pick_method == 'baer':
                    loc, mean_snr = pick_baer(
                        st, picker_config=picker_config, config=config)
                elif pick_method == 'power':
                    loc, mean_snr = pick_power(
                        st, picker_config=picker_config, config=config)
                elif pick_method == 'kalkan':
                    loc, mean_snr = pick_kalkan(st,
                                                picker_config=picker_config,
                                                config=config)
                elif pick_method == 'yeck':
                    loc, mean_snr = pick_kalkan(st)
            except BaseException:
                loc = -1
                mean_snr = np.nan
            row = {
                'Stream': st.get_id(),
                'Method': pick_method,
                'Pick_Time': loc,
                'Mean_SNR': mean_snr
            }
            df = df.append(row, ignore_index=True)

        max_snr = df['Mean_SNR'].max()
        if not np.isnan(max_snr):
            maxrow = df[df['Mean_SNR'] == max_snr].iloc[0]
            tsplit = st[0].stats.starttime + maxrow['Pick_Time']
            preferred_picker = maxrow['Method']
        else:
            tsplit = -1

    # the user may have specified a p_arrival_shift value.
    # this is used to shift the P arrival time (i.e., the break between the
    # noise and signal windows).
    shift = 0.0
    if 'p_arrival_shift' in picker_config:
        shift = picker_config['p_arrival_shift']
        if tsplit + shift >= st[0].stats.starttime:
            tsplit += shift

    if tsplit >= st[0].stats.starttime:
        # Update trace params
        split_params = {
            'split_time': tsplit,
            'method': 'p_arrival',
            'picker_type': preferred_picker
        }
        for tr in st:
            tr.setParameter('signal_split', split_params)

    return st


def signal_end(st, event_time, event_lon, event_lat, event_mag,
               method=None, vmin=None, floor=None,
               model=None, epsilon=2.0):
    """
    Estimate end of signal by using a model of the 5-95% significant
    duration, and adding this value to the "signal_split" time. This probably
    only works well when the split is estimated with a p-wave picker since
    the velocity method often ends up with split times that are well before
    signal actually starts.
    Args:
        st (StationStream):
            Stream of data.
        event_time (UTCDateTime):
            Event origin time.
        event_mag (float):
            Event magnitude.
        event_lon (float):
            Event longitude.
        event_lat (float):
            Event latitude.
        method (str):
            Method for estimating signal end time. Either 'velocity'
            or 'model'.
        vmin (float):
            Velocity (km/s) for estimating end of signal. Only used if
            method="velocity".
        floor (float):
            Minimum duration (sec) applied along with vmin.
        model (str):
            Short name of duration model to use. Must be defined in the
            gmprocess/data/modules.yml file.
        epsilon (float):
            Number of standard deviations; if epsilon is 1.0, then the signal
            window duration is the mean Ds + 1 standard deviation. Only used
            for method="model".
    Returns:
        trace with stats dict updated to include a
        stats['processing_parameters']['signal_end'] dictionary.
    """
    # Load openquake stuff if method="model"
    if method == "model":
        dmodel = load_model(model)

        # Set some "conservative" inputs (in that they will tend to give
        # larger durations).
        rctx = RuptureContext()
        rctx.mag = event_mag
        rctx.rake = -90.0
        rctx.vs30 = np.array([180.0])
        rctx.z1pt0 = np.array([0.51])
        dur_imt = imt.from_string('RSD595')
        stddev_types = [const.StdDev.TOTAL]

    for tr in st:
        if not tr.hasParameter('signal_split'):
            continue
        if method == "velocity":
            if vmin is None:
                raise ValueError('Must specify vmin if method is "velocity".')
            if floor is None:
                raise ValueError('Must specify floor if method is "velocity".')
            epi_dist = gps2dist_azimuth(
                lat1=event_lat,
                lon1=event_lon,
                lat2=tr.stats['coordinates']['latitude'],
                lon2=tr.stats['coordinates']['longitude'])[0] / 1000.0
            end_time = event_time + max(floor, epi_dist / vmin)
        elif method == "model":
            if model is None:
                raise ValueError('Must specify model if method is "model".')
            epi_dist = gps2dist_azimuth(
                lat1=event_lat,
                lon1=event_lon,
                lat2=tr.stats['coordinates']['latitude'],
                lon2=tr.stats['coordinates']['longitude'])[0] / 1000.0
            # Repi >= Rrup, so substitution here should be conservative
            # (leading to larger durations).
            rctx.rrup = np.array([epi_dist])
            rctx.sids = np.array(range(np.size(rctx.rrup)))
            lnmu, lnstd = dmodel.get_mean_and_stddevs(
                rctx, rctx, rctx, dur_imt, stddev_types)
            duration = np.exp(lnmu + epsilon * lnstd[0])
            # Get split time
            split_time = tr.getParameter('signal_split')['split_time']
            end_time = split_time + float(duration)
        else:
            raise ValueError('method must be either "velocity" or "model".')
        # Update trace params
        end_params = {
            'end_time': end_time,
            'method': method,
            'vsplit': vmin,
            'floor': floor,
            'model': model,
            'epsilon': epsilon
        }
        tr.setParameter('signal_end', end_params)

    return st


def check_zero_crossings(st, min_crossings=1.0):
    """
    Check for a large enough density.
    This is intended to screen out instrumental failures or resetting.
    Value determined empirically from observations on the GeoNet network
    by R Lee.
    Args:
        st (StationStream):
            StationStream object.
        min_crossings (float):
            Minimum average number of zero crossings per second for the full
            trace.
    """

    zero_count_tr = []
    delta_t = st[0].stats.delta
    dur = (st[0].stats.npts - 1) * delta_t

    for trace in st:
        # Make a copy of the trace to trim it before counting crossings; we do
        # not want to modify the trace but we only want to count the crossings
        # within the trimmed window
        
        tr = trace.copy()
        df = trace.stats.sampling_rate
        cft = classic_sta_lta(trace.data, int(5 * df), int(10 * df))
        
        print(len(cft))

        # if tr.hasParameter('signal_end'):
        #     etime = tr.getParameter('signal_end')['end_time']
        #     split_time = tr.getParameter('signal_split')['split_time']

        #     sig_start = int((split_time - tr.stats.starttime) / tr.stats.delta)
        #     sig_end = int((etime - tr.stats.starttime) / tr.stats.delta)
        tr_data = tr.data[0:-0]

        zarray = np.multiply(tr_data[0:-1], tr_data[1:])
        zindices = [i for (i, z) in enumerate(zarray) if z < 0]
        zero_count_tr = len(zindices)

        z_rate = zero_count_tr / dur
        
    

        # Fail if zero crossing rate is too low
        if z_rate <= min_crossings:
            return False
        else:
            return True
    
    
    
    