#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 14:25:27 2020

@author: root
"""

import os
from collections import OrderedDict
import numpy as np
import obspy
import scipy.signal as sgn
import matplotlib.pyplot as plt 
import matplotlib.mlab as mlab
from scipy import integrate
from obspy.clients.earthworm import Client
from obspy import UTCDateTime
from obspy.signal.trigger import classic_sta_lta, recursive_sta_lta
from obspy.signal.trigger import plot_trigger, trigger_onset
from obspy import Stream
from obspy import Trace
from numpy import argmax


#kill the code: check saving correctly
#%% open arrays for stations

#inf
st1d = np.zeros(shape=(0,4))
st2d = np.zeros(shape=(0,4))
st3d = np.zeros(shape=(0,4))
st4d = np.zeros(shape=(0,4))
st5d = np.zeros(shape=(0,4))
st6d = np.zeros(shape=(0,4))
st7d = np.zeros(shape=(0,4))
st8d = np.zeros(shape=(0,4))
st9d = np.zeros(shape=(0,4))
st10d = np.zeros(shape=(0,4))
st11d = np.zeros(shape=(0,4))
st12d = np.zeros(shape=(0,4))
st13d = np.zeros(shape=(0,4))
st14d = np.zeros(shape=(0,4))
st15d = np.zeros(shape=(0,4))
st16d = np.zeros(shape=(0,4))
st17d = np.zeros(shape=(0,4))
st18d = np.zeros(shape=(0,4))
st19d = np.zeros(shape=(0,4))
st20d = np.zeros(shape=(0,4))
st21d = np.zeros(shape=(0,4))
st22d = np.zeros(shape=(0,4))
st23d = np.zeros(shape=(0,4))
st24d = np.zeros(shape=(0,4))
st25d = np.zeros(shape=(0,4))
st26d = np.zeros(shape=(0,4))
st27d = np.zeros(shape=(0,4))
st28d = np.zeros(shape=(0,4))
st29d = np.zeros(shape=(0,4))
st30d = np.zeros(shape=(0,4))
st31d = np.zeros(shape=(0,4))
st32d = np.zeros(shape=(0,4))
st33d = np.zeros(shape=(0,4))
#seis
#st34d = np.zeros(shape=(0,4))
#st35d = np.zeros(shape=(0,4))
#st36d = np.zeros(shape=(0,4))
#st37d = np.zeros(shape=(0,4))
#st38d = np.zeros(shape=(0,4))
#st39d = np.zeros(shape=(0,4))
#st40d = np.zeros(shape=(0,4))
#st41d = np.zeros(shape=(0,4))
#%% scan limits
st_num = 41

scan_start = UTCDateTime(2020,9,2,0,0,0).timestamp 
days = 40

Station_info = np.zeros(shape=(days,5,2))


for x in range(0,days,3):    
    
    Day_start = scan_start + (x*60*60*24) 

    
    st1a,st2a,st3a,st4a,st5a,st6a,st7a,st8a,st9a,st10a,st11a,st12a,st13a,st14a,st15a,st16a,st17a,st18a,st19a,st20a,st21a,st22a,st23a,st24a,st25a,st26a,st27a,st28a,st29a,st30a,st31a,st32a,st33a,station_activitya = get_all_Fuego_inf_stations_4d(Day_start)    
   
    all_stations = station_activitya #np.concatenate((station_activitya,seismic_activitya))
    
    print("")
    print('day', x+1, 'of', days)
    print(UTCDateTime(Day_start)+600)
#    print(len(station_activitya)+len(seismic_activitya),"stations running")
    print(len(station_activitya),"Inf stations running")
#    print(len(seismic_activitya),"Seis stations running")
    print("")

#%% loop over stations
    if len(all_stations) > 0 :                   
        for i in range(0,len(all_stations)):
            
            st_ID = int(all_stations[i,0])
         
            if st_ID ==1 :
                tr_day = st1a[0].slice(starttime = st1a[0].stats.starttime + 3, endtime = st1a[0].stats.endtime)
                
            if st_ID == 2:
                tr_day = st2a[0].slice(starttime = st2a[0].stats.starttime + 3, endtime = st2a[0].stats.endtime)
                
            if st_ID == 3:
                tr_day = st3a[0].slice(starttime = st3a[0].stats.starttime + 3, endtime = st3a[0].stats.endtime)
                
            if st_ID == 4:
                tr_day = st4a[0].slice(starttime = st4a[0].stats.starttime + 3, endtime = st4a[0].stats.endtime)
                
            if st_ID == 5:
                tr_day = st5a[0].slice(starttime = st5a[0].stats.starttime + 3, endtime = st5a[0].stats.endtime)
                
            if st_ID == 6:
                tr_day = st6a[0].slice(starttime = st6a[0].stats.starttime + 3, endtime = st6a[0].stats.endtime)
                
            if st_ID == 7:
                tr_day = st7a[0].slice(starttime = st7a[0].stats.starttime + 3, endtime = st7a[0].stats.endtime)
                
            if st_ID == 8:
                tr_day = st8a[0].slice(starttime = st8a[0].stats.starttime + 3, endtime = st8a[0].stats.endtime)
                
            if st_ID == 9:
                tr_day = st9a[0].slice(starttime = st9a[0].stats.starttime + 3, endtime = st9a[0].stats.endtime)
                
            if st_ID == 10:
                tr_day = st10a[0].slice(starttime = st10a[0].stats.starttime + 3, endtime = st10a[0].stats.endtime)
                
            if st_ID == 11:
                tr_day = st11a[0].slice(starttime = st11a[0].stats.starttime + 3, endtime = st11a[0].stats.endtime)
                
            if st_ID == 12:
                tr_day = st12a[0].slice(starttime = st12a[0].stats.starttime + 3, endtime = st12a[0].stats.endtime)
                
            if st_ID == 13:
                tr_day = st13a[0].slice(starttime = st13a[0].stats.starttime + 3, endtime = st13a[0].stats.endtime)
                
            if st_ID == 14:
                tr_day = st14a[0].slice(starttime = st14a[0].stats.starttime + 3, endtime = st14a[0].stats.endtime)
                
            if st_ID == 15:
                tr_day = st15a[0].slice(starttime = st15a[0].stats.starttime + 3, endtime = st15a[0].stats.endtime)
                
            if st_ID == 16:
                tr_day = st16a[0].slice(starttime = st16a[0].stats.starttime + 3, endtime = st16a[0].stats.endtime)
                
            if st_ID == 17:
                tr_day = st17a[0].slice(starttime = st17a[0].stats.starttime + 3, endtime = st17a[0].stats.endtime)
                
            if st_ID == 18:
                tr_day = st18a[0].slice(starttime = st18a[0].stats.starttime + 3, endtime = st18a[0].stats.endtime)
                
            if st_ID == 19:
                tr_day = st19a[0].slice(starttime = st19a[0].stats.starttime + 3, endtime = st19a[0].stats.endtime)
                
            if st_ID == 20:
                tr_day = st20a[0].slice(starttime = st20a[0].stats.starttime + 3, endtime = st20a[0].stats.endtime)
                
            if st_ID == 21:
                tr_day = st21a[0].slice(starttime = st21a[0].stats.starttime + 3, endtime = st21a[0].stats.endtime)
                
            if st_ID == 22:
                tr_day = st22a[0].slice(starttime = st22a[0].stats.starttime + 3, endtime = st22a[0].stats.endtime)
                
            if st_ID == 23:
                tr_day = st23a[0].slice(starttime = st23a[0].stats.starttime + 3, endtime = st23a[0].stats.endtime)
                
            if st_ID == 24:
                tr_day = st24a[0].slice(starttime = st24a[0].stats.starttime + 3, endtime = st24a[0].stats.endtime)
                
            if st_ID == 25:
                tr_day = st25a[0].slice(starttime = st25a[0].stats.starttime + 3, endtime = st25a[0].stats.endtime)
                
            if st_ID == 26:
                tr_day = st26a[0].slice(starttime = st26a[0].stats.starttime + 3, endtime = st26a[0].stats.endtime)
                
            if st_ID == 27:
                tr_day = st27a[0].slice(starttime = st27a[0].stats.starttime + 3, endtime = st27a[0].stats.endtime)
                
            if st_ID == 28:
                tr_day = st28a[0].slice(starttime = st28a[0].stats.starttime + 3, endtime = st28a[0].stats.endtime)
                
            if st_ID == 29:
                tr_day = st29a[0].slice(starttime = st29a[0].stats.starttime + 3, endtime = st29a[0].stats.endtime)
                
            if st_ID == 30:
                tr_day = st30a[0].slice(starttime = st30a[0].stats.starttime + 3, endtime = st30a[0].stats.endtime)
                
            if st_ID == 31:
                tr_day = st31a[0].slice(starttime = st31a[0].stats.starttime + 3, endtime = st31a[0].stats.endtime)
                
            if st_ID == 32:
                tr_day = st32a[0].slice(starttime = st32a[0].stats.starttime + 3, endtime = st32a[0].stats.endtime)
                
            if st_ID == 33:
                tr_day = st33a[0].slice(starttime = st33a[0].stats.starttime + 3, endtime = st33a[0].stats.endtime)

#            
#            if st_ID == 34:
#                tr_day = st34a[0].slice(starttime = st34a[0].stats.starttime + 3, endtime = st34a[0].stats.endtime)
#                
#            if st_ID == 35:
#                tr_day = st35a[0].slice(starttime = st35a[0].stats.starttime + 3, endtime = st35a[0].stats.endtime)
#                
#            if st_ID == 36:
#                tr_day = st36a[0].slice(starttime = st36a[0].stats.starttime + 3, endtime = st36a[0].stats.endtime)
#                
#            if st_ID == 37:
#                tr_day = st37a[0].slice(starttime = st37a[0].stats.starttime + 3, endtime = st37a[0].stats.endtime)
#                
#            if st_ID == 38:
#                tr_day = st38a[0].slice(starttime = st38a[0].stats.starttime + 3, endtime = st38a[0].stats.endtime)
#                
#            if st_ID == 39:
#                tr_day = st39a[0].slice(starttime = st39a[0].stats.starttime + 3, endtime = st39a[0].stats.endtime)
#                
#            if st_ID == 40:
#                tr_day = st40a[0].slice(starttime = st40a[0].stats.starttime + 3, endtime = st40a[0].stats.endtime)
#                
#            if st_ID == 41:
#                tr_day = st41a[0].slice(starttime = st41a[0].stats.starttime + 3, endtime = st41a[0].stats.endtime)

    
            t1 = Day_start  
            t2 = Day_start + (4*24*60*60) 
            
            tr_day.detrend(type='linear')
            tr_day.detrend(type='demean')
            
            if st_ID < 34:
                tr_day.filter(type='bandpass',freqmin=0.5, freqmax=6) # 0.5 to 6 
#            else:
#                tr_day.filter(type='bandpass',freqmin=0.2, freqmax=20) ### find seis freq range
            
#            print(st_ID, tr_day)
            
            sr = tr_day.stats.sampling_rate            
            
            med_trace = np.median(abs(tr_day.data))
            
            tr = tr_day
            tr_a = abs(tr.data)     
            
            pts= int(sr * 60 * 3)
            
            t_step = 10
            step = int(t_step * sr)
            
            lenst = len(tr)-pts
            num_sum = int(lenst/step)
            
            smooth = np.zeros(shape=num_sum)
            
            for u in range(0,len(smooth)):
                q = u*step
                smooth[u] = sum(abs(tr[q:q+pts]))
            
            smooth = smooth - np.percentile(smooth,10)
            smooth=smooth/pts
            
            
            if st_ID < 34:
                on_trig = np.percentile(smooth,85) #   85
                off_trig = np.percentile(smooth,45) #  45
#            else:
#                on_trig = np.percentile(smooth,35)
#                off_trig = np.percentile(smooth,10)
                            
            d_len_trig = 600
            d_trig= d_len_trig/t_step
            
            switch =0
            on_off = np.zeros(shape=(0,5))
            num=0
            
            for m in range(0,len(smooth)):
                
                if switch ==0:
                    if smooth[m] > on_trig:
                        trig_on = m
                        switch = 1
                
                if switch ==1:
                    if smooth[m] < off_trig:
                        switch = 0
                        trig_off = m
                        d_len = (trig_off - trig_on ) 
                        
                        if d_len > d_trig:
                                                    
                            tr_cut = tr_day.slice(starttime=tr_day.stats.starttime+trig_on*t_step, endtime=tr_day.stats.starttime+trig_off*t_step)
                           
                            tr_prev = tr_day.slice(starttime=tr_day.stats.starttime+trig_on*t_step - 10*60, endtime=tr_day.stats.starttime+trig_on*t_step)

                            med_abs_amp = np.median(abs(tr_cut.data))
                            
                            
#                            print("amp compare", med_abs_amp, 1.5*med_trace)
                            
                            if med_abs_amp > 1.5*med_trace and len(tr_prev) > (sr*8*60 ):
                                
                           
                                tr_data = tr_cut.data
                                m=np.mean(tr_data)
                                tr_data = tr_data-m
                                famp = abs(np.fft.fft(tr_data))
                                start = tr_cut.stats.starttime.timestamp
                                end= tr_cut.stats.endtime.timestamp
                                
                                window=end-start
                                
                                ##  range ratio
                                fps=int(len(famp)/(sr/0.5))
                                fpe=int(len(famp)/(sr/1.5))
                                
                                sps=int(len(famp)/(sr/3))
                                spe=int(len(famp)/(sr/4))
                                
                                f_peak_m = np.mean(famp[fps:fpe])/(window*sr)
                                s_peak_m = np.mean(famp[sps:spe])/(window*sr)
                                peak_r = f_peak_m/s_peak_m
                                
                                if peak_r < 2.5: #or st_ID > 33: #currently avoiding this check for seismics
                                    
                                    on_off = np.lib.pad(on_off, ((0,1),(0,0)), 'constant', constant_values=(0))
                                       
                                    on_off[num][0]=trig_on
                                    on_off[num][1]=trig_off
                                    on_off[num][2]=d_len*t_step
                                    
                                    on_off[num][3]=UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
                                    on_off[num][4]=UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
                                    
                                    num+=1
                                    
#                                    if st_ID > 33:
#                                        tr_cut.plot(color='g')
#                                    print(len(tr_prev))
                                    
                                    
                                    if st_ID == 1:
                                        st1d = np.lib.pad(st1d, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        st1d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
                                        st1d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
                                        st1d[-1][2] = d_len*t_step
                                        st1d[-1][3] = len(station_activitya)
                                        
                                        
                                    if st_ID == 2:
                                        st2d = np.lib.pad(st2d, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        st2d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
                                        st2d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
                                        st2d[-1][2] = d_len*t_step
                                        st2d[-1][3] = len(station_activitya)
                                        
                                    if st_ID == 3:
                                        st3d = np.lib.pad(st3d, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        st3d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
                                        st3d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
                                        st3d[-1][2] = d_len*t_step
                                        st3d[-1][3] = len(station_activitya)
                                        
                                    if st_ID == 4:
                                        st4d = np.lib.pad(st4d, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        st4d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
                                        st4d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
                                        st4d[-1][2] = d_len*t_step
                                        st4d[-1][3] = len(station_activitya)
                                        
                                    if st_ID == 5:
                                        st5d = np.lib.pad(st5d, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        st5d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
                                        st5d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
                                        st5d[-1][2] = d_len*t_step
                                        st5d[-1][3] = len(station_activitya)
                                        
                                    if st_ID == 6:
                                        st6d = np.lib.pad(st6d, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        st6d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
                                        st6d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
                                        st6d[-1][2] = d_len*t_step
                                        st6d[-1][3] = len(station_activitya)
                                        
                                    if st_ID == 7:
                                        st7d = np.lib.pad(st7d, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        st7d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
                                        st7d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
                                        st7d[-1][2] = d_len*t_step
                                        st7d[-1][3] = len(station_activitya)
                                        
                                    if st_ID == 8:
                                        st8d = np.lib.pad(st8d, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        st8d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
                                        st8d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
                                        st8d[-1][2] = d_len*t_step
                                        st8d[-1][3] = len(station_activitya)
                                        
                                    if st_ID == 9:
                                        st9d = np.lib.pad(st9d, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        st9d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
                                        st9d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
                                        st9d[-1][2] = d_len*t_step
                                        st9d[-1][3] = len(station_activitya)
                                    
                                    if st_ID == 10:
                                        st10d = np.lib.pad(st10d, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        st10d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
                                        st10d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
                                        st10d[-1][2] = d_len*t_step
                                        st10d[-1][3] = len(station_activitya)
                                        
                                    if st_ID == 11:
                                        st11d = np.lib.pad(st11d, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        st11d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
                                        st11d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
                                        st11d[-1][2] = d_len*t_step
                                        st11d[-1][3] = len(station_activitya)
                                        
                                    if st_ID == 12:
                                        st12d = np.lib.pad(st12d, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        st12d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
                                        st12d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
                                        st12d[-1][2] = d_len*t_step
                                        st12d[-1][3] = len(station_activitya)
                                        
                                    if st_ID == 13:
                                        st13d = np.lib.pad(st13d, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        st13d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
                                        st13d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
                                        st13d[-1][2] = d_len*t_step
                                        st13d[-1][3] = len(station_activitya)
                                        
                                    if st_ID == 14:
                                        st14d = np.lib.pad(st14d, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        st14d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
                                        st14d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
                                        st14d[-1][2] = d_len*t_step
                                        st14d[-1][3] = len(station_activitya)
                                        
                                    if st_ID == 15:
                                        st15d = np.lib.pad(st15d, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        st15d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
                                        st15d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
                                        st15d[-1][2] = d_len*t_step
                                        st15d[-1][3] = len(station_activitya)
                                        
                                    if st_ID == 16:
                                        st16d = np.lib.pad(st16d, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        st16d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
                                        st16d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
                                        st16d[-1][2] = d_len*t_step
                                        st16d[-1][3] = len(station_activitya)
                                        
                                    if st_ID == 17:
                                        st17d = np.lib.pad(st17d, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        st17d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
                                        st17d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
                                        st17d[-1][2] = d_len*t_step
                                        st17d[-1][3] = len(station_activitya)
                                        
                                    if st_ID == 18:
                                        st18d = np.lib.pad(st18d, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        st18d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
                                        st18d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
                                        st18d[-1][2] = d_len*t_step
                                        st18d[-1][3] = len(station_activitya)
                                        
                                    if st_ID == 19:
                                        st19d = np.lib.pad(st19d, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        st19d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
                                        st19d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
                                        st19d[-1][2] = d_len*t_step
                                        st19d[-1][3] = len(station_activitya)
                                        
                                    if st_ID == 20:
                                        st20d = np.lib.pad(st20d, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        st20d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
                                        st20d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
                                        st20d[-1][2] = d_len*t_step
                                        st20d[-1][3] = len(station_activitya)
                                        
                                    if st_ID == 21:
                                        st21d = np.lib.pad(st21d, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        st21d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
                                        st21d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
                                        st21d[-1][2] = d_len*t_step
                                        st21d[-1][3] = len(station_activitya)
                                        
                                    if st_ID == 22:
                                        st22d = np.lib.pad(st22d, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        st22d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
                                        st22d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
                                        st22d[-1][2] = d_len*t_step
                                        st22d[-1][3] = len(station_activitya)
                                        
                                    if st_ID == 23:
                                        st23d = np.lib.pad(st23d, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        st23d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
                                        st23d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
                                        st23d[-1][2] = d_len*t_step
                                        st23d[-1][3] = len(station_activitya)
                                        
                                    if st_ID == 24:
                                        st24d = np.lib.pad(st24d, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        st24d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
                                        st24d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
                                        st24d[-1][2] = d_len*t_step
                                        st24d[-1][3] = len(station_activitya)
                                        
                                    if st_ID == 25:
                                        st25d = np.lib.pad(st25d, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        st25d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
                                        st25d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
                                        st25d[-1][2] = d_len*t_step
                                        st25d[-1][3] = len(station_activitya)
                                        
                                    if st_ID == 26:
                                        st26d = np.lib.pad(st26d, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        st26d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
                                        st26d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
                                        st26d[-1][2] = d_len*t_step
                                        st26d[-1][3] = len(station_activitya)
                                        
                                    if st_ID == 27:
                                        st27d = np.lib.pad(st27d, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        st27d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
                                        st27d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
                                        st27d[-1][2] = d_len*t_step
                                        st27d[-1][3] = len(station_activitya)
                                        
                                    if st_ID == 28:
                                        st28d = np.lib.pad(st28d, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        st28d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
                                        st28d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
                                        st28d[-1][2] = d_len*t_step
                                        st28d[-1][3] = len(station_activitya)
                                        
                                    if st_ID == 29:
                                        st29d = np.lib.pad(st29d, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        st29d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
                                        st29d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
                                        st29d[-1][2] = d_len*t_step
                                        st29d[-1][3] = len(station_activitya)
                                        
                                    if st_ID == 30:
                                        st30d = np.lib.pad(st30d, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        st30d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
                                        st30d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
                                        st30d[-1][2] = d_len*t_step
                                        st30d[-1][3] = len(station_activitya)
                                        
                                    if st_ID == 31:
                                        st31d = np.lib.pad(st31d, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        st31d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
                                        st31d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
                                        st31d[-1][2] = d_len*t_step
                                        st31d[-1][3] = len(station_activitya)
                                        
                                    if st_ID == 32:
                                        st32d = np.lib.pad(st32d, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        st32d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
                                        st32d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
                                        st32d[-1][2] = d_len*t_step
                                        st32d[-1][3] = len(station_activitya)
                                        
                                    if st_ID == 33:
                                        st33d = np.lib.pad(st33d, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        st33d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
                                        st33d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
                                        st33d[-1][2] = d_len*t_step
                                        st33d[-1][3] = len(station_activitya)
                                        
#                                        
#                                    if st_ID == 34:
#                                        st34d = np.lib.pad(st34d, ((0,1),(0,0)), 'constant', constant_values=(0))
#                                        st34d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
#                                        st34d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
#                                        st34d[-1][2] = d_len*t_step
#                                        st34d[-1][3] = len(seismic_activitya)
#                                        
#                                    if st_ID == 35:
#                                        st35d = np.lib.pad(st35d, ((0,1),(0,0)), 'constant', constant_values=(0))
#                                        st35d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
#                                        st35d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
#                                        st35d[-1][2] = d_len*t_step
#                                        st35d[-1][3] = len(seismic_activitya)
#                                        
#                                    if st_ID == 36:
#                                        st36d = np.lib.pad(st36d, ((0,1),(0,0)), 'constant', constant_values=(0))
#                                        st36d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
#                                        st36d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
#                                        st36d[-1][2] = d_len*t_step
#                                        st36d[-1][3] = len(seismic_activitya)
#                                        
#                                    if st_ID == 37:
#                                        st37d = np.lib.pad(st37d, ((0,1),(0,0)), 'constant', constant_values=(0))
#                                        st37d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
#                                        st37d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
#                                        st37d[-1][2] = d_len*t_step
#                                        st37d[-1][3] = len(seismic_activitya)
#                                        
#                                    if st_ID == 38:
#                                        st38d = np.lib.pad(st38d, ((0,1),(0,0)), 'constant', constant_values=(0))
#                                        st38d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
#                                        st38d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
#                                        st38d[-1][2] = d_len*t_step
#                                        st38d[-1][3] = len(seismic_activitya)
#                                        
#                                    if st_ID == 39:
#                                        st39d = np.lib.pad(st39d, ((0,1),(0,0)), 'constant', constant_values=(0))
#                                        st39d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
#                                        st39d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
#                                        st39d[-1][2] = d_len*t_step
#                                        st39d[-1][3] = len(seismic_activitya)
#                                    
#                                    if st_ID == 40:
#                                        st40d = np.lib.pad(st40d, ((0,1),(0,0)), 'constant', constant_values=(0))
#                                        st40d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
#                                        st40d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
#                                        st40d[-1][2] = d_len*t_step
#                                        st40d[-1][3] = len(seismic_activitya)
#                                        
#                                    if st_ID == 41:
#                                        st41d = np.lib.pad(st41d, ((0,1),(0,0)), 'constant', constant_values=(0))
#                                        st41d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
#                                        st41d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
#                                        st41d[-1][2] = d_len*t_step
#                                        st41d[-1][3] = len(seismic_activitya)
#                                        
#                                    

#                                
#                            print(UTCDateTime(t1 + trig_on*t_step),'to',UTCDateTime(t1 + trig_off*t_step), 'length =', d_len*t_step,'s')
#                            print('dom, cf, bwid50 = ',dom,cf, bwid50)
#                            print('mean, med, var = ',mean_abs_amp,med_abs_amp,var_abs_amp)
#            if len(smooth) > 0 and len(on_off) > 0 and st_ID > 33:
#                plt.figure()
#                plt.plot(tr,'b')
#                plt.title(st_ID)
#                plt.figure()
#                plt.plot(smooth,'k')
#                plt.plot([on_off[:,0],on_off[:,0]],[0,max(smooth)],'r')
#                plt.plot([on_off[:,1],on_off[:,1]],[0,max(smooth)],'b')

                
#            
#            print("")
#            print(len(on_off))
#            print("")
#            for e in range(0,len(on_off)):
#                
#                print(UTCDateTime(on_off[e][3]),"{:.2f}".format(on_off[e][2]/60))
                                 
                       
            

    
#%% reorder and remove duplicated events from individual station arrays
      
#reorder          
st1d = st1d[st1d[:,0].argsort()] 
st2d = st2d[st2d[:,0].argsort()]
st3d = st3d[st3d[:,0].argsort()]
st4d = st4d[st4d[:,0].argsort()]
st5d = st5d[st5d[:,0].argsort()]
st6d = st6d[st6d[:,0].argsort()]
st7d = st7d[st7d[:,0].argsort()]
st8d = st8d[st8d[:,0].argsort()]
st9d = st9d[st9d[:,0].argsort()]
st10d = st10d[st10d[:,0].argsort()] 
st11d = st11d[st11d[:,0].argsort()]
st12d = st12d[st12d[:,0].argsort()]
st13d = st13d[st13d[:,0].argsort()]
st14d = st14d[st14d[:,0].argsort()]
st15d = st15d[st15d[:,0].argsort()]
st16d = st16d[st16d[:,0].argsort()]
st17d = st17d[st17d[:,0].argsort()]
st18d = st18d[st18d[:,0].argsort()]
st19d = st19d[st19d[:,0].argsort()] 
st20d = st20d[st20d[:,0].argsort()]
st21d = st21d[st21d[:,0].argsort()]
st22d = st22d[st22d[:,0].argsort()]
st23d = st23d[st23d[:,0].argsort()]
st24d = st24d[st24d[:,0].argsort()]
st25d = st25d[st25d[:,0].argsort()]
st26d = st26d[st26d[:,0].argsort()]
st27d = st27d[st27d[:,0].argsort()]
st28d = st28d[st28d[:,0].argsort()]
st29d = st29d[st29d[:,0].argsort()]
st30d = st30d[st30d[:,0].argsort()]
st31d = st31d[st31d[:,0].argsort()]
st32d = st32d[st32d[:,0].argsort()]
st33d = st33d[st33d[:,0].argsort()]

#st34d = st34d[st34d[:,0].argsort()]
#st35d = st35d[st35d[:,0].argsort()]
#st36d = st36d[st36d[:,0].argsort()]
#st37d = st37d[st37d[:,0].argsort()]
#st38d = st38d[st38d[:,0].argsort()]
#st39d = st39d[st39d[:,0].argsort()]
#st40d = st40d[st40d[:,0].argsort()] 
#st41d = st41d[st41d[:,0].argsort()]

#remove duplicates or zero values
for x in range(len(st1d)-1,0,-1):
    if st1d[x,0]-st1d[x-1,0] < 300:
        if st1d[x,2] > st1d[x-1,2]:
            st1d = np.delete(st1d, x-1, 0)
        else:
            st1d = np.delete(st1d, x, 0)
for x in range(len(st1d)-1,0,-1):
    if st1d[x,0]== 0 or st1d[x,1]== 0 :
        st1d = np.delete(st1d, x, 0)    
            
for x in range(len(st2d)-1,0,-1):
    if st2d[x,0]-st2d[x-1,0] < 300:
        if st2d[x,2] > st2d[x-1,2]:
            st2d = np.delete(st2d, x-1, 0)
        else:
            st2d = np.delete(st2d, x, 0)
for x in range(len(st2d)-1,0,-1):            
    if st2d[x,0]== 0 or st2d[x,1]== 0 :
        st2d = np.delete(st2d, x, 0)    
            
for x in range(len(st3d)-1,0,-1):
    if st3d[x,0]-st3d[x-1,0] < 300:
        if st3d[x,2] > st3d[x-1,2]:
            st3d = np.delete(st3d, x-1, 0)
        else:
            st3d = np.delete(st3d, x, 0)
for x in range(len(st3d)-1,0,-1):
    if st3d[x,0]== 0 or st3d[x,1]== 0 :
        st3d = np.delete(st3d, x, 0)    
                        
for x in range(len(st4d)-1,0,-1):
    if st4d[x,0]-st4d[x-1,0] < 300:
        if st4d[x,2] > st4d[x-1,2]:
            st4d = np.delete(st4d, x-1, 0)
        else:
            st4d = np.delete(st4d, x, 0)
for x in range(len(st4d)-1,0,-1):
    if st4d[x,0]== 0 or st4d[x,1]== 0 :
        st4d = np.delete(st4d, x, 0)    
                        
for x in range(len(st5d)-1,0,-1):
    if st5d[x,0]-st5d[x-1,0] < 300:
        if st5d[x,2] > st5d[x-1,2]:
            st5d = np.delete(st5d, x-1, 0)
        else:
            st5d = np.delete(st5d, x, 0)
for x in range(len(st5d)-1,0,-1):
    if st5d[x,0]== 0 or st5d[x,1]== 0 :
        st5d = np.delete(st5d, x, 0)    
                        
for x in range(len(st6d)-1,0,-1):
    if st6d[x,0]-st6d[x-1,0] < 300:
        if st6d[x,2] > st6d[x-1,2]:
            st6d = np.delete(st6d, x-1, 0)
        else:
            st6d = np.delete(st6d, x, 0)
for x in range(len(st6d)-1,0,-1):
    if st6d[x,0]== 0 or st6d[x,1]== 0 :
        st6d = np.delete(st6d, x, 0)    
            
for x in range(len(st7d)-1,0,-1):
    if st7d[x,0]-st7d[x-1,0] < 300:
        if st7d[x,2] > st7d[x-1,2]:
            st7d = np.delete(st7d, x-1, 0)
        else:
            st7d = np.delete(st7d, x, 0)
for x in range(len(st7d)-1,0,-1):
    if st7d[x,0]== 0 or st7d[x,1]== 0 :
        st7d = np.delete(st7d, x, 0)    
            
for x in range(len(st8d)-1,0,-1):
    if st8d[x,0]-st8d[x-1,0] < 300:
        if st8d[x,2] > st8d[x-1,2]:
            st8d = np.delete(st8d, x-1, 0)
        else:
            st8d = np.delete(st8d, x, 0)
for x in range(len(st8d)-1,0,-1):
    if st8d[x,0]== 0 or st8d[x,1]== 0 :
        st8d = np.delete(st8d, x, 0)    
            
for x in range(len(st9d)-1,0,-1):
    if st9d[x,0]-st9d[x-1,0] < 300:
        if st9d[x,2] > st9d[x-1,2]:
            st9d = np.delete(st9d, x-1, 0)
        else:
            st9d = np.delete(st9d, x, 0)
for x in range(len(st9d)-1,0,-1):
    if st9d[x,0]== 0 or st9d[x,1]== 0 :
        st9d = np.delete(st9d, x, 0)    
            
for x in range(len(st10d)-1,0,-1):
    if st10d[x,0]-st10d[x-1,0] < 300:
        if st10d[x,2] > st10d[x-1,2]:
            st10d = np.delete(st10d, x-1, 0)
        else:
            st10d = np.delete(st10d, x, 0)
for x in range(len(st10d)-1,0,-1):
    if st10d[x,0]== 0 or st10d[x,1]== 0 :
        st10d = np.delete(st10d, x, 0)    
            
for x in range(len(st11d)-1,0,-1):
    if st11d[x,0]-st11d[x-1,0] < 300:
        if st11d[x,2] > st11d[x-1,2]:
            st11d = np.delete(st11d, x-1, 0)
        else:
            st11d = np.delete(st11d, x, 0)
for x in range(len(st11d)-1,0,-1):
    if st11d[x,0]== 0 or st11d[x,1]== 0 :
        st11d = np.delete(st11d, x, 0)    
            
for x in range(len(st12d)-1,0,-1):
    if st12d[x,0]-st12d[x-1,0] < 300:
        if st12d[x,2] > st12d[x-1,2]:
            st12d = np.delete(st12d, x-1, 0)
        else:
            st12d = np.delete(st12d, x, 0)
for x in range(len(st12d)-1,0,-1):
    if st12d[x,0]== 0 or st12d[x,1]== 0 :
        st12d = np.delete(st12d, x, 0)    
            
for x in range(len(st13d)-1,0,-1):
    if st13d[x,0]-st13d[x-1,0] < 300:
        if st13d[x,2] > st13d[x-1,2]:
            st13d = np.delete(st13d, x-1, 0)
        else:
            st13d = np.delete(st13d, x, 0)
for x in range(len(st13d)-1,0,-1):
    if st13d[x,0]== 0 or st13d[x,1]== 0 :
        st13d = np.delete(st13d, x, 0)    
            
for x in range(len(st14d)-1,0,-1):
    if st14d[x,0]-st14d[x-1,0] < 300:
        if st14d[x,2] > st14d[x-1,2]:
            st14d = np.delete(st14d, x-1, 0)
        else:
            st14d = np.delete(st14d, x, 0)
for x in range(len(st14d)-1,0,-1):
    if st14d[x,0]== 0 or st14d[x,1]== 0 :
        st14d = np.delete(st14d, x, 0)    
            
for x in range(len(st15d)-1,0,-1):
    if st15d[x,0]-st15d[x-1,0] < 300:
        if st15d[x,2] > st15d[x-1,2]:
            st15d = np.delete(st15d, x-1, 0)
        else:
            st15d = np.delete(st15d, x, 0)
for x in range(len(st15d)-1,0,-1):
    if st15d[x,0]== 0 or st15d[x,1]== 0 :
        st15d = np.delete(st15d, x, 0)    
            
for x in range(len(st16d)-1,0,-1):
    if st16d[x,0]-st16d[x-1,0] < 300:
        if st16d[x,2] > st16d[x-1,2]:
            st16d = np.delete(st16d, x-1, 0)
        else:
            st16d = np.delete(st16d, x, 0)
for x in range(len(st16d)-1,0,-1):
    if st16d[x,0]== 0 or st16d[x,1]== 0 :
        st16d = np.delete(st16d, x, 0)    
            
for x in range(len(st17d)-1,0,-1):
    if st17d[x,0]-st17d[x-1,0] < 300:
        if st17d[x,2] > st17d[x-1,2]:
            st17d = np.delete(st17d, x-1, 0)
        else:
            st17d = np.delete(st17d, x, 0)
for x in range(len(st17d)-1,0,-1):
    if st17d[x,0]== 0 or st17d[x,1]== 0 :
        st17d = np.delete(st17d, x, 0)    
            
for x in range(len(st18d)-1,0,-1):
    if st18d[x,0]-st18d[x-1,0] < 300:
        if st18d[x,2] > st18d[x-1,2]:
            st18d = np.delete(st18d, x-1, 0)
        else:
            st18d = np.delete(st18d, x, 0)
for x in range(len(st18d)-1,0,-1):
    if st18d[x,0]== 0 or st18d[x,1]== 0 :
        st18d = np.delete(st18d, x, 0)    
            
for x in range(len(st19d)-1,0,-1):
    if st19d[x,0]-st19d[x-1,0] < 300:
        if st19d[x,2] > st19d[x-1,2]:
            st19d = np.delete(st19d, x-1, 0)
        else:
            st19d = np.delete(st19d, x, 0)
for x in range(len(st19d)-1,0,-1):
    if st19d[x,0]== 0 or st19d[x,1]== 0 :
        st19d = np.delete(st19d, x, 0)    
            
for x in range(len(st20d)-1,0,-1):
    if st20d[x,0]-st20d[x-1,0] < 300:
        if st20d[x,2] > st20d[x-1,2]:
            st20d = np.delete(st20d, x-1, 0)
        else:
            st20d = np.delete(st20d, x, 0)
for x in range(len(st20d)-1,0,-1):
    if st20d[x,0]== 0 or st20d[x,1]== 0 :
        st20d = np.delete(st20d, x, 0)    
            
for x in range(len(st21d)-1,0,-1):
    if st21d[x,0]-st21d[x-1,0] < 300:
        if st21d[x,2] > st21d[x-1,2]:
            st21d = np.delete(st21d, x-1, 0)
        else:
            st21d = np.delete(st21d, x, 0)
for x in range(len(st21d)-1,0,-1):
    if st21d[x,0]== 0 or st21d[x,1]== 0 :
        st21d = np.delete(st21d, x, 0)    
            
for x in range(len(st22d)-1,0,-1):
    if st22d[x,0]-st22d[x-1,0] < 300:
        if st22d[x,2] > st22d[x-1,2]:
            st22d = np.delete(st22d, x-1, 0)
        else:
            st22d = np.delete(st22d, x, 0)
for x in range(len(st22d)-1,0,-1):
    if st22d[x,0]== 0 or st22d[x,1]== 0 :
        st22d = np.delete(st22d, x, 0)    
            
for x in range(len(st23d)-1,0,-1):
    if st23d[x,0]-st23d[x-1,0] < 300:
        if st23d[x,2] > st23d[x-1,2]:
            st23d = np.delete(st23d, x-1, 0)
        else:
            st23d = np.delete(st23d, x, 0)
for x in range(len(st23d)-1,0,-1):
    if st23d[x,0]== 0 or st23d[x,1]== 0 :
        st23d = np.delete(st23d, x, 0)    
            
for x in range(len(st24d)-1,0,-1):
    if st24d[x,0]-st24d[x-1,0] < 300:
        if st24d[x,2] > st24d[x-1,2]:
            st24d = np.delete(st24d, x-1, 0)
        else:
            st24d = np.delete(st24d, x, 0)
for x in range(len(st24d)-1,0,-1):
    if st24d[x,0]== 0 or st24d[x,1]== 0 :
        st24d = np.delete(st24d, x, 0)    
            
for x in range(len(st25d)-1,0,-1):
    if st25d[x,0]-st25d[x-1,0] < 300:
        if st25d[x,2] > st25d[x-1,2]:
            st25d = np.delete(st25d, x-1, 0)
        else:
            st25d = np.delete(st25d, x, 0)
for x in range(len(st25d)-1,0,-1):
    if st25d[x,0]== 0 or st25d[x,1]== 0 :
        st25d = np.delete(st25d, x, 0)    
            
for x in range(len(st26d)-1,0,-1):
    if st26d[x,0]-st26d[x-1,0] < 300:
        if st26d[x,2] > st26d[x-1,2]:
            st26d = np.delete(st26d, x-1, 0)
        else:
            st26d = np.delete(st26d, x, 0)
for x in range(len(st26d)-1,0,-1):
    if st26d[x,0]== 0 or st26d[x,1]== 0 :
        st26d = np.delete(st26d, x, 0)    
            
for x in range(len(st27d)-1,0,-1):
    if st27d[x,0]-st27d[x-1,0] < 300:
        if st27d[x,2] > st27d[x-1,2]:
            st27d = np.delete(st27d, x-1, 0)
        else:
            st27d = np.delete(st27d, x, 0)
for x in range(len(st27d)-1,0,-1):
    if st27d[x,0]== 0 or st27d[x,1]== 0 :
        st27d = np.delete(st27d, x, 0)    
            
for x in range(len(st28d)-1,0,-1):
    if st28d[x,0]-st28d[x-1,0] < 300:
        if st28d[x,2] > st28d[x-1,2]:
            st28d = np.delete(st28d, x-1, 0)
        else:
            st28d = np.delete(st28d, x, 0)
for x in range(len(st28d)-1,0,-1):
    if st28d[x,0]== 0 or st28d[x,1]== 0 :
        st28d = np.delete(st28d, x, 0)    
            
for x in range(len(st29d)-1,0,-1):
    if st29d[x,0]-st29d[x-1,0] < 300:
        if st29d[x,2] > st29d[x-1,2]:
            st29d = np.delete(st29d, x-1, 0)
        else:
            st29d = np.delete(st29d, x, 0)
for x in range(len(st29d)-1,0,-1):
    if st29d[x,0]== 0 or st29d[x,1]== 0 :
        st29d = np.delete(st29d, x, 0)    
            
for x in range(len(st30d)-1,0,-1):
    if st30d[x,0]-st30d[x-1,0] < 300:
        if st30d[x,2] > st30d[x-1,2]:
            st30d = np.delete(st30d, x-1, 0)
        else:
            st30d = np.delete(st30d, x, 0)
for x in range(len(st30d)-1,0,-1):
    if st30d[x,0]== 0 or st30d[x,1]== 0 :
        st30d = np.delete(st30d, x, 0)    
            
for x in range(len(st31d)-1,0,-1):
    if st31d[x,0]-st31d[x-1,0] < 300:
        if st31d[x,2] > st31d[x-1,2]:
            st31d = np.delete(st31d, x-1, 0)
        else:
            st31d = np.delete(st31d, x, 0)
for x in range(len(st31d)-1,0,-1):
    if st31d[x,0]== 0 or st31d[x,1]== 0 :
        st31d = np.delete(st31d, x, 0)    
            
for x in range(len(st32d)-1,0,-1):
    if st32d[x,0]-st32d[x-1,0] < 300:
        if st32d[x,2] > st32d[x-1,2]:
            st32d = np.delete(st32d, x-1, 0)
        else:
            st32d = np.delete(st32d, x, 0)
for x in range(len(st32d)-1,0,-1):
    if st32d[x,0]== 0 or st32d[x,1]== 0 :
        st32d = np.delete(st32d, x, 0)    
            
for x in range(len(st33d)-1,0,-1):
    if st33d[x,0]-st33d[x-1,0] < 300:
        if st33d[x,2] > st33d[x-1,2]:
            st33d = np.delete(st33d, x-1, 0)
        else:
            st33d = np.delete(st33d, x, 0)
for x in range(len(st33d)-1,0,-1):
    if st33d[x,0]== 0 or st33d[x,1]== 0 :
        st33d = np.delete(st33d, x, 0)    
            

#for x in range(len(st34d)-1,0,-1):
#    if st34d[x,0]-st34d[x-1,0] < 300:
#        if st34d[x,2] > st34d[x-1,2]:
#            st34d = np.delete(st34d, x-1, 0)
#        else:
#            st34d = np.delete(st34d, x, 0)
#for x in range(len(st34d)-1,0,-1):
#    if st34d[x,0]== 0 or st34d[x,1]== 0 :
#        st34d = np.delete(st34d, x, 0)    
#                        
#for x in range(len(st35d)-1,0,-1):
#    if st35d[x,0]-st35d[x-1,0] < 300:
#        if st35d[x,2] > st35d[x-1,2]:
#            st35d = np.delete(st35d, x-1, 0)
#        else:
#            st35d = np.delete(st35d, x, 0)
#for x in range(len(st35d)-1,0,-1):
#    if st35d[x,0]== 0 or st35d[x,1]== 0 :
#        st35d = np.delete(st35d, x, 0)    
#                        
#for x in range(len(st36d)-1,0,-1):
#    if st36d[x,0]-st36d[x-1,0] < 300:
#        if st36d[x,2] > st36d[x-1,2]:
#            st36d = np.delete(st36d, x-1, 0)
#        else:
#            st36d = np.delete(st36d, x, 0)
#for x in range(len(st36d)-1,0,-1):
#    if st36d[x,0]== 0 or st36d[x,1]== 0 :
#        st36d = np.delete(st36d, x, 0)    
#            
#for x in range(len(st37d)-1,0,-1):
#    if st37d[x,0]-st37d[x-1,0] < 300:
#        if st37d[x,2] > st37d[x-1,2]:
#            st37d = np.delete(st37d, x-1, 0)
#        else:
#            st37d = np.delete(st37d, x, 0)
#for x in range(len(st37d)-1,0,-1):
#    if st37d[x,0]== 0 or st37d[x,1]== 0 :
#        st37d = np.delete(st37d, x, 0)    
#            
#for x in range(len(st38d)-1,0,-1):
#    if st38d[x,0]-st38d[x-1,0] < 300:
#        if st38d[x,2] > st38d[x-1,2]:
#            st38d = np.delete(st38d, x-1, 0)
#        else:
#            st38d = np.delete(st38d, x, 0)
#for x in range(len(st38d)-1,0,-1):
#    if st38d[x,0]== 0 or st38d[x,1]== 0 :
#        st38d = np.delete(st38d, x, 0)    
#            
#for x in range(len(st39d)-1,0,-1):
#    if st39d[x,0]-st39d[x-1,0] < 300:
#        if st39d[x,2] > st39d[x-1,2]:
#            st39d = np.delete(st39d, x-1, 0)
#        else:
#            st39d = np.delete(st39d, x, 0)
#for x in range(len(st39d)-1,0,-1):
#    if st39d[x,0]== 0 or st39d[x,1]== 0 :
#        st39d = np.delete(st39d, x, 0)    
#            
#for x in range(len(st40d)-1,0,-1):
#    if st40d[x,0]-st40d[x-1,0] < 300:
#        if st40d[x,2] > st40d[x-1,2]:
#            st40d = np.delete(st40d, x-1, 0)
#        else:
#            st40d = np.delete(st40d, x, 0)
#for x in range(len(st40d)-1,0,-1):
#    if st40d[x,0]== 0 or st40d[x,1]== 0 :
#        st40d = np.delete(st40d, x, 0)    
#            
#for x in range(len(st41d)-1,0,-1):
#    if st41d[x,0]-st41d[x-1,0] < 300:
#        if st41d[x,2] > st41d[x-1,2]:
#            st41d = np.delete(st41d, x-1, 0)
#        else:
#            st41d = np.delete(st41d, x, 0)
#for x in range(len(st41d)-1,0,-1):
#    if st41d[x,0]== 0 or st41d[x,1]== 0 :
#        st41d = np.delete(st41d, x, 0)               

             
#%% make merged array and populate - Infrasound

#index 0 = number of stations, then all station ticks, index 34 = earliest start, 35 = latest start, 36 = earliest end, 37 = latest end
merged_detec = np.zeros(shape=(len(st1d),39))

merged_detec[:,0] = 1              
merged_detec[:,1] = 1            
merged_detec[:,34] = st1d[:,0]    
merged_detec[:,35] = st1d[:,0]  
merged_detec[:,36] = st1d[:,1] 
merged_detec[:,37] = st1d[:,1]    
merged_detec[:,38] = st1d[:,3]   
                

#add in other stations

#2
if len(merged_detec) > 0:
    for x in range(0,len(st2d)):
        added = 0
        for y in range(0,len(merged_detec)):
            if added == 0:
                if (abs(merged_detec[y,34] - st2d[x,0]) < 60*30 and abs(merged_detec[y,37] - st2d[x,1]) < 60*120) or (abs(merged_detec[y,35] - st2d[x,0]) < 60*10 and abs(merged_detec[y,37] - st2d[x,1]) < 60*120) or (abs(merged_detec[y,37] - st2d[x,1]) < 60*30 and abs(merged_detec[y,34] - st2d[x,0]) < 60*120) or (abs(merged_detec[y,36] - st2d[x,1]) < 60*20 and abs(merged_detec[y,34] - st2d[x,0]) < 60*120):
                    
                    if merged_detec[y,34] < st2d[x,0] < merged_detec[y,37] or merged_detec[y,34] < st2d[x,1] < merged_detec[y,37] or  st2d[x,0] < merged_detec[y,34] < st2d[x,1] :
                        
                        merged_detec[y,0] += 1
                        merged_detec[y,2] = 1
                        
                        if merged_detec[y,34] > st2d[x,0] :
                            merged_detec[y,34] = st2d[x,0]
                        if merged_detec[y,35] < st2d[x,0] :
                            merged_detec[y,35] = st2d[x,0]
                        if merged_detec[y,36] > st2d[x,1] :
                            merged_detec[y,36] = st2d[x,1]
                        if merged_detec[y,37] < st2d[x,1] :
                            merged_detec[y,37] = st2d[x,1]                  
                        
                        added = 1
                    
        if added == 0 :
            merged_detec = np.lib.pad(merged_detec, ((0,1),(0,0)), 'constant', constant_values=(0))
            merged_detec[-1,0] = 1 
            merged_detec[-1,2] = 1    
            merged_detec[-1,34] = st2d[x,0]    
            merged_detec[-1,35] = st2d[x,0]  
            merged_detec[-1,36] = st2d[x,1] 
            merged_detec[-1,37] = st2d[x,1] 
            merged_detec[-1,38] = st2d[x,3]  
                 
else:
    merged_detec = np.lib.pad(merged_detec, ((0,len(st2d)),(0,0)), 'constant', constant_values=(0))
    merged_detec[:,0] = 1              
    merged_detec[:,2] = 1            
    merged_detec[:,34] = st2d[:,0]    
    merged_detec[:,35] = st2d[:,0]  
    merged_detec[:,36] = st2d[:,1] 
    merged_detec[:,37] = st2d[:,1]  
    merged_detec[:,38] = st2d[:,3]  
                  
                
#3
if len(merged_detec) > 0:
    for x in range(0,len(st3d)):
        added = 0
        for y in range(0,len(merged_detec)):
            if added == 0:
                if (abs(merged_detec[y,34] - st3d[x,0]) < 60*30 and abs(merged_detec[y,37] - st3d[x,1]) < 60*120) or (abs(merged_detec[y,35] - st3d[x,0]) < 60*10 and abs(merged_detec[y,37] - st3d[x,1]) < 60*120) or (abs(merged_detec[y,37] - st3d[x,1]) < 60*30 and abs(merged_detec[y,34] - st3d[x,0]) < 60*120) or (abs(merged_detec[y,36] - st3d[x,1]) < 60*20 and abs(merged_detec[y,34] - st3d[x,0]) < 60*120):
                
                    if merged_detec[y,34] < st3d[x,0] < merged_detec[y,37] or merged_detec[y,34] < st3d[x,1] < merged_detec[y,37] or  st3d[x,0] < merged_detec[y,34] < st3d[x,1] :
                        
                        merged_detec[y,0] += 1
                        merged_detec[y,3] = 1
                        
                        if merged_detec[y,34] >  st3d[x,0] :
                            merged_detec[y,34] = st3d[x,0]
                        if merged_detec[y,35] <  st3d[x,0] :
                            merged_detec[y,35] = st3d[x,0]
                        if merged_detec[y,36] >  st3d[x,1] :
                            merged_detec[y,36] = st3d[x,1]
                        if merged_detec[y,37] <  st3d[x,1] :
                            merged_detec[y,37] = st3d[x,1]                  
                        
                        added = 1
                    
        if added == 0 :
            merged_detec = np.lib.pad(merged_detec, ((0,1),(0,0)), 'constant', constant_values=(0))
            merged_detec[-1,0] = 1 
            merged_detec[-1,3] = 1    
            merged_detec[-1,34] = st3d[x,0]    
            merged_detec[-1,35] = st3d[x,0]  
            merged_detec[-1,36] = st3d[x,1] 
            merged_detec[-1,37] = st3d[x,1]  
            merged_detec[-1,38] = st3d[x,3] 
else:
    merged_detec = np.lib.pad(merged_detec, ((0,len(st3d)),(0,0)), 'constant', constant_values=(0))
    merged_detec[:,0] = 1              
    merged_detec[:,3] = 1            
    merged_detec[:,34] = st3d[:,0]    
    merged_detec[:,35] = st3d[:,0]  
    merged_detec[:,36] = st3d[:,1] 
    merged_detec[:,37] = st3d[:,1]
    merged_detec[:,38] = st3d[:,3]  
                
#4
if len(merged_detec) > 0:
    for x in range(0,len(st4d)):
        added = 0
        for y in range(0,len(merged_detec)):
            if added == 0:
                if (abs(merged_detec[y,34] - st4d[x,0]) < 60*30 and abs(merged_detec[y,37] - st4d[x,1]) < 60*120) or (abs(merged_detec[y,35] - st4d[x,0]) < 60*10 and abs(merged_detec[y,37] - st4d[x,1]) < 60*120) or (abs(merged_detec[y,37] - st4d[x,1]) < 60*30 and abs(merged_detec[y,34] - st4d[x,0]) < 60*120) or (abs(merged_detec[y,36] - st4d[x,1]) < 60*20 and abs(merged_detec[y,34] - st4d[x,0]) < 60*120):

                    if merged_detec[y,34] < st4d[x,0] < merged_detec[y,37] or merged_detec[y,34] < st4d[x,1] < merged_detec[y,37] or  st4d[x,0] < merged_detec[y,34] < st4d[x,1] :
                        
                        merged_detec[y,0] += 1
                        merged_detec[y,4] = 1
                        
                        if merged_detec[y,34] >  st4d[x,0] :
                            merged_detec[y,34] = st4d[x,0]
                        if merged_detec[y,35] <  st4d[x,0] :
                            merged_detec[y,35] = st4d[x,0]
                        if merged_detec[y,36] >  st4d[x,1] :
                            merged_detec[y,36] = st4d[x,1]
                        if merged_detec[y,37] <  st4d[x,1] :
                            merged_detec[y,37] = st4d[x,1]                  
                        
                        added = 1
                    
        if added == 0 :
            merged_detec = np.lib.pad(merged_detec, ((0,1),(0,0)), 'constant', constant_values=(0))
            merged_detec[-1,0] = 1 
            merged_detec[-1,4] = 1    
            merged_detec[-1,34] = st4d[x,0]    
            merged_detec[-1,35] = st4d[x,0]  
            merged_detec[-1,36] = st4d[x,1] 
            merged_detec[-1,37] = st4d[x,1] 
            merged_detec[-1,38] = st4d[x,3]
else:
    merged_detec = np.lib.pad(merged_detec, ((0,len(st4d)),(0,0)), 'constant', constant_values=(0))
    merged_detec[:,0] = 1              
    merged_detec[:,4] = 1            
    merged_detec[:,34] = st4d[:,0]    
    merged_detec[:,35] = st4d[:,0]  
    merged_detec[:,36] = st4d[:,1] 
    merged_detec[:,37] = st4d[:,1]
    merged_detec[:,38] = st4d[:,3]  

#5
if len(merged_detec) > 0:
    for x in range(0,len(st5d)):
        added = 0
        for y in range(0,len(merged_detec)):
            if added == 0:
                if (abs(merged_detec[y,34] - st5d[x,0]) < 60*30 and abs(merged_detec[y,37] - st5d[x,1]) < 60*120) or (abs(merged_detec[y,35] - st5d[x,0]) < 60*10 and abs(merged_detec[y,37] - st5d[x,1]) < 60*120) or (abs(merged_detec[y,37] - st5d[x,1]) < 60*30 and abs(merged_detec[y,34] - st5d[x,0]) < 60*120) or (abs(merged_detec[y,36] - st5d[x,1]) < 60*20 and abs(merged_detec[y,34] - st5d[x,0]) < 60*120):

                    if merged_detec[y,34] < st5d[x,0] < merged_detec[y,37] or merged_detec[y,34] < st5d[x,1] < merged_detec[y,37] or  st5d[x,0] < merged_detec[y,34] < st5d[x,1] :
                        
                        merged_detec[y,0] += 1
                        merged_detec[y,5] = 1
                        
                        if merged_detec[y,34] >  st5d[x,0] :
                            merged_detec[y,34] = st5d[x,0]
                        if merged_detec[y,35] <  st5d[x,0] :
                            merged_detec[y,35] = st5d[x,0]
                        if merged_detec[y,36] >  st5d[x,1] :
                            merged_detec[y,36] = st5d[x,1]
                        if merged_detec[y,37] <  st5d[x,1] :
                            merged_detec[y,37] = st5d[x,1]                  
                        
                        added = 1
                        
        if added == 0 :
            merged_detec = np.lib.pad(merged_detec, ((0,1),(0,0)), 'constant', constant_values=(0))
            merged_detec[-1,0] = 1 
            merged_detec[-1,5] = 1    
            merged_detec[-1,34] = st5d[x,0]    
            merged_detec[-1,35] = st5d[x,0]  
            merged_detec[-1,36] = st5d[x,1] 
            merged_detec[-1,37] = st5d[x,1] 
            merged_detec[-1,38] = st5d[x,3]
else:
    merged_detec = np.lib.pad(merged_detec, ((0,len(st5d)),(0,0)), 'constant', constant_values=(0))
    merged_detec[:,0] = 1              
    merged_detec[:,5] = 1            
    merged_detec[:,34] = st5d[:,0]    
    merged_detec[:,35] = st5d[:,0]  
    merged_detec[:,36] = st5d[:,1] 
    merged_detec[:,37] = st5d[:,1] 
    merged_detec[:,38] = st5d[:,3]  

#6
if len(merged_detec) > 0:
    for x in range(0,len(st6d)):
        added = 0
        for y in range(0,len(merged_detec)):
            if added == 0:
                if (abs(merged_detec[y,34] - st6d[x,0]) < 60*30 and abs(merged_detec[y,37] - st6d[x,1]) < 60*120) or (abs(merged_detec[y,35] - st6d[x,0]) < 60*10 and abs(merged_detec[y,37] - st6d[x,1]) < 60*120) or (abs(merged_detec[y,37] - st6d[x,1]) < 60*30 and abs(merged_detec[y,34] - st6d[x,0]) < 60*120) or (abs(merged_detec[y,36] - st6d[x,1]) < 60*20 and abs(merged_detec[y,34] - st6d[x,0]) < 60*120):

                    if merged_detec[y,34] < st6d[x,0] < merged_detec[y,37] or merged_detec[y,34] < st6d[x,1] < merged_detec[y,37] or  st6d[x,0] < merged_detec[y,34] < st6d[x,1] :
                        
                        merged_detec[y,0] += 1
                        merged_detec[y,6] = 1
                        
                        if merged_detec[y,34] >  st6d[x,0] :
                            merged_detec[y,34] = st6d[x,0]
                        if merged_detec[y,35] <  st6d[x,0] :
                            merged_detec[y,35] = st6d[x,0]
                        if merged_detec[y,36] >  st6d[x,1] :
                            merged_detec[y,36] = st6d[x,1]
                        if merged_detec[y,37] <  st6d[x,1] :
                            merged_detec[y,37] = st6d[x,1]                  
                        
                        added = 1
                    
        if added == 0 :
            merged_detec = np.lib.pad(merged_detec, ((0,1),(0,0)), 'constant', constant_values=(0))
            merged_detec[-1,0] = 1 
            merged_detec[-1,6] = 1    
            merged_detec[-1,34] = st6d[x,0]    
            merged_detec[-1,35] = st6d[x,0]  
            merged_detec[-1,36] = st6d[x,1] 
            merged_detec[-1,37] = st6d[x,1] 
            merged_detec[-1,38] = st6d[x,3]                   
else:
    merged_detec = np.lib.pad(merged_detec, ((0,len(st6d)),(0,0)), 'constant', constant_values=(0))
    merged_detec[:,0] = 1              
    merged_detec[:,6] = 1            
    merged_detec[:,34] = st6d[:,0]    
    merged_detec[:,35] = st6d[:,0]  
    merged_detec[:,36] = st6d[:,1] 
    merged_detec[:,37] = st6d[:,1]  
    merged_detec[:,38] = st6d[:,3]  

#7
if len(merged_detec) > 0:
    for x in range(0,len(st7d)):
        added = 0
        for y in range(0,len(merged_detec)):
            if added == 0:
                if (abs(merged_detec[y,34] - st7d[x,0]) < 60*30 and abs(merged_detec[y,37] - st7d[x,1]) < 60*120) or (abs(merged_detec[y,35] - st7d[x,0]) < 60*10 and abs(merged_detec[y,37] - st7d[x,1]) < 60*120) or (abs(merged_detec[y,37] - st7d[x,1]) < 60*30 and abs(merged_detec[y,34] - st7d[x,0]) < 60*120) or (abs(merged_detec[y,36] - st7d[x,1]) < 60*20 and abs(merged_detec[y,34] - st7d[x,0]) < 60*120):

                    if merged_detec[y,34] < st7d[x,0] < merged_detec[y,37] or merged_detec[y,34] < st7d[x,1] < merged_detec[y,37] or  st7d[x,0] < merged_detec[y,34] < st7d[x,1] :
                        
                        merged_detec[y,0] += 1
                        merged_detec[y,7] = 1
                        
                        if merged_detec[y,34] >  st7d[x,0] :
                            merged_detec[y,34] = st7d[x,0]
                        if merged_detec[y,35] <  st7d[x,0] :
                            merged_detec[y,35] = st7d[x,0]
                        if merged_detec[y,36] >  st7d[x,1] :
                            merged_detec[y,36] = st7d[x,1]
                        if merged_detec[y,37] <  st7d[x,1] :
                            merged_detec[y,37] = st7d[x,1]                  
                        
                        added = 1
                    
        if added == 0 :
            merged_detec = np.lib.pad(merged_detec, ((0,1),(0,0)), 'constant', constant_values=(0))
            merged_detec[-1,0] = 1 
            merged_detec[-1,7] = 1    
            merged_detec[-1,34] = st7d[x,0]    
            merged_detec[-1,35] = st7d[x,0]  
            merged_detec[-1,36] = st7d[x,1] 
            merged_detec[-1,37] = st7d[x,1]  
            merged_detec[-1,38] = st7d[x,3]                  
else:
    merged_detec = np.lib.pad(merged_detec, ((0,len(st7d)),(0,0)), 'constant', constant_values=(0))
    merged_detec[:,0] = 1              
    merged_detec[:,7] = 1            
    merged_detec[:,34] = st7d[:,0]    
    merged_detec[:,35] = st7d[:,0]  
    merged_detec[:,36] = st7d[:,1] 
    merged_detec[:,37] = st7d[:,1] 
    merged_detec[:,38] = st7d[:,3]  

#8
if len(merged_detec) > 0:
    for x in range(0,len(st8d)):
        added = 0
        for y in range(0,len(merged_detec)):
            if added == 0:
                if (abs(merged_detec[y,34] - st8d[x,0]) < 60*30 and abs(merged_detec[y,37] - st8d[x,1]) < 60*120) or (abs(merged_detec[y,35] - st8d[x,0]) < 60*10 and abs(merged_detec[y,37] - st8d[x,1]) < 60*120) or (abs(merged_detec[y,37] - st8d[x,1]) < 60*30 and abs(merged_detec[y,34] - st8d[x,0]) < 60*120) or (abs(merged_detec[y,36] - st8d[x,1]) < 60*20 and abs(merged_detec[y,34] - st8d[x,0]) < 60*120):

                    if merged_detec[y,34] < st8d[x,0] < merged_detec[y,37] or merged_detec[y,34] < st8d[x,1] < merged_detec[y,37] or  st8d[x,0] < merged_detec[y,34] < st8d[x,1] :
                        
                        merged_detec[y,0] += 1
                        merged_detec[y,8] = 1
                        
                        if merged_detec[y,34] >  st8d[x,0] :
                            merged_detec[y,34] = st8d[x,0]
                        if merged_detec[y,35] <  st8d[x,0] :
                            merged_detec[y,35] = st8d[x,0]
                        if merged_detec[y,36] >  st8d[x,1] :
                            merged_detec[y,36] = st8d[x,1]
                        if merged_detec[y,37] <  st8d[x,1] :
                            merged_detec[y,37] = st8d[x,1]                  
                        
                        added = 1
                    
        if added == 0 :
            merged_detec = np.lib.pad(merged_detec, ((0,1),(0,0)), 'constant', constant_values=(0))
            merged_detec[-1,0] = 1 
            merged_detec[-1,8] = 1    
            merged_detec[-1,34] = st8d[x,0]    
            merged_detec[-1,35] = st8d[x,0]  
            merged_detec[-1,36] = st8d[x,1] 
            merged_detec[-1,37] = st8d[x,1]  
            merged_detec[-1,38] = st8d[x,3]                  
else:
    merged_detec = np.lib.pad(merged_detec, ((0,len(st8d)),(0,0)), 'constant', constant_values=(0))
    merged_detec[:,0] = 1              
    merged_detec[:,8] = 1            
    merged_detec[:,34] = st8d[:,0]    
    merged_detec[:,35] = st8d[:,0]  
    merged_detec[:,36] = st8d[:,1] 
    merged_detec[:,37] = st8d[:,1] 
    merged_detec[:,38] = st8d[:,3]  

#9
if len(merged_detec) > 0:
    for x in range(0,len(st9d)):
        added = 0
        for y in range(0,len(merged_detec)):
            if added == 0:
                if (abs(merged_detec[y,34] - st9d[x,0]) < 60*30 and abs(merged_detec[y,37] - st9d[x,1]) < 60*120) or (abs(merged_detec[y,35] - st9d[x,0]) < 60*10 and abs(merged_detec[y,37] - st9d[x,1]) < 60*120) or (abs(merged_detec[y,37] - st9d[x,1]) < 60*30 and abs(merged_detec[y,34] - st9d[x,0]) < 60*120) or (abs(merged_detec[y,36] - st9d[x,1]) < 60*20 and abs(merged_detec[y,34] - st9d[x,0]) < 60*120):

                    if merged_detec[y,34] < st9d[x,0] < merged_detec[y,37] or merged_detec[y,34] < st9d[x,1] < merged_detec[y,37] or  st9d[x,0] < merged_detec[y,34] < st9d[x,1] :
                        
                        merged_detec[y,0] += 1
                        merged_detec[y,9] = 1
                        
                        if merged_detec[y,34] >  st9d[x,0] :
                            merged_detec[y,34] = st9d[x,0]
                        if merged_detec[y,35] <  st9d[x,0] :
                            merged_detec[y,35] = st9d[x,0]
                        if merged_detec[y,36] >  st9d[x,1] :
                            merged_detec[y,36] = st9d[x,1]
                        if merged_detec[y,37] <  st9d[x,1] :
                            merged_detec[y,37] = st9d[x,1]                  
                        
                        added = 1
                    
        if added == 0 :
            merged_detec = np.lib.pad(merged_detec, ((0,1),(0,0)), 'constant', constant_values=(0))
            merged_detec[-1,0] = 1 
            merged_detec[-1,9] = 1    
            merged_detec[-1,34] = st9d[x,0]    
            merged_detec[-1,35] = st9d[x,0]  
            merged_detec[-1,36] = st9d[x,1] 
            merged_detec[-1,37] = st9d[x,1]  
            merged_detec[-1,38] = st9d[x,3]                  
else:
    merged_detec = np.lib.pad(merged_detec, ((0,len(st9d)),(0,0)), 'constant', constant_values=(0))
    merged_detec[:,0] = 1              
    merged_detec[:,9] = 1            
    merged_detec[:,34] = st9d[:,0]    
    merged_detec[:,35] = st9d[:,0]  
    merged_detec[:,36] = st9d[:,1] 
    merged_detec[:,37] = st9d[:,1] 
    merged_detec[:,38] = st9d[:,3]  

#10
if len(merged_detec) > 0:
    for x in range(0,len(st10d)):
        added = 0
        for y in range(0,len(merged_detec)):
            if added == 0:
                if (abs(merged_detec[y,34] - st10d[x,0]) < 60*30 and abs(merged_detec[y,37] - st10d[x,1]) < 60*120) or (abs(merged_detec[y,35] - st10d[x,0]) < 60*10 and abs(merged_detec[y,37] - st10d[x,1]) < 60*120) or (abs(merged_detec[y,37] - st10d[x,1]) < 60*30 and abs(merged_detec[y,34] - st10d[x,0]) < 60*120) or (abs(merged_detec[y,36] - st10d[x,1]) < 60*20 and abs(merged_detec[y,34] - st10d[x,0]) < 60*120):

                    if merged_detec[y,34] < st10d[x,0] < merged_detec[y,37] or merged_detec[y,34] < st10d[x,1] < merged_detec[y,37] or  st10d[x,0] < merged_detec[y,34] < st10d[x,1] :
                        
                        merged_detec[y,0] += 1
                        merged_detec[y,10] = 1
                        
                        if merged_detec[y,34] >  st10d[x,0] :
                            merged_detec[y,34] = st10d[x,0]
                        if merged_detec[y,35] <  st10d[x,0] :
                            merged_detec[y,35] = st10d[x,0]
                        if merged_detec[y,36] >  st10d[x,1] :
                            merged_detec[y,36] = st10d[x,1]
                        if merged_detec[y,37] <  st10d[x,1] :
                            merged_detec[y,37] = st10d[x,1]                 
                        
                        added = 1
                    
        if added == 0 :
            merged_detec = np.lib.pad(merged_detec, ((0,1),(0,0)), 'constant', constant_values=(0))
            merged_detec[-1,0] = 1 
            merged_detec[-1,10] = 1    
            merged_detec[-1,34] = st10d[x,0]    
            merged_detec[-1,35] = st10d[x,0]  
            merged_detec[-1,36] = st10d[x,1] 
            merged_detec[-1,37] = st10d[x,1]  
            merged_detec[-1,38] = st10d[x,3]                  
else:
    merged_detec = np.lib.pad(merged_detec, ((0,len(st10d)),(0,0)), 'constant', constant_values=(0))
    merged_detec[:,0] = 1              
    merged_detec[:,10] = 1            
    merged_detec[:,34] = st10d[:,0]    
    merged_detec[:,35] = st10d[:,0]  
    merged_detec[:,36] = st10d[:,1] 
    merged_detec[:,37] = st10d[:,1] 
    merged_detec[:,38] = st10d[:,3]  

#11
if len(merged_detec) > 0:
    for x in range(0,len(st11d)):
        added = 0
        for y in range(0,len(merged_detec)):
            if added == 0:
                if (abs(merged_detec[y,34] - st11d[x,0]) < 60*30 and abs(merged_detec[y,37] - st11d[x,1]) < 60*120) or (abs(merged_detec[y,35] - st11d[x,0]) < 60*10 and abs(merged_detec[y,37] - st11d[x,1]) < 60*120) or (abs(merged_detec[y,37] - st11d[x,1]) < 60*30 and abs(merged_detec[y,34] - st11d[x,0]) < 60*120) or (abs(merged_detec[y,36] - st11d[x,1]) < 60*20 and abs(merged_detec[y,34] - st11d[x,0]) < 60*120):

                    if merged_detec[y,34] < st11d[x,0] < merged_detec[y,37] or merged_detec[y,34] < st11d[x,1] < merged_detec[y,37] or  st11d[x,0] < merged_detec[y,34] < st11d[x,1] :
                        
                        merged_detec[y,0] += 1
                        merged_detec[y,11] = 1
                        
                        if merged_detec[y,34] >  st11d[x,0] :
                            merged_detec[y,34] = st11d[x,0]
                        if merged_detec[y,35] <  st11d[x,0] :
                            merged_detec[y,35] = st11d[x,0]
                        if merged_detec[y,36] >  st11d[x,1] :
                            merged_detec[y,36] = st11d[x,1]
                        if merged_detec[y,37] <  st11d[x,1] :
                            merged_detec[y,37] = st11d[x,1]               
                        
                        added = 1
                    
        if added == 0 :
            merged_detec = np.lib.pad(merged_detec, ((0,1),(0,0)), 'constant', constant_values=(0))
            merged_detec[-1,0] = 1 
            merged_detec[-1,11] = 1    
            merged_detec[-1,34] = st11d[x,0]    
            merged_detec[-1,35] = st11d[x,0]  
            merged_detec[-1,36] = st11d[x,1] 
            merged_detec[-1,37] = st11d[x,1] 
            merged_detec[-1,38] = st11d[x,3]                   
else:
    merged_detec = np.lib.pad(merged_detec, ((0,len(st11d)),(0,0)), 'constant', constant_values=(0))
    merged_detec[:,0] = 1              
    merged_detec[:,11] = 1            
    merged_detec[:,34] = st11d[:,0]    
    merged_detec[:,35] = st11d[:,0]  
    merged_detec[:,36] = st11d[:,1] 
    merged_detec[:,37] = st11d[:,1]
    merged_detec[:,38] = st11d[:,3]  

#12
if len(merged_detec) > 0:
    for x in range(0,len(st12d)):
        added = 0
        for y in range(0,len(merged_detec)):
            if added == 0:
                if (abs(merged_detec[y,34] - st12d[x,0]) < 60*30 and abs(merged_detec[y,37] - st12d[x,1]) < 60*120) or (abs(merged_detec[y,35] - st12d[x,0]) < 60*10 and abs(merged_detec[y,37] - st12d[x,1]) < 60*120) or (abs(merged_detec[y,37] - st12d[x,1]) < 60*30 and abs(merged_detec[y,34] - st12d[x,0]) < 60*120) or (abs(merged_detec[y,36] - st12d[x,1]) < 60*20 and abs(merged_detec[y,34] - st12d[x,0]) < 60*120):

                    if merged_detec[y,34] < st12d[x,0] < merged_detec[y,37] or merged_detec[y,34] < st12d[x,1] < merged_detec[y,37] or  st12d[x,0] < merged_detec[y,34] < st12d[x,1] :
                        
                        merged_detec[y,0] += 1
                        merged_detec[y,12] = 1
                        
                        if merged_detec[y,34] >  st12d[x,0] :
                            merged_detec[y,34] = st12d[x,0]
                        if merged_detec[y,35] <  st12d[x,0] :
                            merged_detec[y,35] = st12d[x,0]
                        if merged_detec[y,36] >  st12d[x,1] :
                            merged_detec[y,36] = st12d[x,1]
                        if merged_detec[y,37] <  st12d[x,1] :
                            merged_detec[y,37] = st12d[x,1]                    
                        
                        added = 1
                    
        if added == 0 :
            merged_detec = np.lib.pad(merged_detec, ((0,1),(0,0)), 'constant', constant_values=(0))
            merged_detec[-1,0] = 1 
            merged_detec[-1,12] = 1    
            merged_detec[-1,34] = st12d[x,0]    
            merged_detec[-1,35] = st12d[x,0]  
            merged_detec[-1,36] = st12d[x,1] 
            merged_detec[-1,37] = st12d[x,1]  
            merged_detec[-1,38] = st12d[x,3]                  
else:
    merged_detec = np.lib.pad(merged_detec, ((0,len(st12d)),(0,0)), 'constant', constant_values=(0))
    merged_detec[:,0] = 1              
    merged_detec[:,12] = 1            
    merged_detec[:,34] = st12d[:,0]    
    merged_detec[:,35] = st12d[:,0]  
    merged_detec[:,36] = st12d[:,1] 
    merged_detec[:,37] = st12d[:,1]
    merged_detec[:,38] = st12d[:,3]  

#13
if len(merged_detec) > 0:
    for x in range(0,len(st13d)):
        added = 0
        for y in range(0,len(merged_detec)):
            if added == 0:
                if (abs(merged_detec[y,34] - st13d[x,0]) < 60*30 and abs(merged_detec[y,37] - st13d[x,1]) < 60*120) or (abs(merged_detec[y,35] - st13d[x,0]) < 60*10 and abs(merged_detec[y,37] - st13d[x,1]) < 60*120) or (abs(merged_detec[y,37] - st13d[x,1]) < 60*30 and abs(merged_detec[y,34] - st13d[x,0]) < 60*120) or (abs(merged_detec[y,36] - st13d[x,1]) < 60*20 and abs(merged_detec[y,34] - st13d[x,0]) < 60*120):

                    if merged_detec[y,34] < st13d[x,0] < merged_detec[y,37] or merged_detec[y,34] < st13d[x,1] < merged_detec[y,37] or  st13d[x,0] < merged_detec[y,34] < st13d[x,1] :
                        
                        merged_detec[y,0] += 1
                        merged_detec[y,13] = 1
                        
                        if merged_detec[y,34] >  st13d[x,0] :
                            merged_detec[y,34] = st13d[x,0]
                        if merged_detec[y,35] <  st13d[x,0] :
                            merged_detec[y,35] = st13d[x,0]
                        if merged_detec[y,36] >  st13d[x,1] :
                            merged_detec[y,36] = st13d[x,1]
                        if merged_detec[y,37] <  st13d[x,1] :
                            merged_detec[y,37] = st13d[x,1]                    
                        
                        added = 1
                    
        if added == 0 :
            merged_detec = np.lib.pad(merged_detec, ((0,1),(0,0)), 'constant', constant_values=(0))
            merged_detec[-1,0] = 1 
            merged_detec[-1,13] = 1    
            merged_detec[-1,34] = st13d[x,0]    
            merged_detec[-1,35] = st13d[x,0]  
            merged_detec[-1,36] = st13d[x,1] 
            merged_detec[-1,37] = st13d[x,1]   
            merged_detec[-1,38] = st13d[x,3]                 
else:
    merged_detec = np.lib.pad(merged_detec, ((0,len(st13d)),(0,0)), 'constant', constant_values=(0))
    merged_detec[:,0] = 1              
    merged_detec[:,13] = 1            
    merged_detec[:,34] = st13d[:,0]    
    merged_detec[:,35] = st13d[:,0]  
    merged_detec[:,36] = st13d[:,1] 
    merged_detec[:,37] = st13d[:,1] 
    merged_detec[:,38] = st13d[:,3]  

#14
if len(merged_detec) > 0:
    for x in range(0,len(st14d)):
        added = 0
        for y in range(0,len(merged_detec)):
            if added == 0:
               if (abs(merged_detec[y,34] - st14d[x,0]) < 60*30 and abs(merged_detec[y,37] - st14d[x,1]) < 60*120) or (abs(merged_detec[y,35] - st14d[x,0]) < 60*10 and abs(merged_detec[y,37] - st14d[x,1]) < 60*120) or (abs(merged_detec[y,37] - st14d[x,1]) < 60*30 and abs(merged_detec[y,34] - st14d[x,0]) < 60*120) or (abs(merged_detec[y,36] - st14d[x,1]) < 60*20 and abs(merged_detec[y,34] - st14d[x,0]) < 60*120):
 
                    if merged_detec[y,34] < st14d[x,0] < merged_detec[y,37] or merged_detec[y,34] < st14d[x,1] < merged_detec[y,37] or  st14d[x,0] < merged_detec[y,34] < st14d[x,1] :
                        
                        merged_detec[y,0] += 1
                        merged_detec[y,14] = 1
                        
                        if merged_detec[y,34] >  st14d[x,0] :
                            merged_detec[y,34] = st14d[x,0]
                        if merged_detec[y,35] <  st14d[x,0] :
                            merged_detec[y,35] = st14d[x,0]
                        if merged_detec[y,36] >  st14d[x,1] :
                            merged_detec[y,36] = st14d[x,1]
                        if merged_detec[y,37] <  st14d[x,1] :
                            merged_detec[y,37] = st14d[x,1]                  
                        
                        added = 1
                    
        if added == 0 :
            merged_detec = np.lib.pad(merged_detec, ((0,1),(0,0)), 'constant', constant_values=(0))
            merged_detec[-1,0] = 1 
            merged_detec[-1,14] = 1    
            merged_detec[-1,34] = st14d[x,0]    
            merged_detec[-1,35] = st14d[x,0]  
            merged_detec[-1,36] = st14d[x,1] 
            merged_detec[-1,37] = st14d[x,1] 
            merged_detec[-1,38] = st14d[x,3]                   
else:
    merged_detec = np.lib.pad(merged_detec, ((0,len(st14d)),(0,0)), 'constant', constant_values=(0))
    merged_detec[:,0] = 1              
    merged_detec[:,14] = 1            
    merged_detec[:,34] = st14d[:,0]    
    merged_detec[:,35] = st14d[:,0]  
    merged_detec[:,36] = st14d[:,1] 
    merged_detec[:,37] = st14d[:,1] 
    merged_detec[:,38] = st14d[:,3]  

#15
if len(merged_detec) > 0:
    for x in range(0,len(st15d)):
        added = 0
        for y in range(0,len(merged_detec)):
            if added == 0:
               if (abs(merged_detec[y,34] - st15d[x,0]) < 60*30 and abs(merged_detec[y,37] - st15d[x,1]) < 60*120) or (abs(merged_detec[y,35] - st15d[x,0]) < 60*10 and abs(merged_detec[y,37] - st15d[x,1]) < 60*120) or (abs(merged_detec[y,37] - st15d[x,1]) < 60*30 and abs(merged_detec[y,34] - st15d[x,0]) < 60*120) or (abs(merged_detec[y,36] - st15d[x,1]) < 60*20 and abs(merged_detec[y,34] - st15d[x,0]) < 60*120):
 
                    if merged_detec[y,34] < st15d[x,0] < merged_detec[y,37] or merged_detec[y,34] < st15d[x,1] < merged_detec[y,37] or  st15d[x,0] < merged_detec[y,34] < st15d[x,1] :
                        
                        merged_detec[y,0] += 1
                        merged_detec[y,15] = 1
                        
                        if merged_detec[y,34] >  st15d[x,0] :
                            merged_detec[y,34] = st15d[x,0]
                        if merged_detec[y,35] <  st15d[x,0] :
                            merged_detec[y,35] = st15d[x,0]
                        if merged_detec[y,36] >  st15d[x,1] :
                            merged_detec[y,36] = st15d[x,1]
                        if merged_detec[y,37] <  st15d[x,1] :
                            merged_detec[y,37] = st15d[x,1]                     
                        
                        added = 1
                    
        if added == 0 :
            merged_detec = np.lib.pad(merged_detec, ((0,1),(0,0)), 'constant', constant_values=(0))
            merged_detec[-1,0] = 1 
            merged_detec[-1,15] = 1    
            merged_detec[-1,34] = st15d[x,0]    
            merged_detec[-1,35] = st15d[x,0]  
            merged_detec[-1,36] = st15d[x,1] 
            merged_detec[-1,37] = st15d[x,1]  
            merged_detec[-1,38] = st15d[x,3]                  
else:
    merged_detec = np.lib.pad(merged_detec, ((0,len(st15d)),(0,0)), 'constant', constant_values=(0))
    merged_detec[:,0] = 1              
    merged_detec[:,15] = 1            
    merged_detec[:,34] = st15d[:,0]    
    merged_detec[:,35] = st15d[:,0]  
    merged_detec[:,36] = st15d[:,1] 
    merged_detec[:,37] = st15d[:,1] 
    merged_detec[:,38] = st15d[:,3]  

#16
if len(merged_detec) > 0:
    for x in range(0,len(st16d)):
        added = 0
        for y in range(0,len(merged_detec)):
            if added == 0:
                if (abs(merged_detec[y,34] - st16d[x,0]) < 60*30 and abs(merged_detec[y,37] - st16d[x,1]) < 60*120) or (abs(merged_detec[y,35] - st16d[x,0]) < 60*10 and abs(merged_detec[y,37] - st16d[x,1]) < 60*120) or (abs(merged_detec[y,37] - st16d[x,1]) < 60*30 and abs(merged_detec[y,34] - st16d[x,0]) < 60*120) or (abs(merged_detec[y,36] - st16d[x,1]) < 60*20 and abs(merged_detec[y,34] - st16d[x,0]) < 60*120):

                    if merged_detec[y,34] < st16d[x,0] < merged_detec[y,37] or merged_detec[y,34] < st16d[x,1] < merged_detec[y,37] or  st16d[x,0] < merged_detec[y,34] < st16d[x,1] :
                        
                        merged_detec[y,0] += 1
                        merged_detec[y,16] = 1
                        
                        if merged_detec[y,34] >  st16d[x,0] :
                            merged_detec[y,34] = st16d[x,0]
                        if merged_detec[y,35] <  st16d[x,0] :
                            merged_detec[y,35] = st16d[x,0]
                        if merged_detec[y,36] >  st16d[x,1] :
                            merged_detec[y,36] = st16d[x,1]
                        if merged_detec[y,37] <  st16d[x,1] :
                            merged_detec[y,37] = st16d[x,1]                     
                        
                        added = 1
                    
        if added == 0 :
            merged_detec = np.lib.pad(merged_detec, ((0,1),(0,0)), 'constant', constant_values=(0))
            merged_detec[-1,0] = 1 
            merged_detec[-1,16] = 1    
            merged_detec[-1,34] = st16d[x,0]    
            merged_detec[-1,35] = st16d[x,0]  
            merged_detec[-1,36] = st16d[x,1] 
            merged_detec[-1,37] = st16d[x,1] 
            merged_detec[-1,38] = st16d[x,3]                   
else:
    merged_detec = np.lib.pad(merged_detec, ((0,len(st16d)),(0,0)), 'constant', constant_values=(0))
    merged_detec[:,0] = 1              
    merged_detec[:,16] = 1            
    merged_detec[:,34] = st16d[:,0]    
    merged_detec[:,35] = st16d[:,0]  
    merged_detec[:,36] = st16d[:,1] 
    merged_detec[:,37] = st16d[:,1] 
    merged_detec[:,38] = st16d[:,3]  

#17
if len(merged_detec) > 0:
    for x in range(0,len(st17d)):
        added = 0
        for y in range(0,len(merged_detec)):
            if added == 0:
                if (abs(merged_detec[y,34] - st17d[x,0]) < 60*30 and abs(merged_detec[y,37] - st17d[x,1]) < 60*120) or (abs(merged_detec[y,35] - st17d[x,0]) < 60*10 and abs(merged_detec[y,37] - st17d[x,1]) < 60*120) or (abs(merged_detec[y,37] - st17d[x,1]) < 60*30 and abs(merged_detec[y,34] - st17d[x,0]) < 60*120) or (abs(merged_detec[y,36] - st17d[x,1]) < 60*20 and abs(merged_detec[y,34] - st17d[x,0]) < 60*120):

                    if merged_detec[y,34] < st17d[x,0] < merged_detec[y,37] or merged_detec[y,34] < st17d[x,1] < merged_detec[y,37] or  st17d[x,0] < merged_detec[y,34] < st17d[x,1] :
                        
                        merged_detec[y,0] += 1
                        merged_detec[y,17] = 1
                        
                        if merged_detec[y,34] >  st17d[x,0] :
                            merged_detec[y,34] = st17d[x,0]
                        if merged_detec[y,35] <  st17d[x,0] :
                            merged_detec[y,35] = st17d[x,0]
                        if merged_detec[y,36] >  st17d[x,1] :
                            merged_detec[y,36] = st17d[x,1]
                        if merged_detec[y,37] <  st17d[x,1] :
                            merged_detec[y,37] = st17d[x,1]              
                        
                        added = 1
                    
        if added == 0 :
            merged_detec = np.lib.pad(merged_detec, ((0,1),(0,0)), 'constant', constant_values=(0))
            merged_detec[-1,0] = 1 
            merged_detec[-1,17] = 1    
            merged_detec[-1,34] = st17d[x,0]    
            merged_detec[-1,35] = st17d[x,0]  
            merged_detec[-1,36] = st17d[x,1] 
            merged_detec[-1,37] = st17d[x,1]  
            merged_detec[-1,38] = st17d[x,3]                  
else:
    merged_detec = np.lib.pad(merged_detec, ((0,len(st17d)),(0,0)), 'constant', constant_values=(0))
    merged_detec[:,0] = 1              
    merged_detec[:,17] = 1            
    merged_detec[:,34] = st17d[:,0]    
    merged_detec[:,35] = st17d[:,0]  
    merged_detec[:,36] = st17d[:,1] 
    merged_detec[:,37] = st17d[:,1] 
    merged_detec[:,38] = st17d[:,3]  

#18
if len(merged_detec) > 0:
    for x in range(0,len(st18d)):
        added = 0
        for y in range(0,len(merged_detec)):
            if added == 0:
                if (abs(merged_detec[y,34] - st18d[x,0]) < 60*30 and abs(merged_detec[y,37] - st18d[x,1]) < 60*120) or (abs(merged_detec[y,35] - st18d[x,0]) < 60*10 and abs(merged_detec[y,37] - st18d[x,1]) < 60*120) or (abs(merged_detec[y,37] - st18d[x,1]) < 60*30 and abs(merged_detec[y,34] - st18d[x,0]) < 60*120) or (abs(merged_detec[y,36] - st18d[x,1]) < 60*20 and abs(merged_detec[y,34] - st18d[x,0]) < 60*120):

                    if merged_detec[y,34] < st18d[x,0] < merged_detec[y,37] or merged_detec[y,34] < st18d[x,1] < merged_detec[y,37] or  st18d[x,0] < merged_detec[y,34] < st18d[x,1] :
                        
                        merged_detec[y,0] += 1
                        merged_detec[y,18] = 1
                        
                        if merged_detec[y,34] >  st18d[x,0] :
                            merged_detec[y,34] = st18d[x,0]
                        if merged_detec[y,35] <  st18d[x,0] :
                            merged_detec[y,35] = st18d[x,0]
                        if merged_detec[y,36] >  st18d[x,1] :
                            merged_detec[y,36] = st18d[x,1]
                        if merged_detec[y,37] <  st18d[x,1] :
                            merged_detec[y,37] = st18d[x,1]                    
                        
                        added = 1
                    
        if added == 0 :
            merged_detec = np.lib.pad(merged_detec, ((0,1),(0,0)), 'constant', constant_values=(0))
            merged_detec[-1,0] = 1 
            merged_detec[-1,18] = 1    
            merged_detec[-1,34] = st18d[x,0]    
            merged_detec[-1,35] = st18d[x,0]  
            merged_detec[-1,36] = st18d[x,1] 
            merged_detec[-1,37] = st18d[x,1]  
            merged_detec[-1,38] = st18d[x,3]                  
else:
    merged_detec = np.lib.pad(merged_detec, ((0,len(st18d)),(0,0)), 'constant', constant_values=(0))
    merged_detec[:,0] = 1              
    merged_detec[:,18] = 1            
    merged_detec[:,34] = st18d[:,0]    
    merged_detec[:,35] = st18d[:,0]  
    merged_detec[:,36] = st18d[:,1] 
    merged_detec[:,37] = st18d[:,1]  
    merged_detec[:,38] = st18d[:,3]  

#19
if len(merged_detec) > 0:
    for x in range(0,len(st19d)):
        added = 0
        for y in range(0,len(merged_detec)):
            if added == 0:
                if (abs(merged_detec[y,34] - st19d[x,0]) < 60*30 and abs(merged_detec[y,37] - st19d[x,1]) < 60*120) or (abs(merged_detec[y,35] - st19d[x,0]) < 60*10 and abs(merged_detec[y,37] - st19d[x,1]) < 60*120) or (abs(merged_detec[y,37] - st19d[x,1]) < 60*30 and abs(merged_detec[y,34] - st19d[x,0]) < 60*120) or (abs(merged_detec[y,36] - st19d[x,1]) < 60*20 and abs(merged_detec[y,34] - st19d[x,0]) < 60*120):

                    if merged_detec[y,34] < st19d[x,0] < merged_detec[y,37] or merged_detec[y,34] < st19d[x,1] < merged_detec[y,37] or  st19d[x,0] < merged_detec[y,34] < st19d[x,1] :
                        
                        merged_detec[y,0] += 1
                        merged_detec[y,19] = 1
                        
                        if merged_detec[y,34] >  st19d[x,0] :
                            merged_detec[y,34] = st19d[x,0]
                        if merged_detec[y,35] <  st19d[x,0] :
                            merged_detec[y,35] = st19d[x,0]
                        if merged_detec[y,36] >  st19d[x,1] :
                            merged_detec[y,36] = st19d[x,1]
                        if merged_detec[y,37] <  st19d[x,1] :
                            merged_detec[y,37] = st19d[x,1]                
                        
                        added = 1
                    
        if added == 0 :
            merged_detec = np.lib.pad(merged_detec, ((0,1),(0,0)), 'constant', constant_values=(0))
            merged_detec[-1,0] = 1 
            merged_detec[-1,19] = 1    
            merged_detec[-1,34] = st19d[x,0]    
            merged_detec[-1,35] = st19d[x,0]  
            merged_detec[-1,36] = st19d[x,1] 
            merged_detec[-1,37] = st19d[x,1]
            merged_detec[-1,38] = st19d[x,3]                    
else:
    merged_detec = np.lib.pad(merged_detec, ((0,len(st19d)),(0,0)), 'constant', constant_values=(0))
    merged_detec[:,0] = 1              
    merged_detec[:,19] = 1            
    merged_detec[:,34] = st19d[:,0]    
    merged_detec[:,35] = st19d[:,0]  
    merged_detec[:,36] = st19d[:,1] 
    merged_detec[:,37] = st19d[:,1] 
    merged_detec[:,38] = st19d[:,3]  

#20
if len(merged_detec) > 0:
    for x in range(0,len(st20d)):
        added = 0
        for y in range(0,len(merged_detec)):
            if added == 0:
                if (abs(merged_detec[y,34] - st20d[x,0]) < 60*30 and abs(merged_detec[y,37] - st20d[x,1]) < 60*120) or (abs(merged_detec[y,35] - st20d[x,0]) < 60*10 and abs(merged_detec[y,37] - st20d[x,1]) < 60*120) or (abs(merged_detec[y,37] - st20d[x,1]) < 60*30 and abs(merged_detec[y,34] - st20d[x,0]) < 60*120) or (abs(merged_detec[y,36] - st20d[x,1]) < 60*20 and abs(merged_detec[y,34] - st20d[x,0]) < 60*120):

                    if merged_detec[y,34] < st20d[x,0] < merged_detec[y,37] or merged_detec[y,34] < st20d[x,1] < merged_detec[y,37] or  st20d[x,0] < merged_detec[y,34] < st20d[x,1] :
                        
                        merged_detec[y,0] += 1
                        merged_detec[y,20] = 1
                        
                        if merged_detec[y,34] >  st20d[x,0] :
                            merged_detec[y,34] = st20d[x,0]
                        if merged_detec[y,35] <  st20d[x,0] :
                            merged_detec[y,35] = st20d[x,0]
                        if merged_detec[y,36] >  st20d[x,1] :
                            merged_detec[y,36] = st20d[x,1]
                        if merged_detec[y,37] <  st20d[x,1] :
                            merged_detec[y,37] = st20d[x,1]                   
                        
                        added = 1
                        
        if added == 0 :
            merged_detec = np.lib.pad(merged_detec, ((0,1),(0,0)), 'constant', constant_values=(0))
            merged_detec[-1,0] = 1 
            merged_detec[-1,20] = 1    
            merged_detec[-1,34] = st20d[x,0]    
            merged_detec[-1,35] = st20d[x,0]  
            merged_detec[-1,36] = st20d[x,1] 
            merged_detec[-1,37] = st20d[x,1] 
            merged_detec[-1,38] = st20d[x,3]                   
else:
    merged_detec = np.lib.pad(merged_detec, ((0,len(st20d)),(0,0)), 'constant', constant_values=(0))
    merged_detec[:,0] = 1              
    merged_detec[:,20] = 1            
    merged_detec[:,34] = st20d[:,0]    
    merged_detec[:,35] = st20d[:,0]  
    merged_detec[:,36] = st20d[:,1] 
    merged_detec[:,37] = st20d[:,1] 
    merged_detec[:,38] = st20d[:,3]  

#21
if len(merged_detec) > 0:
    for x in range(0,len(st21d)):
        added = 0
        for y in range(0,len(merged_detec)):
            if added == 0:
                if (abs(merged_detec[y,34] - st21d[x,0]) < 60*30 and abs(merged_detec[y,37] - st21d[x,1]) < 60*120) or (abs(merged_detec[y,35] - st21d[x,0]) < 60*10 and abs(merged_detec[y,37] - st21d[x,1]) < 60*120) or (abs(merged_detec[y,37] - st21d[x,1]) < 60*30 and abs(merged_detec[y,34] - st21d[x,0]) < 60*120) or (abs(merged_detec[y,36] - st21d[x,1]) < 60*20 and abs(merged_detec[y,34] - st21d[x,0]) < 60*120):

                    if merged_detec[y,34] < st21d[x,0] < merged_detec[y,37] or merged_detec[y,34] < st21d[x,1] < merged_detec[y,37] or  st21d[x,0] < merged_detec[y,34] < st21d[x,1] :
                        
                        merged_detec[y,0] += 1
                        merged_detec[y,21] = 1
                        
                        if merged_detec[y,34] >  st21d[x,0] :
                            merged_detec[y,34] = st21d[x,0]
                        if merged_detec[y,35] <  st21d[x,0] :
                            merged_detec[y,35] = st21d[x,0]
                        if merged_detec[y,36] >  st21d[x,1] :
                            merged_detec[y,36] = st21d[x,1]
                        if merged_detec[y,37] <  st21d[x,1] :
                            merged_detec[y,37] = st21d[x,1]                  
                        
                        added = 1
                    
        if added == 0 :
            merged_detec = np.lib.pad(merged_detec, ((0,1),(0,0)), 'constant', constant_values=(0))
            merged_detec[-1,0] = 1 
            merged_detec[-1,21] = 1    
            merged_detec[-1,34] = st21d[x,0]    
            merged_detec[-1,35] = st21d[x,0]  
            merged_detec[-1,36] = st21d[x,1] 
            merged_detec[-1,37] = st21d[x,1]
            merged_detec[-1,38] = st21d[x,3] 
                  
else:
    merged_detec = np.lib.pad(merged_detec, ((0,len(st21d)),(0,0)), 'constant', constant_values=(0))
    merged_detec[:,0] = 1              
    merged_detec[:,21] = 1            
    merged_detec[:,34] = st21d[:,0]    
    merged_detec[:,35] = st21d[:,0]  
    merged_detec[:,36] = st21d[:,1] 
    merged_detec[:,37] = st21d[:,1] 
    merged_detec[:,38] = st21d[:,3]  

#22
if len(merged_detec) > 0:
    for x in range(0,len(st22d)):
        added = 0
        for y in range(0,len(merged_detec)):
            if added == 0:
                if (abs(merged_detec[y,34] - st22d[x,0]) < 60*30 and abs(merged_detec[y,37] - st22d[x,1]) < 60*120) or (abs(merged_detec[y,35] - st22d[x,0]) < 60*10 and abs(merged_detec[y,37] - st22d[x,1]) < 60*120) or (abs(merged_detec[y,37] - st22d[x,1]) < 60*30 and abs(merged_detec[y,34] - st22d[x,0]) < 60*120) or (abs(merged_detec[y,36] - st22d[x,1]) < 60*20 and abs(merged_detec[y,34] - st22d[x,0]) < 60*120):

                    if merged_detec[y,34] < st22d[x,0] < merged_detec[y,37] or merged_detec[y,34] < st22d[x,1] < merged_detec[y,37] or  st22d[x,0] < merged_detec[y,34] < st22d[x,1] :
                        
                        merged_detec[y,0] += 1
                        merged_detec[y,22] = 1
                        
                        if merged_detec[y,34] >  st22d[x,0] :
                            merged_detec[y,34] = st22d[x,0]
                        if merged_detec[y,35] <  st22d[x,0] :
                            merged_detec[y,35] = st22d[x,0]
                        if merged_detec[y,36] >  st22d[x,1] :
                            merged_detec[y,36] = st22d[x,1]
                        if merged_detec[y,37] <  st22d[x,1] :
                            merged_detec[y,37] = st22d[x,1]                    
                        
                        added = 1
                    
        if added == 0 :
            merged_detec = np.lib.pad(merged_detec, ((0,1),(0,0)), 'constant', constant_values=(0))
            merged_detec[-1,0] = 1 
            merged_detec[-1,22] = 1    
            merged_detec[-1,34] = st22d[x,0]    
            merged_detec[-1,35] = st22d[x,0]  
            merged_detec[-1,36] = st22d[x,1] 
            merged_detec[-1,37] = st22d[x,1] 
            merged_detec[-1,38] = st22d[x,3]                   
else:
    merged_detec = np.lib.pad(merged_detec, ((0,len(st22d)),(0,0)), 'constant', constant_values=(0))
    merged_detec[:,0] = 1              
    merged_detec[:,22] = 1            
    merged_detec[:,34] = st22d[:,0]    
    merged_detec[:,35] = st22d[:,0]  
    merged_detec[:,36] = st22d[:,1] 
    merged_detec[:,37] = st22d[:,1] 
    merged_detec[:,38] = st22d[:,3]  

#23
if len(merged_detec) > 0:
    for x in range(0,len(st23d)):
        added = 0
        for y in range(0,len(merged_detec)):
            if added == 0:
                if (abs(merged_detec[y,34] - st23d[x,0]) < 60*30 and abs(merged_detec[y,37] - st23d[x,1]) < 60*120) or (abs(merged_detec[y,35] - st23d[x,0]) < 60*10 and abs(merged_detec[y,37] - st23d[x,1]) < 60*120) or (abs(merged_detec[y,37] - st23d[x,1]) < 60*30 and abs(merged_detec[y,34] - st23d[x,0]) < 60*120) or (abs(merged_detec[y,36] - st23d[x,1]) < 60*20 and abs(merged_detec[y,34] - st23d[x,0]) < 60*120):

                    if merged_detec[y,34] < st23d[x,0] < merged_detec[y,37] or merged_detec[y,34] < st23d[x,1] < merged_detec[y,37] or  st23d[x,0] < merged_detec[y,34] < st23d[x,1] :
                        
                        merged_detec[y,0] += 1
                        merged_detec[y,23] = 1
                        
                        if merged_detec[y,34] >  st23d[x,0] :
                            merged_detec[y,34] = st23d[x,0]
                        if merged_detec[y,35] <  st23d[x,0] :
                            merged_detec[y,35] = st23d[x,0]
                        if merged_detec[y,36] >  st23d[x,1] :
                            merged_detec[y,36] = st23d[x,1]
                        if merged_detec[y,37] <  st23d[x,1] :
                            merged_detec[y,37] = st23d[x,1]                   
                        
                        added = 1
                    
        if added == 0 :
            merged_detec = np.lib.pad(merged_detec, ((0,1),(0,0)), 'constant', constant_values=(0))
            merged_detec[-1,0] = 1 
            merged_detec[-1,23] = 1    
            merged_detec[-1,34] = st23d[x,0]    
            merged_detec[-1,35] = st23d[x,0]  
            merged_detec[-1,36] = st23d[x,1] 
            merged_detec[-1,37] = st23d[x,1] 
            merged_detec[-1,38] = st23d[x,3]                   
else:
    merged_detec = np.lib.pad(merged_detec, ((0,len(st23d)),(0,0)), 'constant', constant_values=(0))
    merged_detec[:,0] = 1              
    merged_detec[:,23] = 1            
    merged_detec[:,34] = st23d[:,0]    
    merged_detec[:,35] = st23d[:,0]  
    merged_detec[:,36] = st23d[:,1] 
    merged_detec[:,37] = st23d[:,1] 
    merged_detec[:,38] = st23d[:,3]  

#24
if len(merged_detec) > 0:
    for x in range(0,len(st24d)):
        added = 0
        for y in range(0,len(merged_detec)):
            if added == 0:
                if (abs(merged_detec[y,34] - st24d[x,0]) < 60*30 and abs(merged_detec[y,37] - st24d[x,1]) < 60*120) or (abs(merged_detec[y,35] - st24d[x,0]) < 60*10 and abs(merged_detec[y,37] - st24d[x,1]) < 60*120) or (abs(merged_detec[y,37] - st24d[x,1]) < 60*30 and abs(merged_detec[y,34] - st24d[x,0]) < 60*120) or (abs(merged_detec[y,36] - st24d[x,1]) < 60*20 and abs(merged_detec[y,34] - st24d[x,0]) < 60*120):

                    if merged_detec[y,34] < st24d[x,0] < merged_detec[y,37] or merged_detec[y,34] < st24d[x,1] < merged_detec[y,37] or  st24d[x,0] < merged_detec[y,34] < st24d[x,1] :
                        
                        merged_detec[y,0] += 1
                        merged_detec[y,24] = 1
                        
                        if merged_detec[y,34] >  st24d[x,0] :
                            merged_detec[y,34] = st24d[x,0]
                        if merged_detec[y,35] <  st24d[x,0] :
                            merged_detec[y,35] = st24d[x,0]
                        if merged_detec[y,36] >  st24d[x,1] :
                            merged_detec[y,36] = st24d[x,1]
                        if merged_detec[y,37] <  st24d[x,1] :
                            merged_detec[y,37] = st24d[x,1]                  
                        
                        added = 1
                    
        if added == 0 :
            merged_detec = np.lib.pad(merged_detec, ((0,1),(0,0)), 'constant', constant_values=(0))
            merged_detec[-1,0] = 1 
            merged_detec[-1,24] = 1    
            merged_detec[-1,34] = st24d[x,0]    
            merged_detec[-1,35] = st24d[x,0]  
            merged_detec[-1,36] = st24d[x,1] 
            merged_detec[-1,37] = st24d[x,1] 
            merged_detec[-1,38] = st24d[x,3]                   
else:
    merged_detec = np.lib.pad(merged_detec, ((0,len(st24d)),(0,0)), 'constant', constant_values=(0))
    merged_detec[:,0] = 1              
    merged_detec[:,24] = 1            
    merged_detec[:,34] = st24d[:,0]    
    merged_detec[:,35] = st24d[:,0]  
    merged_detec[:,36] = st24d[:,1] 
    merged_detec[:,37] = st24d[:,1] 
    merged_detec[:,38] = st24d[:,3]  

#25
if len(merged_detec) > 0:
    for x in range(0,len(st25d)):
        added = 0
        for y in range(0,len(merged_detec)):
            if added == 0:
                if (abs(merged_detec[y,34] - st25d[x,0]) < 60*30 and abs(merged_detec[y,37] - st25d[x,1]) < 60*120) or (abs(merged_detec[y,35] - st25d[x,0]) < 60*10 and abs(merged_detec[y,37] - st25d[x,1]) < 60*120) or (abs(merged_detec[y,37] - st25d[x,1]) < 60*30 and abs(merged_detec[y,34] - st25d[x,0]) < 60*120) or (abs(merged_detec[y,36] - st25d[x,1]) < 60*20 and abs(merged_detec[y,34] - st25d[x,0]) < 60*120):

                    if merged_detec[y,34] < st25d[x,0] < merged_detec[y,37] or merged_detec[y,34] < st25d[x,1] < merged_detec[y,37] or  st25d[x,0] < merged_detec[y,34] < st25d[x,1] :
                        
                        merged_detec[y,0] += 1
                        merged_detec[y,25] = 1
                        
                        if merged_detec[y,34] >  st25d[x,0] :
                            merged_detec[y,34] = st25d[x,0]
                        if merged_detec[y,35] <  st25d[x,0] :
                            merged_detec[y,35] = st25d[x,0]
                        if merged_detec[y,36] >  st25d[x,1] :
                            merged_detec[y,36] = st25d[x,1]
                        if merged_detec[y,37] <  st25d[x,1] :
                            merged_detec[y,37] = st25d[x,1]                    
                        
                        added = 1
                    
        if added == 0 :
            merged_detec = np.lib.pad(merged_detec, ((0,1),(0,0)), 'constant', constant_values=(0))
            merged_detec[-1,0] = 1 
            merged_detec[-1,25] = 1    
            merged_detec[-1,34] = st25d[x,0]    
            merged_detec[-1,35] = st25d[x,0]  
            merged_detec[-1,36] = st25d[x,1] 
            merged_detec[-1,37] = st25d[x,1] 
            merged_detec[-1,38] = st25d[x,3]                   
else:
    merged_detec = np.lib.pad(merged_detec, ((0,len(st25d)),(0,0)), 'constant', constant_values=(0))
    merged_detec[:,0] = 1              
    merged_detec[:,25] = 1            
    merged_detec[:,34] = st25d[:,0]    
    merged_detec[:,35] = st25d[:,0]  
    merged_detec[:,36] = st25d[:,1] 
    merged_detec[:,37] = st25d[:,1]   
    merged_detec[:,38] = st25d[:,3]  

#26
if len(merged_detec) > 0:
    for x in range(0,len(st26d)):
        added = 0
        for y in range(0,len(merged_detec)):
            if added == 0:
                if (abs(merged_detec[y,34] - st26d[x,0]) < 60*30 and abs(merged_detec[y,37] - st26d[x,1]) < 60*120) or (abs(merged_detec[y,35] - st26d[x,0]) < 60*10 and abs(merged_detec[y,37] - st26d[x,1]) < 60*120) or (abs(merged_detec[y,37] - st26d[x,1]) < 60*30 and abs(merged_detec[y,34] - st26d[x,0]) < 60*120) or (abs(merged_detec[y,36] - st26d[x,1]) < 60*20 and abs(merged_detec[y,34] - st26d[x,0]) < 60*120):

                    if merged_detec[y,34] < st26d[x,0] < merged_detec[y,37] or merged_detec[y,34] < st26d[x,1] < merged_detec[y,37] or  st26d[x,0] < merged_detec[y,34] < st26d[x,1] :
                        
                        merged_detec[y,0] += 1
                        merged_detec[y,26] = 1
                        
                        if merged_detec[y,34] >  st26d[x,0] :
                            merged_detec[y,34] = st26d[x,0]
                        if merged_detec[y,35] <  st26d[x,0] :
                            merged_detec[y,35] = st26d[x,0]
                        if merged_detec[y,36] >  st26d[x,1] :
                            merged_detec[y,36] = st26d[x,1]
                        if merged_detec[y,37] <  st26d[x,1] :
                            merged_detec[y,37] = st26d[x,1]                    
                        
                        added = 1
                        
        if added == 0 :
            merged_detec = np.lib.pad(merged_detec, ((0,1),(0,0)), 'constant', constant_values=(0))
            merged_detec[-1,0] = 1 
            merged_detec[-1,26] = 1    
            merged_detec[-1,34] = st26d[x,0]    
            merged_detec[-1,35] = st26d[x,0]  
            merged_detec[-1,36] = st26d[x,1] 
            merged_detec[-1,37] = st26d[x,1]
            merged_detec[-1,38] = st26d[x,3]                    
else:
    merged_detec = np.lib.pad(merged_detec, ((0,len(st26d)),(0,0)), 'constant', constant_values=(0))
    merged_detec[:,0] = 1              
    merged_detec[:,26] = 1            
    merged_detec[:,34] = st26d[:,0]    
    merged_detec[:,35] = st26d[:,0]  
    merged_detec[:,36] = st26d[:,1] 
    merged_detec[:,37] = st26d[:,1]
    merged_detec[:,38] = st26d[:,3]  

#27
if len(merged_detec) > 0:
    for x in range(0,len(st27d)):
        added = 0
        for y in range(0,len(merged_detec)):
            if added == 0:
                if (abs(merged_detec[y,34] - st27d[x,0]) < 60*30 and abs(merged_detec[y,37] - st27d[x,1]) < 60*120) or (abs(merged_detec[y,35] - st27d[x,0]) < 60*10 and abs(merged_detec[y,37] - st27d[x,1]) < 60*120) or (abs(merged_detec[y,37] - st27d[x,1]) < 60*30 and abs(merged_detec[y,34] - st27d[x,0]) < 60*120) or (abs(merged_detec[y,36] - st27d[x,1]) < 60*20 and abs(merged_detec[y,34] - st27d[x,0]) < 60*120):

                    if merged_detec[y,34] < st27d[x,0] < merged_detec[y,37] or merged_detec[y,34] < st27d[x,1] < merged_detec[y,37] or  st27d[x,0] < merged_detec[y,34] < st27d[x,1] :
                        
                        merged_detec[y,0] += 1
                        merged_detec[y,27] = 1
                        
                        if merged_detec[y,34] >  st27d[x,0] :
                            merged_detec[y,34] = st27d[x,0]
                        if merged_detec[y,35] <  st27d[x,0] :
                            merged_detec[y,35] = st27d[x,0]
                        if merged_detec[y,36] >  st27d[x,1] :
                            merged_detec[y,36] = st27d[x,1]
                        if merged_detec[y,37] <  st27d[x,1] :
                            merged_detec[y,37] = st27d[x,1]                   
                        
                        added = 1
                    
        if added == 0 :
            merged_detec = np.lib.pad(merged_detec, ((0,1),(0,0)), 'constant', constant_values=(0))
            merged_detec[-1,0] = 1 
            merged_detec[-1,27] = 1    
            merged_detec[-1,34] = st27d[x,0]    
            merged_detec[-1,35] = st27d[x,0]  
            merged_detec[-1,36] = st27d[x,1] 
            merged_detec[-1,37] = st27d[x,1] 
            merged_detec[-1,38] = st27d[x,3]                   
else:
    merged_detec = np.lib.pad(merged_detec, ((0,len(st27d)),(0,0)), 'constant', constant_values=(0))
    merged_detec[:,0] = 1              
    merged_detec[:,27] = 1            
    merged_detec[:,34] = st27d[:,0]    
    merged_detec[:,35] = st27d[:,0]  
    merged_detec[:,36] = st27d[:,1] 
    merged_detec[:,37] = st27d[:,1]  
    merged_detec[:,38] = st27d[:,3]  

#28
if len(merged_detec) > 0:
    for x in range(0,len(st28d)):
        added = 0
        for y in range(0,len(merged_detec)):
            if added == 0:
                if (abs(merged_detec[y,34] - st28d[x,0]) < 60*30 and abs(merged_detec[y,37] - st28d[x,1]) < 60*120) or (abs(merged_detec[y,35] - st28d[x,0]) < 60*10 and abs(merged_detec[y,37] - st28d[x,1]) < 60*120) or (abs(merged_detec[y,37] - st28d[x,1]) < 60*30 and abs(merged_detec[y,34] - st28d[x,0]) < 60*120) or (abs(merged_detec[y,36] - st28d[x,1]) < 60*20 and abs(merged_detec[y,34] - st28d[x,0]) < 60*120):

                    if merged_detec[y,34] < st28d[x,0] < merged_detec[y,37] or merged_detec[y,34] < st28d[x,1] < merged_detec[y,37] or  st28d[x,0] < merged_detec[y,34] < st28d[x,1] :
                        
                        merged_detec[y,0] += 1
                        merged_detec[y,28] = 1
                        
                        if merged_detec[y,34] >  st28d[x,0] :
                            merged_detec[y,34] = st28d[x,0]
                        if merged_detec[y,35] <  st28d[x,0] :
                            merged_detec[y,35] = st28d[x,0]
                        if merged_detec[y,36] >  st28d[x,1] :
                            merged_detec[y,36] = st28d[x,1]
                        if merged_detec[y,37] <  st28d[x,1] :
                            merged_detec[y,37] = st28d[x,1]                    
                        
                        added = 1
                    
        if added == 0 :
            merged_detec = np.lib.pad(merged_detec, ((0,1),(0,0)), 'constant', constant_values=(0))
            merged_detec[-1,0] = 1 
            merged_detec[-1,28] = 1    
            merged_detec[-1,34] = st28d[x,0]    
            merged_detec[-1,35] = st28d[x,0]  
            merged_detec[-1,36] = st28d[x,1] 
            merged_detec[-1,37] = st28d[x,1] 
            merged_detec[-1,38] = st28d[x,3]                      
else:
    merged_detec = np.lib.pad(merged_detec, ((0,len(st28d)),(0,0)), 'constant', constant_values=(0))
    merged_detec[:,0] = 1              
    merged_detec[:,28] = 1            
    merged_detec[:,34] = st28d[:,0]    
    merged_detec[:,35] = st28d[:,0]  
    merged_detec[:,36] = st28d[:,1] 
    merged_detec[:,37] = st28d[:,1]
    merged_detec[:,38] = st28d[:,3]  

#29
if len(merged_detec) > 0:
    for x in range(0,len(st29d)):
        added = 0
        for y in range(0,len(merged_detec)):
            if added == 0:
                if (abs(merged_detec[y,34] - st29d[x,0]) < 60*30 and abs(merged_detec[y,37] - st29d[x,1]) < 60*120) or (abs(merged_detec[y,35] - st29d[x,0]) < 60*10 and abs(merged_detec[y,37] - st29d[x,1]) < 60*120) or (abs(merged_detec[y,37] - st29d[x,1]) < 60*30 and abs(merged_detec[y,34] - st29d[x,0]) < 60*120) or (abs(merged_detec[y,36] - st29d[x,1]) < 60*20 and abs(merged_detec[y,34] - st29d[x,0]) < 60*120):

                    if merged_detec[y,34] < st29d[x,0] < merged_detec[y,37] or merged_detec[y,34] < st29d[x,1] < merged_detec[y,37] or  st29d[x,0] < merged_detec[y,34] < st29d[x,1] :
                        
                        merged_detec[y,0] += 1
                        merged_detec[y,29] = 1
                        
                        if merged_detec[y,34] >  st29d[x,0] :
                            merged_detec[y,34] = st29d[x,0]
                        if merged_detec[y,35] <  st29d[x,0] :
                            merged_detec[y,35] = st29d[x,0]
                        if merged_detec[y,36] >  st29d[x,1] :
                            merged_detec[y,36] = st29d[x,1]
                        if merged_detec[y,37] <  st29d[x,1] :
                            merged_detec[y,37] = st29d[x,1]                    
                        
                        added = 1
                    
        if added == 0 :
            merged_detec = np.lib.pad(merged_detec, ((0,1),(0,0)), 'constant', constant_values=(0))
            merged_detec[-1,0] = 1 
            merged_detec[-1,29] = 1    
            merged_detec[-1,34] = st29d[x,0]    
            merged_detec[-1,35] = st29d[x,0]  
            merged_detec[-1,36] = st29d[x,1] 
            merged_detec[-1,37] = st29d[x,1]  
            merged_detec[-1,38] = st29d[x,3]                   
else:
    merged_detec = np.lib.pad(merged_detec, ((0,len(st29d)),(0,0)), 'constant', constant_values=(0))
    merged_detec[:,0] = 1              
    merged_detec[:,29] = 1            
    merged_detec[:,34] = st29d[:,0]    
    merged_detec[:,35] = st29d[:,0]  
    merged_detec[:,36] = st29d[:,1] 
    merged_detec[:,37] = st29d[:,1]  
    merged_detec[:,38] = st29d[:,3]  

#30
if len(merged_detec) > 0:
    for x in range(0,len(st30d)):
        added = 0
        for y in range(0,len(merged_detec)):
            if added == 0:
                if (abs(merged_detec[y,34] - st30d[x,0]) < 60*30 and abs(merged_detec[y,37] - st30d[x,1]) < 60*120) or (abs(merged_detec[y,35] - st30d[x,0]) < 60*10 and abs(merged_detec[y,37] - st30d[x,1]) < 60*120) or (abs(merged_detec[y,37] - st30d[x,1]) < 60*30 and abs(merged_detec[y,34] - st30d[x,0]) < 60*120) or (abs(merged_detec[y,36] - st30d[x,1]) < 60*20 and abs(merged_detec[y,34] - st30d[x,0]) < 60*120):

                    if merged_detec[y,34] < st30d[x,0] < merged_detec[y,37] or merged_detec[y,34] < st30d[x,1] < merged_detec[y,37] or  st30d[x,0] < merged_detec[y,34] < st30d[x,1] :
                        
                        merged_detec[y,0] += 1
                        merged_detec[y,30] = 1
                        
                        if merged_detec[y,34] >  st30d[x,0] :
                            merged_detec[y,34] = st30d[x,0]
                        if merged_detec[y,35] <  st30d[x,0] :
                            merged_detec[y,35] = st30d[x,0]
                        if merged_detec[y,36] >  st30d[x,1] :
                            merged_detec[y,36] = st30d[x,1]
                        if merged_detec[y,37] <  st30d[x,1] :
                            merged_detec[y,37] = st30d[x,1]                  
                        
                        added = 1
                    
        if added == 0 :
            merged_detec = np.lib.pad(merged_detec, ((0,1),(0,0)), 'constant', constant_values=(0))
            merged_detec[-1,0] = 1 
            merged_detec[-1,30] = 1    
            merged_detec[-1,34] = st30d[x,0]    
            merged_detec[-1,35] = st30d[x,0]  
            merged_detec[-1,36] = st30d[x,1] 
            merged_detec[-1,37] = st30d[x,1]  
            merged_detec[-1,38] = st30d[x,3]                   
else:
    merged_detec = np.lib.pad(merged_detec, ((0,len(st30d)),(0,0)), 'constant', constant_values=(0))
    merged_detec[:,0] = 1              
    merged_detec[:,30] = 1            
    merged_detec[:,34] = st30d[:,0]    
    merged_detec[:,35] = st30d[:,0]  
    merged_detec[:,36] = st30d[:,1] 
    merged_detec[:,37] = st30d[:,1]   
    merged_detec[:,38] = st30d[:,3]  

#31
if len(merged_detec) > 0:
    for x in range(0,len(st31d)):
        added = 0
        for y in range(0,len(merged_detec)):
            if added == 0:
                if (abs(merged_detec[y,34] - st31d[x,0]) < 60*30 and abs(merged_detec[y,37] - st31d[x,1]) < 60*120) or (abs(merged_detec[y,35] - st31d[x,0]) < 60*10 and abs(merged_detec[y,37] - st31d[x,1]) < 60*120) or (abs(merged_detec[y,37] - st31d[x,1]) < 60*30 and abs(merged_detec[y,34] - st31d[x,0]) < 60*120) or (abs(merged_detec[y,36] - st31d[x,1]) < 60*20 and abs(merged_detec[y,34] - st31d[x,0]) < 60*120):

                    if merged_detec[y,34] < st31d[x,0] < merged_detec[y,37] or merged_detec[y,34] < st31d[x,1] < merged_detec[y,37] or  st31d[x,0] < merged_detec[y,34] < st31d[x,1] :
                        
                        merged_detec[y,0] += 1
                        merged_detec[y,31] = 1
                        
                        if merged_detec[y,34] >  st31d[x,0] :
                            merged_detec[y,34] = st31d[x,0]
                        if merged_detec[y,35] <  st31d[x,0] :
                            merged_detec[y,35] = st31d[x,0]
                        if merged_detec[y,36] >  st31d[x,1] :
                            merged_detec[y,36] = st31d[x,1]
                        if merged_detec[y,37] <  st31d[x,1] :
                            merged_detec[y,37] = st31d[x,1]                   
                        
                        added = 1
                    
        if added == 0 :
            merged_detec = np.lib.pad(merged_detec, ((0,1),(0,0)), 'constant', constant_values=(0))
            merged_detec[-1,0] = 1 
            merged_detec[-1,31] = 1    
            merged_detec[-1,34] = st31d[x,0]    
            merged_detec[-1,35] = st31d[x,0]  
            merged_detec[-1,36] = st31d[x,1] 
            merged_detec[-1,37] = st31d[x,1] 
            merged_detec[-1,38] = st31d[x,3]                    
else:
    merged_detec = np.lib.pad(merged_detec, ((0,len(st31d)),(0,0)), 'constant', constant_values=(0))
    merged_detec[:,0] = 1              
    merged_detec[:,31] = 1            
    merged_detec[:,34] = st31d[:,0]    
    merged_detec[:,35] = st31d[:,0]  
    merged_detec[:,36] = st31d[:,1] 
    merged_detec[:,37] = st31d[:,1]  
    merged_detec[:,38] = st31d[:,3]  

#32
if len(merged_detec) > 0:
    for x in range(0,len(st32d)):
        added = 0
        for y in range(0,len(merged_detec)):
            if added == 0:
                if (abs(merged_detec[y,34] - st32d[x,0]) < 60*30 and abs(merged_detec[y,37] - st32d[x,1]) < 60*120) or (abs(merged_detec[y,35] - st32d[x,0]) < 60*10 and abs(merged_detec[y,37] - st32d[x,1]) < 60*120) or (abs(merged_detec[y,37] - st32d[x,1]) < 60*30 and abs(merged_detec[y,34] - st32d[x,0]) < 60*120) or (abs(merged_detec[y,36] - st32d[x,1]) < 60*20 and abs(merged_detec[y,34] - st32d[x,0]) < 60*120):

                    if merged_detec[y,34] < st32d[x,0] < merged_detec[y,37] or merged_detec[y,34] < st32d[x,1] < merged_detec[y,37] or  st32d[x,0] < merged_detec[y,34] < st32d[x,1] :
                        
                        merged_detec[y,0] += 1
                        merged_detec[y,32] = 1
                        
                        if merged_detec[y,34] >  st32d[x,0] :
                            merged_detec[y,34] = st32d[x,0]
                        if merged_detec[y,35] <  st32d[x,0] :
                            merged_detec[y,35] = st32d[x,0]
                        if merged_detec[y,36] >  st32d[x,1] :
                            merged_detec[y,36] = st32d[x,1]
                        if merged_detec[y,37] <  st32d[x,1] :
                            merged_detec[y,37] = st32d[x,1]                   
                        
                        added = 1
                    
        if added == 0 :
            merged_detec = np.lib.pad(merged_detec, ((0,1),(0,0)), 'constant', constant_values=(0))
            merged_detec[-1,0] = 1 
            merged_detec[-1,32] = 1    
            merged_detec[-1,34] = st32d[x,0]    
            merged_detec[-1,35] = st32d[x,0]  
            merged_detec[-1,36] = st32d[x,1] 
            merged_detec[-1,37] = st32d[x,1]
            merged_detec[-1,38] = st32d[x,3]                   
else:
    merged_detec = np.lib.pad(merged_detec, ((0,len(st32d)),(0,0)), 'constant', constant_values=(0))
    merged_detec[:,0] = 1              
    merged_detec[:,32] = 1            
    merged_detec[:,34] = st32d[:,0]    
    merged_detec[:,35] = st32d[:,0]  
    merged_detec[:,36] = st32d[:,1] 
    merged_detec[:,37] = st32d[:,1] 
    merged_detec[:,38] = st32d[:,3]  

#33
if len(merged_detec) > 0:
    for x in range(0,len(st33d)):
        added = 0
        for y in range(0,len(merged_detec)):
            if added == 0:
                if (abs(merged_detec[y,34] - st33d[x,0]) < 60*30 and abs(merged_detec[y,37] - st33d[x,1]) < 60*120) or (abs(merged_detec[y,35] - st33d[x,0]) < 60*10 and abs(merged_detec[y,37] - st33d[x,1]) < 60*120) or (abs(merged_detec[y,37] - st33d[x,1]) < 60*30 and abs(merged_detec[y,34] - st33d[x,0]) < 60*120) or (abs(merged_detec[y,36] - st33d[x,1]) < 60*20 and abs(merged_detec[y,34] - st33d[x,0]) < 60*120):

                    if merged_detec[y,34] < st33d[x,0] < merged_detec[y,37] or merged_detec[y,34] < st33d[x,1] < merged_detec[y,37] or  st33d[x,0] < merged_detec[y,34] < st33d[x,1] :
                        
                        merged_detec[y,0] += 1
                        merged_detec[y,33] = 1
                        
                        if merged_detec[y,34] >  st33d[x,0] :
                            merged_detec[y,34] = st33d[x,0]
                        if merged_detec[y,35] <  st33d[x,0] :
                            merged_detec[y,35] = st33d[x,0]
                        if merged_detec[y,36] >  st33d[x,1] :
                            merged_detec[y,36] = st33d[x,1]
                        if merged_detec[y,37] <  st33d[x,1] :
                            merged_detec[y,37] = st33d[x,1]   
                            
                        added = 1
                    
        if added == 0 :
            merged_detec = np.lib.pad(merged_detec, ((0,1),(0,0)), 'constant', constant_values=(0))
            merged_detec[-1,0] = 1 
            merged_detec[-1,33] = 1    
            merged_detec[-1,34] = st33d[x,0]    
            merged_detec[-1,35] = st33d[x,0]  
            merged_detec[-1,36] = st33d[x,1] 
            merged_detec[-1,37] = st33d[x,1]  
            merged_detec[-1,38] = st33d[x,3]                      
else:
    merged_detec = np.lib.pad(merged_detec, ((0,len(st33d)),(0,0)), 'constant', constant_values=(0))
    merged_detec[:,0] = 1              
    merged_detec[:,33] = 1            
    merged_detec[:,34] = st33d[:,0]    
    merged_detec[:,35] = st33d[:,0]  
    merged_detec[:,36] = st33d[:,1] 
    merged_detec[:,37] = st33d[:,1]   
    merged_detec[:,38] = st33d[:,3]  

        
                
#%% reorder by time, remove non-overlapped events from merged array and any duplicates (including missed "chuncked overlapped")

#reorder by earliest starttime
merged_detec= merged_detec[merged_detec[:,34].argsort()] 


#delete 2nd event if two neighboring events have overlapping time // join station ticks 
for p in range(0,10):
    for x in range(len(merged_detec)-1,0,-1):
        if  merged_detec[x,34] < merged_detec[x-1,37] :
#            print("start time new event",UTCDateTime(merged_detec[x,34]),"previous event end",UTCDateTime(merged_detec[x-1,37]))
            for y in range(1,34):
                if merged_detec[x,y] == 1:
                    merged_detec[x-1,y] = 1
            if merged_detec[x,37] > merged_detec[x-1,37]:
                merged_detec[x-1,37] = merged_detec[x,37]
            if merged_detec[x,35] > merged_detec[x-1,35]:
                merged_detec[x-1,35] = merged_detec[x,35]
            if merged_detec[x,36] < merged_detec[x-1,36]:
                merged_detec[x-1,36] = merged_detec[x,36]
            
            merged_detec[x-1,0] = sum(merged_detec[x-1,1:34])
                
            merged_detec = np.delete(merged_detec, x, 0)
        
# delete if less than 3 detecting stations when 4 or more stations active
# or if there is only 1 station active when 3 or more are active
# if 2 or 1 station active, only one station is needed.

for x in range(len(merged_detec),0,-1):
    if merged_detec[x-1,38] > 3:
        if merged_detec[x-1,0] < 3:
            merged_detec = np.delete(merged_detec, x-1, 0)  
for x in range(len(merged_detec),0,-1):
    if merged_detec[x-1,38] == 3:
        if merged_detec[x-1,0] < 2:
            merged_detec = np.delete(merged_detec, x-1, 0)  
   
        
        



#%% print out times and save
        
#print("infrasound")
#sum_dur = 0
print(len(merged_detec), ' detections' )
#for x in range(0,len(merged_detec)):
#    print(UTCDateTime(merged_detec[x,34])," to ", UTCDateTime(merged_detec[x,37]), )
#    sum_dur += merged_detec[x,37] - merged_detec[x,34]
#
#sum_dur = sum_dur/(60*60)
#print('sum of durations (hours) = ', sum_dur)


#%%

np.savetxt("/Users/william/Documents/Fuego_catalogue/Fuego_trem_inf_v8_900_950.csv", merged_detec,delimiter=",",header="Stations,VF01,VF02,VF03,VF04,VF05,VF06,FG8A,FG8B,FG8C,FG10A,FG10B,FG10C,FG11A,FG11B,FG11C,FG11D,FG12A,FG12B,FG12C,FG13A,FG13B,FG13C,FG15A,FG15B,FG15C,FG15D,FG15E,FG15F,FV01,FV02,FV03,FV04,FV08,Early_start,Late_start,Early_end,Late_end,Available_st")

    




























