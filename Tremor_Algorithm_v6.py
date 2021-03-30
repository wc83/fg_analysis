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

#seis
st34d = np.zeros(shape=(0,4))
st35d = np.zeros(shape=(0,4))
st36d = np.zeros(shape=(0,4))
st37d = np.zeros(shape=(0,4))
st38d = np.zeros(shape=(0,4))
st39d = np.zeros(shape=(0,4))
st40d = np.zeros(shape=(0,4))
st41d = np.zeros(shape=(0,4))
#%% scan limits
st_num = 41

scan_start = UTCDateTime(2020,8,30,0,0,0).timestamp 
days = 100

Station_info = np.zeros(shape=(days,5,2))


for x in range(0,days,3):    
    
    Day_start = scan_start + (x*60*60*24) 

    
    st34a,st35a,st36a,st37a,st38a,st39a,st40a,st41a,seismic_activitya = get_seis_Fuego_stations_4d(Day_start)    
   
    all_stations = seismic_activitya
    
    print("")
    print('day', x+1, 'of', days)
    print(UTCDateTime(Day_start)+600)

    print(len(seismic_activitya),"Seis stations running")
    print("")

#%% loop over stations
    if len(all_stations) > 0 :                   
        for i in range(0,len(all_stations)):
            
            st_ID = int(all_stations[i,0])

#            
            if st_ID == 34:
                tr_day = st34a[0].slice(starttime = st34a[0].stats.starttime + 3, endtime = st34a[0].stats.endtime)
                
            if st_ID == 35:
                tr_day = st35a[0].slice(starttime = st35a[0].stats.starttime + 3, endtime = st35a[0].stats.endtime)
                
            if st_ID == 36:
                tr_day = st36a[0].slice(starttime = st36a[0].stats.starttime + 3, endtime = st36a[0].stats.endtime)
                
            if st_ID == 37:
                tr_day = st37a[0].slice(starttime = st37a[0].stats.starttime + 3, endtime = st37a[0].stats.endtime)
                
            if st_ID == 38:
                tr_day = st38a[0].slice(starttime = st38a[0].stats.starttime + 3, endtime = st38a[0].stats.endtime)
                
            if st_ID == 39:
                tr_day = st39a[0].slice(starttime = st39a[0].stats.starttime + 3, endtime = st39a[0].stats.endtime)
                
            if st_ID == 40:
                tr_day = st40a[0].slice(starttime = st40a[0].stats.starttime + 3, endtime = st40a[0].stats.endtime)
                
            if st_ID == 41:
                tr_day = st41a[0].slice(starttime = st41a[0].stats.starttime + 3, endtime = st41a[0].stats.endtime)

    
            t1 = Day_start  
            t2 = Day_start + (4*24*60*60) 
            
            tr_day.detrend(type='linear')
            tr_day.detrend(type='demean')
            
            if st_ID < 34:
                tr_day.filter(type='bandpass',freqmin=0.5, freqmax=6) 
            else:
                tr_day.filter(type='bandpass',freqmin=0.5, freqmax=6) 
            
#            print(st_ID, tr_day)
            
            sr = tr_day.stats.sampling_rate            
            
            med_trace = np.median(abs(tr_day.data))
            
            tr = tr_day
            tr_a = abs(tr.data)     
            
            pts= int(sr * 60 * 2)
            
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
                on_trig = np.percentile(smooth,75) 
                off_trig = np.percentile(smooth,25) # 25
            else:
                on_trig = np.percentile(smooth,65) #50
                off_trig = np.percentile(smooth,25) #10
                            
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
                            
                            if len(tr_prev) > (sr*8*60 ): #not after a break in data
                                
                           
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
                            
                                mps=int(len(famp)/(sr/1.5))
                                mpe=int(len(famp)/(sr/2.5))
                                
                                lps=int(len(famp)/(sr/2.5))
                                lpe=int(len(famp)/(sr/3.5))
                                #
                                f_peak_m = np.mean(famp[fps:fpe])/(window*sr)
                                m_peak_m = np.mean(famp[mps:mpe])/(window*sr)
                                l_peak_m = np.mean(famp[lps:lpe])/(window*sr)
                                peak_r = m_peak_m/((f_peak_m+l_peak_m)/2)
                                
                                if peak_r > 2.2 : 
                                    
                                    on_off = np.lib.pad(on_off, ((0,1),(0,0)), 'constant', constant_values=(0))
                                       
                                    on_off[num][0]=trig_on
                                    on_off[num][1]=trig_off
                                    on_off[num][2]=d_len*t_step
                                    
                                    on_off[num][3]=UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
                                    on_off[num][4]=UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
                                    
                                    num+=1

                                  
                                        
                                    if st_ID == 34:
                                        st34d = np.lib.pad(st34d, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        st34d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
                                        st34d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
                                        st34d[-1][2] = d_len*t_step
                                        st34d[-1][3] = len(seismic_activitya)
                                        
                                    if st_ID == 35:
                                        st35d = np.lib.pad(st35d, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        st35d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
                                        st35d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
                                        st35d[-1][2] = d_len*t_step
                                        st35d[-1][3] = len(seismic_activitya)
                                        
                                    if st_ID == 36:
                                        st36d = np.lib.pad(st36d, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        st36d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
                                        st36d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
                                        st36d[-1][2] = d_len*t_step
                                        st36d[-1][3] = len(seismic_activitya)
                                        
                                    if st_ID == 37:
                                        st37d = np.lib.pad(st37d, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        st37d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
                                        st37d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
                                        st37d[-1][2] = d_len*t_step
                                        st37d[-1][3] = len(seismic_activitya)
                                        
                                    if st_ID == 38:
                                        st38d = np.lib.pad(st38d, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        st38d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
                                        st38d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
                                        st38d[-1][2] = d_len*t_step
                                        st38d[-1][3] = len(seismic_activitya)
                                        
                                    if st_ID == 39:
                                        st39d = np.lib.pad(st39d, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        st39d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
                                        st39d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
                                        st39d[-1][2] = d_len*t_step
                                        st39d[-1][3] = len(seismic_activitya)
                                    
                                    if st_ID == 40:
                                        st40d = np.lib.pad(st40d, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        st40d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
                                        st40d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
                                        st40d[-1][2] = d_len*t_step
                                        st40d[-1][3] = len(seismic_activitya)
                                        
                                    if st_ID == 41:
                                        st41d = np.lib.pad(st41d, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        st41d[-1][0] = UTCDateTime(tr_day.stats.starttime + trig_on*t_step).timestamp
                                        st41d[-1][1] = UTCDateTime(tr_day.stats.starttime + trig_off*t_step).timestamp
                                        st41d[-1][2] = d_len*t_step
                                        st41d[-1][3] = len(seismic_activitya)
                                        
                                    

#                                
#                            print(UTCDateTime(t1 + trig_on*t_step),'to',UTCDateTime(t1 + trig_off*t_step), 'length =', d_len*t_step,'s')
#                            print('dom, cf, bwid50 = ',dom,cf, bwid50)
#                            print('mean, med, var = ',mean_abs_amp,med_abs_amp,var_abs_amp)
#                  #plot 
#            if len(smooth) > 0 and len(on_off) > 0 and st_ID > 33:
#                plt.figure()
#                plt.plot(tr,'b')
#                plt.title(st_ID)
#                plt.figure()
#                plt.plot(smooth,'k')
#                plt.plot([on_off[:,0],on_off[:,0]],[0,max(smooth)],'r')
#                plt.plot([on_off[:,1],on_off[:,1]],[0,max(smooth)],'b')

                
            
          
            

    
#%% reorder and remove duplicated events from individual station arrays
      
#reorder          


st34d = st34d[st34d[:,0].argsort()]
st35d = st35d[st35d[:,0].argsort()]
st36d = st36d[st36d[:,0].argsort()]
st37d = st37d[st37d[:,0].argsort()]
st38d = st38d[st38d[:,0].argsort()]
st39d = st39d[st39d[:,0].argsort()]
st40d = st40d[st40d[:,0].argsort()] 
st41d = st41d[st41d[:,0].argsort()]

#remove duplicates or zero values


for x in range(len(st34d)-1,0,-1):
    if st34d[x,0]-st34d[x-1,0] < 300:
        if st34d[x,2] > st34d[x-1,2]:
            st34d = np.delete(st34d, x-1, 0)
        else:
            st34d = np.delete(st34d, x, 0)
for x in range(len(st34d)-1,0,-1):
    if st34d[x,0]== 0 or st34d[x,1]== 0 :
        st34d = np.delete(st34d, x, 0)    
                        
for x in range(len(st35d)-1,0,-1):
    if st35d[x,0]-st35d[x-1,0] < 300:
        if st35d[x,2] > st35d[x-1,2]:
            st35d = np.delete(st35d, x-1, 0)
        else:
            st35d = np.delete(st35d, x, 0)
for x in range(len(st35d)-1,0,-1):
    if st35d[x,0]== 0 or st35d[x,1]== 0 :
        st35d = np.delete(st35d, x, 0)    
                        
for x in range(len(st36d)-1,0,-1):
    if st36d[x,0]-st36d[x-1,0] < 300:
        if st36d[x,2] > st36d[x-1,2]:
            st36d = np.delete(st36d, x-1, 0)
        else:
            st36d = np.delete(st36d, x, 0)
for x in range(len(st36d)-1,0,-1):
    if st36d[x,0]== 0 or st36d[x,1]== 0 :
        st36d = np.delete(st36d, x, 0)    
            
for x in range(len(st37d)-1,0,-1):
    if st37d[x,0]-st37d[x-1,0] < 300:
        if st37d[x,2] > st37d[x-1,2]:
            st37d = np.delete(st37d, x-1, 0)
        else:
            st37d = np.delete(st37d, x, 0)
for x in range(len(st37d)-1,0,-1):
    if st37d[x,0]== 0 or st37d[x,1]== 0 :
        st37d = np.delete(st37d, x, 0)    
            
for x in range(len(st38d)-1,0,-1):
    if st38d[x,0]-st38d[x-1,0] < 300:
        if st38d[x,2] > st38d[x-1,2]:
            st38d = np.delete(st38d, x-1, 0)
        else:
            st38d = np.delete(st38d, x, 0)
for x in range(len(st38d)-1,0,-1):
    if st38d[x,0]== 0 or st38d[x,1]== 0 :
        st38d = np.delete(st38d, x, 0)    
            
for x in range(len(st39d)-1,0,-1):
    if st39d[x,0]-st39d[x-1,0] < 300:
        if st39d[x,2] > st39d[x-1,2]:
            st39d = np.delete(st39d, x-1, 0)
        else:
            st39d = np.delete(st39d, x, 0)
for x in range(len(st39d)-1,0,-1):
    if st39d[x,0]== 0 or st39d[x,1]== 0 :
        st39d = np.delete(st39d, x, 0)    
            
for x in range(len(st40d)-1,0,-1):
    if st40d[x,0]-st40d[x-1,0] < 300:
        if st40d[x,2] > st40d[x-1,2]:
            st40d = np.delete(st40d, x-1, 0)
        else:
            st40d = np.delete(st40d, x, 0)
for x in range(len(st40d)-1,0,-1):
    if st40d[x,0]== 0 or st40d[x,1]== 0 :
        st40d = np.delete(st40d, x, 0)    
            
for x in range(len(st41d)-1,0,-1):
    if st41d[x,0]-st41d[x-1,0] < 300:
        if st41d[x,2] > st41d[x-1,2]:
            st41d = np.delete(st41d, x-1, 0)
        else:
            st41d = np.delete(st41d, x, 0)
for x in range(len(st41d)-1,0,-1):
    if st41d[x,0]== 0 or st41d[x,1]== 0 :
        st41d = np.delete(st41d, x, 0)               

             


        

#%% make merged array and populate - Seismics

#index 0 = number of stations, then all station ticks, index 34 = earliest start, 35 = latest start, 36 = earliest end, 37 = latest end
merged_detec_s = np.zeros(shape=(len(st34d),14))

merged_detec_s[:,0] = 1              
merged_detec_s[:,1] = 1            
merged_detec_s[:,9] = st34d[:,0]    
merged_detec_s[:,10] = st34d[:,0]  
merged_detec_s[:,11] = st34d[:,1] 
merged_detec_s[:,12] = st34d[:,1] 
merged_detec_s[:,13] = st34d[:,3]        

    
#35
if len(merged_detec_s) > 0:
    for x in range(0,len(st35d)):
        added = 0
        for y in range(0,len(merged_detec_s)):
            if added == 0:
                if (abs(merged_detec_s[y,9] - st35d[x,0]) < 60*30 and abs(merged_detec_s[y,12] - st35d[x,1]) < 60*120) or (abs(merged_detec_s[y,10] - st35d[x,0]) < 60*10 and abs(merged_detec_s[y,12] - st35d[x,1]) < 60*120) or (abs(merged_detec_s[y,12] - st35d[x,1]) < 60*30 and abs(merged_detec_s[y,9] - st35d[x,0]) < 60*120) or (abs(merged_detec_s[y,11] - st35d[x,1]) < 60*20 and abs(merged_detec_s[y,9] - st35d[x,0]) < 60*120):

                    if merged_detec_s[y,9] < st35d[x,0] < merged_detec_s[y,12] or merged_detec_s[y,9] < st35d[x,1] < merged_detec_s[y,12] or  st35d[x,0] < merged_detec_s[y,9] < st35d[x,1] :
                        
                        merged_detec_s[y,0] += 1
                        merged_detec_s[y,2] = 1
                        
                        if merged_detec_s[y,9] >  st35d[x,0] :
                            merged_detec_s[y,9] = st35d[x,0]
                        if merged_detec_s[y,10] <  st35d[x,0] :
                            merged_detec_s[y,10] = st35d[x,0]
                        if merged_detec_s[y,11] >  st35d[x,1] :
                            merged_detec_s[y,11] = st35d[x,1]
                        if merged_detec_s[y,12] <  st35d[x,1] :
                            merged_detec_s[y,12] = st35d[x,1]                  
                        
                        added = 1
                        
        if added == 0 :
            merged_detec_s = np.lib.pad(merged_detec_s, ((0,1),(0,0)), 'constant', constant_values=(0))
            merged_detec_s[-1,0] = 1 
            merged_detec_s[-1,2] = 1    
            merged_detec_s[-1,9] = st35d[x,0]    
            merged_detec_s[-1,10] = st35d[x,0]  
            merged_detec_s[-1,11] = st35d[x,1] 
            merged_detec_s[-1,12] = st35d[x,1] 
            merged_detec_s[-1,13] = st35d[x,3]                       
else:
    merged_detec_s = np.lib.pad(merged_detec_s, ((0,len(st35d)),(0,0)), 'constant', constant_values=(0))
    merged_detec_s[:,0] = 1              
    merged_detec_s[:,2] = 1            
    merged_detec_s[:,9] = st35d[:,0]    
    merged_detec_s[:,10] = st35d[:,0]  
    merged_detec_s[:,11] = st35d[:,1] 
    merged_detec_s[:,12] = st35d[:,1] 
    merged_detec_s[:,13] = st35d[:,3] 

#36
if len(merged_detec_s) > 0:
    for x in range(0,len(st36d)):
        added = 0
        for y in range(0,len(merged_detec_s)):
            if added == 0:
                if (abs(merged_detec_s[y,9] - st36d[x,0]) < 60*30 and abs(merged_detec_s[y,12] - st36d[x,1]) < 60*120) or (abs(merged_detec_s[y,10] - st36d[x,0]) < 60*10 and abs(merged_detec_s[y,12] - st36d[x,1]) < 60*120) or (abs(merged_detec_s[y,12] - st36d[x,1]) < 60*30 and abs(merged_detec_s[y,9] - st36d[x,0]) < 60*120) or (abs(merged_detec_s[y,11] - st36d[x,1]) < 60*20 and abs(merged_detec_s[y,9] - st36d[x,0]) < 60*120):

                    if merged_detec_s[y,9] < st36d[x,0] < merged_detec_s[y,12] or merged_detec_s[y,9] < st36d[x,1] < merged_detec_s[y,12] or  st36d[x,0] < merged_detec_s[y,9] < st36d[x,1] :
                        
                        merged_detec_s[y,0] += 1
                        merged_detec_s[y,3] = 1
                        
                        if merged_detec_s[y,9] >  st36d[x,0] :
                            merged_detec_s[y,9] = st36d[x,0]
                        if merged_detec_s[y,10] <  st36d[x,0] :
                            merged_detec_s[y,10] = st36d[x,0]
                        if merged_detec_s[y,11] >  st36d[x,1] :
                            merged_detec_s[y,11] = st36d[x,1]
                        if merged_detec_s[y,12] <  st36d[x,1] :
                            merged_detec_s[y,12] = st36d[x,1]                  
                        
                        added = 1
                        
        if added == 0 :
            merged_detec_s = np.lib.pad(merged_detec_s, ((0,1),(0,0)), 'constant', constant_values=(0))
            merged_detec_s[-1,0] = 1 
            merged_detec_s[-1,3] = 1    
            merged_detec_s[-1,9] = st36d[x,0]    
            merged_detec_s[-1,10] = st36d[x,0]  
            merged_detec_s[-1,11] = st36d[x,1] 
            merged_detec_s[-1,12] = st36d[x,1]   
            merged_detec_s[-1,13] = st36d[x,3]                     
else:
    merged_detec_s = np.lib.pad(merged_detec_s, ((0,len(st36d)),(0,0)), 'constant', constant_values=(0))
    merged_detec_s[:,0] = 1              
    merged_detec_s[:,3] = 1            
    merged_detec_s[:,9] = st36d[:,0]    
    merged_detec_s[:,10] = st36d[:,0]  
    merged_detec_s[:,11] = st36d[:,1] 
    merged_detec_s[:,12] = st36d[:,1] 
    merged_detec_s[:,13] = st36d[:,3] 

#37
if len(merged_detec_s) > 0:
    for x in range(0,len(st37d)):
        added = 0
        for y in range(0,len(merged_detec_s)):
            if added == 0:
                if (abs(merged_detec_s[y,9] - st37d[x,0]) < 60*30 and abs(merged_detec_s[y,12] - st37d[x,1]) < 60*120) or (abs(merged_detec_s[y,10] - st37d[x,0]) < 60*10 and abs(merged_detec_s[y,12] - st37d[x,1]) < 60*120) or (abs(merged_detec_s[y,12] - st37d[x,1]) < 60*30 and abs(merged_detec_s[y,9] - st37d[x,0]) < 60*120) or (abs(merged_detec_s[y,11] - st37d[x,1]) < 60*20 and abs(merged_detec_s[y,9] - st37d[x,0]) < 60*120):

                    if merged_detec_s[y,9] < st37d[x,0] < merged_detec_s[y,12] or merged_detec_s[y,9] < st37d[x,1] < merged_detec_s[y,12] or  st37d[x,0] < merged_detec_s[y,9] < st37d[x,1] :
                        
                        merged_detec_s[y,0] += 1
                        merged_detec_s[y,4] = 1
                        
                        if merged_detec_s[y,9] >  st37d[x,0] :
                            merged_detec_s[y,9] = st37d[x,0]
                        if merged_detec_s[y,10] <  st37d[x,0] :
                            merged_detec_s[y,10] = st37d[x,0]
                        if merged_detec_s[y,11] >  st37d[x,1] :
                            merged_detec_s[y,11] = st37d[x,1]
                        if merged_detec_s[y,12] <  st37d[x,1] :
                            merged_detec_s[y,12] = st37d[x,1]                  
                        
                        added = 1
                        
        if added == 0 :
            merged_detec_s = np.lib.pad(merged_detec_s, ((0,1),(0,0)), 'constant', constant_values=(0))
            merged_detec_s[-1,0] = 1 
            merged_detec_s[-1,4] = 1    
            merged_detec_s[-1,9] = st37d[x,0]    
            merged_detec_s[-1,10] = st37d[x,0]  
            merged_detec_s[-1,11] = st37d[x,1] 
            merged_detec_s[-1,12] = st37d[x,1] 
            merged_detec_s[-1,13] = st37d[x,3]                       
else:
    merged_detec_s = np.lib.pad(merged_detec_s, ((0,len(st37d)),(0,0)), 'constant', constant_values=(0))
    merged_detec_s[:,0] = 1              
    merged_detec_s[:,4] = 1            
    merged_detec_s[:,9] = st37d[:,0]    
    merged_detec_s[:,10] = st37d[:,0]  
    merged_detec_s[:,11] = st37d[:,1] 
    merged_detec_s[:,12] = st37d[:,1]
    merged_detec_s[:,13] = st37d[:,3] 

#38
if len(merged_detec_s) > 0:
    for x in range(0,len(st38d)):
        added = 0
        for y in range(0,len(merged_detec_s)):
            if added == 0:
                if (abs(merged_detec_s[y,9] - st38d[x,0]) < 60*30 and abs(merged_detec_s[y,12] - st38d[x,1]) < 60*120) or (abs(merged_detec_s[y,10] - st38d[x,0]) < 60*10 and abs(merged_detec_s[y,12] - st38d[x,1]) < 60*120) or (abs(merged_detec_s[y,12] - st38d[x,1]) < 60*30 and abs(merged_detec_s[y,9] - st38d[x,0]) < 60*120) or (abs(merged_detec_s[y,11] - st38d[x,1]) < 60*20 and abs(merged_detec_s[y,9] - st38d[x,0]) < 60*120):

                    if merged_detec_s[y,9] < st38d[x,0] < merged_detec_s[y,12] or merged_detec_s[y,9] < st38d[x,1] < merged_detec_s[y,12] or  st38d[x,0] < merged_detec_s[y,9] < st38d[x,1] :
                        
                        merged_detec_s[y,0] += 1
                        merged_detec_s[y,5] = 1
                        
                        if merged_detec_s[y,9] >  st38d[x,0] :
                            merged_detec_s[y,9] = st38d[x,0]
                        if merged_detec_s[y,10] <  st38d[x,0] :
                            merged_detec_s[y,10] = st38d[x,0]
                        if merged_detec_s[y,11] >  st38d[x,1] :
                            merged_detec_s[y,11] = st38d[x,1]
                        if merged_detec_s[y,12] <  st38d[x,1] :
                            merged_detec_s[y,12] = st38d[x,1]                  
                        
                        added = 1
                        
        if added == 0 :
            merged_detec_s = np.lib.pad(merged_detec_s, ((0,1),(0,0)), 'constant', constant_values=(0))
            merged_detec_s[-1,0] = 1 
            merged_detec_s[-1,5] = 1    
            merged_detec_s[-1,9] = st38d[x,0]    
            merged_detec_s[-1,10] = st38d[x,0]  
            merged_detec_s[-1,11] = st38d[x,1] 
            merged_detec_s[-1,12] = st38d[x,1]
            merged_detec_s[-1,13] = st38d[x,3]                        
else:
    merged_detec_s = np.lib.pad(merged_detec_s, ((0,len(st38d)),(0,0)), 'constant', constant_values=(0))
    merged_detec_s[:,0] = 1              
    merged_detec_s[:,5] = 1            
    merged_detec_s[:,9] = st38d[:,0]    
    merged_detec_s[:,10] = st38d[:,0]  
    merged_detec_s[:,11] = st38d[:,1] 
    merged_detec_s[:,12] = st38d[:,1] 
    merged_detec_s[:,13] = st38d[:,3] 

#39
if len(merged_detec_s) > 0:
    for x in range(0,len(st39d)):
        added = 0
        for y in range(0,len(merged_detec_s)):
            if added == 0:
                if (abs(merged_detec_s[y,9] - st39d[x,0]) < 60*30 and abs(merged_detec_s[y,12] - st39d[x,1]) < 60*120) or (abs(merged_detec_s[y,10] - st39d[x,0]) < 60*10 and abs(merged_detec_s[y,12] - st39d[x,1]) < 60*120) or (abs(merged_detec_s[y,12] - st39d[x,1]) < 60*30 and abs(merged_detec_s[y,9] - st39d[x,0]) < 60*120) or (abs(merged_detec_s[y,11] - st39d[x,1]) < 60*20 and abs(merged_detec_s[y,9] - st39d[x,0]) < 60*120):

                    if merged_detec_s[y,9] < st39d[x,0] < merged_detec_s[y,12] or merged_detec_s[y,9] < st39d[x,1] < merged_detec_s[y,12] or  st39d[x,0] < merged_detec_s[y,9] < st39d[x,1] :
                        
                        merged_detec_s[y,0] += 1
                        merged_detec_s[y,6] = 1
                        
                        if merged_detec_s[y,9] >  st39d[x,0] :
                            merged_detec_s[y,9] = st39d[x,0]
                        if merged_detec_s[y,10] <  st39d[x,0] :
                            merged_detec_s[y,10] = st39d[x,0]
                        if merged_detec_s[y,11] >  st39d[x,1] :
                            merged_detec_s[y,11] = st39d[x,1]
                        if merged_detec_s[y,12] <  st39d[x,1] :
                            merged_detec_s[y,12] = st39d[x,1]                  
                        
                        added = 1
                        
        if added == 0 :
            merged_detec_s = np.lib.pad(merged_detec_s, ((0,1),(0,0)), 'constant', constant_values=(0))
            merged_detec_s[-1,0] = 1 
            merged_detec_s[-1,6] = 1    
            merged_detec_s[-1,9] = st39d[x,0]    
            merged_detec_s[-1,10] = st39d[x,0]  
            merged_detec_s[-1,11] = st39d[x,1] 
            merged_detec_s[-1,12] = st39d[x,1] 
            merged_detec_s[-1,13] = st39d[x,3]                   
else:
    merged_detec_s = np.lib.pad(merged_detec_s, ((0,len(st39d)),(0,0)), 'constant', constant_values=(0))
    merged_detec_s[:,0] = 1              
    merged_detec_s[:,6] = 1            
    merged_detec_s[:,9] = st39d[:,0]    
    merged_detec_s[:,10] = st39d[:,0]  
    merged_detec_s[:,11] = st39d[:,1] 
    merged_detec_s[:,12] = st39d[:,1] 
    merged_detec_s[:,13] = st39d[:,3]    

#40
if len(merged_detec_s) > 0:
    for x in range(0,len(st40d)):
        added = 0
        for y in range(0,len(merged_detec_s)):
            if added == 0:
                if (abs(merged_detec_s[y,9] - st40d[x,0]) < 60*30 and abs(merged_detec_s[y,12] - st40d[x,1]) < 60*120) or (abs(merged_detec_s[y,10] - st40d[x,0]) < 60*10 and abs(merged_detec_s[y,12] - st40d[x,1]) < 60*120) or (abs(merged_detec_s[y,12] - st40d[x,1]) < 60*30 and abs(merged_detec_s[y,9] - st40d[x,0]) < 60*120) or (abs(merged_detec_s[y,11] - st40d[x,1]) < 60*20 and abs(merged_detec_s[y,9] - st40d[x,0]) < 60*120):

                    if merged_detec_s[y,9] < st40d[x,0] < merged_detec_s[y,12] or merged_detec_s[y,9] < st40d[x,1] < merged_detec_s[y,12] or  st40d[x,0] < merged_detec_s[y,9] < st40d[x,1] :
                        
                        merged_detec_s[y,0] += 1
                        merged_detec_s[y,7] = 1
                        
                        if merged_detec_s[y,9] >  st40d[x,0] :
                            merged_detec_s[y,9] = st40d[x,0]
                        if merged_detec_s[y,10] <  st40d[x,0] :
                            merged_detec_s[y,10] = st40d[x,0]
                        if merged_detec_s[y,11] >  st40d[x,1] :
                            merged_detec_s[y,11] = st40d[x,1]
                        if merged_detec_s[y,12] <  st40d[x,1] :
                            merged_detec_s[y,12] = st40d[x,1]                  
                        
                        added = 1
                        
        if added == 0 :
            merged_detec_s = np.lib.pad(merged_detec_s, ((0,1),(0,0)), 'constant', constant_values=(0))
            merged_detec_s[-1,0] = 1 
            merged_detec_s[-1,7] = 1    
            merged_detec_s[-1,9] = st40d[x,0]    
            merged_detec_s[-1,10] = st40d[x,0]  
            merged_detec_s[-1,11] = st40d[x,1] 
            merged_detec_s[-1,12] = st40d[x,1]
            merged_detec_s[-1,13] = st40d[x,3]                    
else:
    merged_detec_s = np.lib.pad(merged_detec_s, ((0,len(st40d)),(0,0)), 'constant', constant_values=(0))
    merged_detec_s[:,0] = 1              
    merged_detec_s[:,7] = 1            
    merged_detec_s[:,9] = st40d[:,0]    
    merged_detec_s[:,10] = st40d[:,0]  
    merged_detec_s[:,11] = st40d[:,1] 
    merged_detec_s[:,12] = st40d[:,1] 
    merged_detec_s[:,13] = st40d[:,3] 

#41
if len(merged_detec_s) > 0:
    for x in range(0,len(st41d)):
        added = 0
        for y in range(0,len(merged_detec_s)):
            if added == 0:
                if (abs(merged_detec_s[y,9] - st41d[x,0]) < 60*30 and abs(merged_detec_s[y,12] - st41d[x,1]) < 60*120) or (abs(merged_detec_s[y,10] - st41d[x,0]) < 60*10 and abs(merged_detec_s[y,12] - st41d[x,1]) < 60*120) or (abs(merged_detec_s[y,12] - st41d[x,1]) < 60*30 and abs(merged_detec_s[y,9] - st41d[x,0]) < 60*120) or (abs(merged_detec_s[y,11] - st41d[x,1]) < 60*20 and abs(merged_detec_s[y,9] - st41d[x,0]) < 60*120):

                    if merged_detec_s[y,9] < st41d[x,0] < merged_detec_s[y,12] or merged_detec_s[y,9] < st41d[x,1] < merged_detec_s[y,12] or  st41d[x,0] < merged_detec_s[y,9] < st41d[x,1] :
                        
                        merged_detec_s[y,0] += 1
                        merged_detec_s[y,8] = 1
                        
                        if merged_detec_s[y,9] >  st41d[x,0] :
                            merged_detec_s[y,9] = st41d[x,0]
                        if merged_detec_s[y,10] <  st41d[x,0] :
                            merged_detec_s[y,10] = st41d[x,0]
                        if merged_detec_s[y,11] >  st41d[x,1] :
                            merged_detec_s[y,11] = st41d[x,1]
                        if merged_detec_s[y,12] <  st41d[x,1] :
                            merged_detec_s[y,12] = st41d[x,1]                  
                        
                        added = 1
                        
        if added == 0 :
            merged_detec_s = np.lib.pad(merged_detec_s, ((0,1),(0,0)), 'constant', constant_values=(0))
            merged_detec_s[-1,0] = 1 
            merged_detec_s[-1,8] = 1    
            merged_detec_s[-1,9] = st41d[x,0]    
            merged_detec_s[-1,10] = st41d[x,0]  
            merged_detec_s[-1,11] = st41d[x,1] 
            merged_detec_s[-1,12] = st41d[x,1] 
            merged_detec_s[-1,13] = st41d[x,3]                   
else:
    merged_detec_s = np.lib.pad(merged_detec_s, ((0,len(st41d)),(0,0)), 'constant', constant_values=(0))
    merged_detec_s[:,0] = 1              
    merged_detec_s[:,8] = 1            
    merged_detec_s[:,9] = st41d[:,0]    
    merged_detec_s[:,10] = st41d[:,0]  
    merged_detec_s[:,11] = st41d[:,1] 
    merged_detec_s[:,12] = st41d[:,1] 
    merged_detec_s[:,13] = st41d[:,3] 

    
#%% reorder by time, remove non-overlapped events from merged array and any duplicates (including missed "chuncked overlapped")

#reorder by earliest starttime
merged_detec_s= merged_detec_s[merged_detec_s[:,9].argsort()] 


#delete 2nd event if two neighboring events have overlapping time // join station ticks 
for p in range(0,10):
    for x in range(len(merged_detec_s)-1,0,-1):
        if  merged_detec_s[x,9] < merged_detec_s[x-1,12] :
#            print("start time new event",UTCDateTime(merged_detec[x,34]),"previous event end",UTCDateTime(merged_detec[x-1,37]))
            for y in range(1,9):
                if merged_detec_s[x,y] == 1:
                    merged_detec_s[x-1,y] = 1
            if merged_detec_s[x,12] > merged_detec_s[x-1,12]:
                merged_detec_s[x-1,12] = merged_detec_s[x,12]
            if merged_detec_s[x,10] > merged_detec_s[x-1,10]:
                merged_detec_s[x-1,10] = merged_detec_s[x,10]
            if merged_detec_s[x,11] < merged_detec_s[x-1,11]:
                merged_detec_s[x-1,11] = merged_detec_s[x,11]
            
            merged_detec_s[x-1,0] = sum(merged_detec_s[x-1,1:9])
                
            merged_detec_s = np.delete(merged_detec_s, x, 0)
        
# delete if less than 2 (i.e. 1) detecting stations when 3 stations or more are active.
            # if 2 or 1 station active, only 1 detection is needed

for x in range(len(merged_detec_s),0,-1):
    if merged_detec_s[x-1,13] > 2:
        if merged_detec_s[x-1,0] < 2:
            merged_detec_s = np.delete(merged_detec_s, x-1, 0)       
        



#%% print out times and save
cat = merged_detec_s
for x in range(len(cat),1,-1):
    if cat[x-1,9] - cat[x-2,12] < 5*60:
        for y in range(1,9):
            if cat[x-1,y] == 1:
                cat[x-2,y] = 1
        cat[x-2,0] = sum(cat[x-2,1:9])
        cat[x-2,12] = cat[x-1,12]
 
        cat = np.delete(cat, x-1, 0)     

sum_dur = 0        
for x in range(0,len(cat)):
    sum_dur += cat[0,12]-cat[0,9]
    


print("")
print(len(cat), 'detections')
print('total duration =', sum_dur/(60*60))

#for x in range(0,len(merged_detec_s)):
#    print(UTCDateTime(merged_detec_s[x,9])," to ", UTCDateTime(merged_detec_s[x,12]), )


#%%

np.savetxt("/Users/william/Documents/Fuego_catalogue/Fuego_seis_trem_v8_900_1000.csv", cat,delimiter=",",header="Stations,FG3S,FG8S,FG10S,FG11S,FG12S,FG13S,FG14S,FG16S,Early_start,Late_start,Early_end,Late_end,Available_st")





























