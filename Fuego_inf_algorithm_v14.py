#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 10:43:32 2020

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
from numpy import argmax

#%% Gas rich reference
# STATION, CHANNEL (DDF --> 400 Hz), NETWWORK AND LOCATION CODES 
sta = 'VF03' # STATION 
cha = 'HDF' # CHANNEL
net = 'XZ'  # 
loc = ''    # location, it depends mostly of which network you are in. 
client = Client('138.253.112.23', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23

t1r = UTCDateTime(2018, 5, 22, 8, 4, 12) #the format is year:day_of_the_year:month
t2r = t1r + 14
gr_r = Stream()
gr_r = client.get_waveforms(net, sta, loc, cha, t1r , t2r)
gr_r.detrend(type='linear')
gr_r.detrend(type='demean')
sr_r = gr_r[0].stats.sampling_rate
gr_r.filter(type='bandpass',freqmin=0.1, freqmax=(sr_r/2)-(sr_r/20))
gr_r.plot(type='relative',color='r',starttime=t1r, endtime=t2r)

#%% Seismic Reference
sta = 'FG8' # STATION 
cha = 'BHZ' # CHANNEL
net = 'GI'  # 
loc = '00'    # location, it depends mostly of which network you are in. 
client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23


t1r2 = UTCDateTime(2018, 5, 20, 2, 58, 20) #the format is year:day_of_the_year:month
t2r2 = t1r2 + 90
st_ref2 = Stream()
st_ref2 = client.get_waveforms(net, sta, loc, cha, t1r2 , t2r2)
st_ref2.detrend(type='linear')
st_ref2.detrend(type='demean')
st_ref2.filter(type='bandpass',freqmin=0.5, freqmax=5)
st_ref2.plot(color='r',starttime=t1r2, endtime=t2r2)

ref2_st = st_ref2[0].data

#%% scan limits
st_num = 41

multi_day_cat = np.zeros(shape=(0,st_num+3))

scan_start = UTCDateTime(2019,5,20,0,0,0).timestamp - 600
days = 14

for x in range(0,days):    
    Day_start = scan_start + (x*60*60*24) 
    
    st1,st2,st3,st4,st5,st6,st7,st8,st9,st10,st11,st12,st13,st14,st15,st16,st17,st18,st19,st20,st21,st22,st23,st24,st25,st26,st27,st28,st29,st30,st31,st32,st33,st34,st35,st36,st37,st38,st39,st40,st41,station_activity,seismic_activity = get_all_Fuego_stations(Day_start)    
    
    #%% loop over stations
    print("")
    print('day', x+1, 'of', days)
    print(UTCDateTime(Day_start)+600)
    print(len(station_activity)+len(seismic_activity),"stations running")
    print(len(station_activity),"Inf stations running")
    print(len(seismic_activity),"Seis stations running")
    print("")
    
    if len(station_activity) > 0:                   
        for i in range(0,len(station_activity)):
            
            st_ID = int(station_activity[i,0])
         
            if st_ID == 1:
                st_day = st1
            if st_ID == 2:
                st_day = st2
            if st_ID == 3:
                st_day = st3
            if st_ID == 4:
                st_day = st4
            if st_ID == 5:
                st_day = st5
            if st_ID == 6:
                st_day = st6
            if st_ID == 7:
                st_day = st7
            if st_ID == 8:
                st_day = st8
            if st_ID == 9:
                st_day = st9
            if st_ID == 10:
                st_day = st10
            if st_ID == 11:
                st_day = st11
            if st_ID == 12:
                st_day = st12
            if st_ID == 13:
                st_day = st13
            if st_ID == 14:
                st_day = st14
            if st_ID == 15:
                st_day = st15
            if st_ID == 16:
                st_day = st16
            if st_ID == 17:
                st_day = st17
            if st_ID == 18:
                st_day = st18
            if st_ID == 19:
                st_day = st19
            if st_ID == 20:
                st_day = st20
            if st_ID == 21:
                st_day = st21
            if st_ID == 22:
                st_day = st22
            if st_ID == 23:
                st_day = st23
            if st_ID == 24:
                st_day = st24
            if st_ID == 25:
                st_day = st25
            if st_ID == 26:
                st_day = st26
            if st_ID == 27:
                st_day = st27
            if st_ID == 28:
                st_day = st28
            if st_ID == 29:
                st_day = st29
            if st_ID == 30:
                st_day = st30
            if st_ID == 31:
                st_day = st31
            if st_ID == 32:
                st_day = st32
            if st_ID == 33:
                st_day = st33
            
    
            t1 = UTCDateTime(station_activity[i,1])
            t2 = UTCDateTime(station_activity[i,2])    
            
            sr = st_day[0].stats.sampling_rate

#            st_day.plot(type='dayplot',starttime=t1, endtime=t2)
            
            trace_t = st_day[0]
            trace_t.filter(type='bandpass',freqmin=0.5, freqmax=20)
            trace_i = st_day[0].slice(starttime = UTCDateTime(Day_start)+600, endtime = (UTCDateTime(Day_start)+24*60*60)+600)
            trace_i.filter(type='bandpass',freqmin=0.2, freqmax=(sr/2)-(sr/20))
            
            med_trace = np.median(abs(trace_t.data))
            #%% Scan for Tremor
            #window endpoints
            
#            print("smoothing")
            tr = trace_t
            tr_a = abs(tr.data)
            
            pts= int(sr * 60 * 2)
            step = 500
            
            t_step = step/sr
            
            lenst = len(tr)-pts
            num_sum = int(lenst/step)
            
            smooth = np.zeros(shape=num_sum)
            
            for u in range(0,len(smooth)):
                q = u*step
                smooth[u] = sum(abs(tr[q:q+pts]))
            
            smooth = smooth - np.percentile(smooth,10)
            smooth=smooth/pts
            
            
            
            on_trig = np.percentile(smooth,45)
            off_trig = np.percentile(smooth,25)
                        
            d_len_trig = 300
            d_trig= d_len_trig/t_step
            
            switch =0
            on_off = np.zeros(shape=(0,5))
            num=0
            
            for m in range(short,len(smooth)-short):
                
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
                            on_off = np.lib.pad(on_off, ((0,1),(0,0)), 'constant', constant_values=(0))
                                   
                            on_off[num][0]=trig_on
                            on_off[num][1]=trig_off
                            on_off[num][2]=d_len*t_step
                            
                            on_off[num][3]=UTCDateTime(t1 + trig_on*t_step).timestamp
                            on_off[num][4]=UTCDateTime(t1 + trig_off*t_step).timestamp
                            
#                            st_day.plot(color='b',starttime=t1 + trig_on*t_step, endtime=t1 + trig_off*t_step)
#                            print(UTCDateTime(t1 + trig_on*t_step),'to',UTCDateTime(t1 + trig_off*t_step), 'length =', d_len*t_step,'s')
                            num+=1
            
            plt.figure()
            plt.plot(tr,'r')
            plt.title(st_ID)
            plt.figure()
            plt.plot(smooth,'k')
            plt.plot([on_off[:,0],on_off[:,0]],[0,max(smooth)],'r')
            plt.plot([on_off[:,1],on_off[:,1]],[0,max(smooth)],'b')
            plt.title(st_ID)
            
            catalogue_t1 = np.zeros(shape=(0,3))
            catalogue_t = np.zeros(shape=(0,3))
            cat_T = np.zeros(shape=(0,3))
            catalogue_a = np.zeros(shape=(0,3))
            num=0
            num1=0
            numa = 0
            
            shift2 = 1000
            for r in range(0,len(on_off)):
                
                if  UTCDateTime(Day_start)+600 < on_off[r,3] < (UTCDateTime(Day_start)+60*60*24)+600:
                

                    t_start = UTCDateTime(on_off[r,3])
                    t_end = UTCDateTime(on_off[r,4])
                    event1 = tr.slice(starttime = t_start , endtime= t_end )
                                        
                    trs = event1
                    
                    dom,cf, bwid50 = freq_info_NqV(trs,t_start,t_end,sr)
                    med_dat = np.median(abs(event1.data))
                    
                    if dom > 1 and cf < 10:
                        
                        tr_amp = np.percentile(abs(event1.data), 50)
                        
                        if tr_amp > 1.5*med_trace:
                            
                            tr.plot(starttime =UTCDateTime(on_off[r,3]),endtime = UTCDateTime(on_off[r,4]) )
                                        
                            catalogue_t1 = np.lib.pad(catalogue_t1, ((0,1),(0,0)), 'constant', constant_values=(0))
                            catalogue_t1[num1][0]= on_off[r,3] #start time
                            catalogue_t1[num1][1]= on_off[r,4] # end time
                            catalogue_t1[num1][2]= on_off[r,2] #duration
                            num1 +=1            
            
            if len(catalogue_t1) > 0:
                catalogue_t = np.lib.pad(catalogue_t, ((0,1),(0,0)), 'constant', constant_values=(0))
                catalogue_t[0][0]= catalogue_t1[0][0]
                catalogue_t[0][1]= catalogue_t1[0][1]
                catalogue_t[0][2]= catalogue_t1[0][2]
                        
                for h in range(1,len(catalogue_t1)):
                    catalogue_t = np.lib.pad(catalogue_t, ((0,1),(0,0)), 'constant', constant_values=(0))
                    catalogue_t[h][1]= catalogue_t1[h][1]
                    
                    t_gap = catalogue_t1[h,0] - catalogue_t1[h-1,1]
                #    T_len = catalogue_t1[h-1,2] + catalogue_t1[h,2]
                    
                    if  t_gap > 120 : 
                        catalogue_t[h][0]= catalogue_t1[h][0]
                        catalogue_t[h][2]= catalogue_t[h][1] - catalogue_t[h][0] 
                    else:
                        
                        catalogue_t[h][0]= catalogue_t[h-1][0]
                        catalogue_t[h][2]= catalogue_t[h][1] - catalogue_t[h][0]  
                                   
                for g in range(0,len(catalogue_t)-1):
                    if catalogue_t[g,0] != catalogue_t[g+1,0]:
                        catalogue_a = np.lib.pad(catalogue_a, ((0,1),(0,0)), 'constant', constant_values=(0))
                        catalogue_a[numa][0]= catalogue_t1[g][0]
                        catalogue_a[numa][1]= catalogue_t1[g][2]
                        catalogue_a[numa][2]= 0
                        
                        cat_T = np.lib.pad(cat_T, ((0,1),(0,0)), 'constant', constant_values=(0))
                        cat_T[numa][0]= catalogue_t1[g][0]
                        cat_T[numa][1]= catalogue_t1[g][1]
                        cat_T[numa][2]= catalogue_t1[g][2]
                        
                        numa +=1
            
                catalogue_a = np.lib.pad(catalogue_a, ((0,1),(0,0)), 'constant', constant_values=(0))
                catalogue_a[numa][0]= catalogue_t1[-1][0]
                catalogue_a[numa][1]= catalogue_t1[-1][2]
                catalogue_a[numa][2]= 0
                
                cat_T = np.lib.pad(cat_T, ((0,1),(0,0)), 'constant', constant_values=(0))
                cat_T[numa][0]= catalogue_t1[-1][0]
                cat_T[numa][1]= catalogue_t1[-1][1]
                cat_T[numa][2]= catalogue_t1[-1][2]       
            
            else:
                do_nothing = 6
                    
            
            #%% Scan for Isolated events
            #window endpoints
            start= t1 + 600 #time window start 
            end= t2
            trs = trace_i.slice(starttime = start  , endtime= end) 
            
            nsta=int(1*sr)
            nlta=int(20*sr)
            stream=trs.data
            cft=recursive_sta_lta(stream, nsta, nlta)
            trig_on=10
            trig_off=0.1
#            plot_trigger(trs, cft, trig_on, trig_off) 
            on_off = trigger_onset(cft,trig_on,trig_off)
           
            
            #%%
            
            catalogue_i = np.zeros(shape=(0,3))
            num_i=0
            
            shift = 200
            if len(on_off) > 0:
                for b in range(0,len(on_off)):
                    elen = len(trace_i[on_off[b,0]:on_off[b,1]])/sr
                    event = trace_i[on_off[b,0]-1000:on_off[b,1]+2000]    
                    event40 = trace_i[on_off[b,0]-400:on_off[b,0]+int(20*sr)]   
                                                    
                    if elen/60 < 2 :
                        
                        e_start = start + on_off[b,0]/sr  
                        e_end = start + on_off[b,1]/sr 
                        event2 = trs.slice(starttime = e_start - 5 , endtime= e_end + 25 )        

                        et = int((start + on_off[b,0]/sr).timestamp)
                        if len(catalogue_t1) > 0:
                            near,ind=find_nearest(catalogue_a[:,0], et )
                        else:
                            near = 0
                        
                        if abs(et-near) > 30:   
                            top_gr,top,corellgr = corel(gr_r[0],event40,shift)
                                        
                            catalogue_i = np.lib.pad(catalogue_i, ((0,1),(0,0)), 'constant', constant_values=(0))
                            catalogue_i[num_i][0]= start + on_off[b,0]/sr 
                            catalogue_i[num_i][1]= elen
                        
                            
#                            plt.figure()
#                            plt.plot(event40,'b')
#                            plt.title(top_gr)
                            
                            
                            if top_gr > 0.5:
                                catalogue_i[num_i][2]= 1
                                
                            else:
                                catalogue_i[num_i][2]= 2
                                
                            num_i += 1
                            
                            catalogue_a = np.lib.pad(catalogue_a, ((0,1),(0,0)), 'constant', constant_values=(0))
                            catalogue_a[numa][0]= start + on_off[b,0]/sr 
                            catalogue_a[numa][1]= elen
                            
                            if top_gr > 0.5:
                                catalogue_a[numa][2]= 1
    #                            event2.plot(color='b')
                            else:
                                catalogue_a[numa][2]= 2
    #                            event2.plot(color='g')
                            numa+=1
            
    #%%
            if len(catalogue_a) > 0:
                catalogue_s = catalogue_a[np.argsort(catalogue_a[0:len(catalogue_a),0])]
            
            
            #%% label events
            
            # 0 = Tremor
            # 1 = Isolated Gas Rich Event
            # 2 = Isolated Ash Rich Event
            # 3 = seismic detection only - with Ground Coupled Airwave
            # 20 = Burst During Tremor
            
                event_num = len(catalogue_s)
                catalogue_d = np.zeros(shape=(event_num,4))
                catalogue_d[:,0]=catalogue_s[:,0]
                catalogue_d[:,1]=catalogue_s[:,2]
                catalogue_d[:,2]=catalogue_s[:,1]
                catalogue_d[:,3]=st_ID
                
                for f in range(0,len(catalogue_d)):
                    if catalogue_d[f,1] == 0:
                        trem_end = catalogue_d[f,0] + catalogue_d[f,2]
                        
                        for v in range(f+1,len(catalogue_s)):
                            if catalogue_d[v,0] <  trem_end :
                                if catalogue_d[v,1] > 0:
                                    catalogue_d[v,1] = catalogue_d[v,1] * 10
                                    
            
            if len(station_activity) == 1:
                Day_catalogue = np.zeros(shape=(len(catalogue_d),st_num+3))
                
                Day_catalogue[:,0] = catalogue_d[:,0]
                Day_catalogue[:,1] = catalogue_d[:,1]
                Day_catalogue[:,st_ID+1] = 1
                Day_catalogue[:,st_num+2] = catalogue_d[:,2]
                
                
            #Combined station events 
            if len(station_activity) > 1:
                if i == 0:
                    li = len(catalogue_i)
                    lt = len(cat_T)
                    
                    cat_i = np.zeros(shape=(li,st_num+3))
                    cat_t = np.zeros(shape=(lt,4))
        
                    cat_i[:,0] = catalogue_i[:,0]
                    cat_i[:,st_num+2] = catalogue_i[:,1]
                    cat_i[:,1] = catalogue_i[:,2]
                    
                    stat_ind = st_ID +1
                    cat_i[:,stat_ind]= 1
    #                print("STATION INDEX i=0 is:",stat_ind)
                    
                    
                    cat_t[:,0] = cat_T[:,0]
                    cat_t[:,1] = cat_T[:,1]
                    cat_t[:,2] = cat_T[:,2]
                    cat_t[:,3]= st_ID
                    
                else:
                     # add new isolated explosions if detected on new stations
                     for e in range(0,len(catalogue_i)):
                        ev = catalogue_i[e,0]
                        near,ind=find_nearest(cat_i[:,0], ev )
                        if abs(ev-near) > 30: 
                            cat_i = np.lib.pad(cat_i, ((0,1),(0,0)), 'constant', constant_values=(0))
                            cat_i[-1][0]=  catalogue_i[e,0]
                            cat_i[-1][st_num+2]=  catalogue_i[e,1]
                            cat_i[-1][1]=  catalogue_i[e,2]
                            stat_ind = st_ID +1
                            cat_i[-1][stat_ind]=1
                        # Add labels of detecting stations
                        else:
                            stat_ind = st_ID +1
                            cat_i[ind][stat_ind]=1
                            
                     # add all tremor from all stations - to be organised later
                     for e in range(0,len(cat_T)):
                        cat_t = np.lib.pad(cat_t, ((0,1),(0,0)), 'constant', constant_values=(0))
                        cat_t[-1][0]=  cat_T[e,0]
                        cat_t[-1][1]=  cat_T[e,1]
                        cat_t[-1][2]=  cat_T[e,2]
                        cat_t[-1][3]=  st_ID
                                        
        #%% order detections in time
                                        
                catalogue_day_iso = cat_i[np.argsort(cat_i[0:len(cat_i),0])] 
                catalogue_day_trem_list = cat_t[np.argsort(cat_t[0:len(cat_t),0])]
                
                if len(catalogue_day_iso) > 0:
                    if catalogue_day_iso[0][0] == 1:
                        catalogue_day_iso = np.delete(catalogue_day_iso, 0, 0) 
                if len(catalogue_day_trem_list) > 0: 
                    if catalogue_day_trem_list[0][0] == 1:
                        catalogue_day_trem_list = np.delete(catalogue_day_trem_list, 0, 0)
                    
        #%% remove iso events detected on a single station only
                for c in range(len(catalogue_day_iso)-1,-1,-1):
                    if sum(catalogue_day_iso[c,2:st_num+2]) == 1:
                        catalogue_day_iso = np.delete(catalogue_day_iso, c, 0)          
                    
            #%% combine tremor detections from multiple stations
        
                    
                catalogue_day_trem = np.zeros(shape=(0,st_num+3))
                
                if len(catalogue_day_trem_list)>0: 
                    # make first entry 
                    catalogue_day_trem = np.lib.pad(catalogue_day_trem, ((0,1),(0,0)), 'constant', constant_values=(0))
                    catalogue_day_trem[0][0]= catalogue_day_trem_list[0][0]
                    # add station tick 
                    st_det = catalogue_day_trem_list[0][3]
                    st_det_row = int(st_det)+1
                    catalogue_day_trem[0][st_det_row] = 1
                    # get list of durations for 1st entry
                    duration_l = np.zeros(shape=(1,1))
                    duration_l[0][0] = catalogue_day_trem_list[0][2]
                
                    # loop over all events and save new events, or add station ticks
                    for e in range(1,len(catalogue_day_trem_list)):
                        #if stating within 120 seconds - add tick
                        if catalogue_day_trem_list[e][0] - catalogue_day_trem[-1][0] < 120:
                            st_det = catalogue_day_trem_list[e,3]
                            st_det_row = int(st_det)+1
                            catalogue_day_trem[-1][st_det_row] = 1
                            #add to list of durations
                            duration_l = np.lib.pad(duration_l, ((0,1),(0,0)), 'constant', constant_values=(0))
                            duration_l[-1][0] = catalogue_day_trem_list[e][2]
                        # get duration for previous event, save new event   
                        else:
                            if len(duration_l)== 2:
                                med_dur = np.min(duration_l)
                            else:
                                med_dur = np.median(duration_l)
        
                            catalogue_day_trem[-1][st_num+2] = med_dur
                            
                            #reset duration list
                            duration_l = np.zeros(shape=(1,1))
                            duration_l[0][0] = catalogue_day_trem_list[e,2]
                            
                            catalogue_day_trem = np.lib.pad(catalogue_day_trem, ((0,1),(0,0)), 'constant', constant_values=(0))
                            catalogue_day_trem[-1][0]= catalogue_day_trem_list[e,0]
                            
                            st_det = catalogue_day_trem_list[e,3]
                            st_det_row = int(st_det)+1
                            catalogue_day_trem[-1][st_det_row] = 1
                       
                if len(catalogue_day_trem)>0:                     
                    if len(duration_l)== 2:
                        med_dur = np.min(duration_l)
                    else:
                        med_dur = np.median(duration_l)
        
                    catalogue_day_trem[-1][st_num+2] = med_dur            
                
                #only start times confirmed by 2 or more
        
                catalogue_day_trem_co = np.zeros(shape=(0,st_num+3))
                
                if len(catalogue_day_trem)>0: 
                    for o in range(0,len(catalogue_day_trem)):
                        if sum(catalogue_day_trem[o,1:st_num+1]) > 1:
                            catalogue_day_trem_co = np.lib.pad(catalogue_day_trem_co, ((0,1),(0,0)), 'constant', constant_values=(0))
                            catalogue_day_trem_co[-1][:]= catalogue_day_trem[o][:]    
                # add stations which catch end half of an event
                        if len(catalogue_day_trem_co) > 0:
                            if catalogue_day_trem_co[-1,0] < catalogue_day_trem[o,0] < catalogue_day_trem_co[-1,0] + catalogue_day_trem_co[-1,8]:
                                # add tremor ID code
                                catalogue_day_trem_co[-1,1] = 0
                                # add station ticks
        
                                if catalogue_day_trem[o,1] == 1:
                                    catalogue_day_trem_co[-1,2] = 1
                                if catalogue_day_trem[o,2] == 1:
                                    catalogue_day_trem_co[-1,3] = 1
                                if catalogue_day_trem[o,3] == 1:
                                    catalogue_day_trem_co[-1,4] = 1
                                if catalogue_day_trem[o,4] == 1:
                                    catalogue_day_trem_co[-1,5] = 1
                                if catalogue_day_trem[o,5] == 1:
                                    catalogue_day_trem_co[-1,6] = 1
                                if catalogue_day_trem[o,6] == 1:
                                    catalogue_day_trem_co[-1,7] = 1
                                if catalogue_day_trem[o,7] == 1:
                                    catalogue_day_trem_co[-1,8] = 1
                                if catalogue_day_trem[o,8] == 1:
                                    catalogue_day_trem_co[-1,9] = 1
                                if catalogue_day_trem[o,9] == 1:
                                    catalogue_day_trem_co[-1,10] = 1
                                if catalogue_day_trem[o,10] == 1:
                                    catalogue_day_trem_co[-1,11] = 1
                                if catalogue_day_trem[o,11] == 1:
                                    catalogue_day_trem_co[-1,12] = 1
                                if catalogue_day_trem[o,12] == 1:
                                    catalogue_day_trem_co[-1,13] = 1
                                if catalogue_day_trem[o,13] == 1:
                                    catalogue_day_trem_co[-1,14] = 1
                                if catalogue_day_trem[o,14] == 1:
                                    catalogue_day_trem_co[-1,15] = 1
                                if catalogue_day_trem[o,15] == 1:
                                    catalogue_day_trem_co[-1,16] = 1
                                if catalogue_day_trem[o,16] == 1:
                                    catalogue_day_trem_co[-1,17] = 1
                                if catalogue_day_trem[o,17] == 1:
                                    catalogue_day_trem_co[-1,18] = 1
                                if catalogue_day_trem[o,18] == 1:
                                    catalogue_day_trem_co[-1,19] = 1
                                if catalogue_day_trem[o,19] == 1:
                                    catalogue_day_trem_co[-1,20] = 1
                                if catalogue_day_trem[o,20] == 1:
                                    catalogue_day_trem_co[-1,21] = 1
                                if catalogue_day_trem[o,21] == 1:
                                    catalogue_day_trem_co[-1,22] = 1
                                if catalogue_day_trem[o,22] == 1:
                                    catalogue_day_trem_co[-1,23] = 1
                                if catalogue_day_trem[o,23] == 1:
                                    catalogue_day_trem_co[-1,24] = 1
                                if catalogue_day_trem[o,24] == 1:
                                    catalogue_day_trem_co[-1,25] = 1
                                if catalogue_day_trem[o,25] == 1:
                                    catalogue_day_trem_co[-1,26] = 1
                                if catalogue_day_trem[o,26] == 1:
                                    catalogue_day_trem_co[-1,27] = 1
                                if catalogue_day_trem[o,27] == 1:
                                    catalogue_day_trem_co[-1,28] = 1
                                if catalogue_day_trem[o,28] == 1:
                                    catalogue_day_trem_co[-1,29] = 1
                                if catalogue_day_trem[o,29] == 1:
                                    catalogue_day_trem_co[-1,30] = 1
                                if catalogue_day_trem[o,30] == 1:
                                    catalogue_day_trem_co[-1,31] = 1
                                if catalogue_day_trem[o,31] == 1:
                                    catalogue_day_trem_co[-1,32] = 1
                                if catalogue_day_trem[o,32] == 1:
                                    catalogue_day_trem_co[-1,33] = 1
                                if catalogue_day_trem[o,33] == 1:
                                    catalogue_day_trem_co[-1,34] = 1
                                    
                                    
                                    
                                    
                # remove joined events
                for o in range(0,len(catalogue_day_trem_co)-1): 
                    try:
                        for y in range(o+1,len(catalogue_day_trem_co)):
                            if catalogue_day_trem_co[y][0] < catalogue_day_trem_co[o][0] + catalogue_day_trem_co[o][8]:
        
        #expand for all stations                                            
                                if catalogue_day_trem_co[y,2] == 1:
                                    catalogue_day_trem_co[o,2] = 1
                                if catalogue_day_trem_co[y,3] == 1:
                                    catalogue_day_trem_co[o,3] = 1
                                if catalogue_day_trem_co[y,4] == 1:
                                    catalogue_day_trem_co[o,4] = 1
                                if catalogue_day_trem_co[y,5] == 1:
                                    catalogue_day_trem_co[o,5] = 1
                                if catalogue_day_trem_co[y,6] == 1:
                                    catalogue_day_trem_co[o,6] = 1
                                if catalogue_day_trem_co[y,7] == 1:
                                    catalogue_day_trem_co[o,7] = 1
                                if catalogue_day_trem_co[y,8] == 1:
                                    catalogue_day_trem_co[o,8] = 1
                                if catalogue_day_trem_co[y,9] == 1:
                                    catalogue_day_trem_co[o,9] = 1
                                if catalogue_day_trem_co[y,10] == 1:
                                    catalogue_day_trem_co[o,10] = 1
                                if catalogue_day_trem_co[y,11] == 1:
                                    catalogue_day_trem_co[o,11] = 1
                                if catalogue_day_trem_co[y,12] == 1:
                                    catalogue_day_trem_co[o,12] = 1
                                if catalogue_day_trem_co[y,13] == 1:
                                    catalogue_day_trem_co[o,13] = 1
                                if catalogue_day_trem_co[y,14] == 1:
                                    catalogue_day_trem_co[o,14] = 1
                                if catalogue_day_trem_co[y,15] == 1:
                                    catalogue_day_trem_co[o,15] = 1                        
                                if catalogue_day_trem_co[y,16] == 1:
                                    catalogue_day_trem_co[o,16] = 1
                                if catalogue_day_trem_co[y,17] == 1:
                                    catalogue_day_trem_co[o,17] = 1
                                if catalogue_day_trem_co[y,18] == 1:
                                    catalogue_day_trem_co[o,18] = 1
                                if catalogue_day_trem_co[y,19] == 1:
                                    catalogue_day_trem_co[o,19] = 1
                                if catalogue_day_trem_co[y,20] == 1:
                                    catalogue_day_trem_co[o,20] = 1
                                if catalogue_day_trem_co[y,21] == 1:
                                    catalogue_day_trem_co[o,21] = 1
                                if catalogue_day_trem_co[y,22] == 1:
                                    catalogue_day_trem_co[o,22] = 1
                                if catalogue_day_trem_co[y,23] == 1:
                                    catalogue_day_trem_co[o,23] = 1
                                if catalogue_day_trem_co[y,24] == 1:
                                    catalogue_day_trem_co[o,24] = 1
                                if catalogue_day_trem_co[y,25] == 1:
                                    catalogue_day_trem_co[o,25] = 1
                                if catalogue_day_trem_co[y,26] == 1:
                                    catalogue_day_trem_co[o,26] = 1
                                if catalogue_day_trem_co[y,27] == 1:
                                    catalogue_day_trem_co[o,27] = 1
                                if catalogue_day_trem_co[y,28] == 1:
                                    catalogue_day_trem_co[o,28] = 1
                                if catalogue_day_trem_co[y,29] == 1:
                                    catalogue_day_trem_co[o,29] = 1
                                if catalogue_day_trem_co[y,30] == 1:
                                    catalogue_day_trem_co[o,30] = 1
                                if catalogue_day_trem_co[y,31] == 1:
                                    catalogue_day_trem_co[o,31] = 1
                                if catalogue_day_trem_co[y,32] == 1:
                                    catalogue_day_trem_co[o,32] = 1
                                if catalogue_day_trem_co[y,33] == 1:
                                    catalogue_day_trem_co[o,33] = 1
                                if catalogue_day_trem_co[y,34] == 1:
                                    catalogue_day_trem_co[o,34] = 1
                                                                
                                    
                                catalogue_day_trem_co = np.delete(catalogue_day_trem_co, y, 0) 
                        
                    except:
                        do_nothing=0
                if len(catalogue_day_trem_co) > 0:
                    if catalogue_day_trem_co[0][0] == 1:
                        catalogue_day_trem_co = np.delete(catalogue_day_trem_co, 0, 0)
                                                
                                                        
            #%% merge tremor and isolation into one catalogue
            
                Day_catalogue_merge = catalogue_day_trem_co
                
                for a in range(0,len(catalogue_day_iso)):
                    ev = catalogue_day_iso[a,0]
                    if len(Day_catalogue_merge) > 0:
                        near,ind=find_nearest(Day_catalogue_merge[:,0], ev )
                    if abs(ev-near) < 10:
                        do_nothing = 1
                    
                    else:
                        Day_catalogue_merge = np.lib.pad(Day_catalogue_merge, ((0,1),(0,0)), 'constant', constant_values=(0))
                        
                        Day_catalogue_merge[-1,:]= catalogue_day_iso[a,:]                       
                
                if len(Day_catalogue_merge) > 0:
                    Day_catalogue = Day_catalogue_merge[np.argsort(Day_catalogue_merge[0:len(Day_catalogue_merge),0])]
        
             
                    for z in range(0,len(Day_catalogue)):        
                        if Day_catalogue[z,1] == 0:
                            trem_end = Day_catalogue[z,0] + Day_catalogue[z,st_num+2]            
                            for y in range(z+1,len(Day_catalogue)):                
                                if (Day_catalogue[y,0] <  trem_end) and (Day_catalogue[y,1] > 0):
                                        Day_catalogue[y,1] = Day_catalogue[y,1] * 10
                    
                    if Day_catalogue[0][0] == 1:
                        Day_catalogue = np.delete(Day_catalogue, 0, 0)
                       


#%% Scan for seismic GCA to add station ticks 
        if len(seismic_activity) > 0:
            for y in range(0,len(seismic_activity)):
                                
                st_ID = int(seismic_activity[y,0])
                
                if st_ID == 34:
                    st_day = st34
                if st_ID == 35:
                    st_day = st35
                if st_ID == 36:
                    st_day = st36
                if st_ID == 37:
                    st_day = st37
                if st_ID == 38:
                    st_day = st38
                if st_ID == 39:
                    st_day = st39
                if st_ID == 40:
                    st_day = st40
                if st_ID == 41:
                    st_day = st41
                    
                t1 = UTCDateTime(seismic_activity[y,1])
                t2 = UTCDateTime(seismic_activity[y,2])    
                
                sr = st_day[0].stats.sampling_rate
                
                trace_l = st_day[0]
                trace_h = st_day[0]
                
                trace_l = st_day[0].slice(starttime = UTCDateTime(Day_start)+600, endtime = (UTCDateTime(Day_start)+24*60*60)+600)
                trace_h = st_day[0].slice(starttime = UTCDateTime(Day_start)+600, endtime = (UTCDateTime(Day_start)+24*60*60)+600)
                
                trace_l.filter(type='bandpass',freqmin=0.5, freqmax=5)
                trace_h.filter(type='bandpass',freqmin=12, freqmax=22)
                
                
                nsta=int(0.5*sr)
                nlta=int(25*sr)
                stream=trace_h.data
                cft=recursive_sta_lta(stream, nsta, nlta)
                trig_on=15
                trig_off=0.5
    #            plot_trigger(trace_h, cft, trig_on, trig_off)             
                on_off = trigger_onset(cft,trig_on,trig_off)
    
                scat = np.zeros(shape=(0,1))
                snum=0
    
                shift3 = 200
                
                for p in range(0,len(on_off)):
                    
                    Sstart = UTCDateTime(trace_l.stats.starttime) + on_off[p,0]/sr 
            
                    scat = np.lib.pad(scat, ((0,1),(0,0)), 'constant', constant_values=(0))
                    scat[snum]= Sstart  
                    snum+=1
                    
                    near,ind=find_nearest(Day_catalogue[:,0], Sstart.timestamp )
                    
                    if abs(near-Sstart.timestamp) > 30:  
                        Day_catalogue[ind,st_ID+1] = 1
                        
#                        print('seismic at',UTCDateTime(Day_catalogue[ind,0]) )
                            
#                            plt.figure()
#                            plt.plot(event_h,color='r')
                        
                

#%% merge days together into one long matrix
            
        multi_day_cat = np.lib.pad(multi_day_cat, ((0,len(Day_catalogue)),(0,0)), 'constant', constant_values=(0))
        row1 = -len(Day_catalogue)
        
        multi_day_cat[row1:,:]=Day_catalogue


#%% Scan for seismic GCA to make new detectiosn when no infrasound available + add to multi day catalogue
    else:
        if len(seismic_activity) > 0:
            
            Day_catalogue = np.zeros(shape=(0,st_num+3))
            
            for y in range(0,len(seismic_activity)):
                
                st_ID = int(seismic_activity[y,0])
                
                if st_ID == 34:
                    st_day = st34
                if st_ID == 35:
                    st_day = st35
                if st_ID == 36:
                    st_day = st36
                if st_ID == 37:
                    st_day = st37
                if st_ID == 38:
                    st_day = st38
                if st_ID == 39:
                    st_day = st39
                if st_ID == 40:
                    st_day = st40
                if st_ID == 41:
                    st_day = st41
                    
                t1 = UTCDateTime(seismic_activity[y,1])
                t2 = UTCDateTime(seismic_activity[y,2])    
                
                sr = st_day[0].stats.sampling_rate
                
                trace_l = st_day[0]
                trace_h = st_day[0]
                
                trace_l = st_day[0].slice(starttime = UTCDateTime(Day_start)+600, endtime = (UTCDateTime(Day_start)+24*60*60)+600)
                trace_h = st_day[0].slice(starttime = UTCDateTime(Day_start)+600, endtime = (UTCDateTime(Day_start)+24*60*60)+600)
                
                trace_l.filter(type='bandpass',freqmin=0.5, freqmax=5)
                trace_h.filter(type='bandpass',freqmin=12, freqmax=22)
                
                
                nsta=int(0.5*sr)
                nlta=int(25*sr)
                stream=trace_h.data
                cft=recursive_sta_lta(stream, nsta, nlta)
                trig_on=15
                trig_off=0.3
    #            plot_trigger(trace_h, cft, trig_on, trig_off)             
                on_off = trigger_onset(cft,trig_on,trig_off)
    
                scat = np.zeros(shape=(0,1))
                snum=0
    
                shift3 = 200
                
                for p in range(0,len(on_off)):
                    
                    Sstart = UTCDateTime(trace_l.stats.starttime) + on_off[p,0]/sr 
                    
                    event_l=trace_l[on_off[p,0]-2000:on_off[p,1]+2000]    
                    event_h=trace_h[on_off[p,0]-2000:on_off[p,1]+2000]   
    
                    top_vl,top3,corell3 = corel(abs(ref2_st),abs(event_l),shift3)
                    
                    if  top_vl > 0.3:
                        
                        Day_catalogue = np.lib.pad(Day_catalogue, ((0,1),(0,0)), 'constant', constant_values=(0))
                        
                        Day_catalogue[snum,0] = Sstart.timestamp
                        Day_catalogue[snum,1] = 3                        
                        Day_catalogue[snum,st_ID+1] = 1
                        dur= (on_off[p,1]-on_off[p,0])/sr
                        Day_catalogue[snum,st_num+2] = dur
                        snum+=1
                        
#                        plt.figure()
#                        plt.plot(event_l,color='k')
#                        plt.plot(event_h,color='r')
#                        plt.title('duration = {}'.format(dur))
                        
                        
        
        multi_day_cat = np.lib.pad(multi_day_cat, ((0,len(Day_catalogue)),(0,0)), 'constant', constant_values=(0))
        row1 = -len(Day_catalogue)
        
        multi_day_cat[row1:,:]=Day_catalogue



#%% print out multi day catalogue
#    
print("")
print("")
print("Multi Day Catalogue:")
print("")
if len(multi_day_cat)>0:       
    for x in range(0,len(multi_day_cat)):
        
        print(UTCDateTime(multi_day_cat[x,0]),multi_day_cat[x,1],sum(multi_day_cat[x,2:st_num+2]),multi_day_cat[x,st_num+2])
##
#print("")
#print('number of detections =', len(multi_day_cat))
#print("")
#%% get all info sorted for save


#expand for all stations
        
        
Final_Catalogue = np.zeros(shape=(len(multi_day_cat),st_num+12))

Final_Catalogue[:,0] = multi_day_cat[:,0] #time

Final_Catalogue[:,7] = multi_day_cat[:,1] # type

Final_Catalogue[:,8] = multi_day_cat[:,st_num+2] #duration

#station ticks
Final_Catalogue[:,12] = multi_day_cat[:,2]
Final_Catalogue[:,13] = multi_day_cat[:,3]
Final_Catalogue[:,14] = multi_day_cat[:,4]
Final_Catalogue[:,15] = multi_day_cat[:,5]
Final_Catalogue[:,16] = multi_day_cat[:,6]
Final_Catalogue[:,17] = multi_day_cat[:,7]
Final_Catalogue[:,18] = multi_day_cat[:,8]
Final_Catalogue[:,19] = multi_day_cat[:,9]
Final_Catalogue[:,20] = multi_day_cat[:,10]
Final_Catalogue[:,21] = multi_day_cat[:,11]
Final_Catalogue[:,22] = multi_day_cat[:,12]
Final_Catalogue[:,23] = multi_day_cat[:,13]
Final_Catalogue[:,24] = multi_day_cat[:,14]
Final_Catalogue[:,25] = multi_day_cat[:,15]
Final_Catalogue[:,26] = multi_day_cat[:,16]
Final_Catalogue[:,27] = multi_day_cat[:,17]
Final_Catalogue[:,28] = multi_day_cat[:,18]
Final_Catalogue[:,29] = multi_day_cat[:,19]
Final_Catalogue[:,30] = multi_day_cat[:,20]
Final_Catalogue[:,31] = multi_day_cat[:,21]
Final_Catalogue[:,32] = multi_day_cat[:,22]
Final_Catalogue[:,33] = multi_day_cat[:,23]
Final_Catalogue[:,34] = multi_day_cat[:,24]
Final_Catalogue[:,35] = multi_day_cat[:,25]
Final_Catalogue[:,36] = multi_day_cat[:,26]
Final_Catalogue[:,37] = multi_day_cat[:,27]
Final_Catalogue[:,38] = multi_day_cat[:,28]
Final_Catalogue[:,39] = multi_day_cat[:,29]
Final_Catalogue[:,40] = multi_day_cat[:,30]
Final_Catalogue[:,41] = multi_day_cat[:,31]
Final_Catalogue[:,42] = multi_day_cat[:,32]
Final_Catalogue[:,43] = multi_day_cat[:,33]
Final_Catalogue[:,44] = multi_day_cat[:,34]

#seismic
Final_Catalogue[:,45] = multi_day_cat[:,35]
Final_Catalogue[:,46] = multi_day_cat[:,36]
Final_Catalogue[:,47] = multi_day_cat[:,37]
Final_Catalogue[:,48] = multi_day_cat[:,38]
Final_Catalogue[:,49] = multi_day_cat[:,39]
Final_Catalogue[:,50] = multi_day_cat[:,40]
Final_Catalogue[:,51] = multi_day_cat[:,41]
Final_Catalogue[:,52] = multi_day_cat[:,42]


for x in range(0,len(multi_day_cat)):
    Final_Catalogue[x,1] = UTCDateTime(multi_day_cat[x,0]).year
    Final_Catalogue[x,2] = UTCDateTime(multi_day_cat[x,0]).month
    Final_Catalogue[x,3] = UTCDateTime(multi_day_cat[x,0]).day
    Final_Catalogue[x,4] = UTCDateTime(multi_day_cat[x,0]).hour
    Final_Catalogue[x,5] = UTCDateTime(multi_day_cat[x,0]).minute
    Final_Catalogue[x,6] = UTCDateTime(multi_day_cat[x,0]).second

    Final_Catalogue[x,9] = sum(multi_day_cat[x,2:35]) # number of inf
    Final_Catalogue[x,10] = sum(multi_day_cat[x,35:43]) # number of seis
    Final_Catalogue[x,11] = sum(multi_day_cat[x,2:43]) # number of stations

if len(Final_Catalogue) > 0:
        final_ordered = Final_Catalogue[np.argsort(Final_Catalogue[0:len(Final_Catalogue),0])]    


#%% save catalogue to .csv file
    
np.savetxt("/Users/william/Documents/Fuego_catalogue/Fuego_scan_v12.csv", final_ordered,delimiter=",",header="Time_stamp,Year,Month,Day,Hour,Min,Sec,Event_type_ID,Duration [s],#inf,#seis,#Stations,VF01,VF02,VF03,VF04,VF05,VF06,FG8A,FG8B,FG8C,FG10A,FG10B,FG10C,FG11A,FG11B,FG11C,FG11D,FG12A,FG12B,FG12C,FG13A,FG13B,FG13C,FG15A,FG15B,FG15C,FG15D,FG15E,FG15F,FV01,FV02,FV03,FV04,FV08,FG3S,FG8S,FG10S,FG11S,FG12S,FG13S,FG14S,FG16S")















