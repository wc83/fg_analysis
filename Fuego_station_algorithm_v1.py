#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 12:01:42 2020

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


# STATION, CHANNEL (DDF --> 400 Hz), NETWWORK AND LOCATION CODES 
sta = 'FG8' # STATION 
cha = 'BHZ' # CHANNEL
net = 'GI'  # 
loc = '00'    # location, it depends mostly of which network you are in. 
client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23

#%% reference waveforms

t1 = UTCDateTime(2018, 5, 20, 17, 13, 45) #the format is year:day_of_the_year:month
t2 = t1 + 40
st_ref = Stream()
st_ref = client.get_waveforms(net, sta, loc, cha, t1 , t2)
st_ref.detrend(type='linear')
st_ref.detrend(type='demean')
st_ref.filter(type='bandpass',freqmin=12, freqmax=22)
st_ref.plot(color='k',starttime=t1, endtime=t2)

ref_st = st_ref[0].data

#
t1 = UTCDateTime(2018, 5, 20, 2, 58, 20) #the format is year:day_of_the_year:month
t2 = t1 + 90
st_ref2 = Stream()
st_ref2 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
st_ref2.detrend(type='linear')
st_ref2.detrend(type='demean')
st_ref2.filter(type='bandpass',freqmin=0.5, freqmax=5)
st_ref2.plot(color='r',starttime=t1, endtime=t2)

ref2_st = st_ref2[0].data


#print("here 0")
#%%

#client = Client('138.253.112.23', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
t1 = UTCDateTime(2019, 10, 1, 0, 0, 0) #the format is year:day_of_the_year:month
t2 = t1 + 24*60*60
sts = Stream()
sts = client.get_waveforms(net, sta, loc, cha, t1 , t2)
#print("here 1")
sts.detrend(type='linear')
sts.detrend(type='demean')

sts1 = sts.copy()
sts2 = sts.copy()
sts3 = sts.copy()



fb1 = sts1.filter(type='bandpass',freqmin=12, freqmax=22)
#fb1.plot(type = 'dayplot',interval=15,starttime=t1, endtime=t2)

fb2 = sts2.filter(type='bandpass',freqmin=0.5, freqmax=22)
fb2.plot(type = 'dayplot',interval=30,starttime=t1, endtime=t2)

fb3=sts3.filter(type='bandpass',freqmin=0.5, freqmax=5)
#fb3.plot(type = 'dayplot',interval=60,starttime=t1, endtime=t2)

#%%



trace_high=fb1[0]
trace_all=fb2[0]
trace_low=fb3[0]

#window endpoints
start= t1 + 10 #time window start 
end= t2
#end=sample[0].stats.endtime 
trs_high = trace_high.slice(starttime = start  , endtime= end) #cut out sample waveform with same window length as chosen event
trs_all = trace_all.slice(starttime = start  , endtime= end)
trs_low = trace_low.slice(starttime = start  , endtime= end)
#trs.filter("bandpass", freqmin=fmin,freqmax=fmax)
#trs_e = obspy.signal.filter.envelope(trs.data)
#print('reference waveform')


trs_high.plot(type='relative',color='b')#, starttime=start , endtime=end)

sr = trace_all.stats.sampling_rate
nsta=int(8*sr)
nlta=int(60*sr)
stream=trs_all.data
cft=recursive_sta_lta(stream, nsta, nlta)
trig_on=2
trig_off=0.5
plot_trigger(trs_all, cft, trig_on, trig_off) 

on_off = trigger_onset(cft,trig_on,trig_off)
#print("here 2")
#%%

catalogue = np.zeros(shape=(0,1))
num=0

shift2 = 200
shift = 10
for x in range(0,len(on_off)):
    event_h = trs_high[on_off[x,0]-2000:on_off[x,1]+2000]
    event_h2 = trs_high[on_off[x,0]-1000:on_off[x,1]+1000]
    event_a=trs_all[on_off[x,0]-2000:on_off[x,1]+2000]
    event_l=trs_low[on_off[x,0]-2000:on_off[x,1]+2000]    
    
    start = t1 + on_off[x,0]/sr 
    end = start + 240
    
    dom,cf, bw = freq_info_NqV(event_l,start,end,sr)
    top_v,top,corell = corel(event_a,event_l,shift)
    top_va,top2,corell2 = corel(abs(ref2_st),abs(event_a),shift2)
    top_vh,top2,corell2 = corel(abs(ref_st),abs(event_h),shift2)
    top_vl,top3,corell3 = corel(abs(ref2_st),abs(event_l),shift2)
    
    
#    plt.figure()
#    plt.plot(event_a,color='b')
#    plt.title(UTCDateTime(start + on_off[x,0]/sr ))
##    
##    plt.figure()
##    plt.plot(event_l,color='k')
##    plt.title(top_vl)
##    
##    plt.figure()
#    plt.plot(event_h,color='r')
#    plt.title(top_vh)
    
    
    
    
    if  top_vl > 0.3:# and top_vh > 0.3:
        
        catalogue = np.lib.pad(catalogue, ((0,1),(0,0)), 'constant', constant_values=(0))
        catalogue[num]= start + on_off[x,0]/sr 
        num+=1
                
#        plt.figure()
#        plt.plot(event_a,color='k')
#        plt.title(UTCDateTime(start + on_off[x,0]/sr ))
##        
##        plt.figure()
#        plt.plot(event_h,color='r')
##        plt.title(top_vh)
##        
##        plt.figure()
##        plt.plot(event_l,color='r')
##        plt.title(top_vl)
#        
#    else:
#        if  top_vl > 0.3:
#            plt.figure()
#            plt.plot(event_a,color='b')
#            plt.title(top_vh)
#        
#        plt.figure()
#        plt.plot(event_a,color='m')
#        plt.title(top_vl)
#        plt.plot(event_l,color='k')
#        
#
#
#
#
##%%
print("")       
print(len(catalogue),"detections from all f detector")
for x in range(0,len(catalogue)):
    
    print(UTCDateTime(catalogue[x]))
#
#




#%%



trace_high2=fb1[0]
trace_all2=fb2[0]
trace_low2=fb3[0]

#window endpoints
start= t1 + 10 #time window start 
end= t2
#end=sample[0].stats.endtime 
trs_high2 = trace_high2.slice(starttime = start  , endtime= end) #cut out sample waveform with same window length as chosen event
trs_all2 = trace_all2.slice(starttime = start  , endtime= end)
trs_low2 = trace_low2.slice(starttime = start  , endtime= end)
#trs.filter("bandpass", freqmin=fmin,freqmax=fmax)
#trs_e = obspy.signal.filter.envelope(trs.data)
#print('reference waveform')


trs_high2.plot(type='relative',color='b')#, starttime=start , endtime=end)

sr = trace_all.stats.sampling_rate
nsta=int(0.5*sr)
nlta=int(25*sr)
stream=trs_high2.data
cft=recursive_sta_lta(stream, nsta, nlta)
trig_on=15
trig_off=0.5
plot_trigger(trs_high2, cft, trig_on, trig_off) 

on_off = trigger_onset(cft,trig_on,trig_off)
#print("here 2")



#%%

catalogue2 = np.zeros(shape=(0,1))
num2=0

shift2 = 200
shift = 10
for x in range(0,len(on_off)):
    event2_h = trs_high2[on_off[x,0]-2000:on_off[x,1]+2000]
    event2_h2 = trs_high2[on_off[x,0]-1000:on_off[x,1]+1000]
    event2_a=trs_all2[on_off[x,0]-2000:on_off[x,1]+2000]
    event2_l=trs_low2[on_off[x,0]-2000:on_off[x,1]+2000]    
    
    start = t1 + on_off[x,0]/sr 
    end = start + 240
    
    dom,cf, bw = freq_info_NqV(event_l,start,end,sr)
    top_v,top,corell = corel(event_a,event_l,shift)
    top_va,top2,corell2 = corel(abs(ref2_st),abs(event_a),shift2)
    top_vh,top2,corell2 = corel(abs(ref_st),abs(event_h),shift2)
    top_vl,top3,corell3 = corel(abs(ref2_st),abs(event_l),shift2)
    
    
#    plt.figure()
#    plt.plot(event_a,color='b')
#    plt.title(UTCDateTime(start + on_off[x,0]/sr ))
##    
##    plt.figure()
##    plt.plot(event_l,color='k')
##    plt.title(top_vl)
##    
##    plt.figure()
#    plt.plot(event_h,color='r')
#    plt.title(top_vh)
    
    
    
    
    if  top_vl > 0.3:# and top_vh > 0.3:
        
        catalogue2 = np.lib.pad(catalogue2, ((0,1),(0,0)), 'constant', constant_values=(0))
        catalogue2[num2]= start + on_off[x,0]/sr 
        num2+=1
                
#        plt.figure()
#        plt.plot(event2_a,color='k')
#        plt.title(UTCDateTime(start + on_off[x,0]/sr ))
##        
###        plt.figure()
#        plt.plot(event2_h,color='r')
##        plt.title(top_vh)
##        
##        plt.figure()
##        plt.plot(event_l,color='r')
##        plt.title(top_vl)
#        
#    else:
#        if  top_vl > 0.3:
#            plt.figure()
#            plt.plot(event_a,color='b')
#            plt.title(top_vh)
#        
#        plt.figure()
#        plt.plot(event_a,color='m')
#        plt.title(top_vl)
#        plt.plot(event_l,color='k')
#        
#
#
#

print("")       
print(len(catalogue2),"detections from high f detector")
for x in range(0,len(catalogue2)):
    
    print(UTCDateTime(catalogue2[x]))

    
#%%
# STATION, CHANNEL (DDF --> 400 Hz), NETWWORK AND LOCATION CODES 
sta = 'FG8' # STATION 
cha = 'BDF' # CHANNEL
net = 'GI'  # 
loc = '02'    # location, it depends mostly of which network you are in. 
client = Client('138.253.112.23', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23

        
        
st_inf = Stream()
st_inf = client.get_waveforms(net, sta, loc, cha, t1 , t2)
st_inf.detrend(type='linear')
st_inf.detrend(type='demean')
print('inf in')
st_inf1 = st_inf.copy()

fb = st_inf1.filter(type='bandpass',freqmin=0.1, freqmax=22)
trace=fb[0]

cha2 = 'BHZ' # CHANNEL
loc2 = '00'  

sts = Stream()
sts = client.get_waveforms(net, sta, loc2, cha2, t1 , t2)
sts.detrend(type='linear')
sts.detrend(type='demean')
print('seis in')
sts1 = sts.copy()
sts2 = sts.copy()
sts0 = sts.copy()
fb1 = sts0.filter(type='bandpass',freqmin=0.5, freqmax=22)
fb2 = sts1.filter(type='bandpass',freqmin=12, freqmax=22)
fb3 = sts2.filter(type='bandpass',freqmin=0.5, freqmax=5)
trace1=fb1[0]
trace2=fb2[0]
trace3=fb3[0]
#%%
#window endpoints
start= t1 + 10 #time window start 
end= t2
trs = trace.slice(starttime = start  , endtime= end) 

trs.plot(type='relative',color='b')#, starttime=start , endtime=end)

sr = trace.stats.sampling_rate
nsta=int(0.5*sr)
nlta=int(25*sr)
stream=trs.data
cft=recursive_sta_lta(stream, nsta, nlta)
trig_on=20
trig_off=1
plot_trigger(trs, cft, trig_on, trig_off) 
on_off = trigger_onset(cft,trig_on,trig_off)




#%%

catalogue3 = np.zeros(shape=(0,1))
num3=0

shift2 = 1000
for x in range(0,len(on_off)):
    event = trs[on_off[x,0]-500:on_off[x,0]+1000]
    eventsh = trace2[on_off[x,0]-2500:on_off[x,0]+4500]
    eventsl = trace3[on_off[x,0]-2500:on_off[x,0]+4500]
    eventsa = trace1[on_off[x,0]-2500:on_off[x,0]+4500]
    
    start = t1 + on_off[x,0]/sr -20
    end = start + 40
    
    trace_fft = trace.slice(starttime = start  , endtime= end) 
    
    window=end-start
    tr_data=trace_fft.data
    m=np.mean(tr_data)
    tr_data = tr_data-m
    famp = abs(np.fft.fft(tr_data))
    
    
    # double peak ratio
    fps=int(len(famp)/(sr/0.5))
    midp=int(len(famp)/(sr/2))
    spe=int(len(famp)/(sr/4.5))
    
    f_peak_m = np.mean(famp[fps:midp])/(window*sr)
    s_peak_m = np.mean(famp[midp:spe])/(window*sr)
    peak_r = f_peak_m/s_peak_m
    
    # 2Hz drop ratio
    
    fpbs=int(len(famp)/(sr/0.5))
    fpbe=int(len(famp)/(sr/1.8))
    
    spbs=int(len(famp)/(sr/2.2))
    spbe=int(len(famp)/(sr/3.5))
    
    peaks_ma = np.mean(famp[fpbs:fpbe])/(window*sr)
    peaks_mb = np.mean(famp[spbs:spbe])/(window*sr)
    peaks_m = (peaks_ma + peaks_mb)/2
    
    trough_m = np.mean(famp[fpbe:spbs])/(window*sr)    
    trough_r = peaks_m/trough_m
    
    
    
    top_vl,top3,corell3 = corel(abs(ref2_st),abs(eventsl),shift2)
    
   
    if  top_vl > 0.3 and trough_r > 1:
        catalogue3 = np.lib.pad(catalogue3, ((0,1),(0,0)), 'constant', constant_values=(0))
        catalogue3[num3]= start + on_off[x,0]/sr 
        num3+=1
        
        
#        plt.figure()
#        plt.plot(event,color='b')
#        plt.title(UTCDateTime(start + on_off[x,0]/sr ))
#        
#        plt.figure()
#        plt.plot(eventsl,color='b')
#        plt.plot(eventsh,color='r')
#        
#        plt.title(UTCDateTime(start + on_off[x,0]/sr ))
#    
#    else:
#        plt.figure()
#        plt.plot(event,color='k')
#        plt.title(trough_r)
#        plt.figure()
#        plt.plot(eventsa,color='r')
#        plt.plot(eventsh,color='k')
#        plt.title(top_vl)
    
        
#%%
print('number of detections from inf detector =',num3)
print('number of rejections =',len(on_off)-num3)

for x in range(0,len(catalogue3)):
    
    print(UTCDateTime(catalogue3[x]))




 
















