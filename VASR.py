#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 14:07:45 2021

@author: root
"""

import os
from collections import OrderedDict
import numpy as np
import obspy
import obspy.signal
import scipy.signal as sgn
import matplotlib.pyplot as plt 
import matplotlib.mlab as mlab
from scipy import integrate
from obspy.clients.earthworm import Client
from obspy import UTCDateTime
from obspy.signal.trigger import plot_trigger, trigger_onset
from obspy.signal.trigger import classic_sta_lta, recursive_sta_lta
from obspy import Stream
from numpy import argmax
import matplotlib.dates as md
import datetime as dt
#%%

#PARAMETERS BELOW - THESE ARE OK FOR FUEGO
rhoA = 1 # Atmospheric density
va = 310 # sound velocity
rhoS = 2700 # rock density
vs = 2500 # seismic velocity
winSec = 10 # window length for intensity calculation in seconds (palacios uses 60 seconds)
noverlapSec = 5 # amount of overlap bnetween windows in seconds

fs_default = 50 # default frequency

calSeis = 1 # calibration given in counts/m/s
calAcous = 1 # calibration given in counts/Pa

#COORDINATES OF VENT AND STATION
srcCoords = [14.4743, -90.8811] # Fuego vent
staCoords = [14.43651, -90.83606] # FG12 location
#staCoords = [14.40677 -90.81859]; % FG13 location

sta = 'FG12' # STATION 
cha1 = 'BHZ' # CHANNEL
cha2 = 'BDF' # CHANNEL
net = 'GI'  # 
loc1 = '00'    # location, it depends mostly of which network you are in. 
loc2 = '02' 
client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23

t1 = UTCDateTime(2018,11,2,8,0,0)
t2 = UTCDateTime(2018,11,2,10,0,0)


st1 = Stream()
st1 = client.get_waveforms(net, sta, loc1, cha1, t1 , t2 )
tr1 = st1[0]

print("stream 1 in")

sr1 = st1[0].stats.sampling_rate

tr1.detrend(type='linear')
tr1.detrend(type='demean')
tr1.filter(type='bandpass',freqmin=0.5, freqmax=15) 


r_FG12 = 1#6962

tAcous = r_FG12/va;
tSeis = r_FG12/vs;

ttDiff = tAcous-tSeis


st2 = Stream()
st2 = client.get_waveforms(net, sta, loc2, cha2, t1+ttDiff , t2+ttDiff )
tr2 = st2[0]
tr2.detrend(type='linear')
tr2.detrend(type='demean')
tr2.filter(type='bandpass',freqmin=0.5, freqmax=15)

plt.figure()
plt.plot(tr1,'k-')
plt.plot(tr2,'b--')


sr2=st2[0].stats.sampling_rate

shiftSamples = ttDiff*sr2;

seis_dat = tr1.data
ac_dat = tr2.data

len_diff = len(seis_dat) - len(ac_dat)

if len_diff > 0:
    
    for x in range(0,len_diff):
        ac_dat = np.append(ac_dat,1)
        
if len_diff < 0:
    
    for x in range(0,-len_diff):
        seis_dat = np.append(seis_dat,1)
        
el_start = 0
el_end = int(winSec*sr1)
num=0
VASR_a = np.zeros(shape=(0,2))


while el_end < len(seis_dat):
    
    
    y1 = seis_dat[el_start:el_end]
    y2 = ac_dat[el_start:el_end]
    
    
    
    Is = ((rhoS*vs)/winSec)*sum((y1**2)/sr1);
    Ia = (1/(rhoA*va*winSec))*sum((y2**2)/sr2);
    VASR = Ia/Is;
    
    VASR_a = np.lib.pad(VASR_a, ((0,1),(0,0)), 'constant', constant_values=(0))
    
    VASR_a[num,0] = t1 + ((el_start+el_end)/2)/sr1
    VASR_a[num,1] = VASR
    
    num+=1
    el_start += int((winSec - noverlapSec) *sr1)
    el_end = int(el_start + winSec*sr1)
    
   
    
    
    
   #%% 
plt.figure(figsize=(12,3))  

dates=[dt.datetime.fromtimestamp(ts) for ts in VASR_a[:,0]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y:%H')
plt.xticks( rotation=0 )
ax.xaxis.set_major_formatter(xfmt)
plt.semilogy(datenums, VASR_a[:,1],'b-')
plt.xlabel('Date and Time')
plt.ylabel('VASR [a.u.]')


















































