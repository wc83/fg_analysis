#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 12:41:31 2021

@author: root
"""

import obspy
from obspy import read
import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
from obspy import Stream
from obspy.signal.trigger import classic_sta_lta, recursive_sta_lta
from obspy.signal.trigger import plot_trigger, trigger_onset
from obspy import UTCDateTime
from matplotlib import pyplot
import scipy
from obspy import UTCDateTime
import datetime as dt
import matplotlib.dates as md

#%% Read in catalogue from csv file

tt = genfromtxt("/Users/william/Documents/Fuego_catalogue/analysis/tt2_all_FG12_02b.csv", delimiter=',',skip_header=0)
vasr = genfromtxt("/Users/william/Documents/Fuego_catalogue/analysis/VASR_all_FG12_02b.csv", delimiter=',',skip_header=0)


tt_conv = np.zeros(shape=(len(tt),1))
for x in range(0,len(tt)):
    tt_conv[x] = UTCDateTime((UTCDateTime(dt.datetime.fromordinal(int(tt[x])-366)).timestamp)+6 + int((24*(tt[x]-int(tt[x])))*60*60)).timestamp

#vasr_conv =  np.zeros(shape=(len(vasr),1))
#for x in range(0,len(vasr)):
#    if tt_conv[x] < 1571184000:
#        vasr_conv[x] = vasr[x]
#    else:
#        vasr_conv[x] = vasr[x]*(10**6.5)

xl=min(tt)
xr=max(tt)

plt.figure(figsize=(12,3))
dates=[dt.datetime.fromtimestamp(ts) for ts in  tt_conv[:,0]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=0 )
ax.xaxis.set_major_formatter(xfmt)
plt.semilogy(datenums, vasr,'r',linewidth=0.25)
#plt.plot(datenums, mod[:,2],'bo',markersize=5)

plt.xlabel('Date')
plt.ylabel('VASR [a.u.]')
plt.title("VASR: FG12")
ax.set_xlim([736695,737791])
#ax.set_ylim([10e-11,10e-5])

#ax.set_xlim([737200,737500])

#%% smoothed
tt_conv_s = np.zeros(shape=(len(tt),1))
sum_tt = 0

for x in range(0,len(tt)-12):
    for y in range(0,12):
        sum_tt += UTCDateTime((UTCDateTime(dt.datetime.fromordinal(int(tt[x+y])-366)).timestamp)+6 + int((24*(tt[x+y]-int(tt[x+y])))*60*60)).timestamp
                          
    tt_conv_s[x] = sum_tt/12
    sum_tt = 0
    
    

vasr_conv_s =  np.zeros(shape=(len(vasr),1))
sum_vasr = 0
for x in range(0,len(vasr)-12):
    for y in range(0,12):
        sum_vasr += vasr[x+y] * 7096786000
    
    vasr_conv_s[x] = sum_vasr/12
    sum_vasr = 0
    
    
#    if tt_conv[x] < 1571184000:
#        vasr_conv[x] = vasr[x]
#    else:
#        vasr_conv[x] = vasr[x]*(10**6.5)

xl=min(tt)
xr=max(tt)

plt.figure(figsize=(12,3))
dates=[dt.datetime.fromtimestamp(ts) for ts in  tt_conv_s[:,0]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=0 )
ax.xaxis.set_major_formatter(xfmt)
plt.semilogy(datenums[0:-12], vasr_conv_s[0:-12],'b',linewidth=0.5)
#plt.plot(datenums, mod[:,2],'bo',markersize=5)

plt.xlabel('Date')
plt.ylabel('VASR [a.u.]')
plt.title("VASR: FG12")
ax.set_xlim([736695,737791])
#ax.set_ylim([10e-11,10e-5])

#ax.set_xlim([737200,737500])


