#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 16:18:05 2019

@author: william
"""

import io
import os
from collections import OrderedDict
import numpy as np
import obspy
import scipy.signal as sgn
import matplotlib.pyplot as plt 
import matplotlib.mlab as mlab
from obspy.core import read
from obspy.clients.earthworm import Client
from obspy import UTCDateTime
from obspy.signal.trigger import trigger_onset
from numpy import genfromtxt
#from scipy.signal import welch
from obspy import Stream
from scipy import integrate
from scipy.stats import norm

cat = genfromtxt("/Users/william/Documents/scanner/all_stations/Final_Catalogue_2014_2018.csv", delimiter=',',skip_header=1)


Duration = cat[:,24]
Energy = cat[:,23]
evt = cat[:,0]


#%%
day_no=cat[:,1]
epd=np.zeros(shape=(0,4)) 
exp_no = 0


for x in range(0,1500): # loop over the number of days
    exp_pd = 0
    en_pd =0
    
    epd = np.lib.pad(epd, ((0,1),(0,0)), 'constant', constant_values=(0))
    epd[x][0] = x+1 #day number
    epd[x][2] = 1416787200.0 + x*24*60*60 # time stamp for day
    

    for e in range(0,len(cat)): # loop over catalogue to find number of eplosions for given day
        
        if day_no[e] == x:
            exp_pd +=1
            en_pd += cat[e,23]
            
    
    epd[x][1]=exp_pd #number of exposions for given day
    epd[x][3]=en_pd   
        
#%%

epw=np.zeros(shape=(0,5))
week=1
event_c =0
day_c =0

tot_e =0

for x in range(0,len(epd)): # loop over days in "events per day"
    day_c += 1 # count number of days
    event_c += epd[x][1] # count number of events in week (1 day at a time)
    tot_e  += epd[x][3] # sum energy (1 day at a time)
    
    if day_c == 7: # save number counted once 7 days is reached and reset counters for next week
        epw = np.lib.pad(epw, ((0,1),(0,0)), 'constant', constant_values=(0))
        epw[week-1][0]=week
        epw[week-1][1] = 1416787200.0 + (week-1)*7*24*60*60
        epw[week-1][2]= event_c
        epw[week-1][3] = tot_e
        epw[week-1][4] = tot_e/event_c
        week += 1
        event_c = 0
        day_c =0       
        tot_e =0
        
        
plt.figure()
plt.semilogy(epw[:,0],epw[:,3],'rx')       
plt.figure()
plt.semilogy(epw[:,0],epw[:,4],'rx')            
        
#np.savetxt("/Users/william/Documents/scanner/output_data/Energy_per_week.csv", epw ,delimiter=",",header="week,time,detections,total_E,average_E")      
        
        
        
        
        
        
        
        
        
        
        
        
        
    