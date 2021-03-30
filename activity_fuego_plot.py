#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 13:30:44 2019

@author: root
"""

import os
from collections import OrderedDict
import numpy as np
import obspy
import scipy.signal as sgn
import matplotlib.pyplot as plt 
import matplotlib.mlab as mlab

from obspy.clients.earthworm import Client
from obspy import UTCDateTime
from obspy.signal.trigger import trigger_onset
from numpy import genfromtxt
#from scipy.signal import welch
from obspy import Stream




cat= genfromtxt("/Users/william/Documents/scanner/all_stations/active_stations_fuego.csv", delimiter=',',skip_header=1)

day = cat[:,0]
num_all= cat[:,1]
nums= cat[:,2]
FG3s= cat[:,3]
FG8s= cat[:,4]
FG10s= cat[:,5]
FG11s= cat[:,6]
FG12s= cat[:,7]
FG13s= cat[:,8]
FG14s= cat[:,9]
FG16s= cat[:,10]
    
numa= cat[:,11]
FG8a= cat[:,12]
FG10a= cat[:,13]
FG11a= cat[:,14]
FG12a= cat[:,15]
FG13a= cat[:,16]
FG15a= cat[:,17]
FV01a= cat[:,18]
FV02a= cat[:,19]
FV03a= cat[:,20]
FV04a= cat[:,21]
VF01a= cat[:,22]
VF02a= cat[:,23]
VF03a= cat[:,24]
VF04a= cat[:,25]
VF05a= cat[:,26]
VF06a= cat[:,27]
    
    
    
    
    
#%% Seismic Station Activity


for y in range(3,11):
    switch =0
    count =0
    off_count =0
    num=0
    turned=0
    station = np.zeros(shape=(0,1))
    for x in range(0,len(cat[:,y])):
        if switch == 0:
            if int(cat[:,y][x])==0:
                off_count += 1
            if int(cat[:,y][x])==1:
                station = np.lib.pad(station, ((0,1),(0,0)), 'constant', constant_values=(0))
                station[num]=x
                num+=1
                count = 1 
                off_count = 0
                switch = 1
                turned = 1
                
        if switch ==1 and turned == 0:
            if int(cat[:,y][x])==1:
                count += 1
            if int(cat[:,y][x])==0:
                station = np.lib.pad(station, ((0,1),(0,0)), 'constant', constant_values=(0))
                station[num] = count
                count=0
                off_count = 1
                num +=1
                switch = 0
                
        if turned == 1:
            turned =0
            switch =1
        
    if ((len(station)) % 2) == 0:
        even=1
    else:
        station = np.lib.pad(station, ((0,1),(0,0)), 'constant', constant_values=(0))
        station[num]=len(cat[:,y])-station[num-1]
    
    if y == 3:
        FG3s_ac = station[:]
    if y == 4:
        FG8s_ac = station[:]
    if y == 5:
        FG10s_ac = station[:]
    if y == 6:
        FG11s_ac = station[:]
    if y == 7:
        FG12s_ac = station[:]
    if y == 8:
        FG13s_ac = station[:]
    if y == 9:
        FG14s_ac = station[:]
    if y == 10:
        FG16s_ac = station[:]
    
#    numb=0
#    station_b = np.zeros(shape=(2,1))
#    station_b[0]=0
#    station_b[1]=station[0]
#    for x in range(2,len(station),2):
#        station_b = np.lib.pad(station_b, ((0,1),(0,0)), 'constant', constant_values=(0))
#        station_b[x]=station[x-2]+station[x-1]
#        station_b = np.lib.pad(station_b, ((0,1),(0,0)), 'constant', constant_values=(0))
#        station_b[x+1]= station[x]-station_b[x]
#                
#    if y == 3:
#        FG3s_acb = station_b[:]
#    if y == 4:
#        FG8s_acb = station_b[:]
#    if y == 5:
#        FG10s_acb = station_b[:]
#    if y == 6:
#        FG11s_acb = station_b[:]
#    if y == 7:
#        FG12s_acb = station_b[:]
#    if y == 8:
#        FG13s_acb = station_b[:]
#    if y == 9:
#        FG14s_acb = station_b[:]
#    if y == 10:
#        FG16s_acb = station_b[:]
        
        
        
        
        
    
#%%    
    
    
    
fig, ax = plt.subplots()
ax.set_xlabel('Day [14 Mar2018 - 4 Oct2019]')
ax.set_ylabel('Station')
ax.set_title('Seismic Station Activity')
ax.set_xlim(0, 570)
#ax.set_xlim(0, 45)
ax.set_yticks([3, 8, 13, 18, 23, 28,33,38])
ax.set_yticklabels(['FG3s', 'FG8s','FG10s','FG11s','FG12s','FG13s','FG14s', 'FG16s'])
for x in range(0,len(FG3s_ac),2):
        ax.broken_barh([(FG3s_ac[x], FG3s_ac[x+1])], (1, 4), facecolors='k')
for x in range(0,len(FG8s_ac),2):
        ax.broken_barh([(FG8s_ac[x], FG8s_ac[x+1])], (6, 4), facecolors='k')
for x in range(0,len(FG10s_ac),2):
        ax.broken_barh([(FG10s_ac[x], FG10s_ac[x+1])], (11,4), facecolors='k')
for x in range(0,len(FG11s_ac),2):
        ax.broken_barh([(FG11s_ac[x], FG11s_ac[x+1])], (16,4), facecolors='k')
for x in range(0,len(FG12s_ac),2):
        ax.broken_barh([(FG12s_ac[x], FG12s_ac[x+1])], (21,4), facecolors='k')
for x in range(0,len(FG13s_ac),2):
        ax.broken_barh([(FG13s_ac[x], FG13s_ac[x+1])], (26,4), facecolors='k')
for x in range(0,len(FG14s_ac),2):
        ax.broken_barh([(FG14s_ac[x], FG14s_ac[x+1])], (31,4), facecolors='k') 
for x in range(0,len(FG14s_ac),2):
        ax.broken_barh([(FG14s_ac[x], FG14s_ac[x+1])], (36,4), facecolors='k') 
    
    
    
ax2 = ax.twinx()        
ax2.set_ylabel('Total Active Seismic Stations [#]')  
ax2.plot(day,nums,'r--')            
    
    
    
    
    
#%% Seismic Station Activity


for y in range(12,28):
    switch =0
    count =0
    off_count =0
    num=0
    turned=0
    station = np.zeros(shape=(0,1))
    for x in range(0,len(cat[:,y])):
        if switch == 0:
            if int(cat[:,y][x])==0:
                off_count += 1
            if int(cat[:,y][x])==1:
                station = np.lib.pad(station, ((0,1),(0,0)), 'constant', constant_values=(0))
                station[num]=x
                num+=1
                count = 1 
                off_count = 0
                switch = 1
                turned = 1
                
        if switch ==1 and turned == 0:
            if int(cat[:,y][x])==1:
                count += 1
            if int(cat[:,y][x])==0:
                station = np.lib.pad(station, ((0,1),(0,0)), 'constant', constant_values=(0))
                station[num] = count
                count=0
                off_count = 1
                num +=1
                switch = 0
                
        if turned == 1:
            turned =0
            switch =1
        
    if ((len(station)) % 2) == 0:
        even=1
    else:
        station = np.lib.pad(station, ((0,1),(0,0)), 'constant', constant_values=(0))
        station[num]=len(cat[:,y])-station[num-1]
    
    if y == 12:
        FG8a_ac = station[:]
    if y == 13:
        FG10a_ac = station[:]
    if y == 14:
        FG11a_ac = station[:]
    if y == 15:
        FG12a_ac = station[:]
    if y == 16:
        FG13a_ac = station[:]
    if y == 17:
        FG15a_ac = station[:]
    if y == 18:
        FV01a_ac = station[:]
    if y == 19:
        FV02a_ac = station[:]    
    if y == 20:
        FV03a_ac = station[:]    
    if y == 21:
        FV04a_ac = station[:]    
    if y == 22:
        VF01a_ac = station[:]    
    if y == 23:
        VF02a_ac = station[:]    
    if y == 24:
        VF03a_ac = station[:]    
    if y == 25:
        VF04a_ac = station[:]    
    if y == 26:
        VF05a_ac = station[:]    
    if y == 27:
        VF06a_ac = station[:]    
  
 
    
    #%%
    
    
        
fig, ax = plt.subplots()
ax.set_xlabel('Day [14 Mar2018 - 4 Oct2019]')
ax.set_ylabel('Station')
ax.set_title('Infrasound Station Activity')
ax.set_xlim(0, 570)
ax.set_yticks([3, 8, 13, 18, 23, 28 , 33, 38, 43, 48, 53, 58, 63, 68, 73, 78])
ax.set_yticklabels(['FG8a','FG10a','FG11a','FG12a','FG13a','FG15a', 'FV01a', 'FV02a', 'FV03a', 'FV04a', 'VF01a', 'VF02a', 'VF03a', 'VF04a', 'VF05a', 'VF06a'])
for x in range(0,len(FG8a_ac),2):
        ax.broken_barh([(FG8a_ac[x], FG8a_ac[x+1])], (1, 4), facecolors='k')
for x in range(0,len(FG10a_ac),2):
        ax.broken_barh([(FG10a_ac[x], FG10a_ac[x+1])], (6, 4), facecolors='k')
for x in range(0,len(FG11a_ac),2):
        ax.broken_barh([(FG11a_ac[x], FG11a_ac[x+1])], (11,4), facecolors='k')
for x in range(0,len(FG12a_ac),2):
        ax.broken_barh([(FG12a_ac[x], FG12a_ac[x+1])], (16,4), facecolors='k')
for x in range(0,len(FG13a_ac),2):
        ax.broken_barh([(FG13a_ac[x], FG13a_ac[x+1])], (21,4), facecolors='k')
for x in range(0,len(FG15a_ac),2):
        ax.broken_barh([(FG15a_ac[x], FG15a_ac[x+1])], (26,4), facecolors='k')
for x in range(0,len(FV01a_ac),2):
        ax.broken_barh([(FV01a_ac[x], FV01a_ac[x+1])], (31,4), facecolors='k') 
for x in range(0,len(FV02a_ac),2):
        ax.broken_barh([(FV02a_ac[x], FV02a_ac[x+1])], (36,4), facecolors='k') 
for x in range(0,len(FV03a_ac),2):
        ax.broken_barh([(FV03a_ac[x], FV03a_ac[x+1])], (41,4), facecolors='k') 
for x in range(0,len(FV04a_ac),2):
        ax.broken_barh([(FV04a_ac[x], FV04a_ac[x+1])], (46,4), facecolors='k') 
for x in range(0,len(VF01a_ac),2):
        ax.broken_barh([(VF01a_ac[x], VF01a_ac[x+1])], (51,4), facecolors='k')      
for x in range(0,len(VF02a_ac),2):
        ax.broken_barh([(VF02a_ac[x], VF02a_ac[x+1])], (56,4), facecolors='k')        
for x in range(0,len(VF03a_ac),2):
        ax.broken_barh([(VF03a_ac[x], VF03a_ac[x+1])], (61,4), facecolors='k')
for x in range(0,len(VF04a_ac),2):
        ax.broken_barh([(VF04a_ac[x], VF04a_ac[x+1])], (66,4), facecolors='k')
for x in range(0,len(VF05a_ac),2):
        ax.broken_barh([(VF05a_ac[x], VF05a_ac[x+1])], (71,4), facecolors='k')
for x in range(0,len(VF06a_ac),2):
        ax.broken_barh([(VF06a_ac[x], VF06a_ac[x+1])], (76,4), facecolors='k')        
        
        
        
ax2 = ax.twinx()        
ax2.set_ylabel('Total Active Infrasound Stations [#]')  
ax2.plot(day,numa,'r--')            
        
#%%        
        
fig, ax = plt.subplots()    
ax.plot(day,numa,'k--',label='Inf')        
ax.plot(day,nums,'b--',label='Seis')        
ax.plot(day,num_all,'r-',label='Both')         
ax.set_ylim(0, 14)      
ax.set_xlim(0, 570)  
ax.legend()   
ax.set_ylabel('Number of Stations [#]')
ax.set_xlabel('Day [14 Mar2018 - 4 Oct2019]')       
ax.set_title('Station Activity')      
        
        
        
        
    
  