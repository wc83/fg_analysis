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


cat= genfromtxt("/Users/william/Documents/Fuego_catalogue/active_stations_fuego_all_instruments_v2.csv", delimiter=',',skip_header=1)

day = cat[:,0]

date_t=np.zeros(shape=(len(day),1))

start_date = 1520640000 # (10th march 2018)

for x in range(0,len(date_t)):
    date_t[x,0]= (x*60*60*24) + start_date


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

FG8_1= cat[:,12]
FG8_2= cat[:,13]
FG8_3= cat[:,14]

FG10_1 = cat[:,15]
FG10_2 = cat[:,16]
FG10_3 = cat[:,17]

FG11_1= cat[:,18]
FG11_2= cat[:,19]
FG11_4= cat[:,20]
FG11_5= cat[:,21]

FG12_1= cat[:,22]
FG12_2= cat[:,23]
FG12_3= cat[:,24]

FG13_1= cat[:,25]
FG13_2= cat[:,26]
FG13_3= cat[:,27]

FG15_1= cat[:,28]
FG15_2= cat[:,29]
FG15_3= cat[:,30]
FG15_4= cat[:,31]
FG15_5= cat[:,32]
FG15_6= cat[:,33]

FV01a= cat[:,34]
FV02a= cat[:,35]
FV03a= cat[:,36]
FV04a= cat[:,37]
VF01a= cat[:,38]
VF02a= cat[:,39]
VF03a= cat[:,40]
VF04a= cat[:,41]
VF05a= cat[:,42]
VF06a= cat[:,43]
    
    
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
    
    
    
fig, ax = plt.subplots(figsize=(8,5))
ax.set_xlabel('Day ') 
ax.set_ylabel('Station')
ax.set_title('Seismic Station Activity')
ax.set_xlim(0, 820)
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


for y in range(12,44):
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
        FG8_1ac = station[:]
    if y == 13:
        FG8_2ac = station[:]
    if y == 14:
        FG8_3ac = station[:]
        
    if y == 15:
        FG10_1ac = station[:]
    if y == 16:
        FG10_2ac = station[:]
    if y == 17:
        FG10_3ac = station[:]
        
    if y == 18:
        FG11_1ac = station[:]
    if y == 19:
        FG11_2ac = station[:]
    if y == 20:
        FG11_4ac = station[:]
    if y == 21:
        FG11_5ac = station[:]
        
    if y == 22:
        FG12_1ac = station[:]
    if y == 23:
        FG12_2ac = station[:]
    if y == 24:
        FG12_3ac = station[:]
        
    if y == 25:
        FG13_1ac = station[:]
    if y == 26:
        FG13_2ac = station[:]
    if y == 27:
        FG13_3ac = station[:]
        
    if y == 28:
        FG15_1ac = station[:]
    if y == 29:
        FG15_2ac = station[:]
    if y == 30:
        FG15_3ac = station[:]
    if y == 31:
        FG15_4ac = station[:]
    if y == 32:
        FG15_5ac = station[:]
    if y == 33:
        FG15_6ac = station[:]
        
    if y == 34:
        FV01a_ac = station[:]
    if y == 35:
        FV02a_ac = station[:]    
    if y == 36:
        FV03a_ac = station[:]    
    if y == 37:
        FV04a_ac = station[:]    
    if y == 38:
        VF01a_ac = station[:]    
    if y == 39:
        VF02a_ac = station[:]    
    if y == 40:
        VF03a_ac = station[:]    
    if y == 41:
        VF04a_ac = station[:]    
    if y == 42:
        VF05a_ac = station[:]    
    if y == 43:
        VF06a_ac = station[:]    
  
 
    
    #%%
    
    
        
fig, ax = plt.subplots(figsize=(8,10))
ax.set_xlabel('Day ') 
ax.set_ylabel('Station')
ax.set_title('Infrasound Station Activity')
ax.set_xlim(0, 820)
ax.set_ylim(0, 141)
ax.set_yticks([3, 8, 13, 18, 23, 28 , 33, 38, 43, 48, 53, 58, 63, 68, 73, 78,83,88,93,98,103,108,113,118,123,128,133,138])
ax.set_yticklabels(['FG8_01','FG8_02','FG8_03','FG10_01','FG10_02','FG10_03', 'FG12_01','FG12_02','FG12_03','FG13_01','FG13_02','FG13_03','FG15_01','FG15_02','FG15_03','FG15_04','FG15_05','FG15_06', 'FV01', 'FV02', 'FV03', 'FV04', 'VF01', 'VF02', 'VF03', 'VF04', 'VF05', 'VF06'])
for x in range(0,len(FG8_1ac),2):
        ax.broken_barh([(FG8_1ac[x], FG8_1ac[x+1])], (1, 4), facecolors='k')
for x in range(0,len(FG8_2ac),2):
        ax.broken_barh([(FG8_2ac[x], FG8_2ac[x+1])], (6, 4), facecolors='k')
for x in range(0,len(FG8_3ac),2):
        ax.broken_barh([(FG8_3ac[x], FG8_3ac[x+1])], (11, 4), facecolors='k')
        
for x in range(0,len(FG10_1ac),2):
        ax.broken_barh([(FG10_1ac[x], FG10_1ac[x+1])], (16, 4), facecolors='k')
for x in range(0,len(FG10_2ac),2):
        ax.broken_barh([(FG10_2ac[x], FG10_2ac[x+1])], (21, 4), facecolors='k')
for x in range(0,len(FG10_3ac),2):
        ax.broken_barh([(FG10_3ac[x], FG10_3ac[x+1])], (26, 4), facecolors='k')
        
for x in range(0,len(FG12_1ac),2):
        ax.broken_barh([(FG12_1ac[x], FG12_1ac[x+1])], (31,4), facecolors='k')
for x in range(0,len(FG12_2ac),2):
        ax.broken_barh([(FG12_2ac[x], FG12_2ac[x+1])], (36,4), facecolors='k')
for x in range(0,len(FG12_3ac),2):
        ax.broken_barh([(FG12_3ac[x], FG12_3ac[x+1])], (41,4), facecolors='k')
        
for x in range(0,len(FG13_1ac),2):
        ax.broken_barh([(FG13_1ac[x], FG13_1ac[x+1])], (46,4), facecolors='k')
for x in range(0,len(FG13_2ac),2):
        ax.broken_barh([(FG13_2ac[x], FG13_2ac[x+1])], (51,4), facecolors='k')
for x in range(0,len(FG13_3ac),2):
        ax.broken_barh([(FG13_3ac[x], FG13_3ac[x+1])], (56,4), facecolors='k')
        
for x in range(0,len(FG15_1ac),2):
        ax.broken_barh([(FG15_1ac[x], FG15_1ac[x+1])], (61,4), facecolors='k')
for x in range(0,len(FG15_2ac),2):
        ax.broken_barh([(FG15_2ac[x], FG15_2ac[x+1])], (66,4), facecolors='k')
for x in range(0,len(FG15_3ac),2):
        ax.broken_barh([(FG15_3ac[x], FG15_3ac[x+1])], (71,4), facecolors='k')
for x in range(0,len(FG15_4ac),2):
        ax.broken_barh([(FG15_4ac[x], FG15_4ac[x+1])], (76,4), facecolors='k')
for x in range(0,len(FG15_5ac),2):
        ax.broken_barh([(FG15_5ac[x], FG15_5ac[x+1])], (81,4), facecolors='k')
for x in range(0,len(FG15_6ac),2):
        ax.broken_barh([(FG15_6ac[x], FG15_6ac[x+1])], (86,4), facecolors='k')
        
for x in range(0,len(FV01a_ac),2):
        ax.broken_barh([(FV01a_ac[x], FV01a_ac[x+1])], (91,4), facecolors='k') 
for x in range(0,len(FV02a_ac),2):
        ax.broken_barh([(FV02a_ac[x], FV02a_ac[x+1])], (96,4), facecolors='k') 
for x in range(0,len(FV03a_ac),2):
        ax.broken_barh([(FV03a_ac[x], FV03a_ac[x+1])], (101,4), facecolors='k') 
for x in range(0,len(FV04a_ac),2):
        ax.broken_barh([(FV04a_ac[x], FV04a_ac[x+1])], (106,4), facecolors='k') 
for x in range(0,len(VF01a_ac),2):
        ax.broken_barh([(VF01a_ac[x], VF01a_ac[x+1])], (111,4), facecolors='k')      
for x in range(0,len(VF02a_ac),2):
        ax.broken_barh([(VF02a_ac[x], VF02a_ac[x+1])], (116,4), facecolors='k')        
for x in range(0,len(VF03a_ac),2):
        ax.broken_barh([(VF03a_ac[x], VF03a_ac[x+1])], (121,4), facecolors='k')
for x in range(0,len(VF04a_ac),2):
        ax.broken_barh([(VF04a_ac[x], VF04a_ac[x+1])], (126,4), facecolors='k')
for x in range(0,len(VF05a_ac),2):
        ax.broken_barh([(VF05a_ac[x], VF05a_ac[x+1])], (131,4), facecolors='k')
for x in range(0,len(VF06a_ac),2):
        ax.broken_barh([(VF06a_ac[x], VF06a_ac[x+1])], (136,4), facecolors='k')        
        
        
        
ax2 = ax.twinx()        
ax2.set_ylabel('Total Active Infrasound Stations [#]')  
ax2.set_ylim(0, 20)
ax2.plot(day,numa,'r--')            
        
#%%        
        
fig, ax = plt.subplots(figsize=(8,5))    
ax.plot(day,numa,'k--',label='Inf')        
ax.plot(day,nums,'b--',label='Seis')        
ax.plot(day,num_all,'r-',label='Both')         
#ax.set_ylim(0, 14)      
ax.set_xlim(0, 820)  
ax.legend()   
ax.set_ylabel('Number of Stations [#]')
ax.set_xlabel('Day ')       
ax.set_title('Station Activity')      
        
        
        
        
    
  