#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 12:29:01 2020

@author: root
"""



#%%
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




#%% read catalogue of events

cat = genfromtxt("/Users/william/Documents/Fuego_catalogue/Fuego_Final_scan.csv", delimiter=',',skip_header=1)
ev_time=cat[:,0]
ev_type=cat[:,7]
dur=cat[:,8]
#%% percentage joint and list of joint events

solo_ac = 0
solo_s = 0
joint_ac = 0
joint_s =0

duration_ac = np.zeros(shape=(0,1))
duration_s = np.zeros(shape=(0,1))

seis_joint_list = np.zeros(shape=(0,1))
ac_joint_list = np.zeros(shape=(0,1))

for x in range(0,len(cat)):
    if ev_type[x] == 1:
        
        duration_ac = np.lib.pad(duration_ac, ((0,1),(0,0)), 'constant', constant_values=(0))
        duration_ac[-1] = dur[x]
        
        solo_ac += 1
    if ev_type[x] == 10:
        joint_ac += 1
        ac_joint_list = np.lib.pad(ac_joint_list, ((0,1),(0,0)), 'constant', constant_values=(0))
        ac_joint_list[-1] = ev_time[x]
        
        duration_ac = np.lib.pad(duration_ac, ((0,1),(0,0)), 'constant', constant_values=(0))
        duration_ac[-1] = dur[x]
    
    if ev_type[x] == 2:
        solo_s += 1
        
        duration_s = np.lib.pad(duration_s, ((0,1),(0,0)), 'constant', constant_values=(0))
        duration_s[-1] = dur[x]
        
    if ev_type[x] == 20:
        joint_s += 1
        seis_joint_list = np.lib.pad(seis_joint_list, ((0,1),(0,0)), 'constant', constant_values=(0))
        seis_joint_list[-1] = ev_time[x]
        
        duration_s = np.lib.pad(duration_s, ((0,1),(0,0)), 'constant', constant_values=(0))
        duration_s[-1] = dur[x]
        
#%%
        
ac_percent_joint = joint_ac / (joint_ac + solo_ac) *100
s_percent_joint = joint_s / (joint_s + solo_s) *100

print('')
print('percent acoustic trem with seis = {:.1f}%'.format(ac_percent_joint) )
print('percent seis trem with acoustic = {:.1f}%'.format(s_percent_joint) )
print('')

#%%

trem_diff = np.zeros(shape=(0,3))


for x in range(0,len(ac_joint_list)):
    
    near_time,idx = find_nearest(seis_joint_list,ac_joint_list[x])

    t_diff = near_time[0] - ac_joint_list[x] 
    
    trem_diff = np.lib.pad(trem_diff, ((0,1),(0,0)), 'constant', constant_values=(0))
    trem_diff[-1][0] = t_diff/60
    trem_diff[-1][1] = ac_joint_list[x] #Acoustic time
    trem_diff[-1][2] = near_time #seismic time
    
plt.figure()
plt.hist(trem_diff[:,0],bins=200)
plt.xlim([-60,60])
plt.xlabel('Time difference [mins]')
plt.ylabel('Count')
plt.title('Seismic Tremor start compared to Acoustic Tremor')

med_time = np.median(trem_diff)
mean_time = np.mean(trem_diff)

print('meadian time difference = {:.2f} mins'.format(med_time))
print('mean time difference = {:.2f} mins'.format(mean_time))



#%% durations


mean_ac_dur = np.mean(duration_ac)/60

mean_s_dur = np.mean(duration_s)/60
print("")
print("mean acoustic tremor duration = {:.1f} min".format(mean_ac_dur))
print("mean seismic tremor duration = {:.1f} min".format(mean_s_dur))




#%%


for x in range(200,400):
    
    if -10 < trem_diff[x,0] < -7:
        
        print("Acoustic:",UTCDateTime(trem_diff[x,1]))
        print("Seismic:",UTCDateTime(trem_diff[x,2]))
        print(trem_diff[x,2])
        print("Delay: {:.2f}".format(trem_diff[x,0]))
        print("")





























