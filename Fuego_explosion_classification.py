#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  5 16:14:21 2021

@author: root
"""
# import packages

import numpy as np
import matplotlib.pyplot as plt 
from obspy.core import read
from obspy.clients.earthworm import Client
from obspy import UTCDateTime
from obspy.signal.trigger import trigger_onset
from numpy import genfromtxt
from obspy import Stream
from scipy import integrate





#%% import catalogue

cat = genfromtxt("/Users/william/Documents/Fuego_catalogue/Final_Catalogues/Fuego_explosions.csv", delimiter=',',skip_header=1)

#%% open matrix

Event_Type = np.zeros(shape=(len(cat),1))

dur = cat[:,1]
#%% loop over events


for x in range(90000,len(cat)):
    print(x," of ",len(cat))
    num_seis = cat[x,10]
    num_inf = cat[x,9]
    
    start = cat[x,0]
    duration = cat[x,1]

# if sum inf detections = 0, event type = 5
    if num_seis > 0 and num_inf == 0:
        e_type = 5

# if sum inf > 0 and sum seis > 0, go over tests - assign GCA code 
    if num_seis > 0 and num_inf > 0:
        event_class = np.zeros(shape=(1,int(num_inf)))
        num=0
            
        for y in range(12,45):
            if cat[x,y] == 1:
                
                s = y - 7
                
                try:
                    class_exp = classify_exps(s,start,duration)
                except:
                    print(x, "failed")
                    class_exp = 0.5
                
                event_class[0,num] = class_exp
                num += 1
        
        if sum(event_class[0,:]) >= len(event_class[0])/2:
            e_type = 2
#            print("e_type = gas" )
        else:
            e_type = 4
#            print("e_type = ash" )


# if sum inf > 0 and sum seis = 0, go over tests - assign non-GCA code 
    if num_seis == 0 and num_inf > 0:
        event_class = np.zeros(shape=(1,int(num_inf)))
        num=0
            
        for y in range(12,45):
            if cat[x,y] == 1:
                
                s = y - 7
                
                try:
                    class_exp = classify_exps(s,start,duration)
                except:
                    class_exp = 0.5
                    print(x, "failed")
                
                event_class[0,num] = class_exp
                num += 1
        
        if sum(event_class[0,:]) >= len(event_class[0])/2:
            e_type = 1
#            print("e_type = gas" )
        else:
            e_type = 3
#            print("e_type = ash" )
        
# write class to matrix
    Event_Type[x] = e_type
    
   
    
#%%

np.savetxt("/Users/william/Documents/Fuego_catalogue/tests/Fuego_exp_class_v2c.csv", Event_Type,delimiter=",",header="class_ID")









