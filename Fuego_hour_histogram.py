#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 09:19:11 2021

@author: root
"""

import numpy as np
import matplotlib.pyplot as plt 
from obspy.core import read
from obspy.clients.earthworm import Client
from obspy import UTCDateTime
from obspy.signal.trigger import trigger_onset
from numpy import genfromtxt
from obspy import Stream
from scipy import integrate
import datetime as dt
import matplotlib.dates as md


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
#%%

cat = genfromtxt("/Users/william/Documents/Fuego_catalogue/Final_catalogues/Fuego_explosions.csv", delimiter=',',skip_header=1)
Energy = genfromtxt("/Users/william/Documents/Fuego_catalogue/output_data/explosion_inf_energy.csv", delimiter=',',skip_header=1)
s_Energy=  genfromtxt("/Users/william/Documents/Fuego_catalogue/output_data/explosion_seis_energy.csv", delimiter=',',skip_header=1)

#%%



E_time=Energy[:,0]
En = Energy[:,1]

h0=0
h1=0
h2=0
h3=0
h4=0
h5=0
h6=0
h7=0
h8=0
h9=0
h10=0
h11=0
h12=0
h13=0
h14=0
h15=0
h16=0
h17=0
h18=0
h19=0
h20=0
h21=0
h22=0
h23=0

n0=0
n1=0
n2=0
n3=0
n4=0
n5=0
n6=0
n7=0
n8=0
n9=0
n10=0
n11=0
n12=0
n13=0
n14=0
n15=0
n16=0
n17=0
n18=0
n19=0
n20=0
n21=0
n22=0
n23=0

for x in range(0,len(E_time)):
    if UTCDateTime(E_time[x]).hour == 0:
        h0 += En[x]
        n0 += 1
    if UTCDateTime(E_time[x]).hour == 1:
        h1 += En[x]
        n1 += 1
    if UTCDateTime(E_time[x]).hour == 2:
        h2 += En[x]
        n2 += 1
    if UTCDateTime(E_time[x]).hour == 3:
        h3 += En[x]
        n3 += 1
    if UTCDateTime(E_time[x]).hour == 4:
        h4 += En[x]
        n4 += 1
    if UTCDateTime(E_time[x]).hour == 5:
        h5 += En[x]
        n5 += 1
    if UTCDateTime(E_time[x]).hour == 6:
        h6 += En[x]
        n6 += 1
    if UTCDateTime(E_time[x]).hour == 7:
        h7 += En[x]
        n7 += 1
    if UTCDateTime(E_time[x]).hour == 8:
        h8 += En[x]
        n8 += 1
    if UTCDateTime(E_time[x]).hour == 9:
        h9 += En[x]
        n9 += 1
    if UTCDateTime(E_time[x]).hour == 10:
        h10 += En[x]
        n10 += 1
    if UTCDateTime(E_time[x]).hour == 11:
        h11 += En[x]
        n11 += 1
    if UTCDateTime(E_time[x]).hour == 12:
        h12 += En[x]
        n12 += 1
    if UTCDateTime(E_time[x]).hour == 13:
        h13 += En[x]
        n13 += 1
    if UTCDateTime(E_time[x]).hour == 14:
        h14 += En[x]
        n14 += 1
    if UTCDateTime(E_time[x]).hour == 15:
        h15 += En[x]
        n15 += 1
    if UTCDateTime(E_time[x]).hour == 16:
        h16 += En[x]
        n16 += 1
    if UTCDateTime(E_time[x]).hour == 17:
        h17 += En[x]
        n17 += 1
    if UTCDateTime(E_time[x]).hour == 18:
        h18 += En[x]
        n18 += 1
    if UTCDateTime(E_time[x]).hour == 19:
        h19 += En[x]
        n19 += 1
    if UTCDateTime(E_time[x]).hour == 20:
        h20 += En[x]
        n20 += 1
    if UTCDateTime(E_time[x]).hour == 21:
        h21 += En[x]
        n21 += 1
    if UTCDateTime(E_time[x]).hour == 22:
        h22 += En[x]
        n22 += 1
    if UTCDateTime(E_time[x]).hour == 23:
        h23 += En[x]
        n23 += 1


h0m = h0/n0
h1m = h1/n1
h2m = h2/n2
h3m = h3/n3
h4m = h4/n4
h5m = h5/n5
h6m = h6/n6
h7m = h7/n7
h8m = h8/n8
h9m = h9/n9
h10m = h10/n10
h11m = h11/n11
h12m = h12/n12
h13m = h13/n13
h14m = h14/n14
h15m = h15/n15
h16m = h16/n16
h17m = h17/n17
h18m = h18/n18
h19m = h19/n19
h20m = h20/n20
h21m = h21/n21
h22m = h22/n22
h23m = h23/n23



En_time = np.zeros(shape =(24,2))

En_time[:,0]=[0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5,17.5,18.5,19.5,20.5,21.5,22.5,23.5]

En_time[0,1] = h0m
En_time[1,1] = h1m
En_time[2,1] = h2m
En_time[3,1] = h3m
En_time[4,1] = h4m
En_time[5,1] = h5m
En_time[6,1] = h6m
En_time[7,1] = h7m
En_time[8,1] = h8m
En_time[9,1] = h9m
En_time[10,1] = h10m
En_time[11,1] = h11m
En_time[12,1] = h12m
En_time[13,1] = h13m
En_time[14,1] = h14m
En_time[15,1] = h15m
En_time[16,1] = h16m
En_time[17,1] = h17m
En_time[18,1] = h18m
En_time[19,1] = h19m
En_time[20,1] = h20m
En_time[21,1] = h21m
En_time[22,1] = h22m
En_time[23,1] = h23m


#%%

plt.figure(figsize=(12,3))
plt.bar(En_time[:,0],En_time[:,1])
plt.xlabel('Hour (UTC time)')
plt.ylabel('Explosion Acoustic Energy (J)')
plt.title('Mean Energy vs Time of Day [Full Catalogue]')
plt.xlim([0,24])

#%% partial time


E_time=Energy[:,0]
En = Energy[:,1]

h0=0
h1=0
h2=0
h3=0
h4=0
h5=0
h6=0
h7=0
h8=0
h9=0
h10=0
h11=0
h12=0
h13=0
h14=0
h15=0
h16=0
h17=0
h18=0
h19=0
h20=0
h21=0
h22=0
h23=0

n0=0
n1=0
n2=0
n3=0
n4=0
n5=0
n6=0
n7=0
n8=0
n9=0
n10=0
n11=0
n12=0
n13=0
n14=0
n15=0
n16=0
n17=0
n18=0
n19=0
n20=0
n21=0
n22=0
n23=0

for x in range(1716,8349):
    if UTCDateTime(E_time[x]).hour == 0:
        h0 += En[x]
        n0 += 1
    if UTCDateTime(E_time[x]).hour == 1:
        h1 += En[x]
        n1 += 1
    if UTCDateTime(E_time[x]).hour == 2:
        h2 += En[x]
        n2 += 1
    if UTCDateTime(E_time[x]).hour == 3:
        h3 += En[x]
        n3 += 1
    if UTCDateTime(E_time[x]).hour == 4:
        h4 += En[x]
        n4 += 1
    if UTCDateTime(E_time[x]).hour == 5:
        h5 += En[x]
        n5 += 1
    if UTCDateTime(E_time[x]).hour == 6:
        h6 += En[x]
        n6 += 1
    if UTCDateTime(E_time[x]).hour == 7:
        h7 += En[x]
        n7 += 1
    if UTCDateTime(E_time[x]).hour == 8:
        h8 += En[x]
        n8 += 1
    if UTCDateTime(E_time[x]).hour == 9:
        h9 += En[x]
        n9 += 1
    if UTCDateTime(E_time[x]).hour == 10:
        h10 += En[x]
        n10 += 1
    if UTCDateTime(E_time[x]).hour == 11:
        h11 += En[x]
        n11 += 1
    if UTCDateTime(E_time[x]).hour == 12:
        h12 += En[x]
        n12 += 1
    if UTCDateTime(E_time[x]).hour == 13:
        h13 += En[x]
        n13 += 1
    if UTCDateTime(E_time[x]).hour == 14:
        h14 += En[x]
        n14 += 1
    if UTCDateTime(E_time[x]).hour == 15:
        h15 += En[x]
        n15 += 1
    if UTCDateTime(E_time[x]).hour == 16:
        h16 += En[x]
        n16 += 1
    if UTCDateTime(E_time[x]).hour == 17:
        h17 += En[x]
        n17 += 1
    if UTCDateTime(E_time[x]).hour == 18:
        h18 += En[x]
        n18 += 1
    if UTCDateTime(E_time[x]).hour == 19:
        h19 += En[x]
        n19 += 1
    if UTCDateTime(E_time[x]).hour == 20:
        h20 += En[x]
        n20 += 1
    if UTCDateTime(E_time[x]).hour == 21:
        h21 += En[x]
        n21 += 1
    if UTCDateTime(E_time[x]).hour == 22:
        h22 += En[x]
        n22 += 1
    if UTCDateTime(E_time[x]).hour == 23:
        h23 += En[x]
        n23 += 1


h0m = h0/n0
h1m = h1/n1
h2m = h2/n2
h3m = h3/n3
h4m = h4/n4
h5m = h5/n5
h6m = h6/n6
h7m = h7/n7
h8m = h8/n8
h9m = h9/n9
h10m = h10/n10
h11m = h11/n11
h12m = h12/n12
h13m = h13/n13
h14m = h14/n14
h15m = h15/n15
h16m = h16/n16
h17m = h17/n17
h18m = h18/n18
h19m = h19/n19
h20m = h20/n20
h21m = h21/n21
h22m = h22/n22
h23m = h23/n23



En_time = np.zeros(shape =(24,2))

En_time[:,0]=[0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5,17.5,18.5,19.5,20.5,21.5,22.5,23.5]

En_time[0,1] = h0m
En_time[1,1] = h1m
En_time[2,1] = h2m
En_time[3,1] = h3m
En_time[4,1] = h4m
En_time[5,1] = h5m
En_time[6,1] = h6m
En_time[7,1] = h7m
En_time[8,1] = h8m
En_time[9,1] = h9m
En_time[10,1] = h10m
En_time[11,1] = h11m
En_time[12,1] = h12m
En_time[13,1] = h13m
En_time[14,1] = h14m
En_time[15,1] = h15m
En_time[16,1] = h16m
En_time[17,1] = h17m
En_time[18,1] = h18m
En_time[19,1] = h19m
En_time[20,1] = h20m
En_time[21,1] = h21m
En_time[22,1] = h22m
En_time[23,1] = h23m


#%%

plt.figure(figsize=(12,3))
plt.bar(En_time[:,0],En_time[:,1])
plt.xlabel('Hour (UTC time)')
plt.ylabel('Explosion Acoustic Energy (J)')
plt.title('Mean Energy vs Time of Day [Oct 1 - Nov 15 2018]')
plt.xlim([0,24])




#%% full time - seis


E_time=s_Energy[:,0]
En = s_Energy[:,1]

h0=0
h1=0
h2=0
h3=0
h4=0
h5=0
h6=0
h7=0
h8=0
h9=0
h10=0
h11=0
h12=0
h13=0
h14=0
h15=0
h16=0
h17=0
h18=0
h19=0
h20=0
h21=0
h22=0
h23=0

n0=0
n1=0
n2=0
n3=0
n4=0
n5=0
n6=0
n7=0
n8=0
n9=0
n10=0
n11=0
n12=0
n13=0
n14=0
n15=0
n16=0
n17=0
n18=0
n19=0
n20=0
n21=0
n22=0
n23=0

for x in range(0,len(E_time)):
    if UTCDateTime(E_time[x]).hour == 0:
        h0 += En[x]
        n0 += 1
    if UTCDateTime(E_time[x]).hour == 1:
        h1 += En[x]
        n1 += 1
    if UTCDateTime(E_time[x]).hour == 2:
        h2 += En[x]
        n2 += 1
    if UTCDateTime(E_time[x]).hour == 3:
        h3 += En[x]
        n3 += 1
    if UTCDateTime(E_time[x]).hour == 4:
        h4 += En[x]
        n4 += 1
    if UTCDateTime(E_time[x]).hour == 5:
        h5 += En[x]
        n5 += 1
    if UTCDateTime(E_time[x]).hour == 6:
        h6 += En[x]
        n6 += 1
    if UTCDateTime(E_time[x]).hour == 7:
        h7 += En[x]
        n7 += 1
    if UTCDateTime(E_time[x]).hour == 8:
        h8 += En[x]
        n8 += 1
    if UTCDateTime(E_time[x]).hour == 9:
        h9 += En[x]
        n9 += 1
    if UTCDateTime(E_time[x]).hour == 10:
        h10 += En[x]
        n10 += 1
    if UTCDateTime(E_time[x]).hour == 11:
        h11 += En[x]
        n11 += 1
    if UTCDateTime(E_time[x]).hour == 12:
        h12 += En[x]
        n12 += 1
    if UTCDateTime(E_time[x]).hour == 13:
        h13 += En[x]
        n13 += 1
    if UTCDateTime(E_time[x]).hour == 14:
        h14 += En[x]
        n14 += 1
    if UTCDateTime(E_time[x]).hour == 15:
        h15 += En[x]
        n15 += 1
    if UTCDateTime(E_time[x]).hour == 16:
        h16 += En[x]
        n16 += 1
    if UTCDateTime(E_time[x]).hour == 17:
        h17 += En[x]
        n17 += 1
    if UTCDateTime(E_time[x]).hour == 18:
        h18 += En[x]
        n18 += 1
    if UTCDateTime(E_time[x]).hour == 19:
        h19 += En[x]
        n19 += 1
    if UTCDateTime(E_time[x]).hour == 20:
        h20 += En[x]
        n20 += 1
    if UTCDateTime(E_time[x]).hour == 21:
        h21 += En[x]
        n21 += 1
    if UTCDateTime(E_time[x]).hour == 22:
        h22 += En[x]
        n22 += 1
    if UTCDateTime(E_time[x]).hour == 23:
        h23 += En[x]
        n23 += 1


h0m = h0/n0
h1m = h1/n1
h2m = h2/n2
h3m = h3/n3
h4m = h4/n4
h5m = h5/n5
h6m = h6/n6
h7m = h7/n7
h8m = h8/n8
h9m = h9/n9
h10m = h10/n10
h11m = h11/n11
h12m = h12/n12
h13m = h13/n13
h14m = h14/n14
h15m = h15/n15
h16m = h16/n16
h17m = h17/n17
h18m = h18/n18
h19m = h19/n19
h20m = h20/n20
h21m = h21/n21
h22m = h22/n22
h23m = h23/n23



En_time = np.zeros(shape =(24,2))

En_time[:,0]=[0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5,17.5,18.5,19.5,20.5,21.5,22.5,23.5]

En_time[0,1] = h0m
En_time[1,1] = h1m
En_time[2,1] = h2m
En_time[3,1] = h3m
En_time[4,1] = h4m
En_time[5,1] = h5m
En_time[6,1] = h6m
En_time[7,1] = h7m
En_time[8,1] = h8m
En_time[9,1] = h9m
En_time[10,1] = h10m
En_time[11,1] = h11m
En_time[12,1] = h12m
En_time[13,1] = h13m
En_time[14,1] = h14m
En_time[15,1] = h15m
En_time[16,1] = h16m
En_time[17,1] = h17m
En_time[18,1] = h18m
En_time[19,1] = h19m
En_time[20,1] = h20m
En_time[21,1] = h21m
En_time[22,1] = h22m
En_time[23,1] = h23m


#%%

plt.figure(figsize=(12,3))
plt.bar(En_time[:,0],En_time[:,1])
plt.xlabel('Hour (UTC time)')
plt.ylabel('Explosion Seismic Energy (J)')
plt.title('Mean Sismic Energy vs Time of Day [Full Catalogue]')
plt.xlim([0,24])

#%% partial time - seis


E_time= s_Energy[:,0]
En = s_Energy[:,1]

h0=0
h1=0
h2=0
h3=0
h4=0
h5=0
h6=0
h7=0
h8=0
h9=0
h10=0
h11=0
h12=0
h13=0
h14=0
h15=0
h16=0
h17=0
h18=0
h19=0
h20=0
h21=0
h22=0
h23=0

n0=0
n1=0
n2=0
n3=0
n4=0
n5=0
n6=0
n7=0
n8=0
n9=0
n10=0
n11=0
n12=0
n13=0
n14=0
n15=0
n16=0
n17=0
n18=0
n19=0
n20=0
n21=0
n22=0
n23=0

for x in range(1716,8349):
    if UTCDateTime(E_time[x]).hour == 0:
        h0 += En[x]
        n0 += 1
    if UTCDateTime(E_time[x]).hour == 1:
        h1 += En[x]
        n1 += 1
    if UTCDateTime(E_time[x]).hour == 2:
        h2 += En[x]
        n2 += 1
    if UTCDateTime(E_time[x]).hour == 3:
        h3 += En[x]
        n3 += 1
    if UTCDateTime(E_time[x]).hour == 4:
        h4 += En[x]
        n4 += 1
    if UTCDateTime(E_time[x]).hour == 5:
        h5 += En[x]
        n5 += 1
    if UTCDateTime(E_time[x]).hour == 6:
        h6 += En[x]
        n6 += 1
    if UTCDateTime(E_time[x]).hour == 7:
        h7 += En[x]
        n7 += 1
    if UTCDateTime(E_time[x]).hour == 8:
        h8 += En[x]
        n8 += 1
    if UTCDateTime(E_time[x]).hour == 9:
        h9 += En[x]
        n9 += 1
    if UTCDateTime(E_time[x]).hour == 10:
        h10 += En[x]
        n10 += 1
    if UTCDateTime(E_time[x]).hour == 11:
        h11 += En[x]
        n11 += 1
    if UTCDateTime(E_time[x]).hour == 12:
        h12 += En[x]
        n12 += 1
    if UTCDateTime(E_time[x]).hour == 13:
        h13 += En[x]
        n13 += 1
    if UTCDateTime(E_time[x]).hour == 14:
        h14 += En[x]
        n14 += 1
    if UTCDateTime(E_time[x]).hour == 15:
        h15 += En[x]
        n15 += 1
    if UTCDateTime(E_time[x]).hour == 16:
        h16 += En[x]
        n16 += 1
    if UTCDateTime(E_time[x]).hour == 17:
        h17 += En[x]
        n17 += 1
    if UTCDateTime(E_time[x]).hour == 18:
        h18 += En[x]
        n18 += 1
    if UTCDateTime(E_time[x]).hour == 19:
        h19 += En[x]
        n19 += 1
    if UTCDateTime(E_time[x]).hour == 20:
        h20 += En[x]
        n20 += 1
    if UTCDateTime(E_time[x]).hour == 21:
        h21 += En[x]
        n21 += 1
    if UTCDateTime(E_time[x]).hour == 22:
        h22 += En[x]
        n22 += 1
    if UTCDateTime(E_time[x]).hour == 23:
        h23 += En[x]
        n23 += 1


h0m = h0/n0
h1m = h1/n1
h2m = h2/n2
h3m = h3/n3
h4m = h4/n4
h5m = h5/n5
h6m = h6/n6
h7m = h7/n7
h8m = h8/n8
h9m = h9/n9
h10m = h10/n10
h11m = h11/n11
h12m = h12/n12
h13m = h13/n13
h14m = h14/n14
h15m = h15/n15
h16m = h16/n16
h17m = h17/n17
h18m = h18/n18
h19m = h19/n19
h20m = h20/n20
h21m = h21/n21
h22m = h22/n22
h23m = h23/n23



En_time = np.zeros(shape =(24,2))

En_time[:,0]=[0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5,17.5,18.5,19.5,20.5,21.5,22.5,23.5]

En_time[0,1] = h0m
En_time[1,1] = h1m
En_time[2,1] = h2m
En_time[3,1] = h3m
En_time[4,1] = h4m
En_time[5,1] = h5m
En_time[6,1] = h6m
En_time[7,1] = h7m
En_time[8,1] = h8m
En_time[9,1] = h9m
En_time[10,1] = h10m
En_time[11,1] = h11m
En_time[12,1] = h12m
En_time[13,1] = h13m
En_time[14,1] = h14m
En_time[15,1] = h15m
En_time[16,1] = h16m
En_time[17,1] = h17m
En_time[18,1] = h18m
En_time[19,1] = h19m
En_time[20,1] = h20m
En_time[21,1] = h21m
En_time[22,1] = h22m
En_time[23,1] = h23m


#%%

plt.figure(figsize=(12,3))
plt.bar(En_time[:,0],En_time[:,1])
plt.xlabel('Hour (UTC time)')
plt.ylabel('Explosion Seismic Energy (J)')
plt.title('Mean Seismic Energy vs Time of Day [Oct 1 - Nov 15 2018]')
plt.xlim([0,24])

#%% dayplot

# inf calibs
c1 = 1.22E-04
c2 = 6.10E-05
c3 = 2.78E-04
c4 = 6.25E-06
c5 = 2.50E-06
#for seismics
c6 = 3.30E-09
c7 = 1.33E-09


sta = 'FG12' # STATION 
cha = 'BDF' # CHANNEL
net = 'GI'  # 
loc = '02'    # location, it depends mostly of which network you are in. 
client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23

t1 = UTCDateTime(2018,11,12,0,0,0) #the format is year:day_of_the_year:month
t2 = t1 +  24*60*60 

st = Stream()
st = client.get_waveforms(net, sta, loc, cha, t1 , t2 + 2)
tr = st[0] 
tr.data = tr.data * c3

sr = st[0].stats.sampling_rate
tr.detrend(type='linear')
tr.detrend(type='demean')
tr.filter(type='bandpass',freqmin=0.1, freqmax=20) 
tr.plot(type='dayplot', starttime=t1, endtime=t2)
tr.plot(starttime=t1, endtime=t2,color='b')


sta = 'FG12' # STATION 
cha = 'BHZ' # CHANNEL
net = 'GI'  # 
loc = '00'    # location, it depends mostly of which network you are in. 

st = Stream()
st = client.get_waveforms(net, sta, loc, cha, t1 , t2 + 2)
tr = st[0] 
tr.data = tr.data * c6

sr = st[0].stats.sampling_rate
tr.detrend(type='linear')
tr.detrend(type='demean')
tr.filter(type='bandpass',freqmin=0.1, freqmax=20) 
tr.plot(type='dayplot',  starttime=t1, endtime=t2)
tr.plot(starttime=t1, endtime=t2,color='r')


