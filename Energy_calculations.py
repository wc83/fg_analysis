#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 16:42:44 2019

@author: william
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

#%% constants

A=1
rhoE=2500
cE=2000
pi=3.14159
r1=4630
r2=3370
r3=2310
r4=1300
r5=810
r6=7660
rs1=5710
rs2=5490
rs3=3900
rs4=5520
rs5=4290
rs6=2610


#%% read catalogue of events

cat = genfromtxt("/Users/william/Documents/scanner/all_stations/Final_Catalogue_2014_2018.csv", delimiter=',',skip_header=1)
ev_t = cat[:,0]
day_num= cat[:,1]
doy=cat[:,2]
year= cat[:,3]
month = cat[:,4]
day = cat[:,5]
hour = cat[:,6]
minute = cat[:,7]
second = cat[:,8]
milisec = cat[:,9]
lb01 = cat[:,10]
lb02 = cat[:,11]
lb03 = cat[:,12]
lb04 = cat[:,13]
lb05 = cat[:,14]
lb06 = cat[:,15]
ls01 = cat[:,16]
ls02 = cat[:,17]
ls03 = cat[:,18]
ls04 = cat[:,19]
ls05 = cat[:,20]
ls06 = cat[:,21]
trust = cat[:,22]

dur = cat[:,24]

    

net = 'Z4' 
st_n = Stream()
st_e = Stream() 
client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23



enz = np.zeros(shape=(0,2))
numz=0

s=200
e=s+1
for x in range(s,e):
    for y in range(0,11):
        


        t1=UTCDateTime(ev_t[x] -10)
        t2=t1+80
        try:  
            # chose station 
            if y == 0:
                sta ='LB01' 
            if y == 1:
                sta ='LB02' 
            if y == 2:
                sta ='LB03' 
            if y == 3:
                sta ='LB04' 
            if y == 4:
                sta ='LB05' 
            if y == 5:
                sta ='LS01' 
            if y == 6:
                sta ='LS02' 
            if y == 7:
                sta ='LS03' 
            if y == 8:
                sta ='LS04' 
            if y == 9:
                sta ='LS05' 
            if y == 10:
                sta ='LS06' 

            
            
            # converstion factors to LB03 
            if sta == 'LB03':
                conf = 1
                cha1='HHZ'
                cha2='HHE'
                cha3='HHN'
                r=r3
            if sta == 'LB01' and ev_t[x] < 1451606400.0:
                conf = 0.4394913
                cha1='HHZ'
                cha2='HHE'
                cha3='HHN'
                r=r1
            if sta == 'LB02' and ev_t[x] < 1451606400.0:
                conf = 0.17515
                cha1='HHZ'
                cha2='HHE'
                cha3='HHN'
                r=r2
            if sta == 'LS05' and ev_t[x] < 1451606400.0:
                conf = 0.6900589
                cha1='EHZ'
                cha2='EHE'
                cha3='EHN'
                r=rs5
            if sta == 'LS06' and ev_t[x] < 1451606400.0:
                conf = 0.1898134 
                cha1='EHZ'
                cha2='EHE'
                cha3='EHN'
                r=rs6
            if sta == 'LB01' and ev_t[x] > 1451606400.0:
                conf = 0.308727
                cha1='HHZ'
                cha2='HHE'
                cha3='HHN'
                r=r1
            if sta == 'LB02' and ev_t[x] > 1451606400.0:
                conf = 0.13611627
                cha1='HHZ'
                cha2='HHE'
                cha3='HHN'
                r=r2
            if sta == 'LS05' and ev_t[x] > 1451606400.0:
                conf = 1.78750 
                cha1='EHZ'
                cha2='EHE'
                cha3='EHN'
                r=rs5
            if sta == 'LS06' and ev_t[x] > 1451606400.0:
                conf = 0.1501310
                cha1='EHZ'
                cha2='EHE'
                cha3='EHN'
                r=rs6
                
            if sta == 'LB04':
                conf = 0.0413037
                cha1='HHZ'
                cha2='HHE'
                cha3='HHN'
                r=r4
            if sta == 'LB05':
                conf = 0.3399529
                cha1='HHZ'
                cha2='HHE'
                cha3='HHN'
                r=r5
            if sta == 'LB06':
                conf = 0.00254462
                cha1='HHZ'
                cha2='HHE'
                cha3='HHN'
                r=r6
            if sta == 'LS01':
                conf = 0.10720354
                cha1='EHZ'
                cha2='EHE'
                cha3='EHN'
                r=rs1
            if sta == 'LS02':
                conf = 1.221094559
                cha1='EHZ'
                cha2='EHE'
                cha3='EHN'
                r=rs2
            if sta == 'LS03':
                conf = 0.03363487
                cha1='EHZ'
                cha2='EHE'
                cha3='EHN'
                r=rs3
            if sta == 'LS04':
                conf = 0.0537372
                cha1='EHZ'
                cha2='EHE'
                cha3='EHN'
                r=rs4
            
            #Get each chanel data
            st_z = client.get_waveforms(net, sta, '', cha1, t1 , t2)
            trsz= st_z[0]
            trsz.detrend(type='linear')
            trsz.detrend(type='demean')
            trsz.filter(type='bandpass',freqmin=0.1, freqmax=10)
            st_acz = calibrate1(trsz)
                
            st_e = client.get_waveforms(net, sta, '', cha2, t1 , t2)
            trse= st_e[0]
            trse.detrend(type='linear')
            trse.detrend(type='demean')
            trse.filter(type='bandpass',freqmin=0.1, freqmax=10)
            st_ace = calibrate1(trse)
            
            st_n = client.get_waveforms(net, sta, '', cha3, t1 , t2)
            trsn= st_n[0]
            trsn.detrend(type='linear')
            trsn.detrend(type='demean')
            trsn.filter(type='bandpass',freqmin=0.1, freqmax=10)
            st_acn = calibrate1(trsn)
            
            #energy calculation
            B=2*pi*rhoE*cE*(1/A)
            dl=len(st_acz[0].data)
            p = np.linspace(0,dl/100, num=dl)
            
            y= np.sqrt(st_acz[0].data**2 + st_ace[0].data**2 + st_acn[0].data**2)
            y2= y**2
            
            
            y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
            EI = y_int2[-1]
            EE1= B*(r*r)*EI*conf
            
            
            y_int3 = B*(r*r)*conf*integrate.cumtrapz(y2, p, initial=0)
            
            M1= ((2/3)*np.log10(EE1) )- 2.9
            
            
            
            plt.figure()
            plt.plot(st_acz[0].data)
            plt.title(EE1/1000000)
            plt.xlabel('{}'.format(sta))
            
            enz  = np.lib.pad(enz, ((0,1),(0,0)), 'constant', constant_values=(0))
            enz[numz][0]=EE1/1000000
            enz[numz][1]=M1 
            
            
            
            numz+=1
            
        except:
            print('something broke')

#%%
print(enz)

std_e=np.std(enz[:,0])
std_m=np.std(enz[:,1])

print("std_e =", std_e)
print("std_m =", std_m)     
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      

