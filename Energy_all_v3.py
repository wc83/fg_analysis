#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 16:42:44 2019

@author: william
"""
kill the code 
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




#%% get all 2018 data

stream1 = read("/Volumes/Wills_Disc/2018_seismic_data/STG8/HHZ.D/GI.STG8..HHZ.D.2018.*")
stream2 = read("/Volumes/Wills_Disc/2018_seismic_data/STG8/HHE.D/GI.STG8..HHE.D.2018.*")
stream3 = read("/Volumes/Wills_Disc/2018_seismic_data/STG8/HHN.D/GI.STG8..HHN.D.2018.*")
print('2018 read in') 

#%%
Energy_trip = np.zeros(shape=( len(cat),3))
Energy_trip[:,0]=ev_t

#%% 2018 data

for x in range(0,len(stream1)):
    tr_z = stream1[x]
    tr_n = stream2[x]
    tr_e = stream3[x]
    
    s_start = stream1[x].stats.starttime.timestamp
    s_end = stream1[x].stats.endtime.timestamp
    
    for x in range(0,len(cat)):
        if s_start +10 < ev_t[x] < s_end -70 :
            t1=UTCDateTime(ev_t[x] - 10)
            t2=t1 + 80
            trs_z=tr_z.slice(starttime=t1, endtime=t2)
            trs_n=tr_n.slice(starttime=t1, endtime=t2)
            trs_e=tr_e.slice(starttime=t1, endtime=t2)
            
            
            trs_z.detrend(type='linear')
            trs_z.detrend(type='demean')
            trs_z.filter(type='bandpass',freqmin=0.1, freqmax=10)
            
            trs_n.detrend(type='linear')
            trs_n.detrend(type='demean')
            trs_n.filter(type='bandpass',freqmin=0.1, freqmax=10)
            
            trs_e.detrend(type='linear')
            trs_e.detrend(type='demean')
            trs_e.filter(type='bandpass',freqmin=0.1, freqmax=10)

            r=r3
            st_ac_z = calibrate1(trs_z)
            st_ac_n = calibrate1(trs_n)
            st_ac_e = calibrate1(trs_e)
            
            B=2*pi*rhoE*cE*(1/A)
            dl=len(st_ac_z[0].data)
            p = np.linspace(0,dl/100, num=dl)
            
            y= np.sqrt(st_ac_z[0].data**2 + st_ac_n[0].data**2 + st_ac_e[0].data**2 )
                        
            y2= np.square(y)
            
            
            y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
            EI = y_int2[-1]
            EE1= B*(r*r)*EI
            
            y_int3 = B*(r*r)*integrate.cumtrapz(y2, p, initial=0)
            
            E25=EE1*0.025
            E95=EE1*0.95
            
            n25, ind25 = find_nearest(y_int3, E25)
            n95, ind95 = find_nearest(y_int3, E95)
            duration = (ind95/100) - (ind25/100)
            
            Energy_trip[x,1]=EE1
            Energy_trip[x,2]=duration
            print(duration)
            if EE1 < 500:
                print(t1, EE1)
            
print('2018 done')            
        

#%% get all 2014-2017 data
net = 'Z4' 
st_n = Stream()
st_e = Stream() 
client = Client('138.253.112.23', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23

for x in range(0,len(cat)):
#    print(x)
    if    year[x] < 2018 :
#        print(x)
        if  Energy_trip[x,1] ==0 :
#        if 1444215750 < ev_t[x] < 1445353000:

    
#            print(sta)
            t1=UTCDateTime(ev_t[x] -10)
            t2=t1+80
            try:
                # chose station 
                sta ='LS06' 
                cha1='EHZ'
                cha2='EHE'
                cha3='EHN'
                r=rs6
                
                # converstion factors to LB03 
                if sta == 'LB03':
                    conf = 1
                if sta == 'LB01' and ev_t[x] < 1451606400.0:
                    conf = 0.4394913
                if sta == 'LB02' and ev_t[x] < 1451606400.0:
                    conf = 0.17515
                if sta == 'LS05' and ev_t[x] < 1451606400.0:
                    conf = 0.6900589
                if sta == 'LS06' and ev_t[x] < 1451606400.0:
                    conf = 0.1898134 
                    
                if sta == 'LB01' and ev_t[x] > 1451606400.0:
                    conf = 0.308727
                if sta == 'LB02' and ev_t[x] > 1451606400.0:
                    conf = 0.13611627
                if sta == 'LS05' and ev_t[x] > 1451606400.0:
                    conf = 1.78750 
                if sta == 'LS06' and ev_t[x] > 1451606400.0:
                    conf = 0.1501310
                    
                if sta == 'LB04':
                    conf = 0.0413037
                if sta == 'LB05':
                    conf = 0.3399529
                if sta == 'LB06':
                    conf = 0.00254462
                if sta == 'LS01':
                    conf = 0.10720354
                if sta == 'LS02':
                    conf = 1.221094559
                if sta == 'LS03':
                    conf = 0.03363487
                if sta == 'LS04':
                    conf = 0.0537372
                
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
                
                E25=EE1*0.025
                E95=EE1*0.95
                
                n25, ind25 = find_nearest(y_int3, E25)
                n95, ind95 = find_nearest(y_int3, E95)
                duration = (ind95/100) - (ind25/100)
                
                Energy_trip[x,1]=EE1
                Energy_trip[x,2]=duration
#                print(duration)
                
                
                
                if EE1 < 250:
                    print(UTCDateTime(ev_t[x]),EE1)
                    break
                 
            except:
                print('could not get data for ', t1, sta)





#%%  
#E1 = sum(Energy_trip[:,1])
#EE = sum(Energy_trip[:,4])
#Er=EE/E1
#print(EE/E1) 
#for x in range(0,len(Energy_trip)):
#     if Energy_trip[x,2] == 0:
#         Energy_trip[x,1] = cat[x,23]
#         Energy_trip[x,4] = cat[x,23]*Er


#%%
  
for x in range(0, len(cat)):
    if year[x]==2018:
        Energy_trip[x,1] = cat[x,23]
              
 #%%               

for x in range(0,len(Energy_trip)):
    if Energy_trip[x,2] ==0:# and Energy_trip[x,0] < 1451606400.0:
        print(UTCDateTime(Energy_trip[x,0]))
#        Energy_trip[x,2] = cat[x,24]
        
        
#%%
        

        
new_cat=np.zeros(shape=(len(cat),25))
new_cat[:,0]=cat[:,0]
new_cat[:,1]=cat[:,1]
new_cat[:,2]=cat[:,2]
new_cat[:,3]=cat[:,3]
new_cat[:,4]=cat[:,4]
new_cat[:,5]=cat[:,5]
new_cat[:,6]=cat[:,6]
new_cat[:,7]=cat[:,7]
new_cat[:,8]=cat[:,8]
new_cat[:,9]=cat[:,9]
new_cat[:,10]=cat[:,10]
new_cat[:,11]=cat[:,11]
new_cat[:,12]=cat[:,12]
new_cat[:,13]=cat[:,13]
new_cat[:,14]=cat[:,14]
new_cat[:,15]=cat[:,15]
new_cat[:,16]=cat[:,16]
new_cat[:,17]=cat[:,17]
new_cat[:,18]=cat[:,18]
new_cat[:,19]=cat[:,19]
new_cat[:,20]=cat[:,20]
new_cat[:,21]=cat[:,21]
new_cat[:,22]=cat[:,22]
#new_cat[:,23]=Energy_trip[:,1]
new_cat[:,23]=cat[:,23]
#new_cat[:,24]=cat[:,24]
new_cat[:,24]=Energy_trip[:,2]
        
 



#%%
        
#np.savetxt("/Users/william/Documents/scanner/all_stations/Final_Catalogue_2014_2018_2.csv", new_cat,delimiter=",",header="Time_stamp,Day_number,Day_of_year,Year,Month,Day,Hour,Min,Sec,Milisec,LBO1,LB02,LB03,LB04,LB05,LB06,LSO1,LS02,LS03,LS04,LS05,LS06,Trust,Energy,Duration")

#%%  
#for x in range(0,len(new_cat)):
#    if new_cat[x,23] ==0 or new_cat[x,24] ==0:
#        print(new_cat[x,0])
#  
  
  
  
  
for x in range(0,len(new_cat)):
    if new_cat[x,23] < 1000:
        print((new_cat[x,0]),new_cat[x,23])
  
  
  
  
  
  
  