#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 16:42:44 2019

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


client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
#%% read catalogue of events

cat = genfromtxt("/Users/william/Documents/scanner/all_stations/Full_Catalogue_full_data.csv", delimiter=',',skip_header=1)
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


#%% get all 2018 data
reclen = 512
chunksize = 100000 * reclen # Around 50 MB
st_a = Stream()

# ******THIS still has 31 noise waveforms not removed from manual cleaning********
with io.open("/Users/william/Documents/scanner/output_data/EXP_LB03_2018_V4.mseed", "rb") as fh: 
    while True:
        with io.BytesIO() as buf:
            c = fh.read(chunksize)
            if not c:
                break
            buf.write(c)
            buf.seek(0, 0)
            st = obspy.read(buf)

            for x in range (0,len(st)):
                st[x].detrend(type='linear')
                st[x].detrend(type='demean')
                st[x].filter(type='bandpass',freqmin=0.1, freqmax=10)
                st_a.append(st[x])
print('2018 read in') 
#%% get all 2014-2017 data
st_f = Stream()
saved=np.zeros(shape=(1,1))
n=1

t1=UTCDateTime(1416791642)
t2=t1+60
first = Stream()
first = client.get_waveforms('Z4', 'LB01', '', 'HHZ', t1-20 , t2)

st_f.append(first[0])
saved[0][0]=1416791642


with io.open("/Users/william/Documents/scanner/output_data/m32.mseed", "rb") as fh:
        # just month 2
#with io.open("/Users/william/Documents/scanner/output_data/EXP_all_data_stream_2_month_2.mseed", "rb") as fh:
    while True:
        with io.BytesIO() as buf:
            c = fh.read(chunksize)
            if not c:
                break
            buf.write(c)
            buf.seek(0, 0)
            st = obspy.read(buf)
            
            
#            for x in range (0,len(st)):
#                st_f.append(st[x])
            for x in range(0,len(st)):
                etim=st[x].stats.starttime.timestamp
                near_save,ind=find_nearest(saved,etim )  
                near_cat,indx=find_nearest(cat[:,0],etim )  
                
                if abs(etim-near_save[0]) >  30 and abs(etim-near_cat) < 90 :
                    saved = np.lib.pad(saved, ((0,1),(0,0)), 'constant', constant_values=(0))
                    saved[n][0]=etim
                    n +=1
                    st[x].detrend(type='linear')
                    st[x].detrend(type='demean')
                    st[x].filter(type='bandpass',freqmin=0.1, freqmax=10)
                    st_f.append(st[x])
print('all read in')                
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


#%% make new catalogue


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

#%% get 2018 energy
nu=0
e2018=np.zeros(shape=(2378,3))
for x in range(0,len(st_a)):
    ET = st_a[x].stats.starttime.timestamp
    near_event,ind=find_nearest(cat[:,0],ET )  #find how near an event in Olivers catalogue is to mine          
    if abs(ET-near_event) <  20:
        r=r3
        st_ac = calibrate1(st_a[x])
        
        B=2*pi*rhoE*cE*(1/A)
        dl=len(st_ac[0].data)
        p = np.linspace(0,dl/100, num=dl)
        
        y= st_ac[0].data
        y2= np.square(st_ac[0].data)
        
        
        y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
        EI = y_int2[-1]
        EE= B*(r*r)*EI
        E25=EE*0.025
        E975=EE*0.975
        
        y_int3 = B*(r*r)*integrate.cumtrapz(y2, p, initial=0) # cumulative energy
        
        n25, ind25 = find_nearest(y_int3, E25)
        n975, ind975 = find_nearest(y_int3, E975)
        duration = (ind975/100) - (ind25/100)
        
        
        e2018[nu][0]=EE
        e2018[nu][1]=duration
        e2018[nu][2]=ET
        nu+=1
        
        

#%% get 2014-2017 energy

nur=0
e2014=np.zeros(shape=(16524,3))
for x in range(0,len(st_f)):
    ET1 = st_f[x].stats.starttime.timestamp
    near_event1,ind=find_nearest(cat[:,0],ET1 )  #find how near an event in Olivers catalogue is to mine          
    if abs(ET1-near_event1) <  20:
        
        if st_f[x].stats.station == "LB01":
            r=r1
        if st_f[x].stats.station == "LB02":
            r=r2
        if st_f[x].stats.station == "LB03":
            r=r3
        if st_f[x].stats.station == "LB04":
            r=r4
        if st_f[x].stats.station == "LB05":
            r=r5
        if st_f[x].stats.station == "LB06":
            r=r6
        if st_f[x].stats.station == "LS01":
            r=rs1
        if st_f[x].stats.station == "LS02":
            r=rs2
        if st_f[x].stats.station == "LS03":
            r=rs3
        if st_f[x].stats.station == "LS04":
            r=rs4
        if st_f[x].stats.station == "LS05":
            r=rs5
        if st_f[x].stats.station == "LS06":
            r=rs6


        st_fc = calibrate1(st_f[x])
        
        B=2*pi*rhoE*cE*(1/A)
        dl=len(st_fc[0].data)
        p = np.linspace(0,dl/100, num=dl)
        
        y= st_fc[0].data
        y2= np.square(st_fc[0].data)
        
        
        y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
        EI = y_int2[-1]
        EE= B*(r*r)*EI
        E25=EE*0.025
        E975=EE*0.975
        
        y_int3 = B*(r*r)*integrate.cumtrapz(y2, p, initial=0) # cumulative energy
        
        n25, ind25 = find_nearest(y_int3, E25)
        n975, ind975 = find_nearest(y_int3, E975)
        duration = (ind975/100) - (ind25/100)

        e2014[nur][0]=EE
        e2014[nur][1]=duration
        e2014[nur][2]=ET1
        nur+=1
        
        
        

#%% make new catalogue energy column - match event with energy

for x in range(0,len(new_cat)):
    if new_cat[x,3] < 2018:
        ET2 = new_cat[x,0]
        
        near_2014,ind=find_nearest(e2014[:,2],ET2 ) 
        
        if abs(near_2014 - ET2) < 20:
            new_cat[x,23] = e2014[ind,0]
            new_cat[x,24] = e2014[ind,1]
        else:
#            print(UTCDateTime(ET2), 'no energy found')
            k=4235


    if new_cat[x,3] == 2018:
        ET3 = new_cat[x,0]
        
        near_2018,indx=find_nearest(e2018[:,2],ET3 ) 
        
        if abs(near_2018 - ET3) < 20:
            new_cat[x,23] = e2018[indx,0]
            new_cat[x,24] = e2018[indx,1]
        else:
#            print(UTCDateTime(ET3), 'no energy found')
            k=234143



#%%g et 2014-2017 data for missing events from previous step
print('adding missed')  
for x in range(0,len(new_cat)):
    if new_cat[x,23]==0:
        print(x)
        if new_cat[x,10] == 1:
            r=r1
            sta='LB01'
            cha='HHZ'
        else:
            if new_cat[x,11] == 1:
                r=r2
                sta='LB02'
                cha='HHZ'
            else:
                if new_cat[x,12] == 1:
                    r=r3
                    sta='LB03'
                    cha='HHZ'
                else:
                    if new_cat[x,13] == 1:
                        r=r4
                        sta='LB04'
                        cha='HHZ'
                    else:
                        if new_cat[x,14] == 1:
                           r=r5
                           sta='LB05'
                           cha='HHZ'
                        else:
                            if new_cat[x,15] == 1:
                               r=r6
                               sta='LB06'
                               cha='HHZ'
                            else:
                                if new_cat[x,16] == 1:
                                    r=rs1
                                    sta='LS01'
                                    cha='EHZ'
                                else:
                                    if new_cat[x,17] == 1:
                                        r=rs2
                                        sta='LS02'
                                        cha='EHZ'
                                    else:
                                        if new_cat[x,18] == 1:
                                           r=rs3
                                           sta='LS03'
                                           cha='EHZ'
                                        else:
                                            if new_cat[x,19] == 1:
                                               r=rs4
                                               sta='LS04'
                                               cha='EHZ'
                                            else:
                                                if new_cat[x,20] == 1:
                                                   r=rs5
                                                   sta='LS05'
                                                   cha='EHZ'
                                                else:
                                                    if new_cat[x,21] == 1:
                                                       r=rs6
                                                       sta='LS06'
                                                       cha='EHZ'
                                                    else:
                                                        print(new_cat[x,0])
            
            
        t1 = UTCDateTime(new_cat[x,0]) #the format is year:day_of_the_year:month
        t2 = t1 + 80
        st = Stream()
        st = client.get_waveforms('Z4', sta, '', cha, t1-20 , t2) 
        st.detrend(type='linear')
        st.detrend(type='demean')
        st.filter(type='bandpass',freqmin=0.1, freqmax=10)
        
        tr=st[0].slice(starttime=t1, endtime=t2)
        st_c = calibrate1(tr)
    
        B=2*pi*rhoE*cE*(1/A)
        
        dl=len(st_c[0].data)
        p = np.linspace(0,dl/100, num=dl)
        y= st_c[0].data
        y2= np.square(st_c[0].data)
        
        
        y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
        EI = y_int2[-1]
        EE= B*(r*r)*EI
        E25=EE*0.025
        E975=EE*0.975
        
        y_int3 = B*(r*r)*integrate.cumtrapz(y2, p, initial=0) # cumulative energy
        
        n25, ind25 = find_nearest(y_int3, E25)
        n975, ind975 = find_nearest(y_int3, E975)
        duration = (ind975/100) - (ind25/100)

        new_cat[x,23]= EE
        new_cat[x,24]=duration
        
        
        
        
        if EE > 1e11:
            print('2015 ',EE)




#%%
        
#np.savetxt("/Users/william/Documents/scanner/all_stations/Final_Catalogue_2014_2018.csv", new_cat,delimiter=",",header="Time_stamp,Day_number,Day_of_year,Year,Month,Day,Hour,Min,Sec,Milisec,LBO1,LB02,LB03,LB04,LB05,LB06,LSO1,LS02,LS03,LS04,LS05,LS06,Trust,Energy,Duration")

#%%  
#for x in range(0,len(new_cat)):
#    if new_cat[x,23] ==0 or new_cat[x,24] ==0:
#        print(new_cat[x,0])
#  
  
for x in range(0,len(new_cat)):
    if new_cat[x,23] < 500:
        print(new_cat[x,0], new_cat[x,23])
  
  
  
  
  
  
  
  
  
  