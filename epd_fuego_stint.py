#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 16 10:49:37 2018

@author: william
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

#%% Read in catalogue from csv file

cat = genfromtxt("/Users/william/Documents/Fuego_catalogue/Fuego_Final_scan.csv", delimiter=',',skip_header=1)


stations = genfromtxt("/Users/william/Documents/Fuego_catalogue/active_stations_fuego_all_instruments.csv", delimiter=',',skip_header=1)

#%% extract columns of information from catalogue

evt=cat[:,0]
etype=cat[:,7]

days_c = stations[:,0]

inf_num=stations[:,2]
seis_num=stations[:,3]
st_num=stations[:,4]

dlen = 60*60*24

catalogue_start = 1520985600

start_day = 1538352000
days = 92


#%% events per day (epd) for catalogue


epd=np.zeros(shape=(0,11)) 

ev_n = 0
exp_n = 0
gas_n =0
ash_n = 0

inf_t = 0
seis_t = 0

day = 0


day0 = start_day


for y in range(0,days):
    start = day0 + y*dlen
    end = day0 + (y+1)*dlen
    
    epd = np.lib.pad(epd, ((0,1),(0,0)), 'constant', constant_values=(0))
    epd[y][0] = y
    epd[y][1] = start
    
    epd[y][2] = st_num[y]
    epd[y][3] = inf_num[y]
    epd[y][4] = seis_num[y]
    
    
    for x in range(0,len(cat)):
        
        if start < evt[x] < end:
            ev_n += 1
            
            if etype[x] == 3 or etype[x] == 4 or etype[x] == 5 or etype[x] == 6 or etype[x] == 7 or etype[x] == 30 or etype[x] == 40 or etype[x] == 50 or etype[x] == 60:    
                exp_n +=1
            if etype[x] == 3 or etype[x] == 30 or etype[x] == 5 or etype[x] == 50:    
                gas_n +=1
            if etype[x] == 4 or etype[x] == 40 or etype[x] == 6 or etype[x] == 60:    
                ash_n +=1
            if etype[x] == 1 or  etype[x] == 10 :    
                inf_t +=1
            if etype[x] == 2  or etype[x] == 20:    
                seis_t +=1
            
    epd[y][5] = ev_n
    epd[y][6] = exp_n
    epd[y][7] = gas_n
    epd[y][8] = ash_n
    epd[y][9] = inf_t
    epd[y][10] = seis_t
    
    ev_n = 0
    exp_n =0
    gas_n =0
    ash_n = 0
    
    inf_t = 0
    seis_t = 0
    
    

    
#np.savetxt("/Users/william/Documents/Fuego_catalogue/output_data/events_per_day.csv", epd ,delimiter=",",header="day,time,staions,inf_S,seis_S,detections,exp,gas,ash,inf_t,seis_t")



#%%    plot events per day and save csv file with data
    



plt.figure()
plt.plot(epd[:,0], epd[:,5],'r')
plt.xlim([0,days])
#plt.ylim([0,500])
plt.xlabel('Days')
plt.ylabel('Events')
plt.title("Total Events")

plt.figure()
plt.plot(epd[:,0], epd[:,6],'k')
plt.xlim([0,days])
#plt.ylim([0,500])
plt.xlabel('Days')
plt.ylabel('Explosions')
plt.title("Oct 1st '18 - Jan 1st '19 Explosions")

#plt.figure()
#plt.plot(epd[:,0], epd[:,4],'b')
#plt.plot(epd[:,0], epd[:,3],'g')
#plt.plot(epd[:,0], epd[:,2],'k')
#plt.xlim([0,days])
##plt.ylim([0,30])
#plt.xlabel('Days')
#plt.ylabel('Active stations')
#plt.title("Sation Availability")
#plt.legend(['Seismic','Infrasound','Total'])
#
#
#fig, ax1 = plt.subplots()
#ax1.set_title('Events Detected per Day - Starting 15th March 2018')
#ax1.plot(epd[:,0],epd[:,2],'k--')
#ax1.plot(epd[:,0], epd[:,3],'g--')
#ax1.plot(epd[:,0], epd[:,4],'b--')
#ax1.set_ylabel('Stations')
#ax1.set_xlabel('Time (Days)')
##ax1.set_ylim([0,30])
#ax1.set_xlim([0,days])
#ax2 = ax1.twinx()
#ax2.plot(epd[:,0],epd[:,5],'r')
#ax2.set_ylabel('Number of Explosions')
#ax2.set_ylim([0,500])



fig, ax1 = plt.subplots()
ax1.set_title('Gas vs Ash')
ax1.plot(epd[:,0],epd[:,7],'k')
ax1.set_ylabel('Gas')
ax1.set_xlabel('Time (Days)')
#ax1.set_ylim([0,350])
ax1.set_xlim([0,days])
ax2 = ax1.twinx()
ax2.plot(epd[:,0],epd[:,8],'r')
ax2.set_ylabel('Ash')
#ax2.set_ylim([0,350])
ax1.legend(['Gas'],bbox_to_anchor=(0.8, 1), loc='upper left')
ax2.legend(['Ash'],bbox_to_anchor=(0.8, 0.9), loc='upper left')


plt.figure()
plt.plot(epd[:,0], epd[:,9],'k')
plt.xlim([0,days])
#plt.ylim([0,500])
plt.xlabel('Days')
plt.ylabel('Inf trem')
plt.title("Total Inf tremor events")


plt.figure()
plt.plot(epd[:,0], epd[:,10],'k')
plt.xlim([0,days])
#plt.ylim([0,500])
plt.xlabel('Days')
plt.ylabel('Seis trem')
plt.title("Total Seis tremor events")







#%% Events per week (epw)

epw=np.zeros(shape=(0,11))
week=1
event_c =0
exp_c = 0
gas_c = 0
ash_c = 0

st_c = 0
inf_c = 0
seis_c = 0

inft_c =0
seist_c =0


day_c =0
w_en = 0

for x in range(0,len(epd)): # loop over days in "events per day"
    day_c += 1 # count number of days
    event_c += epd[x][5] # count number of events in week (1 day at a time)
    exp_c += epd[x][6]
    gas_c += epd[x][7]
    ash_c += epd[x][8]
    inft_c += epd[x][9]
    seist_c += epd[x][10]
    
    st_c += st_num[x]
    inf_c += inf_num[x]
    seis_c += seis_num[x]
    
    
    if day_c == 7: # save number counted once 7 days is reached and reset counters for next week
        epw = np.lib.pad(epw, ((0,1),(0,0)), 'constant', constant_values=(0))
        epw[week-1][0]=week
        epw[week-1][1] = 1521072000 + (week-1)*7*24*60*60
        
        epw[week-1][2]=st_c/7
        epw[week-1][3]=inf_c/7
        epw[week-1][4]=seis_c/7
        
        epw[week-1][5]=event_c
        epw[week-1][6]=exp_c
        epw[week-1][7]=gas_c
        epw[week-1][8]=ash_c
        epw[week-1][9]=inft_c
        epw[week-1][10]=seist_c
        
        
        week += 1
        event_c = 0
        exp_c = 0
        gas_c = 0
        ash_c = 0
        day_c =0
        inft_c =0
        seist_c =0
        
        st_c = 0
        inf_c = 0
        seis_c = 0
    


    
#        
#%% plot events per week and save csv file with data




plt.figure()
plt.plot(epw[:,0], epw[:,5],'r')
plt.xlim([0,111])
plt.ylim([0,2600])
plt.xlabel('Weeks')
plt.ylabel('Events')
plt.title("Total Events")

plt.figure()
plt.plot(epw[:,0], epw[:,6],'k')
plt.xlim([0,111])
plt.ylim([0,2600])
plt.xlabel('Weeks')
plt.ylabel('Explosions')
plt.title("Total Explosions")


plt.figure()
plt.plot(epw[:,0], epw[:,4],'b')
plt.plot(epw[:,0], epw[:,3],'g')
plt.plot(epw[:,0], epw[:,2],'k')
plt.xlim([0,111])
plt.ylim([0,30])
plt.xlabel('Weeks')
plt.ylabel('Active stations')
plt.title("Sation Availability")
plt.legend(['Seismic','Infrasound','Total'])

fig, ax1 = plt.subplots()
ax1.set_title('Events Detected per week - Starting 15th March 2018')
ax1.plot(epw[:,0], epw[:,2],'k--')
ax1.plot(epw[:,0], epw[:,3],'g--')
ax1.plot(epw[:,0], epw[:,4],'b--')
ax1.set_xlabel('Time (Weeks)')
ax1.set_ylabel('Stations')
ax1.set_ylim([0,30])
ax1.set_xlim([0,111])
ax2 = ax1.twinx()
ax2.plot(epw[:,0],epw[:,5],'r')
ax2.set_ylabel('Number of Events')
ax2.set_ylim([0,2500])


fig, ax1 = plt.subplots()
ax1.set_title('Ash vs Gas')
ax1.plot(epw[:,0], epw[:,7],'k')
ax1.set_xlabel('week')
ax1.set_ylabel('Gas')
ax1.set_ylim([0,1500])
ax1.set_xlim([0,111])
ax2 = ax1.twinx()
ax2.plot(epw[:,0],epw[:,8],'r')
ax2.set_ylabel('Ash')
ax2.set_ylim([0,1500])
ax1.legend(['Gas'],bbox_to_anchor=(0.8, 1), loc='upper left')
ax2.legend(['Ash'],bbox_to_anchor=(0.8, 0.9), loc='upper left')

plt.figure()
plt.plot(epw[:,0], epw[:,9],'k')
plt.xlim([0,111])
#plt.ylim([0,500])
plt.xlabel('Week')
plt.ylabel('Inf trem')
plt.title("Total Inf tremor events")


plt.figure()
plt.plot(epw[:,0], epw[:,10],'k')
plt.xlim([0,111])
#plt.ylim([0,500])
plt.xlabel('Week')
plt.ylabel('Seis trem')
plt.title("Total Seis tremor events")


plt.figure()
plt.plot(epw[:,0], epw[:,10],'k')
plt.plot(epw[:,0], epw[:,9],'r')
plt.xlim([0,111])
#plt.ylim([0,500])
plt.xlabel('Week')
plt.ylabel('Tremor events')
plt.title("Total Tremor events")
plt.legend(['Seis','Inf'])

#np.savetxt("/Users/william/Documents/Fuego_catalogue/output_data/events_per_week.csv", epw ,delimiter=",",header="week,time,staions,inf_S,seis_S,detections,explosions,gas,ash,inf_t,seis_t")













