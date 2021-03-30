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
import datetime as dt
import matplotlib.dates as md

#%% Read in catalogue from csv file

cat = genfromtxt("/Users/william/Documents/Fuego_catalogue/Final_catalogues/Fuego_explosions.csv", delimiter=',',skip_header=1)

i_trem = genfromtxt("/Users/william/Documents/Fuego_catalogue/Final_catalogues/Fuego_trem_Inf_final.csv", delimiter=',',skip_header=1)
s_trem = genfromtxt("/Users/william/Documents/Fuego_catalogue/Final_catalogues/Fuego_trem_Seis_final.csv", delimiter=',',skip_header=1)


stations = genfromtxt("/Users/william/Documents/Fuego_catalogue/Final_catalogues/Fuego_Final_scan_stations.csv", delimiter=',',skip_header=1)

#%% extract columns of information from catalogue

evt=cat[:,0]

evit=i_trem[:,0]
evst=s_trem[:,0]

evid=i_trem[:,1]
evsd=s_trem[:,1]

etype=cat[:,8]

days_c = stations[:,0]

inf_num=stations[:,2]
seis_num=stations[:,3]
st_num=stations[:,4]

dlen = 60*60*24

#%% events per day (epd) for catalogue


epd=np.zeros(shape=(0,15)) 

ev_n = 0
exp_n = 0
gas_n =0
ash_n = 0
gca_n = 0
no_gca_n = 0

inf_t = 0
seis_t = 0

inf_d = 0
seis_d = 0

day = 0


day0 = 1520985600

for y in range(0,len(days_c)):
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
            
            if etype[x] == 1 or etype[x] == 2 or etype[x] == 3 or etype[x] == 4 or etype[x] == 5:    
                exp_n +=1
            if etype[x] == 1 or etype[x] == 2 :    
                gas_n +=1
            if etype[x] == 3 or etype[x] == 4 :    
                ash_n +=1
            if etype[x] == 2 or etype[x] == 4 or etype[x] == 5 :    
                gca_n +=1
            if etype[x] == 1 or etype[x] == 3 :    
                no_gca_n +=1
                
    for x in range(0,len(i_trem)):
        if start < evit[x] < end:
            ev_n += 1
            inf_t += 1
            inf_d += evid[x]/(60*60)
            
    for x in range(0,len(s_trem)):
        if start < evst[x] < end:
            ev_n += 1
            seis_t += 1
            seis_d += evsd[x]/(60*60)
                


            
    epd[y][5] = ev_n
    epd[y][6] = exp_n
    epd[y][7] = gas_n
    epd[y][8] = ash_n
    epd[y][9] = gca_n
    epd[y][10] = no_gca_n
    epd[y][11] = inf_t
    epd[y][12] = seis_t
    epd[y][13] = inf_d
    epd[y][14] = seis_d
    
    ev_n = 0
    exp_n =0
    gas_n =0
    ash_n = 0
    gca_n = 0
    no_gca_n = 0
    inf_d = 0
    seis_d = 0
    
    inf_t = 0
    seis_t = 0
    
    

    
#np.savetxt("/Users/william/Documents/Fuego_catalogue/output_data/events_per_day.csv", epd ,delimiter=",",header="day,time,staions,inf_S,seis_S,detections,exp,gas,ash,inf_t,seis_t")



#%%    plot events per day and save csv file with data
    



#plt.figure()
#dates=[dt.datetime.fromtimestamp(ts) for ts in  epd[:,1]]
#datenums=md.date2num(dates)
#ax=plt.gca()
#xfmt = md.DateFormatter('%d-%m-%Y')
#plt.xticks( rotation=45 )
#ax.xaxis.set_major_formatter(xfmt)
#plt.plot(datenums, epd[:,5],'r')
##plt.xlim([0,782])
##plt.ylim([0,500])
#plt.xlabel('Date')
#plt.ylabel('Events')
#plt.title("Total Events")

plt.figure()
dates=[dt.datetime.fromtimestamp(ts) for ts in  epd[:,1]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums, epd[:,6],'k')
#plt.xlim([0,782])
#plt.ylim([0,500])
plt.xlabel('Date')
plt.ylabel('Explosions')
plt.title("Total Explosions")


plt.figure()
dates=[dt.datetime.fromtimestamp(ts) for ts in  epd[:,1]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums, epd[:,4],'b')
plt.plot(datenums, epd[:,3],'g')
plt.plot(datenums, epd[:,2],'k')
#plt.xlim([0,782])
plt.ylim([0,30])
plt.xlabel('Date')
plt.ylabel('Active stations')
plt.title("Sation Availability")
plt.legend(['Seismic','Infrasound','Total'])


fig, ax1 = plt.subplots()
dates=[dt.datetime.fromtimestamp(ts) for ts in  epd[:,1]]
datenums=md.date2num(dates)
ax1=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax1.xaxis.set_major_formatter(xfmt)

ax1.set_title('Events Detected per Day - Starting 15th March 2018')
ax1.plot(datenums,epd[:,2],'k--')
ax1.plot(datenums, epd[:,3],'g--')
ax1.plot(datenums, epd[:,4],'b--')
ax1.set_ylabel('Stations')
ax1.set_xlabel('Date')
ax1.set_ylim([0,30])
#ax1.set_xlim([0,800])
ax2 = ax1.twinx()
ax2.plot(datenums,epd[:,5],'r')
ax2.set_ylabel('Number of Explosions')
#ax2.set_ylim([0,500])



fig, ax1 = plt.subplots()
dates=[dt.datetime.fromtimestamp(ts) for ts in  epd[:,1]]
datenums=md.date2num(dates)
ax1=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax1.xaxis.set_major_formatter(xfmt)

ax1.set_title('Gas vs Ash')
ax1.plot(datenums,epd[:,7],'k')
ax1.set_ylabel('Gas')
ax1.set_xlabel('Date')
ax1.set_ylim([0,350])
#ax1.set_xlim([0,800])
ax2 = ax1.twinx()
ax2.plot(datenums,epd[:,8],'r')
ax2.set_ylabel('Ash')
#ax2.set_ylim([0,350])
ax1.legend(['Gas'],bbox_to_anchor=(0.8, 1), loc='upper left')
ax2.legend(['Ash'],bbox_to_anchor=(0.8, 0.9), loc='upper left')


fig, ax1 = plt.subplots()
dates=[dt.datetime.fromtimestamp(ts) for ts in  epd[:,1]]
datenums=md.date2num(dates)
ax1=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax1.xaxis.set_major_formatter(xfmt)

ax1.set_title('GCA vs NO GCA')
ax1.plot(datenums,epd[:,9],'k')
ax1.set_ylabel('GCA')
ax1.set_xlabel('Date')
#ax1.set_ylim([0,350])
#ax1.set_xlim([0,800])
ax2 = ax1.twinx()
ax2.plot(datenums,epd[:,10],'r')
ax2.set_ylabel('NO GCA')
#ax2.set_ylim([0,350])
ax1.legend(['GCA'],bbox_to_anchor=(0.75, 1), loc='upper left')
ax2.legend(['NO GCA'],bbox_to_anchor=(0.75, 0.9), loc='upper left')

plt.figure()
dates=[dt.datetime.fromtimestamp(ts) for ts in  epd[:,1]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums, epd[:,11],'k')
#plt.xlim([0,782])
#plt.ylim([0,500])
plt.xlabel('Date')
plt.ylabel('Inf trem')
plt.title("Total Inf tremor events")


plt.figure()
dates=[dt.datetime.fromtimestamp(ts) for ts in  epd[:,1]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums, epd[:,12],'k')
#plt.xlim([0,782])
#plt.ylim([0,500])
plt.xlabel('Date')
plt.ylabel('Seis trem')
plt.title("Total Seis tremor events")



plt.figure()
dates=[dt.datetime.fromtimestamp(ts) for ts in  epd[:,1]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums, epd[:,13],'k')
#plt.xlim([0,782])
#plt.ylim([0,500])
plt.xlabel('Date')
plt.ylabel('Inf trem [hours]')
plt.title("Daily Inf tremor duration")


plt.figure()
dates=[dt.datetime.fromtimestamp(ts) for ts in  epd[:,1]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums, epd[:,14],'k')
#plt.xlim([0,782])
#plt.ylim([0,500])
plt.xlabel('Date')
plt.ylabel('Seis trem [hours]')
plt.title("Daily Seis tremor duration")



#%% Events per week (epw)

epw=np.zeros(shape=(0,15))
week=1
event_c =0
exp_c = 0
gas_c = 0
ash_c = 0
gca_c = 0
no_gca_c = 0
inft_c =0
seist_c =0
inft_cd =0
seist_cd =0

st_c = 0
inf_c = 0
seis_c = 0




day_c =0
w_en = 0

for x in range(0,len(epd)): # loop over days in "events per day"
    day_c += 1 # count number of days
    event_c += epd[x][5] # count number of events in week (1 day at a time)
    exp_c += epd[x][6]
    gas_c += epd[x][7]
    ash_c += epd[x][8]
    gca_c += epd[x][9]
    no_gca_c += epd[x][10]
    
    inft_c += epd[x][11]
    seist_c += epd[x][12]
    inft_cd += epd[x][13]
    seist_cd += epd[x][14]
    
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
        epw[week-1][9]=gca_c
        epw[week-1][10]=no_gca_c
        
        epw[week-1][11]=inft_c
        epw[week-1][12]=seist_c
        epw[week-1][13]=inft_cd
        epw[week-1][14]=seist_cd
        
        
        week += 1
        event_c = 0
        exp_c = 0
        gas_c = 0
        ash_c = 0
        day_c =0
        inft_c =0
        seist_c =0
        gca_c = 0
        no_gca_c = 0
        inft_cd =0
        seist_cd =0
        
        
        st_c = 0
        inf_c = 0
        seis_c = 0
    


    
#        
#%% plot events per week and save csv file with data




#plt.figure()
#dates=[dt.datetime.fromtimestamp(ts) for ts in  epw[:,1]]
#datenums=md.date2num(dates)
#ax=plt.gca()
#xfmt = md.DateFormatter('%d-%m-%Y')
#plt.xticks( rotation=45 )
#ax.xaxis.set_major_formatter(xfmt)
#plt.plot(datenums, epw[:,5],'k')
##plt.xlim([0,111])
##plt.ylim([0,2600])
#plt.xlabel('Date')
#plt.ylabel('Events')
#plt.title("Weekly Events")

plt.figure(figsize=(12,3))
dates=[dt.datetime.fromtimestamp(ts) for ts in  epw[:,1]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums, epw[:,6],'k')
#plt.xlim([0,111])
#plt.ylim([0,2600])
plt.xlabel('Date')
plt.ylabel('Explosions')
plt.title("Weekly Explosions")


plt.figure()
dates=[dt.datetime.fromtimestamp(ts) for ts in  epw[:,1]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums, epw[:,4],'r')
plt.plot(datenums, epw[:,3],'b')
plt.plot(datenums, epw[:,2],'k')
#plt.xlim([0,111])
plt.ylim([0,30])
plt.xlabel('Date')
plt.ylabel('Active stations')
plt.title("Weekly Sation Availability")
plt.legend(['Seismic','Infrasound','Total'])

fig, ax1 = plt.subplots()
dates=[dt.datetime.fromtimestamp(ts) for ts in  epw[:,1]]
datenums=md.date2num(dates)
ax1=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax1.xaxis.set_major_formatter(xfmt)

ax1.set_title('Events Detected per week - Starting 15th March 2018')
ax1.plot(datenums, epw[:,2],'k--')
ax1.plot(datenums, epw[:,3],'g--')
ax1.plot(datenums, epw[:,4],'b--')
ax1.set_xlabel('Date')
ax1.set_ylabel('Stations')
ax1.set_ylim([0,30])
#ax1.set_xlim([0,111])
ax2 = ax1.twinx()
ax2.plot(datenums,epw[:,5],'r')
ax2.set_ylabel('Number of Events')
#ax2.set_ylim([0,2500])


fig, ax1 = plt.subplots()
dates=[dt.datetime.fromtimestamp(ts) for ts in  epw[:,1]]
datenums=md.date2num(dates)
ax1=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax1.xaxis.set_major_formatter(xfmt)

ax1.set_title('Weekly Ash vs Gas')
ax1.plot(datenums, epw[:,7],'k')
ax1.set_xlabel('Date')
ax1.set_ylabel('Gas')
ax1.set_ylim([0,1500])
#ax1.set_xlim([0,111])
ax2 = ax1.twinx()
ax2.plot(datenums,epw[:,8],'r')
ax2.set_ylabel('Ash')
#ax2.set_ylim([0,1500])
ax1.legend(['Gas'],bbox_to_anchor=(0.8, 1), loc='upper left')
ax2.legend(['Ash'],bbox_to_anchor=(0.8, 0.9), loc='upper left')



fig, ax1 = plt.subplots()
dates=[dt.datetime.fromtimestamp(ts) for ts in  epw[:,1]]
datenums=md.date2num(dates)
ax1=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax1.xaxis.set_major_formatter(xfmt)

ax1.set_title('Weekly GCA vs NO GCA')
ax1.plot(datenums, epw[:,9],'k')
ax1.set_xlabel('Date')
ax1.set_ylabel('GCA')
ax1.set_ylim([0,1500])
#ax1.set_xlim([0,111])
ax2 = ax1.twinx()
ax2.plot(datenums,epw[:,10],'r')
ax2.set_ylabel('NO GCA')
#ax2.set_ylim([0,1500])
ax1.legend(['GCA'],bbox_to_anchor=(0.75, 1), loc='upper left')
ax2.legend(['NO GCA'],bbox_to_anchor=(0.75, 0.9), loc='upper left')

plt.figure()
dates=[dt.datetime.fromtimestamp(ts) for ts in  epw[:,1]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums, epw[:,11],'b')
#plt.xlim([0,111])
#plt.ylim([0,500])
plt.xlabel('Date')
plt.ylabel('Inf trem')
plt.title("Weekly Total Inf tremor events")


plt.figure()
dates=[dt.datetime.fromtimestamp(ts) for ts in  epw[:,1]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums, epw[:,12],'r')
#plt.xlim([0,111])
#plt.ylim([0,500])
plt.xlabel('Date')
plt.ylabel('Seis trem')
plt.title("Weekly Total Seis tremor events")


plt.figure()
dates=[dt.datetime.fromtimestamp(ts) for ts in  epw[:,1]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums, epw[:,12],'r')
plt.plot(datenums, epw[:,11],'b')
#plt.xlim([0,111])
#plt.ylim([0,500])
plt.xlabel('Date')
plt.ylabel('Tremor events')
plt.title("Weekly Total Tremor events")
plt.legend(['Seis','Inf'])

plt.figure()
dates=[dt.datetime.fromtimestamp(ts) for ts in  epw[:,1]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums, epw[:,13],'b')
#plt.xlim([0,111])
#plt.ylim([0,500])
plt.xlabel('Date')
plt.ylabel('Inf trem [hours]')
plt.title("Weekly Total Inf tremor duration")

plt.figure()
dates=[dt.datetime.fromtimestamp(ts) for ts in  epw[:,1]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums, epw[:,14],'r')
#plt.xlim([0,111])
#plt.ylim([0,500])
plt.xlabel('Date')
plt.ylabel('Seis trem [hours]')
plt.title("Weekly Total Seis tremor duration")


plt.figure()
dates=[dt.datetime.fromtimestamp(ts) for ts in  epw[:,1]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums, epw[:,14],'r')
plt.plot(datenums, epw[:,13],'b')
#plt.xlim([0,111])
#plt.ylim([0,500])
plt.xlabel('Date')
plt.ylabel('Tremor duration [hours]')
plt.title("Weekly Total Tremor duration")
plt.legend(['Seis','Inf'])


#np.savetxt("/Users/william/Documents/Fuego_catalogue/output_data/events_per_week.csv", epw ,delimiter=",",header="week,time,staions,inf_S,seis_S,detections,explosions,gas,ash,inf_t,seis_t")













