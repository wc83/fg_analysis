#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 16:02:37 2020

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
#%%

cat = genfromtxt("/Users/william/Documents/Fuego_catalogue/Fuego_Final_scan.csv", delimiter=',',skip_header=1)

Energy = genfromtxt("/Users/william/Documents/Fuego_catalogue/Fuego_Energy_Matrix.csv", delimiter=',',skip_header=1)
#%%


En_t = Energy[:,0]
e_type = cat[:,7]

En_inf = Energy[:,1:30]
En_seis = Energy[:,30:38]

med_en = np.zeros(shape =(0,4))
med_ens = np.zeros(shape =(0,4))
#med_en[:,0]=En_t
#med_en[:,5]=e_type


for x in range(0,len(En_inf)):
    if 1538352000 < UTCDateTime(En_t[x]).timestamp < 1546300800:
               
   
        nzi = En_inf[x,:][np.nonzero(En_inf[x,:])]
        
        
        if len(nzi) > 0:
            med_en_ev = np.median(nzi)
            if med_en_ev < 7500000000:
                med_en = np.lib.pad(med_en, ((0,1),(0,0)), 'constant', constant_values=(0))
                med_en[-1,1]=med_en_ev
                med_en[-1,0]=En_t[x]
                med_en[-1,3]=e_type[x]
            else:
                med_en = np.lib.pad(med_en, ((0,1),(0,0)), 'constant', constant_values=(0))
                med_en[-1,1]=10000000
                med_en[-1,0]=En_t[x]
                med_en[-1,3]=e_type[x]
        else:
            med_en = np.lib.pad(med_en, ((0,1),(0,0)), 'constant', constant_values=(0))
            med_en[-1,0]=En_t[x]
            med_en[-1,3]=e_type[x]
            med_en[-1,1]=0
            
            
for x in range(0,len(En_seis)):
    if 1538352000 < UTCDateTime(En_t[x]).timestamp < 1546300800:

        
        nzs = En_seis[x,:][np.nonzero(En_seis[x,:])]
        
        if len(nzs) > 0:
            med_en_ev = np.median(nzs)
            if med_en_ev < 2000000000:
                med_ens = np.lib.pad(med_ens, ((0,1),(0,0)), 'constant', constant_values=(0))
                med_ens[-1,1]=med_en_ev
                med_ens[-1,0]=En_t[x]
                med_ens[-1,3]=e_type[x]
            else:
                med_ens = np.lib.pad(med_ens, ((0,1),(0,0)), 'constant', constant_values=(0))
                med_ens[-1,1]=2000000  
                med_ens[-1,0]=En_t[x]
                med_ens[-1,3]=e_type[x]
        else:
            med_ens = np.lib.pad(med_ens, ((0,1),(0,0)), 'constant', constant_values=(0))
            med_ens[-1,0]=En_t[x]
            med_ens[-1,3]=e_type[x]
            med_ens[-1,1]=0            
    

cum_inf = 0
cum_seis = 0

for b in range(0,len(med_en)):
    cum_inf += med_en[b,1]
    
    med_en[b,2] = cum_inf

    
for b in range(0,len(med_ens)):
    
    cum_seis += med_ens[b,1]
    
    med_ens[b,2] = cum_seis

    
#%%
plt.figure()    
dates=[dt.datetime.fromtimestamp(ts) for ts in med_ens[:,0]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=30 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums,med_ens[:,1],'bx')
plt.title('seismic energy')
plt.xlabel('Date')
plt.ylabel('Event seismic energy [J]')

#
plt.figure()
dates=[dt.datetime.fromtimestamp(ts) for ts in med_en[:,0]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=30 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums,med_en[:,1],'rx')
plt.title('acoustic energy')
plt.xlabel('Date')
plt.ylabel('Event acoustic energy [J]')
#
plt.figure() 
dates=[dt.datetime.fromtimestamp(ts) for ts in med_ens[:,0]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=30 )
ax.xaxis.set_major_formatter(xfmt)   
plt.plot(datenums,med_ens[:,2],'b')
plt.title('cumulative seismic energy')
plt.xlabel('Date')
plt.ylabel('Cumulative seismic energy [J]')
#
plt.figure()
dates=[dt.datetime.fromtimestamp(ts) for ts in med_en[:,0]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=30 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums,med_en[:,2],'r')
plt.title('cumualtive acoustic energy')
plt.xlabel('Date')
plt.ylabel('Cumulative acoustic energy [J]')
#


#%%

e_type_i = med_en[:,3]
e_type_s = med_ens[:,3]

all_exp_en_inf = np.zeros(shape =(0,3))
all_exp_en_seis = np.zeros(shape =(0,3))
gas_exp_en_inf = np.zeros(shape =(0,3))
ash_exp_en_inf = np.zeros(shape =(0,3))
all_exp_en_inf = np.zeros(shape =(0,3))
all_inf_trem_en = np.zeros(shape =(0,3))
all_seis_trem_en = np.zeros(shape =(0,3))
all_SA_trem_en_inf = np.zeros(shape =(0,3))
all_SA_trem_en_seis = np.zeros(shape =(0,3))
all_SA_exp_en_inf = np.zeros(shape =(0,3))
all_SA_exp_en_seis = np.zeros(shape =(0,3))


for x in range(0,len(med_en)):
    
    if e_type_i[x] == 3 or e_type_i[x] == 4 or e_type_i[x] == 5 or e_type_i[x] == 6 or e_type_i[x] == 30 or e_type_i[x] == 40 or e_type_i[x] == 50 or e_type_i[x] == 60:
        all_exp_en_inf = np.lib.pad(all_exp_en_inf, ((0,1),(0,0)), 'constant', constant_values=(0))
        all_exp_en_inf[-1][0]= med_en[x,0]
        all_exp_en_inf[-1][1]= med_en[x,1]
        if len(all_exp_en_inf) == 1:
            all_exp_en_inf[-1][2] = med_en[x,1]
        else:
            all_exp_en_inf[-1][2] = all_exp_en_inf[-2][2] + med_en[x,1]
        
        
    if e_type_i[x] == 3 or e_type_i[x] == 5  or e_type_i[x] == 30  or e_type_i[x] == 50 :
        gas_exp_en_inf = np.lib.pad(gas_exp_en_inf, ((0,1),(0,0)), 'constant', constant_values=(0))
        gas_exp_en_inf[-1][0]= med_en[x,0]
        gas_exp_en_inf[-1][1]= med_en[x,1]
        if len(gas_exp_en_inf) == 1:
            gas_exp_en_inf[-1][2] = med_en[x,1]
        else:
            gas_exp_en_inf[-1][2] = gas_exp_en_inf[-2][2] + med_en[x,1]
        
    if e_type_i[x] == 4 or e_type_i[x] == 6  or e_type_i[x] == 40  or e_type_i[x] == 60 :
        ash_exp_en_inf = np.lib.pad(ash_exp_en_inf, ((0,1),(0,0)), 'constant', constant_values=(0))
        ash_exp_en_inf[-1][0]= med_en[x,0]
        ash_exp_en_inf[-1][1]= med_en[x,1]
        if len(ash_exp_en_inf) == 1:
            ash_exp_en_inf[-1][2] = med_en[x,1]
        else:
            ash_exp_en_inf[-1][2] = ash_exp_en_inf[-2][2] + med_en[x,1]
               
    if e_type_i[x] == 1 or e_type_i[x] == 10  :
        all_inf_trem_en = np.lib.pad(all_inf_trem_en, ((0,1),(0,0)), 'constant', constant_values=(0))
        all_inf_trem_en[-1][0]= med_en[x,0]
        all_inf_trem_en[-1][1]= med_en[x,1]
        if len(all_inf_trem_en) == 1:
            all_inf_trem_en[-1][2] = med_en[x,1]
        else:
            all_inf_trem_en[-1][2] = all_inf_trem_en[-2][2] + med_en[x,1]
        
    if e_type_i[x] == 10 :
        all_SA_trem_en_inf = np.lib.pad(all_SA_trem_en_inf, ((0,1),(0,0)), 'constant', constant_values=(0))
        all_SA_trem_en_inf[-1][0]= med_en[x,0]
        all_SA_trem_en_inf[-1][1]= med_en[x,1]
        if len(all_SA_trem_en_inf) == 1:
            all_SA_trem_en_inf[-1][2] = med_en[x,1]
        else:
            all_SA_trem_en_inf[-1][2] = all_SA_trem_en_inf[-2][2] + med_en[x,1]
        
    if e_type_i[x] == 50 or  e_type_i[x] == 60:
        all_SA_exp_en_inf = np.lib.pad(all_SA_exp_en_inf, ((0,1),(0,0)), 'constant', constant_values=(0))
        all_SA_exp_en_inf[-1][0]= med_en[x,0]
        all_SA_exp_en_inf[-1][1]= med_en[x,1]
        if len(all_SA_exp_en_inf) == 1:
            all_SA_exp_en_inf[-1][2] = med_en[x,1]
        else:
            all_SA_exp_en_inf[-1][2] = all_SA_exp_en_inf[-2][2] + med_en[x,1]
        
   
  
for x in range(0,len(med_ens)):

    if e_type_s[x] == 7  or e_type_s[x] == 50  or e_type_s[x] == 60 :
        all_exp_en_seis = np.lib.pad(all_exp_en_seis, ((0,1),(0,0)), 'constant', constant_values=(0))
        all_exp_en_seis[-1][0]= med_ens[x,0]
        all_exp_en_seis[-1][1]= med_ens[x,2]
        if len(all_exp_en_seis) == 1:
            all_exp_en_seis[-1][2] = med_ens[x,2]
        else:
            all_exp_en_seis[-1][2] = all_exp_en_seis[-2][2] + med_ens[x,2]
    
    if e_type_s[x] == 2 or e_type_s[x] == 20  :
        all_seis_trem_en = np.lib.pad(all_seis_trem_en, ((0,1),(0,0)), 'constant', constant_values=(0))
        all_seis_trem_en[-1][0]= med_ens[x,0]
        all_seis_trem_en[-1][1]= med_ens[x,2]
        if len(all_seis_trem_en) == 1:
            all_seis_trem_en[-1][2] = med_ens[x,2]
        else:
            all_seis_trem_en[-1][2] = all_seis_trem_en[-2][2] + med_ens[x,2]  
            
    if e_type_s[x] == 20 :
        all_SA_trem_en_seis = np.lib.pad(all_SA_trem_en_seis, ((0,1),(0,0)), 'constant', constant_values=(0))
        all_SA_trem_en_seis[-1][0]= med_ens[x,0]
        all_SA_trem_en_seis[-1][1]= med_ens[x,2]
        if len(all_SA_trem_en_seis) == 1:
            all_SA_trem_en_seis[-1][2] = med_ens[x,2]
        else:
            all_SA_trem_en_seis[-1][2] = all_SA_trem_en_seis[-2][2] + med_ens[x,2]
            
    if e_type_s[x] == 50 or e_type_s[x] == 60:
        all_SA_exp_en_seis = np.lib.pad(all_SA_exp_en_seis, ((0,1),(0,0)), 'constant', constant_values=(0))
        all_SA_exp_en_seis[-1][0]= med_ens[x,0]
        all_SA_exp_en_seis[-1][1]= med_ens[x,2]
        if len(all_SA_exp_en_seis) == 1:
            all_SA_exp_en_seis[-1][2] = med_ens[x,2]
        else:
            all_SA_exp_en_seis[-1][2] = all_SA_exp_en_seis[-2][2] + med_ens[x,2]
        
    
 #%% plot separate events   
    
            

plt.figure()
dates=[dt.datetime.fromtimestamp(ts) for ts in all_exp_en_inf[:,0]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums,all_exp_en_inf[:,2],'b')
plt.title('Cumulative Acoustic Explosion Energy')
plt.xlabel('Date')
plt.ylabel('Cumulative Energy [J]')


plt.figure()
dates=[dt.datetime.fromtimestamp(ts) for ts in gas_exp_en_inf[:,0]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums,gas_exp_en_inf[:,2],'b')
plt.title('Cumulative Acoustic Gas-rich Explosion Energy')
plt.xlabel('Date')
plt.ylabel('Cumulative Energy [J]')


plt.figure()
dates=[dt.datetime.fromtimestamp(ts) for ts in ash_exp_en_inf[:,0]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums,ash_exp_en_inf[:,2],'b')
plt.title('Cumulative Acoustic Ash-rich Explosion Energy')
plt.xlabel('Date')
plt.ylabel('Cumulative Energy [J]')

plt.figure()
dates1=[dt.datetime.fromtimestamp(ts) for ts in all_exp_en_inf[:,0]]
datenums1=md.date2num(dates1)
dates2=[dt.datetime.fromtimestamp(ts) for ts in gas_exp_en_inf[:,0]]
datenums2=md.date2num(dates2)
dates3=[dt.datetime.fromtimestamp(ts) for ts in ash_exp_en_inf[:,0]]
datenums3=md.date2num(dates3)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums1,all_exp_en_inf[:,2],'k')
plt.plot(datenums2,gas_exp_en_inf[:,2],'b')
plt.plot(datenums3,ash_exp_en_inf[:,2],'r')
plt.title('Cumulative Acoustic Explosion Energy')
plt.xlabel('Date')
plt.ylabel('Cumulative Energy [J]')
plt.legend(['All explosions','gas rich','ash rich'])



plt.figure()
dates=[dt.datetime.fromtimestamp(ts) for ts in all_exp_en_seis[:,0]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums,all_exp_en_seis[:,2],'b')
plt.title('Cumulative Seismic Explosion Energy')
plt.xlabel('Date')
plt.ylabel('Cumulative Energy [J]')



plt.figure()
dates=[dt.datetime.fromtimestamp(ts) for ts in all_seis_trem_en[:,0]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums,all_seis_trem_en[:,2],'b')
plt.title('Cumulative Seismic Tremor Energy')
plt.xlabel('Date')
plt.ylabel('Cumulative Energy [J]')


plt.figure()
dates=[dt.datetime.fromtimestamp(ts) for ts in all_inf_trem_en[:,0]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums,all_inf_trem_en[:,2],'b')
plt.title('Cumulative Acoustic Tremor Energy')
plt.xlabel('Date')
plt.ylabel('Cumulative Energy [J]')




plt.figure()
dates=[dt.datetime.fromtimestamp(ts) for ts in all_SA_trem_en_inf[:,0]]
datenums=md.date2num(dates)
dates2=[dt.datetime.fromtimestamp(ts) for ts in all_SA_trem_en_seis[:,0]]
datenums2=md.date2num(dates2)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums,all_SA_trem_en_inf[:,2],'b')
plt.plot(datenums2,all_SA_trem_en_seis[:,2],'r')
plt.title('Cumulative Seismo-Acoustic Tremor Energy')
plt.xlabel('Date')
plt.ylabel('Cumulative Energy [J]')
plt.legend(['Acoustic','Seismic'])


plt.figure()
dates=[dt.datetime.fromtimestamp(ts) for ts in  all_SA_exp_en_inf[:,0]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums, all_SA_exp_en_inf[:,2],'b')
plt.title('Cumulative SeismoAcoustic Explosion Energy')
plt.xlabel('Date')
plt.ylabel('Cumulative Energy [J]')
plt.legend(['Acoustic'])

plt.figure()
dates=[dt.datetime.fromtimestamp(ts) for ts in  all_SA_exp_en_seis[:,0]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums, all_SA_exp_en_seis[:,2],'r')
plt.title('Cumulative SeismoAcoustic Explosion Energy')
plt.xlabel('Date')
plt.ylabel('Cumulative Energy [J]')
plt.legend(['Seismic'])



































