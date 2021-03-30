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

cat = genfromtxt("/Users/william/Documents/Fuego_catalogue/Final_catalogues/Fuego_explosions.csv", delimiter=',',skip_header=1)
IT_cat =genfromtxt("/Users/william/Documents/Fuego_catalogue/Final_catalogues/Fuego_trem_Inf_final.csv", delimiter=',',skip_header=1)
ST_cat =genfromtxt("/Users/william/Documents/Fuego_catalogue/Final_catalogues/Fuego_trem_Seis_final.csv", delimiter=',',skip_header=1)

Energy = genfromtxt("/Users/william/Documents/Fuego_catalogue/analysis/Fuego_exp_Energy_Matrix.csv", delimiter=',',skip_header=1)
ITE = genfromtxt("/Users/william/Documents/Fuego_catalogue/analysis/Fuego_inf_tremor_Energy_Matrix_v8.csv", delimiter=',',skip_header=1)
STE = genfromtxt("/Users/william/Documents/Fuego_catalogue/analysis/Fuego_seis_tremor_Energy_Matrix_v8.csv", delimiter=',',skip_header=1)


i_trem = genfromtxt("/Users/william/Documents/Fuego_catalogue/Final_catalogues/Fuego_trem_Inf_final.csv", delimiter=',',skip_header=1)
s_trem = genfromtxt("/Users/william/Documents/Fuego_catalogue/Final_catalogues/Fuego_trem_Seis_final.csv", delimiter=',',skip_header=1)


stations = genfromtxt("/Users/william/Documents/Fuego_catalogue/Final_catalogues/Fuego_Final_scan_stations.csv", delimiter=',',skip_header=1)

#%% explosions


En_t = Energy[:,0]
e_type = cat[:,8]

En_inf = Energy[:,1:30]
En_seis = Energy[:,30:38]

med_en = np.zeros(shape =(len(En_inf),6))
med_en[:,0]=En_t
med_en[:,5]=e_type


for x in range(0,len(En_inf)):
   
    nz = En_inf[x,:][np.nonzero(En_inf[x,:])]
    
    if len(nz) > 0:
        med_en_ev = np.median(nz)
        if med_en_ev < 7500000000:
            med_en[x,1]=med_en_ev
        else:
            med_en[x,1]=10000000
    else:
        med_en[x,1]=0
        
            
for x in range(0,len(En_seis)):
   
    nz = En_seis[x,:][np.nonzero(En_seis[x,:])]
    
    if len(nz) > 0:
        med_en_ev = np.median(nz)
        if med_en_ev < 7500000000:
            med_en[x,2]=med_en_ev
        else:
            med_en[x,2]=10000000    
    else:
        med_en[x,2]=0            
    

cum_inf = 0
cum_seis = 0

for x in range(0,len(med_en)):
    cum_inf += med_en[x,1]
    cum_seis += med_en[x,2]
    
    med_en[x,3] = cum_inf
    med_en[x,4] = cum_seis

    
#%%
    
print("")
print("")
print("Exploions")
print("")

plt.figure()    
dates=[dt.datetime.fromtimestamp(ts) for ts in En_t]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums,med_en[:,2],'rx')
plt.title('seismic energy')
plt.xlabel('Date')
plt.ylabel('Explosion seismic energy [J]')

plt.figure()
dates=[dt.datetime.fromtimestamp(ts) for ts in En_t]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums,med_en[:,1],'bx')
plt.title('acoustic energy')
plt.xlabel('Date')
plt.ylabel('Explosion acoustic energy [J]')

plt.figure() 
dates=[dt.datetime.fromtimestamp(ts) for ts in En_t]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)   
plt.plot(datenums,med_en[:,4],'r')
plt.title('Cumulative explosion seismic energy')
plt.xlabel('Date')
plt.ylabel('Cumulative explosion seismic energy [J]')

plt.figure()
dates=[dt.datetime.fromtimestamp(ts) for ts in En_t]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums,med_en[:,3],'b')
plt.title('cumualtive explosion acoustic energy')
plt.xlabel('Date')
plt.ylabel('Cumulative explosion acoustic energy [J]')



#%%



all_exp = np.zeros(shape =(0,5))
all_ash_poor = np.zeros(shape =(0,5))
all_ash_rich = np.zeros(shape =(0,5))
ash_rich_gca = np.zeros(shape =(0,5))
ash_rich_no_gca = np.zeros(shape =(0,3))
ash_poor_gca = np.zeros(shape =(0,5))
ash_poor_no_gca = np.zeros(shape =(0,3))
all_gca = np.zeros(shape =(0,5))
all_no_gca = np.zeros(shape =(0,3))



all_inf_trem_en = np.zeros(shape =(0,3))
all_seis_trem_en = np.zeros(shape =(0,3))
all_SA_trem_en_inf = np.zeros(shape =(0,3))
all_SA_trem_en_seis = np.zeros(shape =(0,3))
all_SA_exp_en_inf = np.zeros(shape =(0,3))
all_SA_exp_en_seis = np.zeros(shape =(0,3))


for x in range(0,len(cat)):
    
    if e_type[x] == 1 or e_type[x] == 2 or e_type[x] == 3 or e_type[x] == 4 or e_type[x] == 5:# all exp 
        all_exp = np.lib.pad(all_exp, ((0,1),(0,0)), 'constant', constant_values=(0))
        all_exp[-1][0]= cat[x,0]
        all_exp[-1][1]= med_en[x,1]
        all_exp[-1][3]= med_en[x,2]
        
        if len(all_exp) == 1:
            all_exp[-1][2] = med_en[x,1]
            all_exp[-1][4] = med_en[x,2]
        else:
            all_exp[-1][2] = all_exp[-2][2] + med_en[x,1]
            all_exp[-1][4] = all_exp[-2][4] + med_en[x,2]
        
            
    if e_type[x] == 1 or e_type[x] == 2 : # Ash poor
        all_ash_poor = np.lib.pad(all_ash_poor, ((0,1),(0,0)), 'constant', constant_values=(0))
        all_ash_poor[-1][0]= cat[x,0]
        all_ash_poor[-1][1]= med_en[x,1]
        all_ash_poor[-1][3] = med_en[x,2]
        if len(all_ash_poor) == 1:
            all_ash_poor[-1][2] = med_en[x,1]
            all_ash_poor[-1][4] = med_en[x,2]
        else:
            all_ash_poor[-1][2] = all_ash_poor[-2][2] + med_en[x,1]
            all_ash_poor[-1][4] = all_ash_poor[-2][4] + med_en[x,2]
        
        
    if e_type[x] == 3 or e_type[x] == 4 : # Ash rich
        all_ash_rich = np.lib.pad(all_ash_rich, ((0,1),(0,0)), 'constant', constant_values=(0))
        all_ash_rich[-1][0]= cat[x,0]
        all_ash_rich[-1][1]= med_en[x,1]
        all_ash_rich[-1][3]= med_en[x,2]
        
        if len(all_ash_rich) == 1:
            all_ash_rich[-1][2] = med_en[x,1]
            all_ash_rich[-1][4] = med_en[x,2]
        else:
            all_ash_rich[-1][2] = all_ash_rich[-2][2] + med_en[x,1]
            all_ash_rich[-1][4] = all_ash_rich[-2][4] + med_en[x,2]
    
        
    if e_type[x] == 1 : # Ash poor no gca
        ash_poor_no_gca = np.lib.pad(ash_poor_no_gca, ((0,1),(0,0)), 'constant', constant_values=(0))
        ash_poor_no_gca[-1][0]= cat[x,0]
        ash_poor_no_gca[-1][1]= med_en[x,1]
        if len(ash_poor_no_gca) == 1:
            ash_poor_no_gca[-1][2] = med_en[x,1]
        else:
            ash_poor_no_gca[-1][2] = ash_poor_no_gca[-2][2] + med_en[x,1]
            
    if e_type[x] == 2 : # Ash poor gca
        ash_poor_gca = np.lib.pad(ash_poor_gca, ((0,1),(0,0)), 'constant', constant_values=(0))
        ash_poor_gca[-1][0]= cat[x,0]
        ash_poor_gca[-1][1]= med_en[x,1]
        ash_poor_gca[-1][3]= med_en[x,2]
        
        if len(ash_poor_gca) == 1:
            ash_poor_gca[-1][2] = med_en[x,1]
            ash_poor_gca[-1][4] = med_en[x,2]
            
        else:
            ash_poor_gca[-1][2] = ash_poor_gca[-2][2] + med_en[x,1]
            ash_poor_gca[-1][4] = ash_poor_gca[-2][4] + med_en[x,2]
            
        
    if e_type[x] == 3 : # Ash rich no gca
        ash_rich_no_gca = np.lib.pad(ash_rich_no_gca, ((0,1),(0,0)), 'constant', constant_values=(0))
        ash_rich_no_gca[-1][0]= cat[x,0]
        ash_rich_no_gca[-1][1]= med_en[x,1]
        if len(ash_rich_no_gca) == 1:
            ash_rich_no_gca[-1][2] = med_en[x,1]
        else:
            ash_rich_no_gca[-1][2] = ash_rich_no_gca[-2][2] + med_en[x,1]
        
    if e_type[x] == 4 : # Ash rich gca
        ash_rich_gca = np.lib.pad(ash_rich_gca, ((0,1),(0,0)), 'constant', constant_values=(0))
        ash_rich_gca[-1][0]= cat[x,0]
        ash_rich_gca[-1][1]= med_en[x,1]
        ash_rich_gca[-1][3]= med_en[x,2]
        
        if len(ash_rich_gca) == 1:
            ash_rich_gca[-1][2] = med_en[x,1]
            ash_rich_gca[-1][4] = med_en[x,2]
        else:
            ash_rich_gca[-1][2] = ash_rich_gca[-2][2] + med_en[x,1]
            ash_rich_gca[-1][4] = ash_rich_gca[-2][4] + med_en[x,2]
            
    if e_type[x] == 2  or e_type[x] == 4  or e_type[x] == 5 : # All gca
        all_gca = np.lib.pad(all_gca, ((0,1),(0,0)), 'constant', constant_values=(0))
        all_gca[-1][0]= cat[x,0]
        all_gca[-1][1]= med_en[x,1]
        all_gca[-1][3]= med_en[x,2]
        
        if len(all_gca) == 1:
            all_gca[-1][2] = med_en[x,1]
            all_gca[-1][4] = med_en[x,2]
        else:
            all_gca[-1][2] = all_gca[-2][2] + med_en[x,1]
            all_gca[-1][4] = all_gca[-2][4] + med_en[x,2]
    
    if e_type[x] == 1 or e_type[x] == 3  : # All no gca 
        all_no_gca = np.lib.pad(all_no_gca, ((0,1),(0,0)), 'constant', constant_values=(0))
        all_no_gca[-1][0]= cat[x,0]
        all_no_gca[-1][1]= med_en[x,1]
        if len(all_no_gca) == 1:
            all_no_gca[-1][2] = med_en[x,1]
        else:
            all_no_gca[-1][2] = all_no_gca[-2][2] + med_en[x,1]
        

    
 #%% plot separate events   
    
            

plt.figure()
dates=[dt.datetime.fromtimestamp(ts) for ts in all_exp[:,0]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums,all_exp[:,2],'b')
plt.title('Cumulative Acoustic Explosion Energy')
plt.xlabel('Date')
plt.ylabel('Cumulative Energy [J]')

plt.figure()
dates=[dt.datetime.fromtimestamp(ts) for ts in all_exp[:,0]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums,all_exp[:,4],'r')
plt.title('Cumulative Seismic Explosion Energy')
plt.xlabel('Date')
plt.ylabel('Cumulative Energy [J]')


plt.figure()
dates=[dt.datetime.fromtimestamp(ts) for ts in all_ash_poor[:,0]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums,all_ash_poor[:,2],'b')
plt.title('Cumulative Ash-poor Acoustic Explosion Energy')
plt.xlabel('Date')
plt.ylabel('Cumulative Energy [J]')

plt.figure()
dates=[dt.datetime.fromtimestamp(ts) for ts in all_ash_rich[:,0]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums,all_ash_rich[:,2],'b')
plt.title('Cumulative Acoustic Ash-rich Explosion Energy')
plt.xlabel('Date')
plt.ylabel('Cumulative Energy [J]')


plt.figure()
dates1=[dt.datetime.fromtimestamp(ts) for ts in all_exp[:,0]]
datenums1=md.date2num(dates1)
dates2=[dt.datetime.fromtimestamp(ts) for ts in all_ash_poor[:,0]]
datenums2=md.date2num(dates2)
dates3=[dt.datetime.fromtimestamp(ts) for ts in all_ash_rich[:,0]]
datenums3=md.date2num(dates3)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums1,all_exp[:,2],'k')
plt.plot(datenums2,all_ash_poor[:,2],'b')
plt.plot(datenums3,all_ash_rich[:,2],'r')
plt.title('Cumulative Acoustic Explosion Energy')
plt.xlabel('Date')
plt.ylabel('Cumulative Energy [J]')
plt.legend(['All explosions','Ash poor','Ash rich'])


plt.figure()
dates=[dt.datetime.fromtimestamp(ts) for ts in ash_poor_no_gca[:,0]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums,ash_poor_no_gca[:,2],'b')
plt.title('Cumulative Acoustic Ash-poor, no GCA Explosion Energy')
plt.xlabel('Date')
plt.ylabel('Cumulative Energy [J]')


plt.figure()
dates=[dt.datetime.fromtimestamp(ts) for ts in ash_poor_gca[:,0]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums,ash_poor_gca[:,2],'b')
plt.title('Cumulative Acoustic Ash-poor with GCA Explosion Energy')
plt.xlabel('Date')
plt.ylabel('Cumulative Energy [J]')

plt.figure()
dates=[dt.datetime.fromtimestamp(ts) for ts in ash_poor_gca[:,0]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums,ash_poor_gca[:,4],'r')
plt.title('Cumulative Seismic Ash-poor with GCA Explosion Energy')
plt.xlabel('Date')
plt.ylabel('Cumulative Energy [J]')


plt.figure()
dates=[dt.datetime.fromtimestamp(ts) for ts in ash_rich_no_gca[:,0]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums,ash_rich_no_gca[:,2],'b')
plt.title('Cumulative Acoustic Ash-rich, no GCA Explosion Energy')
plt.xlabel('Date')
plt.ylabel('Cumulative Energy [J]')


plt.figure()
dates=[dt.datetime.fromtimestamp(ts) for ts in ash_rich_gca[:,0]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums,ash_rich_gca[:,2],'b')
plt.title('Cumulative Acoustic Ash-rich with GCA Explosion Energy')
plt.xlabel('Date')
plt.ylabel('Cumulative Energy [J]')

plt.figure()
dates=[dt.datetime.fromtimestamp(ts) for ts in ash_rich_gca[:,0]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums,ash_rich_gca[:,4],'r')
plt.title('Cumulative Seismic Ash-rich with GCA Explosion Energy')
plt.xlabel('Date')
plt.ylabel('Cumulative Energy [J]')


plt.figure()
dates1=[dt.datetime.fromtimestamp(ts) for ts in all_exp[:,0]]
datenums1=md.date2num(dates1)
dates2=[dt.datetime.fromtimestamp(ts) for ts in ash_poor_no_gca[:,0]]
datenums2=md.date2num(dates2)
dates3=[dt.datetime.fromtimestamp(ts) for ts in ash_poor_gca[:,0]]
datenums3=md.date2num(dates3)
dates4=[dt.datetime.fromtimestamp(ts) for ts in ash_rich_no_gca[:,0]]
datenums4=md.date2num(dates4)
dates5=[dt.datetime.fromtimestamp(ts) for ts in ash_rich_gca[:,0]]
datenums5=md.date2num(dates5)

ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums1,all_exp[:,2],'k')
plt.plot(datenums2,ash_poor_no_gca[:,2],'b--')
plt.plot(datenums3,ash_poor_gca[:,2],'b')
plt.plot(datenums4,ash_rich_no_gca[:,2],'r--')
plt.plot(datenums5,ash_rich_gca[:,2],'r')
plt.title('Cumulative Acoustic Explosion Energy')
plt.xlabel('Date')
plt.ylabel('Cumulative Energy [J]')
plt.legend(['All explosions','Ash poor, no gca','Ash poor, gca','Ash rich, no gca','Ash rich, gca'])


plt.figure()
dates=[dt.datetime.fromtimestamp(ts) for ts in all_gca[:,0]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums,all_gca[:,2],'b')
plt.title('Cumulative Acoustic GCA Energy')
plt.xlabel('Date')
plt.ylabel('Cumulative Energy [J]')

plt.figure()
dates=[dt.datetime.fromtimestamp(ts) for ts in all_gca[:,0]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums,all_gca[:,4],'r')
plt.title('Cumulative Seismic GCA Energy')
plt.xlabel('Date')
plt.ylabel('Cumulative Energy [J]')


plt.figure()
dates1=[dt.datetime.fromtimestamp(ts) for ts in all_exp[:,0]]
datenums1=md.date2num(dates1)
dates2=[dt.datetime.fromtimestamp(ts) for ts in all_gca[:,0]]
datenums2=md.date2num(dates2)
dates3=[dt.datetime.fromtimestamp(ts) for ts in all_no_gca[:,0]]
datenums3=md.date2num(dates3)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums1,all_exp[:,2],'k')
plt.plot(datenums2,all_gca[:,2],'b')
plt.plot(datenums3,all_no_gca[:,2],'r')
plt.title('Cumulative Acoustic Explosion Energy')
plt.xlabel('Date')
plt.ylabel('Cumulative Energy [J]')
plt.legend(['All explosions','GCA','No GCA'])

#%%




#%%



En_it = ITE[:,0]
En_st = STE[:,0]


En_inf_t = ITE[:,1:34]
En_seis_t = STE[:,1:9]

med_en_it = np.zeros(shape =(len(En_inf_t),3))
med_en_it[:,0]=En_it


med_en_st = np.zeros(shape =(len(En_seis_t),3))
med_en_st[:,0]=En_st


for x in range(0,len(En_inf_t)):
   
    nz = En_inf_t[x,:][np.nonzero(En_inf_t[x,:])]
    
    if len(nz) > 0:
        med_en_ev = np.median(nz)
        if med_en_ev < 4.5e11:
            med_en_it[x,1]=med_en_ev
        else:
            med_en_it[x,1]=1e11
    else:
        med_en_it[x,1]=0
        
            
for x in range(0,len(En_seis_t)):
   
    nz = En_seis_t[x,:][np.nonzero(En_seis_t[x,:])]
    
    if len(nz) > 0:
        med_en_ev = np.median(nz)
        if med_en_ev < 5e11:
            med_en_st[x,1]=med_en_ev
        else:
            med_en_st[x,1]=1e11  
    else:
        med_en_st[x,1]=0            
    

cum_inf_t = 0
cum_seis_t = 0

for x in range(0,len(med_en_it)):
    cum_inf_t += med_en_it[x,1]
    med_en_it[x,2] = cum_inf_t
    
for x in range(0,len(med_en_st)):    
    cum_seis_t += med_en_st[x,1]
    med_en_st[x,2] = cum_seis_t


#%%
print("")
print("")
print("TREMOR")
print("")

plt.figure() 
dates=[dt.datetime.fromtimestamp(ts) for ts in En_it]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)   
plt.plot(datenums,med_en_it[:,1],'b')
plt.title('Acoustic Tremor Energy')
plt.xlabel('Date')
plt.ylabel('Event energy [J]')

plt.figure()
dates=[dt.datetime.fromtimestamp(ts) for ts in En_st]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums,med_en_st[:,1],'r')
plt.title('Seismic Tremor Energy')
plt.xlabel('Date')
plt.ylabel('Event energy [J]')


plt.figure() 
dates=[dt.datetime.fromtimestamp(ts) for ts in En_it]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)   
plt.plot(datenums,med_en_it[:,2],'b')
plt.title('Cumulative Acoustic Tremor Energy')
plt.xlabel('Date')
plt.ylabel('Cumulative energy [J]')

plt.figure()
dates=[dt.datetime.fromtimestamp(ts) for ts in En_st]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums,med_en_st[:,2],'r')
plt.title('Cumualtive Seismic Tremor Energy')
plt.xlabel('Date')
plt.ylabel('Cumulative energy [J]')



plt.figure()
dates1=[dt.datetime.fromtimestamp(ts) for ts in En_it]
datenums1=md.date2num(dates1)
dates2=[dt.datetime.fromtimestamp(ts) for ts in En_st]
datenums2=md.date2num(dates2)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums1,med_en_it[:,2],'b')
plt.plot(datenums2,med_en_st[:,2],'r')
plt.title('Cumulative Tremor Energy')
plt.xlabel('Date')
plt.ylabel('Cumulative Energy [J]')
plt.legend(['Acoustic','Seismic'])




plt.figure()
dates1=[dt.datetime.fromtimestamp(ts) for ts in all_exp[:,0]]
datenums1=md.date2num(dates1)
dates2=[dt.datetime.fromtimestamp(ts) for ts in En_it]
datenums2=md.date2num(dates2)
dates3=[dt.datetime.fromtimestamp(ts) for ts in En_st]
datenums3=md.date2num(dates3)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums1,all_exp[:,2],'k')
plt.plot(datenums2,med_en_it[:,2],'b')
plt.plot(datenums3,med_en_st[:,2],'r')
plt.title('Cumulative Energies')
plt.xlabel('Date')
plt.ylabel('Cumulative Energy [J]')
plt.legend(['Explosions','Acoustic Trem','Seismic Trem'])













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


epd=np.zeros(shape=(0,18)) 

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

ITD = 0
STD = 0
ETD = 0

day0 = 1538352000
last_day = 1543622400

for y in range(0,len(days_c)):
    start = day0 + y*dlen
    
    if start < last_day:
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
        
        for x in range(0,len(med_en_it)):
            if start < med_en_it[x,0] < end:
                ITD += med_en_it[x,1]
                
        for x in range(0,len(med_en_st)):
            if start < med_en_st[x,0] < end:
                STD += med_en_st[x,1]
                
        for x in range(0,len(all_exp)):
            if start < all_exp[x,0] < end:
                ETD += all_exp[x,1]
                
    
    
                
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
        epd[y][15] = ITD
        epd[y][16] = STD
        epd[y][17] = ETD
        
        
        
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
        
        ITD = 0
        STD = 0
        ETD = 0
        
    

    
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

epw=np.zeros(shape=(0,21))
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

ITE_tot = 0
STE_tot = 0
ETE_tot = 0

ITE_cum = 0
STE_cum = 0
ETE_cum = 0


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
    
    ITE_tot += epd[x][15]
    STE_tot += epd[x][16]
    ETE_tot += epd[x][17]
    
    ITE_cum += epd[x][15]
    STE_cum += epd[x][16]
    ETE_cum += epd[x][17]
    
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
        
        epw[week-1][15]= ITE_tot
        epw[week-1][16]= STE_tot
        epw[week-1][17]= ETE_tot
        
        epw[week-1][18]= ITE_cum
        epw[week-1][19]= STE_cum
        epw[week-1][20]= ETE_cum
        
        
        
        
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
        
        ITE_tot = 0
        STE_tot = 0
        ETE_tot = 0
    


    
#        
#%% plot events per week and save csv file with data




plt.figure()
dates=[dt.datetime.fromtimestamp(ts) for ts in  epw[:,1]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums, epw[:,5],'k')
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


fig, ax1 = plt.subplots(figsize=(12,3))
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










#%%  Plots for poster

#Explosion number per week and cumulative explosion Energy 

plt.figure(figsize=(12,3))
dates=[dt.datetime.fromtimestamp(ts) for ts in  epw[28:38,1]]
datenums=md.date2num(dates)
ax1=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax1.xaxis.set_major_formatter(xfmt)
ax1.plot(datenums, epw[28:38,6],'b', label = 'Weekly Detections')
ax1.set_ylabel('Number of Explosions [#]')
ax2 = ax1.twinx()
ax2.plot(datenums, epw[28:38,20],'r', label = 'Cumulative Energy')
ax2.set_ylabel('Explosions Energy [J]')
ax1.legend(loc='upper left', bbox_to_anchor=(0, 1))
ax2.legend(loc='upper left', bbox_to_anchor=(0, 0.85))
plt.title("Explosions")
plt.savefig('/Users/william/Desktop/exp.png', dpi=600, format='png')

#Inf trem duration per week and cumulative Energy 


plt.figure(figsize=(12,3))
dates=[dt.datetime.fromtimestamp(ts) for ts in  epw[28:38,1]]
datenums=md.date2num(dates)
ax1=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax1.xaxis.set_major_formatter(xfmt)
ax1.plot(datenums, epw[28:38,13],'b',label = 'Duration per Week')
ax1.set_ylabel('Tremor Duration [Hours]')
ax2 = ax1.twinx()
ax2.plot(datenums, epw[28:38,18],'r',label = 'Cumulative Energy')
ax2.set_ylabel('Tremor Energy [J]')
ax1.legend(loc='upper left', bbox_to_anchor=(0, 1))
ax2.legend(loc='upper left', bbox_to_anchor=(0, 0.85))
plt.title("Infrasound tremor")
plt.savefig('/Users/william/Desktop/Inf_trem.png', dpi=600, format='png')


#plt.figure(figsize=(12,3))
#dates=[dt.datetime.fromtimestamp(ts) for ts in  epw[:,1]]
#datenums=md.date2num(dates)
#ax=plt.gca()
#xfmt = md.DateFormatter('%d-%m-%Y')
#plt.xticks( rotation=45 )
#ax.xaxis.set_major_formatter(xfmt)
#ax.plot(datenums, epw[:,13],'b')
#ax.set_xlim([datenums[9], datenums[-5]])
#plt.title("Inf trem duration per week")
#
#plt.figure(figsize=(12,3))
#dates1=[dt.datetime.fromtimestamp(ts) for ts in En_it]
#datenums1=md.date2num(dates1)
#ax=plt.gca()
#xfmt = md.DateFormatter('%d-%m-%Y')
#plt.xticks( rotation=45 )
#ax.xaxis.set_major_formatter(xfmt)
#ax.plot(datenums1,med_en_it[:,2],'r')
#ax.set_xlim([datenums1[0], datenums1[-1]])
#plt.title("Inf trem cum energy")




#Seis trem duration per week and cumulative Energy 

plt.figure(figsize=(12,3))
dates=[dt.datetime.fromtimestamp(ts) for ts in  epw[28:38,1]]
datenums=md.date2num(dates)
ax1=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax1.xaxis.set_major_formatter(xfmt)
ax1.plot(datenums, epw[28:38,14],'b',label = 'Duration per Week')
ax1.set_ylabel('Tremor Duration [Hours]')
ax2 = ax1.twinx()
ax2.plot(datenums, epw[28:38,19],'r',label = 'Cumulative Energy')
ax2.set_ylabel('Tremor Energy [J]')
ax1.legend(loc='upper left', bbox_to_anchor=(0, 1))
ax2.legend(loc='upper left', bbox_to_anchor=(0, 0.85))
plt.title("Seismic Tremor")
plt.savefig('/Users/william/Desktop/Seis_trem.png', dpi=600, format='png')



#plt.figure(figsize=(12,3))
#dates=[dt.datetime.fromtimestamp(ts) for ts in  epw[:,1]]
#datenums=md.date2num(dates)
#ax=plt.gca()
#xfmt = md.DateFormatter('%d-%m-%Y')
#plt.xticks( rotation=45 )
#ax.xaxis.set_major_formatter(xfmt)
#ax.plot(datenums, epw[:,14],'b')
#ax.set_xlim([datenums[9], datenums[-5]])
#plt.title("Seis trem duration per week")
#
#
#plt.figure(figsize=(12,3))
#dates2=[dt.datetime.fromtimestamp(ts) for ts in En_st]
#datenums2=md.date2num(dates2)
#ax=plt.gca()
#xfmt = md.DateFormatter('%d-%m-%Y')
#plt.xticks( rotation=45 )
#ax.xaxis.set_major_formatter(xfmt)
#ax.plot(datenums2,med_en_st[:,2],'r')
#ax.set_xlim([datenums2[0], datenums2[-1]])
#plt.title("Seis trem cum energy")

