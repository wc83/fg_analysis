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

Energy = genfromtxt("/Users/william/Documents/Fuego_catalogue/analysis/Fuego_exp_Energy_Matrix_final.csv", delimiter=',',skip_header=1)
ITE = genfromtxt("/Users/william/Documents/Fuego_catalogue/analysis/Fuego_inf_tremor_Energy_Matrix_v8.csv", delimiter=',',skip_header=1)
STE = genfromtxt("/Users/william/Documents/Fuego_catalogue/analysis/Fuego_seis_tremor_Energy_Matrix_v8.csv", delimiter=',',skip_header=1)

ITE12 = genfromtxt("/Users/william/Documents/Fuego_catalogue/analysis/Fuego_inf_trem_fg12.csv", delimiter=',',skip_header=1)
STE12 = genfromtxt("/Users/william/Documents/Fuego_catalogue/analysis/Fuego_seis_trem_fg12.csv", delimiter=',',skip_header=1)


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

#%%
plt.figure(figsize=(12,3))
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
plt.title("Sation Availability")
plt.legend(['Seismic','Infrasound','Total'])
#plt.savefig('/Users/william/Desktop/Stations.png', dpi=600, format='png')

#%%
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

#%%

fig, ax1 = plt.subplots(figsize=(12,3))
dates=[dt.datetime.fromtimestamp(ts) for ts in  epw[:,1]]
datenums=md.date2num(dates)
ax1=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax1.xaxis.set_major_formatter(xfmt)


#ax1.set_title('Weekly Ash vs Gas')
ax1.plot(datenums, epw[:,7],'b')
ax1.set_xlabel('Date')
ax1.set_ylabel('Gas Exposions per week [#]')
ax1.set_ylim([0,1500])
#ax1.set_xlim([0,111])
ax2 = ax1.twinx()
ax2.plot(datenums,epw[:,8],'r')
ax2.set_ylabel('Ash Exposions per week [#]')
ax2.set_ylim([0,1500])
ax1.legend(['Gas'],bbox_to_anchor=(0.8, 1), loc='upper left')
ax2.legend(['Ash'],bbox_to_anchor=(0.8, 0.9), loc='upper left')
#%%


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


plt.figure(figsize=(12,3))
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
dates=[dt.datetime.fromtimestamp(ts) for ts in  epw[5:len(epw),1]]
datenums=md.date2num(dates)
ax1=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax1.xaxis.set_major_formatter(xfmt)
ax1.plot(datenums, epw[5:len(epw),6],'b', label = 'Weekly Detections')
ax1.set_ylabel('Number of Explosions [#]')
#ax2 = ax1.twinx()
#ax2.plot(datenums, epw[5:len(epw),20],'r', label = 'Cumulative Energy')
#ax2.set_ylabel('Explosions Energy [J]')
#ax1.legend(loc='upper left', bbox_to_anchor=(0, 1))
#ax2.legend(loc='upper left', bbox_to_anchor=(0, 0.85))
plt.title("Explosions")
#plt.savefig('/Users/william/Desktop/exp.png', dpi=600, format='png')

#Inf trem duration per week and cumulative Energy 


plt.figure(figsize=(12,3))
dates=[dt.datetime.fromtimestamp(ts) for ts in  epw[5:len(epw),1]]
datenums=md.date2num(dates)
ax1=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax1.xaxis.set_major_formatter(xfmt)
ax1.plot(datenums, epw[5:len(epw),13],'b',label = 'Duration per Week')
ax1.set_ylabel('Tremor Duration [Hours]')
#ax2 = ax1.twinx()
#ax2.plot(datenums, epw[5:len(epw),18],'r',label = 'Cumulative Energy')
#ax2.set_ylabel('Tremor Energy [J]')
#ax1.legend(loc='upper left', bbox_to_anchor=(0, 1))
#ax2.legend(loc='upper left', bbox_to_anchor=(0, 0.85))
plt.title("Infrasound tremor")
#plt.savefig('/Users/william/Desktop/Inf_trem.png', dpi=600, format='png')


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
dates=[dt.datetime.fromtimestamp(ts) for ts in  epw[5:len(epw),1]]
datenums=md.date2num(dates)
ax1=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=45 )
ax1.xaxis.set_major_formatter(xfmt)
ax1.plot(datenums, epw[5:len(epw),14],'b',label = 'Duration per Week')
ax1.set_ylabel('Tremor Duration [Hours]')
#ax2 = ax1.twinx()
#ax2.plot(datenums, epw[5:len(epw),19],'r',label = 'Cumulative Energy')
#ax2.set_ylabel('Tremor Energy [J]')
#ax1.legend(loc='upper left', bbox_to_anchor=(0, 1))
#ax2.legend(loc='upper left', bbox_to_anchor=(0, 0.85))
plt.title("Seismic Tremor")
#plt.savefig('/Users/william/Desktop/Seis_trem.png', dpi=600, format='png')



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






#%% explosions


En_t = Energy[:,0]
e_type = cat[:,8]

En_inf = Energy[:,1:34]
En_seis = Energy[:,34:42]

med_en = np.zeros(shape =(len(En_inf),6))
med_en[:,0]=En_t
med_en[:,5]=e_type

inf_ex = np.zeros(shape =(0,2))

seis_ex = np.zeros(shape =(0,2))


for x in range(0,len(En_inf)):
    nz = En_inf[x]
    for p in range(len(nz)-1,-1,-1):
        if 1000 > nz[p] or  nz[p] > 1E9:
            nz = np.delete(nz,p)
    
#    nz = En_inf[x,:][np.nonzero(En_inf[x,:])]
    
    if len(nz) > 0:
        med_en_ev = np.percentile(nz,50)

        med_en[x,1]=med_en_ev
        
        inf_ex = np.lib.pad(inf_ex, ((0,1),(0,0)), 'constant', constant_values=(0))
        inf_ex[-1,0] = med_en[x,0]
        inf_ex[-1,1] = med_en[x,1]
        



    else:
        med_en[x,1]=0
        
            
for x in range(0,len(En_seis)):
    nz = En_seis[x]
    for p in range(len(nz)-1,-1,-1):
        if 1000 > nz[p] or  nz[p] > 1E8:
            nz = np.delete(nz,p)
#    nz = En_seis[x,:][np.nonzero(En_seis[x,:])]
    
    if len(nz) > 0:
        med_en_ev = np.percentile(nz,50)
#        med_en_ev = max(nz)
       
        med_en[x,2]=med_en_ev

        seis_ex = np.lib.pad(seis_ex, ((0,1),(0,0)), 'constant', constant_values=(0))
        seis_ex[-1,0] = med_en[x,0]
        seis_ex[-1,1] = med_en[x,2]  
        
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

plt.figure(figsize=(12,3))    
dates=[dt.datetime.fromtimestamp(ts) for ts in En_t]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=0 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums,med_en[:,2],'rx')
plt.title('seismic energy')
#plt.ylim([0,1E10])
plt.xlabel('Date')
plt.ylabel('Explosion seismic energy [J]')

plt.figure(figsize=(12,3))
dates=[dt.datetime.fromtimestamp(ts) for ts in En_t]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=0 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums,med_en[:,1],'bx')
#plt.ylim([0,1E10])
plt.title('acoustic energy')
plt.xlabel('Date')
plt.ylabel('Explosion acoustic energy [J]')

plt.figure(figsize=(12,3)) 
dates=[dt.datetime.fromtimestamp(ts) for ts in En_t]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=0 )
ax.xaxis.set_major_formatter(xfmt)   
plt.plot(datenums,med_en[:,4],'r')
plt.title('Cumulative explosion seismic energy')
plt.xlabel('Date')
plt.ylabel('Cumulative explosion seismic energy [J]')

plt.figure(figsize=(12,3))
dates=[dt.datetime.fromtimestamp(ts) for ts in En_t]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=0 )
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
    
            

plt.figure(figsize=(12,3))
dates=[dt.datetime.fromtimestamp(ts) for ts in all_exp[:,0]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=0 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums,all_exp[:,2],'b')
plt.title('Cumulative Acoustic Explosion Energy')
plt.xlabel('Date')
plt.ylabel('Cumulative Energy [J]')

plt.figure(figsize=(12,3))
dates=[dt.datetime.fromtimestamp(ts) for ts in all_exp[:,0]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=0 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums,all_exp[:,4],'r')
plt.title('Cumulative Seismic Explosion Energy')
plt.xlabel('Date')
plt.ylabel('Cumulative Energy [J]')


plt.figure(figsize=(12,3))
dates=[dt.datetime.fromtimestamp(ts) for ts in all_ash_poor[:,0]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=0 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums,all_ash_poor[:,2],'b')
plt.title('Cumulative Ash-poor Acoustic Explosion Energy')
plt.xlabel('Date')
plt.ylabel('Cumulative Energy [J]')

plt.figure(figsize=(12,3))
dates=[dt.datetime.fromtimestamp(ts) for ts in all_ash_rich[:,0]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=0 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums,all_ash_rich[:,2],'b')
plt.title('Cumulative Acoustic Ash-rich Explosion Energy')
plt.xlabel('Date')
plt.ylabel('Cumulative Energy [J]')


plt.figure(figsize=(12,3))
dates1=[dt.datetime.fromtimestamp(ts) for ts in all_exp[:,0]]
datenums1=md.date2num(dates1)
dates2=[dt.datetime.fromtimestamp(ts) for ts in all_ash_poor[:,0]]
datenums2=md.date2num(dates2)
dates3=[dt.datetime.fromtimestamp(ts) for ts in all_ash_rich[:,0]]
datenums3=md.date2num(dates3)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=0 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums1,all_exp[:,2],'k')
plt.plot(datenums2,all_ash_poor[:,2],'b')
plt.plot(datenums3,all_ash_rich[:,2],'r')
plt.title('Cumulative Acoustic Explosion Energy')
plt.xlabel('Date')
plt.ylabel('Cumulative Energy [J]')
plt.legend(['All explosions','Ash poor','Ash rich'])


plt.figure(figsize=(12,3))
dates=[dt.datetime.fromtimestamp(ts) for ts in ash_poor_no_gca[:,0]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=0 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums,ash_poor_no_gca[:,2],'b')
plt.title('Cumulative Acoustic Ash-poor, no GCA Explosion Energy')
plt.xlabel('Date')
plt.ylabel('Cumulative Energy [J]')


plt.figure(figsize=(12,3))
dates=[dt.datetime.fromtimestamp(ts) for ts in ash_poor_gca[:,0]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=0 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums,ash_poor_gca[:,2],'b')
plt.title('Cumulative Acoustic Ash-poor with GCA Explosion Energy')
plt.xlabel('Date')
plt.ylabel('Cumulative Energy [J]')

plt.figure(figsize=(12,3))
dates=[dt.datetime.fromtimestamp(ts) for ts in ash_poor_gca[:,0]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=0 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums,ash_poor_gca[:,4],'r')
plt.title('Cumulative Seismic Ash-poor with GCA Explosion Energy')
plt.xlabel('Date')
plt.ylabel('Cumulative Energy [J]')


plt.figure(figsize=(12,3))
dates=[dt.datetime.fromtimestamp(ts) for ts in ash_rich_no_gca[:,0]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=0 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums,ash_rich_no_gca[:,2],'b')
plt.title('Cumulative Acoustic Ash-rich, no GCA Explosion Energy')
plt.xlabel('Date')
plt.ylabel('Cumulative Energy [J]')


plt.figure(figsize=(12,3))
dates=[dt.datetime.fromtimestamp(ts) for ts in ash_rich_gca[:,0]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=0 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums,ash_rich_gca[:,2],'b')
plt.title('Cumulative Acoustic Ash-rich with GCA Explosion Energy')
plt.xlabel('Date')
plt.ylabel('Cumulative Energy [J]')

plt.figure(figsize=(12,3))
dates=[dt.datetime.fromtimestamp(ts) for ts in ash_rich_gca[:,0]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=0 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums,ash_rich_gca[:,4],'r')
plt.title('Cumulative Seismic Ash-rich with GCA Explosion Energy')
plt.xlabel('Date')
plt.ylabel('Cumulative Energy [J]')


plt.figure(figsize=(12,3))
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
plt.xticks( rotation=0 )
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


plt.figure(figsize=(12,3))
dates=[dt.datetime.fromtimestamp(ts) for ts in all_gca[:,0]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=0 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums,all_gca[:,2],'b')
plt.title('Cumulative Acoustic GCA Energy')
plt.xlabel('Date')
plt.ylabel('Cumulative Energy [J]')

plt.figure(figsize=(12,3))
dates=[dt.datetime.fromtimestamp(ts) for ts in all_gca[:,0]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=0 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums,all_gca[:,4],'r')
plt.title('Cumulative Seismic GCA Energy')
plt.xlabel('Date')
plt.ylabel('Cumulative Energy [J]')


plt.figure(figsize=(12,3))
dates1=[dt.datetime.fromtimestamp(ts) for ts in all_exp[:,0]]
datenums1=md.date2num(dates1)
dates2=[dt.datetime.fromtimestamp(ts) for ts in all_gca[:,0]]
datenums2=md.date2num(dates2)
dates3=[dt.datetime.fromtimestamp(ts) for ts in all_no_gca[:,0]]
datenums3=md.date2num(dates3)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=0 )
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
    
    nz = En_inf_t[x]
    for p in range(len(nz)-1,-1,-1):
        if 1E4 > nz[p] or  nz[p] > 1E12:
            nz = np.delete(nz,p)
            
#    nz = En_inf_t[x,:][np.nonzero(En_inf_t[x,:])]
    
    if len(nz) > 0:
        med_en_ev = np.median(nz)
#        med_en_ev = max(nz)
        
        if med_en_ev < 5e11:
            med_en_it[x,1]=med_en_ev
        else:
            med_en[x,1]=1e11
    else:
        med_en[x,1]=0
        
            
for x in range(0,len(En_seis_t)):
   
#    nz = En_seis_t[x,:][np.nonzero(En_seis_t[x,:])]
    
    nz = En_seis_t[x]
    for p in range(len(nz)-1,-1,-1):
        if 1E4 > nz[p] or  nz[p] > 1E12:
            nz = np.delete(nz,p)
    
    if len(nz) > 0:
        med_en_ev = np.median(nz)
#        med_en_ev = max(nz)
        
        if med_en_ev < 5e11:
            med_en_st[x,1]=med_en_ev
        else:
            med_en[x,2]=1e11  
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

plt.figure(figsize=(12,3)) 
dates=[dt.datetime.fromtimestamp(ts) for ts in En_it]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=0 )
ax.xaxis.set_major_formatter(xfmt)   
plt.plot(datenums,med_en_it[:,1],'b')
plt.ylim([0,3E11])
plt.title('Acoustic Tremor Energy')
plt.xlabel('Date')
plt.ylabel('Event energy [J]')

plt.figure(figsize=(12,3))
dates=[dt.datetime.fromtimestamp(ts) for ts in En_st]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=0 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums,med_en_st[:,1],'r')
plt.ylim([0,2E9])
plt.title('Seismic Tremor Energy')
plt.xlabel('Date')
plt.ylabel('Event energy [J]')


plt.figure(figsize=(12,3)) 
dates=[dt.datetime.fromtimestamp(ts) for ts in En_it]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=0 )
ax.xaxis.set_major_formatter(xfmt)   
plt.plot(datenums,med_en_it[:,2],'b')
plt.title('Cumulative Acoustic Tremor Energy')
plt.xlabel('Date')
plt.ylabel('Cumulative energy [J]')

plt.figure(figsize=(12,3))
dates=[dt.datetime.fromtimestamp(ts) for ts in En_st]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=0 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums,med_en_st[:,2],'r')
plt.title('Cumualtive Seismic Tremor Energy')
plt.xlabel('Date')
plt.ylabel('Cumulative energy [J]')



plt.figure(figsize=(12,3))
dates1=[dt.datetime.fromtimestamp(ts) for ts in En_it]
datenums1=md.date2num(dates1)
dates2=[dt.datetime.fromtimestamp(ts) for ts in En_st]
datenums2=md.date2num(dates2)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=0 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums1,med_en_it[:,2],'b')
plt.plot(datenums2,med_en_st[:,2],'r')
plt.title('Cumulative Tremor Energy')
plt.xlabel('Date')
plt.ylabel('Cumulative Energy [J]')
plt.legend(['Acoustic','Seismic'])




plt.figure(figsize=(12,3))
dates1=[dt.datetime.fromtimestamp(ts) for ts in all_exp[:,0]]
datenums1=md.date2num(dates1)
dates2=[dt.datetime.fromtimestamp(ts) for ts in En_it]
datenums2=md.date2num(dates2)
dates3=[dt.datetime.fromtimestamp(ts) for ts in En_st]
datenums3=md.date2num(dates3)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=0 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums1,all_exp[:,2],'k')
plt.plot(datenums2,med_en_it[:,2],'b')
plt.plot(datenums3,med_en_st[:,2],'r')
plt.title('Cumulative Energies')
plt.xlabel('Date')
plt.ylabel('Cumulative Energy [J]')
plt.legend(['Explosions','Acoustic Trem','Seismic Trem'])







#%% average energy per week

av_inf_ex = np.zeros(shape =(0,2))

av_seis_ex = np.zeros(shape =(0,2))

sum_inf = 0
sum_seis = 0

num_inf = 0
num_seis = 0

s_time = UTCDateTime(2018,3,15).timestamp
w_len = 7*24*60*60
e_time = s_time + w_len

for x in range(0,len(inf_ex)):
    if s_time <= inf_ex[x,0] < e_time:
        
        sum_inf += inf_ex[x,1]
        num_inf += 1
    else:
        if num_inf > 0:
            av_inf_en = sum_inf/num_inf
            
        
            av_inf_ex = np.lib.pad(av_inf_ex, ((0,1),(0,0)), 'constant', constant_values=(0))
            
            av_inf_ex[-1,0] = s_time 
            av_inf_ex[-1,1] = av_inf_en
            
            num_inf = 1
            sum_inf = inf_ex[x,1]
        
        s_time += 7*24*60*60
        e_time += 7*24*60*60
        
dates_inf = av_inf_ex[:,0]

plt.figure(figsize=(12,3))
dates1=[dt.datetime.fromtimestamp(ts) for ts in dates_inf]
datenums1=md.date2num(dates1)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=0 )    
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums1, av_inf_ex[:,1],'b')
plt.title('Mean Inf Explosion Energy per week')
plt.xlabel('Date')
plt.ylabel('Mean Energy [J]')       
        
        
s_time = UTCDateTime(2018,3,15).timestamp
w_len = 7*24*60*60
e_time = s_time + w_len



for x in range(0,len(seis_ex)):
    if s_time <= seis_ex[x,0] < e_time:
        
        sum_seis += seis_ex[x,1]
        num_seis += 1
    else:
        if num_seis > 0:
            av_seis_en = sum_seis/num_seis

        
            av_seis_ex = np.lib.pad(av_seis_ex, ((0,1),(0,0)), 'constant', constant_values=(0))
            
            av_seis_ex[-1,0] = s_time 
            av_seis_ex[-1,1] = av_seis_en
            
            num_seis = 1
            sum_seis = seis_ex[x,1]
        
        s_time += 7*24*60*60
        e_time += 7*24*60*60
        
        
dates_seis = av_seis_ex[:,0]

plt.figure(figsize=(12,3))
dates2=[dt.datetime.fromtimestamp(ts) for ts in dates_seis ]
datenums2=md.date2num(dates2)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=0 )  
ax.xaxis.set_major_formatter(xfmt)  
plt.plot(datenums2, av_seis_ex[:,1],'r')
plt.title('Mean Seis Explosion Energy per week')
plt.xlabel('Date')
plt.ylabel('Mean Energy [J]')   




#%%


fig, ax1 = plt.subplots(figsize=(12,3))
dates=[dt.datetime.fromtimestamp(ts) for ts in  epw[:,1]]
datenums=md.date2num(dates)
dates1=[dt.datetime.fromtimestamp(ts) for ts in dates_inf]
datenums1=md.date2num(dates1)
dates2=[dt.datetime.fromtimestamp(ts) for ts in dates_seis ]
datenums2=md.date2num(dates2)


ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=0 )    
ax1.xaxis.set_major_formatter(xfmt)
ax1.plot(datenums, epw[:,6],'k',  label='Explosion Count')
ax2 = ax1.twinx()
ax2.plot(datenums2, av_seis_ex[:,1],'r',  label='Seismic Energy')
ax3 = ax1.twinx()
ax3.spines["right"].set_position(("axes", 1.07))
ax3.plot(datenums1, av_inf_ex[:,1],'b',  label='Acoustic Energy')

ax1.legend(bbox_to_anchor=(0.99, 0.95))
ax2.legend(bbox_to_anchor=(0.99, 0.85))
ax3.legend(bbox_to_anchor=(0.99, 0.75))

plt.title('Mean Explosion Energy per week')
ax1.set_xlabel('Date')
ax3.set_ylabel('Mean Acoustic Energy [J]')  
ax2.set_ylabel('Mean Seismic Energy [J]')   
ax1.set_ylabel('Weekly explosions')   



#%%

ash_r = np.zeros(shape =(0,2))
gas_r = np.zeros(shape =(0,2))


for x in range(0,len(med_en)):
    if med_en[x,1] > 0 and med_en[x,2] > 0 and med_en[x,1] != 10000000 and med_en[x,2] != 10000000:
        if med_en[x,5] == 1 or med_en[x,5] == 2:
            gas_r = np.lib.pad(gas_r, ((0,1),(0,0)), 'constant', constant_values=(0))
            gas_r[-1,0] = med_en[x,1] #inf
            gas_r[-1,1] = med_en[x,2] #seis
            
        if med_en[x,5] == 3 or med_en[x,5] == 4:
            ash_r = np.lib.pad(ash_r, ((0,1),(0,0)), 'constant', constant_values=(0))
            ash_r[-1,0] = med_en[x,1] #inf
            ash_r[-1,1] = med_en[x,2] #seis
            

plt.figure(figsize=(6,6))
plt.loglog(gas_r[:,1],gas_r[:,0],'b.')
plt.loglog(ash_r[:,1],ash_r[:,0],'r.')
plt.legend(['Gas-rich','Ash-rich'])
plt.xlabel('Seismic Energy')            
plt.ylabel('Acoustic Energy')  
plt.title("Median Energies")  




#%%

cat_st = genfromtxt("/Users/william/Documents/Fuego_catalogue/Final_catalogues/Fuego_trem_Seis_final.csv", delimiter=',',skip_header=1)

SAT = np.zeros(shape =(0,9))
#AT = np.zeros(shape =(0,3))
#ST = np.zeros(shape =(0,3))




for x in range(0,len(ITE12)): #go over all fg12 inf trem
    print(x,'of',len(ITE12),'for SAT / AT')
        
    seis_start,idx = find_nearest(cat_st[:,0],ITE12[x,0]) #find closest seis trem detected
    
        
    if ((cat_st[idx,0] < ITE12[x,0] < cat_st[idx,0] + cat_st[idx,1]) or (cat_st[idx-1,0] < ITE12[x,0] < cat_st[idx-1,0] + cat_st[idx-1,1])): # if overlap

        #Save SAT 
        SAT = np.lib.pad(SAT, ((0,1),(0,0)), 'constant', constant_values=(0))  
        SAT[-1,0] = ITE12[x,0] #time
        SAT[-1,2] = ITE12[x,1] #inf 
        SAT[-1,1] = ITE12[x,2] #seis
            
#    else:  # no overlap
#            
#        #Save AT 
#        AT = np.lib.pad(AT, ((0,1),(0,0)), 'constant', constant_values=(0)) 
#        AT[-1,0] = ITE12[x,0] #time
#        AT[-1,1] = ITE12[x,1] #inf 
#        AT[-1,2] = ITE12[x,2] #seis


#for x in range(0,len(STE12)):
#    print(x,'of',len(STE12),'for ST')
#    
#    sat_start,idx = find_nearest(SAT[:,0],STE12[x,0])
#    
#    if abs(sat_start - STE12[x,0]) > 1200:
#        
#        #Save ST 
#        ST = np.lib.pad(ST, ((0,1),(0,0)), 'constant', constant_values=(0)) 
#        ST[-1,0] = STE12[x,0] #time
#        ST[-1,1] = STE12[x,1] #inf 
#        ST[-1,2] = STE12[x,2] #seis
 
#%%       

# inf calibs
c3 = 2.78E-04
#for seismics
c6 = 3.30E-09

sta = 'FG12' # STATION 
cha1 = 'BHZ' # CHANNEL
cha2 = 'BDF'
net = 'GI'  # 
loc1 = '00'    # location, it depends mostly of which network you are in. 
loc2 = '02'
client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23

for x in range(0,len(SAT)):
    try:
        print("SAT", x, "of", len(SAT))
        t1 = UTCDateTime(int(SAT[x,0]) - 45)
        t2 = UTCDateTime(int(SAT[x,0]) + 90)
        
        
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc1, cha1, t1 , t2 )
        tr1 = st1[0] 
        tr1.data = tr1.data * c6
        sf1 = tr1.stats.sampling_rate
        tr1.detrend(type='linear')
        tr1.detrend(type='demean')
        tr1.filter(type='bandpass',freqmin=0.1, freqmax=10)
        
        st2 = Stream()
        st2 = client.get_waveforms(net, sta, loc2, cha2, t1 , t2 )
        tr2 = st2[0] 
        tr2.data = tr2.data * c3
        sf2 = tr2.stats.sampling_rate
        tr2.detrend(type='linear')
        tr2.detrend(type='demean')
        tr2.filter(type='bandpass',freqmin=0.1, freqmax=10)
        
        peak1,cf1, bwid501 = freq_info_NqV(tr1,t1,t2,sf1) 
        peak2,cf2, bwid502 = freq_info_NqV(tr2,t1,t2,sf2)
        
        max_s_amp = np.max(abs(tr1.data))
        max_i_amp = np.max(abs(tr2.data))
        
        
        SAT[x,3] = peak1 #seis dominant f
        SAT[x,4] = peak2 #inf dominant f
        SAT[x,5] = cf1 #seis central f
        SAT[x,6] = cf2 #inf central f
        SAT[x,7] = max_s_amp # seismic max amplitude
        SAT[x,8] = max_i_amp # inf max amplitude  
    
    except:
        print("Missing data")
        SAT[x,3] = 0
        SAT[x,4] = 0
        SAT[x,5] = 0
        SAT[x,6] = 0
        SAT[x,7] = 0
        SAT[x,8] = 0   
        


#%% get seismic vs acoustic information


#exp_only = genfromtxt("/Users/william/Documents/Fuego_catalogue/analysis/Fuego_fg12_exp_noGCA.csv", delimiter=',',skip_header=1)




#energy

ash_r2 = np.zeros(shape =(0,9))
gas_r2 = np.zeros(shape =(0,9))

#ash_r2_o = np.zeros(shape =(0,3))
#gas_r2_o = np.zeros(shape =(0,3))



#for x in range(0,len(exp_only)):
#        if exp_only[x,1]== 1 or exp_only[x,1]== 2 : 
#        
#            gas_r2_o = np.lib.pad(gas_r2_o, ((0,1),(0,0)), 'constant', constant_values=(0))
#            gas_r2_o[-1,0] = exp_only[x,2] #inf
#            gas_r2_o[-1,1] = exp_only[x,3] #seis
#            gas_r2_o[-1,2] = exp_only[x,0] #time
#            
#        
#        if exp_only[x,1] == 3 or exp_only[x,1] == 4 : 
#        
#            ash_r2_o = np.lib.pad(ash_r2_o, ((0,1),(0,0)), 'constant', constant_values=(0))
#            ash_r2_o[-1,0] = exp_only[x,2] #inf
#            ash_r2_o[-1,1] = exp_only[x,3] #seis
#            ash_r2_o[-1,2] = exp_only[x,0] #time

#FG12 - SA Events
for x in range(0,len(cat)):
    if Energy[x,23] > 0 and Energy[x,38] > 0:
        if cat[x,8]== 1 or cat[x,8]== 2 : 
        
            gas_r2 = np.lib.pad(gas_r2, ((0,1),(0,0)), 'constant', constant_values=(0))
            gas_r2[-1,0] = Energy[x,23] #inf
            gas_r2[-1,1] = Energy[x,38] #seis
            gas_r2[-1,2] = Energy[x,0] #time
        
        if cat[x,8] == 3 or cat[x,8] == 4 : 
        
            ash_r2 = np.lib.pad(ash_r2, ((0,1),(0,0)), 'constant', constant_values=(0))
            ash_r2[-1,0] = Energy[x,23] #inf
            ash_r2[-1,1] = Energy[x,38] #seis
            ash_r2[-1,2] = Energy[x,0] #time





# Frequency, amplitude, duration,



for x in range(0,len(gas_r2)):
    try:
        print("gas", x, "of", len(gas_r2))
        t1 = UTCDateTime(gas_r2[x,2] - 45)
        t2 = UTCDateTime(gas_r2[x,2] + 90)
        
        
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc1, cha1, t1 , t2 )
        tr1 = st1[0] 
        tr1.data = tr1.data * c6
        sf1 = tr1.stats.sampling_rate
        tr1.detrend(type='linear')
        tr1.detrend(type='demean')
        tr1.filter(type='bandpass',freqmin=0.1, freqmax=10)
        
        st2 = Stream()
        st2 = client.get_waveforms(net, sta, loc2, cha2, t1 , t2 )
        tr2 = st2[0] 
        tr2.data = tr2.data * c3
        sf2 = tr2.stats.sampling_rate
        tr2.detrend(type='linear')
        tr2.detrend(type='demean')
        tr2.filter(type='bandpass',freqmin=0.1, freqmax=10)
        
        peak1,cf1, bwid501 = freq_info_NqV(tr1,t1,t2,sf1) 
        peak2,cf2, bwid502 = freq_info_NqV(tr2,t1,t2,sf2)
        
        max_s_amp = np.max(abs(tr1.data))
        max_i_amp = np.max(abs(tr2.data))
        
        
        gas_r2[x,3] = peak1 #seis dominant f
        gas_r2[x,4] = peak2 #inf dominant f
        gas_r2[x,5] = cf1 #seis central f
        gas_r2[x,6] = cf2 #inf central f
        gas_r2[x,7] = max_s_amp # seismic max amplitude
        gas_r2[x,8] = max_i_amp # inf max amplitude
    except:
        print("Missing data")
        gas_r2[x,3] = 0
        gas_r2[x,4] = 0
        gas_r2[x,5] = 0
        gas_r2[x,6] = 0
        gas_r2[x,7] = 0
        gas_r2[x,8] = 0
    
    
for x in range(0,len(ash_r2)):
    try:
        print("ash", x, "of", len(ash_r2))
        t1 = UTCDateTime(ash_r2[x,2] - 45)
        t2 = UTCDateTime(ash_r2[x,2] + 90)
        
        
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc1, cha1, t1 , t2 )
        tr1 = st1[0] 
        tr1.data = tr1.data * c6
        sf1 = tr1.stats.sampling_rate
        tr1.detrend(type='linear')
        tr1.detrend(type='demean')
        tr1.filter(type='bandpass',freqmin=0.1, freqmax=10)
        
        st2 = Stream()
        st2 = client.get_waveforms(net, sta, loc2, cha2, t1 , t2 )
        tr2 = st2[0] 
        tr2.data = tr2.data * c3
        sf2 = tr2.stats.sampling_rate
        tr2.detrend(type='linear')
        tr2.detrend(type='demean')
        tr2.filter(type='bandpass',freqmin=0.1, freqmax=10)
        
        peak1,cf1, bwid501 = freq_info_NqV(tr1,t1,t2,sf1) 
        peak2,cf2, bwid502 = freq_info_NqV(tr2,t1,t2,sf2)
        
        max_s_amp = np.max(abs(tr1.data))
        max_i_amp = np.max(abs(tr2.data))
        
        
        ash_r2[x,3] = peak1 #seis dominant f
        ash_r2[x,4] = peak2 #inf dominant f
        ash_r2[x,5] = cf1 #seis central f
        ash_r2[x,6] = cf2 #inf central f
        ash_r2[x,7] = max_s_amp # seismic max amplitude
        ash_r2[x,8] = max_i_amp # inf max amplitude    
        
    except:
        print("Missing data")
        ash_r2[x,3] = 0
        ash_r2[x,4] = 0
        ash_r2[x,5] = 0
        ash_r2[x,6] = 0
        ash_r2[x,7] = 0
        ash_r2[x,8] = 0
    




#%% plot event comparisons


v_1000 =  [[1,1e20], [1e3,1e23]]
v_100 =  [[1,1e20], [1e2,1e22]]
v_10 =  [[1,1e20], [1e1,1e21]]
v_1 =  [[1,1e20], [1,1e20]]
v_01 =  [[1,1e20], [1e-1,1e19]]
v_001 =  [[1,1e20], [1e-2,1e18]]
v_0001 =  [[1,1e20], [1e-3,1e17]]
        
        
plt.figure(figsize=(8,8))

plt.loglog(SAT[:,1],SAT[:,2],'k.')
plt.loglog(gas_r2[:,1],gas_r2[:,0],'b.')
plt.loglog(ash_r2[:,1],ash_r2[:,0],'r.')

plt.loglog(v_1000[0],v_1000[1],'--',color='grey')
plt.loglog(v_100[0],v_100[1],'--',color='grey')
plt.loglog(v_10[0],v_10[1],'--',color='grey')
plt.loglog(v_1[0],v_1[1],'--',color='grey')
plt.loglog(v_01[0],v_01[1],'--',color='grey')
plt.loglog(v_001[0],v_001[1],'--',color='grey')
plt.loglog(v_0001[0],v_0001[1],'--',color='grey')

plt.text(7e11,1e12,'VASR=1',rotation=45,color='grey')
plt.text(7e10,1e12,'VASR=10',rotation=45,color='grey')
plt.text(7e9,1e12,'VASR=100',rotation=45,color='grey')
plt.text(7e8,1e12,'VASR=1000',rotation=45,color='grey')
plt.text(7e11,1e11,'VASR=0.1',rotation=45,color='grey')
plt.text(7e11,1e10,'VASR=0.01',rotation=45,color='grey')
plt.text(7e11,1e9,'VASR=0.001',rotation=45,color='grey')

plt.legend(['Seismo-Acoustic Trem','Gas-rich Exp w.GCA','Ash-rich Exp w.GCA'])
plt.xlabel('Seismic Energy [J]')            
plt.ylabel('Acoustic Energy [J]')  
plt.title("FG12 VASR: Seismo-Acoustic Events")  
plt.xlim([1e4,1e13])
plt.ylim([1e4,1e13])


#plt.figure(figsize=(8,8))
##plt.loglog(STE12[:,2],STE12[:,1],'c.')
##plt.loglog(ITE12[:,2],ITE12[:,1],'g.')
#
#plt.loglog(AT[:,2],AT[:,1],'g.')
#plt.loglog(ST[:,2],ST[:,1],'c.')
#
#plt.loglog(ash_r2_o[:,1],ash_r2_o[:,0],'k.')
#plt.loglog(gas_r2_o[:,1],gas_r2_o[:,0],'m.')
#
#plt.loglog(v_1000[0],v_1000[1],'--',color='grey')
#plt.loglog(v_100[0],v_100[1],'--',color='grey')
#plt.loglog(v_10[0],v_10[1],'--',color='grey')
#plt.loglog(v_1[0],v_1[1],'--',color='grey')
#plt.loglog(v_01[0],v_01[1],'--',color='grey')
#plt.loglog(v_001[0],v_001[1],'--',color='grey')
#plt.loglog(v_0001[0],v_0001[1],'--',color='grey')
#
#plt.text(7e11,1e12,'VASR=1',rotation=45,color='grey')
#plt.text(7e10,1e12,'VASR=10',rotation=45,color='grey')
#plt.text(7e9,1e12,'VASR=100',rotation=45,color='grey')
#plt.text(7e8,1e12,'VASR=1000',rotation=45,color='grey')
#plt.text(7e11,1e11,'VASR=0.1',rotation=45,color='grey')
#plt.text(7e11,1e10,'VASR=0.01',rotation=45,color='grey')
#plt.text(7e11,1e9,'VASR=0.001',rotation=45,color='grey')
#
#plt.legend(['Acoustic Trem','Seismic Trem','Ash-rich Exp wo.GCA','Gas-rich Exp wo.GCA'])
#plt.xlabel('Seismic Energy [J]')            
#plt.ylabel('Acoustic Energy [J]')  
#plt.title("FG12 VASR")  
#plt.xlim([1e4,1e13])
#plt.ylim([1e4,1e13])




        
plt.figure(figsize=(8,8))
plt.loglog(SAT[:,3],SAT[:,4],'k.')
plt.loglog(gas_r2[:,3],gas_r2[:,4],'b.')
plt.loglog(ash_r2[:,3],ash_r2[:,4],'r.')

plt.legend(['Seismo-Acoustic Trem','Gas-rich Exp w.GCA','Ash-rich Exp w.GCA'])
plt.xlabel('Seismic Dom f [Hz]')            
plt.ylabel('Acoustic Dom f [Hz]')  
plt.title("FG12 frequency: Seismo-Acoustic Events")  
#plt.xlim([1e4,1e13])
#plt.ylim([1e4,1e13])




plt.figure(figsize=(8,8))
plt.loglog(SAT[:,5],SAT[:,6],'k.')
plt.loglog(gas_r2[:,5],gas_r2[:,6],'b.')
plt.loglog(ash_r2[:,5],ash_r2[:,6],'r.')

plt.legend(['Seismo-Acoustic Trem','Gas-rich Exp w.GCA','Ash-rich Exp w.GCA'])
plt.xlabel('Seismic central f [Hz]')            
plt.ylabel('Acoustic central f [Hz]')  
plt.title("FG12 frequency: Seismo-Acoustic Events")  
#plt.xlim([1e4,1e13])
#plt.ylim([1e4,1e13])




plt.figure(figsize=(8,8))
plt.loglog(SAT[:,7],SAT[:,8],'k.')
plt.loglog(gas_r2[:,7],gas_r2[:,8],'b.')
plt.loglog(ash_r2[:,7],ash_r2[:,8],'r.')

plt.legend(['Seismo-Acoustic Trem','Gas-rich Exp w.GCA','Ash-rich Exp w.GCA'])
plt.xlabel('Seismic Amp [m/s]]')            
plt.ylabel('Acoustic Amp [Pa]]')  
plt.title("FG12 maximum Amplitude: Seismo-Acoustic Events")  
#plt.xlim([1e4,1e13])
#plt.ylim([1e4,1e13])
















