#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 14:35:17 2020

@author: root
"""

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



#%%



Fill = np.zeros(shape =(77,1))
dates = np.zeros(shape =(77,1))



dates[0]=UTCDateTime(2018,1,1).timestamp
dates[1]=UTCDateTime(2018,1,15).timestamp
dates[2]=UTCDateTime(2018,1,20).timestamp
dates[3]=UTCDateTime(2018,1,25).timestamp
dates[4]=UTCDateTime(2018,2,4).timestamp
dates[5]=UTCDateTime(2018,2,9).timestamp
dates[6]=UTCDateTime(2018,3,1).timestamp
dates[7]=UTCDateTime(2018,3,6).timestamp
dates[8]=UTCDateTime(2018,3,11).timestamp
dates[9]=UTCDateTime(2018,3,16).timestamp
dates[10]=UTCDateTime(2018,3,21).timestamp
dates[11]=UTCDateTime(2018,4,20).timestamp
dates[12]=UTCDateTime(2018,7,4).timestamp
dates[13]=UTCDateTime(2018,7,24).timestamp
dates[14]=UTCDateTime(2018,9,17).timestamp
dates[15]=UTCDateTime(2018,10,2).timestamp
dates[16]=UTCDateTime(2018,11,11).timestamp
dates[17]=UTCDateTime(2018,11,16).timestamp
dates[18]=UTCDateTime(2018,12,6).timestamp
dates[19]=UTCDateTime(2018,12,16).timestamp
dates[20]=UTCDateTime(2018,12,26).timestamp
dates[21]=UTCDateTime(2018,12,31).timestamp
dates[22]=UTCDateTime(2019,1,5).timestamp
dates[23]=UTCDateTime(2019,1,20).timestamp
dates[24]=UTCDateTime(2019,1,25).timestamp
dates[25]=UTCDateTime(2019,2,14).timestamp
dates[26]=UTCDateTime(2019,2,19).timestamp
dates[27]=UTCDateTime(2019,3,6).timestamp
dates[28]=UTCDateTime(2019,3,11).timestamp
dates[29]=UTCDateTime(2019,3,26).timestamp
dates[30]=UTCDateTime(2019,3,31).timestamp
dates[31]=UTCDateTime(2019,4,15).timestamp
dates[32]=UTCDateTime(2019,4,25).timestamp
dates[33]=UTCDateTime(2019,6,19).timestamp
dates[34]=UTCDateTime(2019,7,4).timestamp
dates[35]=UTCDateTime(2019,7,19).timestamp
dates[36]=UTCDateTime(2019,7,29).timestamp
dates[37]=UTCDateTime(2019,11,1).timestamp
dates[38]=UTCDateTime(2019,11,6).timestamp
dates[39]=UTCDateTime(2019,11,16).timestamp
dates[40]=UTCDateTime(2019,11,26).timestamp
dates[41]=UTCDateTime(2019,12,1).timestamp
dates[42]=UTCDateTime(2019,12,11).timestamp
dates[43]=UTCDateTime(2019,12,16).timestamp
dates[44]=UTCDateTime(2020,1,5).timestamp
dates[45]=UTCDateTime(2020,1,15).timestamp
dates[46]=UTCDateTime(2020,2,14).timestamp
dates[47]=UTCDateTime(2020,2,19).timestamp
dates[48]=UTCDateTime(2020,2,29).timestamp
dates[49]=UTCDateTime(2020,3,5).timestamp
dates[50]=UTCDateTime(2020,3,10).timestamp
dates[51]=UTCDateTime(2020,3,15).timestamp
dates[52]=UTCDateTime(2020,3,20).timestamp
dates[53]=UTCDateTime(2020,3,25).timestamp
dates[54]=UTCDateTime(2020,4,9).timestamp
dates[55]=UTCDateTime(2020,4,14).timestamp
dates[56]=UTCDateTime(2020,4,19).timestamp
dates[57]=UTCDateTime(2020,4,24).timestamp
dates[58]=UTCDateTime(2020,5,9).timestamp
dates[59]=UTCDateTime(2020,6,8).timestamp
dates[60]=UTCDateTime(2020,6,23).timestamp
dates[61]=UTCDateTime(2020,7,8).timestamp
dates[62]=UTCDateTime(2020,8,12).timestamp
dates[63]=UTCDateTime(2020,9,1).timestamp
dates[64]=UTCDateTime(2020,9,11).timestamp
dates[65]=UTCDateTime(2020,9,21).timestamp
dates[66]=UTCDateTime(2020,10,11).timestamp
dates[67]=UTCDateTime(2020,10,31).timestamp
dates[68]=UTCDateTime(2020,11,10).timestamp
dates[69]=UTCDateTime(2020,11,15).timestamp
dates[70]=UTCDateTime(2020,11,20).timestamp
dates[71]=UTCDateTime(2020,11,25).timestamp
dates[72]=UTCDateTime(2020,11,30).timestamp
dates[73]=UTCDateTime(2020,12,10).timestamp
dates[74]=UTCDateTime(2020,12,20).timestamp
dates[75]=UTCDateTime(2020,12,25).timestamp
dates[76]=UTCDateTime(2020,12,30).timestamp


Fill[0]=25
Fill[1]=20
Fill[2]=0
Fill[3]=0
Fill[4]=50
Fill[5]=50
Fill[6]=20
Fill[7]=0
Fill[8]=50
Fill[9]=60
Fill[10]=75
Fill[11]=90
Fill[12]=0
Fill[13]=10
Fill[14]=0
Fill[15]=0
Fill[16]=20
Fill[17]=30
Fill[18]=10
Fill[19]=10
Fill[20]=0
Fill[21]=0
Fill[22]=20
Fill[23]=30
Fill[24]=40
Fill[25]=50
Fill[26]=60
Fill[27]=70
Fill[28]=80
Fill[29]=90
Fill[30]=90
Fill[31]=80
Fill[32]=80
Fill[33]=75
Fill[34]=70
Fill[35]=90
Fill[36]=90
Fill[37]=95
Fill[38]=95
Fill[39]=100
Fill[40]=100
Fill[41]=100
Fill[42]=100
Fill[43]=100
Fill[44]=100
Fill[45]=100
Fill[46]=100
Fill[47]=100
Fill[48]=100
Fill[49]=100
Fill[50]=100
Fill[51]=100
Fill[52]=100
Fill[53]=100
Fill[54]=100
Fill[55]=100
Fill[56]=100
Fill[57]=100
Fill[58]=100
Fill[59]=100
Fill[60]=100
Fill[61]=100
Fill[62]=95
Fill[63]=100
Fill[64]=100
Fill[65]=100
Fill[66]=100
Fill[67]=100
Fill[68]=100
Fill[69]=100
Fill[70]=100
Fill[71]=90
Fill[72]=80
Fill[73]=90
Fill[74]=95
Fill[75]=80
Fill[76]=80






#%%


plt.figure(figsize=(12,3))
Dates=[dt.datetime.fromtimestamp(ts) for ts in  dates[:]]
datenums=md.date2num(Dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=0)
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums, Fill[:],'kx-')

plt.xlabel('Date')
plt.ylabel('Fill [%]')
plt.title("Crater Fill Estimate")
ax.set_xlim([736695,737791])
#plt.savefig('/Users/william/Desktop/crater_fill.png', dpi=600, format='png')




#%%

filling = np.zeros(shape =(8,2))
emptying = np.zeros(shape =(9,2))
LA = np.zeros(shape =(8,2))

filling[0] = [UTCDateTime(2018,1,25).timestamp,(UTCDateTime(2018,2,9).timestamp-UTCDateTime(2018,1,25).timestamp)]
filling[1] = [UTCDateTime(2018,3,6).timestamp,(UTCDateTime(2018,4,20).timestamp-UTCDateTime(2018,3,6).timestamp)]
filling[2] = [UTCDateTime(2018,7,4).timestamp,(UTCDateTime(2018,7,24).timestamp-UTCDateTime(2018,7,4).timestamp)]
filling[3] = [UTCDateTime(2018,10,2).timestamp,(UTCDateTime(2018,11,16).timestamp-UTCDateTime(2018,10,2).timestamp)]
filling[4] = [UTCDateTime(2018,12,31).timestamp,(UTCDateTime(2019,3,31).timestamp-UTCDateTime(2018,12,31).timestamp)]
filling[5] = [UTCDateTime(2019,7,4).timestamp,(UTCDateTime(2020,7,8).timestamp-UTCDateTime(2019,7,4).timestamp)]
filling[6] = [UTCDateTime(2020,8,12).timestamp,(UTCDateTime(2020,11,20).timestamp-UTCDateTime(2020,8,12).timestamp)]
filling[7] = [UTCDateTime(2020,11,30).timestamp,(UTCDateTime(2020,12,20).timestamp-UTCDateTime(2020,11,30).timestamp)]



emptying[0] = [UTCDateTime(2018,1,1).timestamp,(UTCDateTime(2018,1,25).timestamp-UTCDateTime(2018,1,1).timestamp)]
emptying[1] = [UTCDateTime(2018,2,9).timestamp,(UTCDateTime(2018,3,6).timestamp-UTCDateTime(2018,2,9).timestamp)]
emptying[2] = [UTCDateTime(2018,5,27).timestamp,(UTCDateTime(2018,7,4).timestamp-UTCDateTime(2018,5,27).timestamp)]
emptying[3] = [UTCDateTime(2018,7,24).timestamp,(UTCDateTime(2018,10,2).timestamp-UTCDateTime(2018,7,24).timestamp)]
emptying[4] = [UTCDateTime(2018,11,27).timestamp,(UTCDateTime(2018,12,31).timestamp-UTCDateTime(2018,11,27).timestamp)]
emptying[5] = [UTCDateTime(2019,5,18).timestamp,(UTCDateTime(2019,7,4).timestamp-UTCDateTime(2019,5,18).timestamp)]
emptying[6] = [UTCDateTime(2020,7,8).timestamp,(UTCDateTime(2020,8,12).timestamp-UTCDateTime(2020,7,8).timestamp)]
emptying[7] = [UTCDateTime(2020,11,20).timestamp,(UTCDateTime(2020,11,30).timestamp-UTCDateTime(2020,11,20).timestamp)]
emptying[8] = [UTCDateTime(2020,12,20).timestamp,(UTCDateTime(2021,1,11).timestamp-UTCDateTime(2020,12,20).timestamp)]



LA[0] = [UTCDateTime(2018,1,31).timestamp,(UTCDateTime(2018,2,5).timestamp-UTCDateTime(2018,1,31).timestamp)]
LA[1] = [UTCDateTime(2018,4,16).timestamp,(UTCDateTime(2018,5,27).timestamp-UTCDateTime(2018,4,16).timestamp)]
LA[2] = [UTCDateTime(2018,10,12).timestamp,(UTCDateTime(2018,11,27).timestamp-UTCDateTime(2018,10,12).timestamp)]
LA[3] = [UTCDateTime(2019,3,24).timestamp,(UTCDateTime(2019,5,18).timestamp-UTCDateTime(2019,3,24).timestamp)]
LA[4] = [UTCDateTime(2019,10,2).timestamp,(UTCDateTime(2019,10,7).timestamp-UTCDateTime(2019,10,2).timestamp)]
LA[5] = [UTCDateTime(2019,11,6).timestamp,(UTCDateTime(2019,12,28).timestamp-UTCDateTime(2019,11,6).timestamp)]
LA[6] = [UTCDateTime(2020,2,12).timestamp,(UTCDateTime(2020,6,16).timestamp-UTCDateTime(2020,2,12).timestamp)]
LA[7] = [UTCDateTime(2020,8,31).timestamp,(UTCDateTime(2020,9,25).timestamp-UTCDateTime(2020,8,31).timestamp)]





fig, ax = plt.subplots(figsize=(12,3))
ax.broken_barh([filling[0], filling[1],filling[2],filling[3],filling[4],filling[5],filling[6],filling[7]],(0.1,0.2),facecolors='tab:blue')
ax.broken_barh([emptying[0], emptying[1],emptying[2],emptying[3],emptying[4],emptying[5],emptying[6],emptying[7],emptying[8]],(0.4,0.2),facecolors='tab:red')
ax.broken_barh([LA[0], LA[1],LA[2],LA[3],LA[4],LA[5], LA[6],LA[7]],(0.1,0.2),facecolors='tab:green')

ax.set_yticks([0.2, 0.5 ])
ax.set_yticklabels(['Crater Filling/Full/Overspill', 'Crater Empty/Emptying'])

ax.set_xlim([1514764800,1609459200])

#plt.savefig('/Users/william/Desktop/crater_fill_broken_bar.png', dpi=600, format='png')




#%%

#
#
#filling_z = np.zeros(shape =(6,2))
#emptying_z = np.zeros(shape =(5,2))
#LA_z = np.zeros(shape =(8,2))
#
#filling_z[0] = [UTCDateTime(2018,10,1).timestamp,(UTCDateTime(2018,11,16).timestamp-UTCDateTime(2018,10,1).timestamp)]
#
#emptying_z[0] = [UTCDateTime(2018,11,16).timestamp,(UTCDateTime(2018,12,1).timestamp-UTCDateTime(2018,11,16).timestamp)]
#
#LA_z[0] = [UTCDateTime(2018,10,12).timestamp,(UTCDateTime(2018,10,15).timestamp-UTCDateTime(2018,10,12).timestamp)]
#LA_z[1] = [UTCDateTime(2018,11,6).timestamp,(UTCDateTime(2018,11,26).timestamp-UTCDateTime(2018,11,6).timestamp)]
#
#
#
#
#fig, ax = plt.subplots(figsize=(12,3))
#ax.broken_barh([filling_z[0]],(0.1,0.2),facecolors='tab:blue')
#ax.broken_barh([emptying_z[0]],(0.4,0.2),facecolors='tab:red')
#ax.broken_barh([LA_z[0], LA_z[1]],(0.7,0.2),facecolors='tab:green')
#ax.set_yticks([0.2, 0.5, 0.8])
#ax.set_yticklabels(['Crater Full/Filling', 'Crater Empty/Emptying', 'Lava Flow Advancement'])
#
#ax.set_xlim([1538352000,1543622400])
##plt.savefig('/Users/william/Desktop/crater_fill_nov_broken_bar.png', dpi=600, format='png')
#
#




