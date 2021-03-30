#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 11:22:46 2021

@author: root
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

mod = genfromtxt("/Users/william/Documents/Fuego_catalogue/analysis/Modis_records.csv", delimiter=',',skip_header=1)



plt.figure(figsize=(12,3))
dates=[dt.datetime.fromtimestamp(ts) for ts in  mod[:,0]]
datenums=md.date2num(dates)
ax=plt.gca()
xfmt = md.DateFormatter('%d-%m-%Y')
plt.xticks( rotation=0 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums, mod[:,2],'b-')
plt.plot(datenums, mod[:,2],'bo',markersize=5)

plt.xlabel('Date')
plt.ylabel('Radient Flux [MW]')
plt.title("Volcanogenic Radient Flux")
ax.set_xlim([736695,737791])
























