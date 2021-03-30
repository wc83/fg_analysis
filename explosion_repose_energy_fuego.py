#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 11:36:41 2020

@author: root
"""

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
#IT_cat =genfromtxt("/Users/william/Documents/Fuego_catalogue/Final_catalogues/Fuego_trem_Inf_final.csv", delimiter=',',skip_header=1)
#ST_cat =genfromtxt("/Users/william/Documents/Fuego_catalogue/Final_catalogues/Fuego_trem_Seis_final.csv", delimiter=',',skip_header=1)

Energy = genfromtxt("/Users/william/Documents/Fuego_catalogue/analysis/Fuego_exp_Energy_Matrix_final.csv", delimiter=',',skip_header=1)
#ITE = genfromtxt("/Users/william/Documents/Fuego_catalogue/analysis/Fuego_inf_tremor_Energy_Matrix_v8.csv", delimiter=',',skip_header=1)
#STE = genfromtxt("/Users/william/Documents/Fuego_catalogue/analysis/Fuego_seis_tremor_Energy_Matrix_v8.csv", delimiter=',',skip_header=1)


start_timestamp = UTCDateTime(2018,10,1).timestamp
end_timestamp = UTCDateTime(2019,1,1).timestamp

#%% explosions



En_t = Energy[:,0]
e_type = cat[:,8]

En_inf = Energy[:,1:30]
En_seis = Energy[:,30:38]

med_en = np.zeros(shape =(0,5))
med_ens = np.zeros(shape =(0,5))



for x in range(1,len(En_inf)-1):
    if start_timestamp < UTCDateTime(En_t[x]).timestamp < end_timestamp:
        
        nz = En_inf[x,:][np.nonzero(En_inf[x,:])]

        if len(nz) > 0:
            med_en_ev = np.max(nz)
            if med_en_ev < 1e10 and abs(En_t[x] - En_t[x-1]) < 60*60*2 and abs(En_t[x] - En_t[x+1]) < 60*60*2:
                med_en = np.lib.pad(med_en, ((0,1),(0,0)), 'constant', constant_values=(0))
                med_en[-1,0]=En_t[x]/1000000 # event time
                med_en[-1,1]=med_en_ev # event energy
                med_en[-1,2]=e_type[x] # event type
                med_en[-1,3] = abs(En_t[x] - En_t[x-1])/60 # before repose time
                med_en[-1,4] = abs(En_t[x] - En_t[x+1])/60 # after repose time
            
            
        

                
plt.figure()
plt.title("Before Repose and size")
plt.plot(med_en[:,3],med_en[:,1],'bx')
plt.xlabel("Repose time before event [min]")
plt.ylabel("Energy [MJ]")
         
plt.figure()
plt.title("after Repose and size")
plt.plot(med_en[:,4],med_en[:,1],'rx')            
plt.xlabel("Repose time after event [min]")            
plt.ylabel("Energy [MJ]")            
            
            
            
            
            
            
            
    