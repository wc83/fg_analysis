#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 13:43:35 2020

@author: root
"""


#kill the code 
#%%

import numpy as np
import matplotlib.pyplot as plt 
from obspy.core import read
from obspy.clients.earthworm import Client
from obspy import UTCDateTime
from obspy.signal.trigger import trigger_onset
from numpy import genfromtxt
from obspy import Stream
from scipy import integrate

#%% constants

pi=3.14159

#seismic
A=1
rhoE=2500
cE=2000

consts=2*pi*rhoE*cE*(1/A)

#acoustic
rhoA = 1 #air density (kg/m3) - may need changing slightly
cA = 330 #sound speed (m/s)
norm_P = 0 # 8*1E4 #normal air_pressure (N/m2)


consta = (2*pi)/(rhoA*cA)

# ip, port - ip's 138.253.113.19 or 138.253.112.23
client = Client('138.253.113.19', 16022) 

#%% Infrasound Station Distances

r_FV01 = 7865
r_FV02 = 6965
r_FV03 = 4593
r_FV04 = 2606
r_FV08 = 8017
r_FV11 = 9253

r_VF05 = 9308
r_VF04 = 5762
r_VF01 = 8041
r_VF06 = 7866
r_VF03 = 1088
r_VF02 = 8103

r_FG8_1 = 7979
r_FG8_2 = 8052
r_FG8_3 = 8016

r_FG10_1 = 8045
r_FG10_2 = 8045
r_FG10_3 = 8045

r_FG11_1 = 9235
r_FG11_2 = 9235
r_FG11_4 = 9235
r_FG11_5 = 9235

r_FG12_1 = 6931
r_FG12_2 = 6962
r_FG12_3 = 6980

r_FG13_1 = 10534
r_FG13_2 = 10488
r_FG13_3 = 10524

r_FG15_1 = 11737
r_FG15_2 = 11709
r_FG15_3 = 11659
r_FG15_4 = 11645
r_FG15_5 = 11692
r_FG15_6 = 11740

#%% Seismic distances

r_FG3 = 5817
r_FG8 = 8050
r_FG10 = 8045
r_FG11 = 9235
r_FG12 = 6962
r_FG13 = 10475
r_FG14 = 11712
r_FG16 = 4419

#%% Station calibrations - check network spreadsheet for each station

#for inf
c1 = 1.22E-04
c2 = 6.10E-05
c3 = 2.78E-04
c4 = 6.25E-06
c5 = 2.50E-06

#for seismics
c6 = 3.30E-09
c7 = 1.33E-09


#%% read catalogue of events

cat = genfromtxt("/Users/william/Documents/Fuego_catalogue/Final_catalogues/Fuego_explosions.csv", delimiter=',',skip_header=1)

#%%

Energies = np.zeros(shape=(len(cat),43))
Energies[:,0] = cat[:,0]
#%%

#open matrix





# go over each event and each station
s1= 0
e1= 0

#for s in range(s1,41): #loop over 32 infrasound stations (no FG11)
for s in range(0,1):
    print("station", s)
    # if statements to select correct station, channel, calibrations, catalogue_column
    if s == 0:
        r = r_FV01 
        sta = 'FV01' 
        cal=c1
        net='GI'
        cha='HDF'
        loc=''
        col=12
        
    if s == 1:
        r= r_FV02 
        sta = 'FV02'
        cal=c1
        net='GI'
        cha='HDF'
        loc=''
        col=13
        
    if s == 2:
        r= r_FV03 
        sta = 'FV03'
        cal=c1
        net='GI'
        cha='HDF'
        loc=''
        col=14
        
    if s == 3:
        r= r_FV04 
        sta='FV04'
        cal=c2
        net='GI'
        cha='HDF'
        loc=''
        col=15
        
    if s == 4:
        r= r_FV08 
        sta = 'FV08' 
        cal=c1
        net='GI'
        cha='HDF'
        loc=''
        col=16
        
    if s == 5:        
        r= r_VF01 
        sta='VF01'
        cal=c1
        net='XZ'
        cha='HDF'
        loc=''
        col=17
    
    if s == 6:
        r= r_VF02 
        sta='VF02'
        cal=c1
        net='XZ'
        cha='HDF'
        loc=''
        col=18
        
    if s == 7:
        r= r_VF03
        sta='VF03'
        cal=c1
        net='XZ'
        cha='HDF'
        loc=''
        col=19

    if s == 8:
        r= r_VF04 
        sta='VF04'
        cal=c1
        net='XZ'
        cha='HDF'
        loc=''
        col=20

    if s == 9:
        r= r_VF05
        sta='VF05'
        cal=c1
        net='XZ'
        cha='HDF'
        loc=''
        col=21
        
    if s == 10:
        r= r_VF06
        sta='VF06'
        cal=c1
        net='XZ'
        cha='HDF'
        loc=''
        col=22

    if s == 11:
        r= r_FG8_1 
        sta='FG8'
        cal=c3
        net='GI'
        cha='BDF'
        loc='01'
        col=23

    if s == 12:
        r= r_FG8_2 
        sta='FG8'
        cal=c4
        net='GI'
        cha='BDF'
        loc='02'
        col=24
        
    if s == 13:
        r= r_FG8_3 
        sta='FG8'
        cal=c3
        net='GI'
        cha='BDF'
        loc='03'
        col=25

    if s == 14:
        r= r_FG10_1
        sta='FG10'
        cal=c5
        net='GI'
        cha='BDF'
        loc='01'
        col=26
    
    if s == 15:
        r= r_FG10_2 
        sta='FG10'
        cal=c5
        net='GI'
        cha='BDF'
        loc='02'
        col=27
        
    if s == 16:
        r= r_FG10_3 
        sta='FG10'
        cal=c5
        net='GI'
        cha='BDF'
        loc='03'
        col=28
        
    if s == 17:
        r= r_FG11_1
        sta='FG11'
        cal=c5
        net='GI'
        cha='BDF'
        loc='01'
        col=29
    
    if s == 18:
        r= r_FG11_2 
        sta='FG11'
        cal=c5
        net='GI'
        cha='BDF'
        loc='02'
        col=30
        
    if s == 19:
        r= r_FG11_4 
        sta='FG11'
        cal=c5
        net='GI'
        cha='BDF'
        loc='04'
        col=31
        
    if s == 20:
        r= r_FG11_5 
        sta='FG11'
        cal=c5
        net='GI'
        cha='BDF'
        loc='05'
        col=32

    if s == 21:
        r=r_FG12_1 
        sta='FG12'
        cal=c3
        net='GI'
        cha='BDF'
        loc='01'
        col=33
        
    if s == 22:
        r=r_FG12_2 
        sta='FG12'
        cal=c3
        net='GI'
        cha='BDF'
        loc='02'
        col=34
        
    if s == 23:
        r= r_FG12_3
        sta='FG12'
        cal=c3
        net='GI'
        cha='BDF'
        loc='03'
        col=35
        
    if s == 24:
        r= r_FG13_1
        sta='FG13'
        cal=c3
        net='GI'
        cha='BDF'
        loc='01'
        col=36
        
    if s == 25:
        r=r_FG13_2 
        sta='FG13'
        cal=c3
        net='GI'
        cha='BDF'
        loc='02'
        col=37
        
    if s == 26:
        r=r_FG13_3
        sta='FG13'
        cal=c3
        net='GI'
        cha='BDF'
        loc='03'
        col=38

    if s == 27:
        r=r_FG15_1
        sta='FG15'
        cal=c3
        net='GI'
        cha='BDF'
        loc='01'
        col=39
        
    if s == 28:
        r=r_FG15_2 
        sta='FG15'
        cal=c3
        net='GI'
        cha='BDF'
        loc='02'
        col=40
        
    if s == 29:
        r=r_FG15_3 
        sta='FG15'
        cal=c3
        net='GI'
        cha='BDF'
        loc='03'
        col=41
    
    if s == 30:
        r= r_FG15_4 
        sta='FG15'
        cal=c3
        net='GI'
        cha='BDF'
        loc='04'
        col=42
        
    if s == 31:
        r=r_FG15_5 
        sta='FG15'
        cal=c3
        net='GI'
        cha='BDF'
        loc='05'
        col=43
        
    if s == 32:
        r=r_FG15_6 
        sta='FG15'
        cal=c3
        net='GI'
        cha='BDF'
        loc='06'  
        col=44
        
#    
    
    
    if s == 33:
        r=r_FG3 
        sta='FG3'
        cal=c7        # need to find true value
        net='GI'
        cha1='SHZ'
        cha2='SHE'
        cha3='SHN'
        loc='01'  
        col=45
        
    if s == 34:
        r=r_FG8
        sta='FG8'
        cal=c6
        net='GI'
        cha1='BHZ'
        cha2='BHE'
        cha3='BHN'
        loc='00'  
        col=46
    
    if s == 35:
        r=r_FG10 
        sta='FG10'
        cal=c7
        net='GI'
        cha1='BHZ'
        cha2='BHE'
        cha3='BHN'
        loc=''  
        col=47
    
    if s == 36:
        r=r_FG11 
        sta='FG11'
        cal=c7 
        net='GI'
        cha1='BHZ'
        cha2='BHE'
        cha3='BHN'
        loc=''  
        col=48
        
    if s == 37:
        r=r_FG12 
        sta='FG12'
        cal=c6
        net='GI'
        cha1='BHZ'
        cha2='BHE'
        cha3='BHN'
        loc='00'  
        col=49
        
    if s == 38:
        r=r_FG13 
        sta='FG13'
        cal=c7
        net='GI'
        cha1='BHZ'
        cha2='BHE'
        cha3='BHN'
        loc='00'  
        col=50
    
    if s == 39:
        r=r_FG14 
        sta='FG14'
        cal=c6
        net='GI'
        cha1='BHZ'
        cha2='BHE'
        cha3='BHN'
        loc=''  
        col=51
        
    if s == 40:
        r=r_FG16 
        sta='FG16'
        cal=c7 
        net='GI'
        cha1='BHZ'
        cha2='BHE'
        cha3='BHN'
        loc='00'  
        col=52
        
    
#    for e in range(0,len(cat)):
    for e in range(0,1):
        
#        print(sta,loc,e)

        ev_tick =cat[:,col]
        
#        #acoustic
        if s < 33:
            #if recorded by station...
            if ev_tick[e] == 1:
                    print("event found",e)
                    
#                try:
    #                print(s,e)
          
                # import waveform, trim, caibrate, filter
                
                    dur = cat[e,1]
                    t1 = UTCDateTime(cat[e,0]) - 45  # !! may need to make wider time buffer !!
                    t2 = t1 + dur + 90
                    
                    t_len = t2-t1
                    
                    
                    st_ac = Stream()
                    st_ac = client.get_waveforms(net, sta, loc, cha, t1 , t2)
                    
                    sr = st_ac[0].stats.sampling_rate
                    st_ac.detrend(type='linear')
                    st_ac.detrend(type='demean')
                    st_ac.filter(type='bandpass',freqmin=0.2, freqmax= 10)
                    
                    st_cal = cal*st_ac[0].data[int(5*sr):int(-5*sr)]
                    
                    #plot waveform
#                    plt.figure()
#                    points = np.linspace(0,t_len-10,len(st_cal))
#                    plt.plot(points,st_cal,'b')
#                    plt.title("{}: ".format(st_ac[0].stats.station)+" {}".format(st_ac[0].stats.starttime))
#                    plt.xlabel('Time [s]')
#                    plt.ylabel('Acoustic units')   
                
                # calculate Acoustic energy               
                    dl=len(st_cal)
                    p = np.linspace(0, dl/sr, num=dl) 
                
                    Delta_P = st_cal - norm_P #pressure difference 
                    P2 = np.square(Delta_P) #sqared pressure difference
                    
                    int_sq_press = integrate.cumtrapz(P2, p, initial=0) #inegral of squared pressure differences
                    
                    E_Atmos = (consta*(r**2)) * int_sq_press[-1]   #Energy Equation
                    
                    #save to matrix
                    
                    Energies[e,s+1] = E_Atmos
                    
                    
            #plot cumulative energy
    #                cum_E_Atmos = (consta*(r**2)) * int_sq_press / 1000000
    #                plt.figure()
    #                plt.plot(points, cum_E_Atmos,'r')
    #                plt.title("Explosion Energy = {:.2f} MJ".format(E_Atmos/1000000))
    #                plt.xlabel('Time [s]')
    #                plt.ylabel('Cumulative Energy [MJ]')
#                except:
#                    Energies[e,s+1] = 0
#                    print(sta, loc, cha, " Cant find event at: ", t1)
#            
        #seismics
        if s > 32:
            
                #if recorded by station...
            if ev_tick[e] == 1:
                try:
      
            # import waveform, trim, caibrate, filter
            
                    dur = cat[e,1]
                    t1 = UTCDateTime(cat[e,0]) - 45  # !! may beed to make wider time buffer !!
#                    print(t1)
                    
                    t2 = t1 + dur + 90
                    
#                    print(dur)
                    
                    t_len = t2-t1
                    
                    
                    st_z = Stream()
                    st_e = Stream()
                    st_n = Stream()
                    st_z = client.get_waveforms(net, sta, loc, cha1, t1 , t2)
#                    print("Vert in")
                    st_e = client.get_waveforms(net, sta, loc, cha2, t1 , t2)
#                    print("East in")
                    st_n = client.get_waveforms(net, sta, loc, cha3, t1 , t2)
#                    print("North in")
                    
                    sr = st_z[0].stats.sampling_rate
                    
                    st_z.detrend(type='linear')
                    st_z.detrend(type='demean')
                    st_z.filter(type='bandpass',freqmin=0.2, freqmax= 15)
                    
                    st_e.detrend(type='linear')
                    st_e.detrend(type='demean')
                    st_e.filter(type='bandpass',freqmin=0.2, freqmax= 15)
                    
                    st_n.detrend(type='linear')
                    st_n.detrend(type='demean')
                    st_n.filter(type='bandpass',freqmin=0.2, freqmax= 15)
                    
                    #ensure all combined components are same length
                    lens = [len(st_z[0].data),len(st_n[0].data),len(st_e[0].data)]
                    min_len = np.min(lens)
                    s_point = int(5*sr)
                    e_point = int(min_len - 5*sr)
                    
                    #combine Z,N,E and calibrate
                    st_cal= cal*np.sqrt(st_z[0].data[s_point:e_point]**2 + st_n[0].data[s_point:e_point]**2 + st_e[0].data[s_point:e_point]**2 )
                    
##                    plot seismic waveform
#                    points = np.linspace(0,t_len-10,len(st_cal))
#                    plt.figure()
#                    plt.plot(points,st_cal,'b')
#                    plt.title("{}: ".format(st_z[0].stats.station)+" {}".format(st_z[0].stats.starttime))
#                    plt.xlabel('Time [s]')
#                    plt.ylabel('Seismic units')
                                   
                    # calculate Seismic energy
                    
                    dl=len(st_cal)
                    p = np.linspace(0, dl/sr, num=dl) 
                            
                    y2= np.square(st_cal)
                                   
                    y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
                    EI = y_int2[-1]
    
                    E_seis= consts*(r*r)*EI                               
                
            # save to matrix (in J)
                    Energies[e,s+1] = E_seis
#                    print(E_seis)
                
                #plot cumulative energy
    #                cum_E_seis = consts*(r*r)*y_int2
    #                plt.figure()
    #                plt.plot(points, cum_E_seis,'r')
    #                plt.title("Explosion Energy = {:.2f} MJ".format(E_seis/1000000))
    #                plt.xlabel('Time [s]')
    #                plt.ylabel('Cumulative Energy [MJ]')
                except:
                    Energies[e,s+1] = 0
                    print(sta, "seis can't find event at: ", t1)



        


#np.savetxt("/Users/william/Documents/Fuego_catalogue/analysis/Fuego_exp_corrected_b.csv", Energies,delimiter=",",header="Time,FV01,FV02,FV03,FV04,FV08,VF01,VF02,VF03,VF04,VF05,VF06,FG8A,FG8B,FG8C,FG10A,FG10B,FG10C,FG11A,FG11B,FG11C,FG11D,FG12A,FG12B,FG12C,FG13A,FG13B,FG13C,FG15A,FG15B,FG15C,FG15D,FG15E,FG15F,FG3S,FG8S,FG10S,FG11S,FG12S,FG13S,FG14S,FG16S")


#%%
                    
                    
#%%

Energy = genfromtxt("/Users/william/Documents/Fuego_catalogue/analysis/Fuego_exp_corrected.csv", delimiter=',',skip_header=1)
                
inf_only = np.zeros(shape=(0,4))

fg12i = cat[:,34]
fg12s = cat[:,49]

r=r_FG12 
sta='FG12'
cal=c6
net='GI'
cha1='BHZ'
cha2='BHE'
cha3='BHN'
loc='00'  

num_io = 0
for x in range(0,len(cat)):
    if fg12i[x] == 1 and fg12s[x] == 0:
        try:
            inf_only = np.lib.pad(inf_only, ((0,1),(0,0)), 'constant', constant_values=(0))
            inf_only[num_io,0] = cat[x,0]
            inf_only[num_io,1] = cat[x,8]
            
            dur = cat[x,1]
            
            t1 = UTCDateTime(cat[x,0]) - 45  # !! may beed to make wider time buffer !!            
            t2 = t1 + dur + 90
                        
            t_len = t2-t1
                        
                        
            st_z = Stream()
            st_e = Stream()
            st_n = Stream()
            st_z = client.get_waveforms(net, sta, loc, cha1, t1 , t2)
    #                    print("Vert in")
            st_e = client.get_waveforms(net, sta, loc, cha2, t1 , t2)
    #                    print("East in")
            st_n = client.get_waveforms(net, sta, loc, cha3, t1 , t2)
    #                    print("North in")
            
            sr = st_z[0].stats.sampling_rate
            
            st_z.detrend(type='linear')
            st_z.detrend(type='demean')
            st_z.filter(type='bandpass',freqmin=0.2, freqmax= 15)
            
            st_e.detrend(type='linear')
            st_e.detrend(type='demean')
            st_e.filter(type='bandpass',freqmin=0.2, freqmax= 15)
            
            st_n.detrend(type='linear')
            st_n.detrend(type='demean')
            st_n.filter(type='bandpass',freqmin=0.2, freqmax= 15)
            
            #ensure all combined components are same length
            lens = [len(st_z[0].data),len(st_n[0].data),len(st_e[0].data)]
            min_len = np.min(lens)
            s_point = int(5*sr)
            e_point = int(min_len - 5*sr)
            
            #combine Z,N,E and calibrate
            st_cal= cal*np.sqrt(st_z[0].data[s_point:e_point]**2 + st_n[0].data[s_point:e_point]**2 + st_e[0].data[s_point:e_point]**2 )
            
                           
            # calculate Seismic energy
            
            dl=len(st_cal)
            p = np.linspace(0, dl/sr, num=dl) 
                    
            y2= np.square(st_cal)
                           
            y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
            EI = y_int2[-1]
    
            E_seis= consts*(r*r)*EI  
    
            inf_only[num_io,2] = Energy[x,23]
            inf_only[num_io,3] = E_seis 
            
            num_io += 1
            print("seis found", (100*x)/len(cat),"%")
        except:
            print("seis not found", (100*x)/len(cat),"%")
            


        
np.savetxt("/Users/william/Documents/Fuego_catalogue/analysis/Fuego_fg12_exp_noGCA.csv", inf_only,delimiter=",",header="Time,type,Inf,Seis")
        




#%% read catalogue of events

cat_st = genfromtxt("/Users/william/Documents/Fuego_catalogue/Final_catalogues/Fuego_trem_Seis_final.csv", delimiter=',',skip_header=1)
cat_it = genfromtxt("/Users/william/Documents/Fuego_catalogue/Final_catalogues/Fuego_trem_Inf_final.csv", delimiter=',',skip_header=1)

print("cats in")

Energies_st = np.zeros(shape=(len(cat_st),3))
Energies_st[:,0] = cat_st[:,0]

Energies_it = np.zeros(shape=(len(cat_it),3))
Energies_it[:,0] = cat_it[:,0]





r=r_FG12_2 
sta='FG12'

cal_i=c3
net='GI'
cha='BDF'
loc_i='02'
col_i=26

cal_s=c6
cha1='BHZ'
cha2='BHE'
cha3='BHN'
loc_s='00'  
col_s=13


# seismic trem

for e in range(0,len(cat_st)):
        print(e)

        ev_tick = cat_st[:,col_s]

        if ev_tick[e] == 1:
                    
                try:
                    #acoustics
                
                    dur = cat_st[e,1]
                    t1 = UTCDateTime(cat_st[e,0]) - 45  # !! may need to make wider time buffer !!
                    t2 = t1 + dur + 90
                    
                    t_len = t2-t1
                    
                    
                    st_ac = Stream()
                    st_ac = client.get_waveforms(net, sta, loc_i, cha, t1 , t2)
                    
                    sr = st_ac[0].stats.sampling_rate
                    st_ac.detrend(type='linear')
                    st_ac.detrend(type='demean')
                    st_ac.filter(type='bandpass',freqmin=0.2, freqmax= 5)
                    
                    st_cal = cal_i*st_ac[0].data[int(5*sr):int(-5*sr)]
                    
                
                # calculate Acoustic energy               
                    dl=len(st_cal)
                    p = np.linspace(0, dl/sr, num=dl) 
                
                    Delta_P = st_cal - norm_P #pressure difference 
                    P2 = np.square(Delta_P) #sqared pressure difference
                    
                    int_sq_press = integrate.cumtrapz(P2, p, initial=0) #inegral of squared pressure differences
                    
                    E_Atmos = (consta*(r**2)) * int_sq_press[-1]   #Energy Equation
                    
                    #save to matrix
                    
                    Energies_st[e,1] = E_Atmos
                    
                    #seismics
                    
                    st_z = Stream()
                    st_e = Stream()
                    st_n = Stream()
                    st_z = client.get_waveforms(net, sta, loc_s, cha1, t1 , t2)
#                    print("Vert in")
                    st_e = client.get_waveforms(net, sta, loc_s, cha2, t1 , t2)
#                    print("East in")
                    st_n = client.get_waveforms(net, sta, loc_s, cha3, t1 , t2)
#                    print("North in")
                    
                    sr = st_z[0].stats.sampling_rate
                    
                    st_z.detrend(type='linear')
                    st_z.detrend(type='demean')
                    st_z.filter(type='bandpass',freqmin=0.2, freqmax= 10)
                    
                    st_e.detrend(type='linear')
                    st_e.detrend(type='demean')
                    st_e.filter(type='bandpass',freqmin=0.2, freqmax= 10)
                    
                    st_n.detrend(type='linear')
                    st_n.detrend(type='demean')
                    st_n.filter(type='bandpass',freqmin=0.2, freqmax= 10)
                    
                    #ensure all combined components are same length
                    lens = [len(st_z[0].data),len(st_n[0].data),len(st_e[0].data)]
                    min_len = np.min(lens)
                    s_point = int(5*sr)
                    e_point = int(min_len - 5*sr)
                    
                    #combine Z,N,E and calibrate
                    st_cal= cal_s*np.sqrt(st_z[0].data[s_point:e_point]**2 + st_n[0].data[s_point:e_point]**2 + st_e[0].data[s_point:e_point]**2 )
                    
                                   
                    # calculate Seismic energy
                    
                    dl=len(st_cal)
                    p = np.linspace(0, dl/sr, num=dl) 
                            
                    y2= np.square(st_cal)
                                   
                    y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
                    EI = y_int2[-1]
    
                    E_seis= consts*(r*r)*EI                               
                
            # save to matrix (in J)
                    Energies_st[e,2] = E_seis

                except:
                    Energies[e,1] = 0
                    Energies[e,2] = 0
                    print(" Cant find event at: ", t1)


np.savetxt("/Users/william/Documents/Fuego_catalogue/analysis/Fuego_seis_trem_corrected.csv", Energies_st,delimiter=",",header="Time,Inf,Seis")

# Acoustic trem

for e in range(0,len(cat_it)):
        print(e)

        ev_tick = cat_it[:,col_i]

        if ev_tick[e] == 1:
                    
                try:
                    #acoustics
                
                    dur = cat_it[e,1]
                    t1 = UTCDateTime(cat_it[e,0]) - 45  # !! may need to make wider time buffer !!
                    t2 = t1 + dur + 90
                    
                    t_len = t2-t1
                    
                    
                    st_ac = Stream()
                    st_ac = client.get_waveforms(net, sta, loc_i, cha, t1 , t2)
                    
                    sr = st_ac[0].stats.sampling_rate
                    st_ac.detrend(type='linear')
                    st_ac.detrend(type='demean')
                    st_ac.filter(type='bandpass',freqmin=0.2, freqmax= 5)
                    
                    st_cal = cal_i*st_ac[0].data[int(5*sr):int(-5*sr)]
                    
                
                # calculate Acoustic energy               
                    dl=len(st_cal)
                    p = np.linspace(0, dl/sr, num=dl) 
                
                    Delta_P = st_cal - norm_P #pressure difference 
                    P2 = np.square(Delta_P) #sqared pressure difference
                    
                    int_sq_press = integrate.cumtrapz(P2, p, initial=0) #inegral of squared pressure differences
                    
                    E_Atmos = (consta*(r**2)) * int_sq_press[-1]   #Energy Equation
                    
                    #save to matrix
                    
                    Energies_it[e,1] = E_Atmos
                    
                    #seismics
                    
                    st_z = Stream()
                    st_e = Stream()
                    st_n = Stream()
                    st_z = client.get_waveforms(net, sta, loc_s, cha1, t1 , t2)
#                    print("Vert in")
                    st_e = client.get_waveforms(net, sta, loc_s, cha2, t1 , t2)
#                    print("East in")
                    st_n = client.get_waveforms(net, sta, loc_s, cha3, t1 , t2)
#                    print("North in")
                    
                    sr = st_z[0].stats.sampling_rate
                    
                    st_z.detrend(type='linear')
                    st_z.detrend(type='demean')
                    st_z.filter(type='bandpass',freqmin=0.2, freqmax= 10)
                    
                    st_e.detrend(type='linear')
                    st_e.detrend(type='demean')
                    st_e.filter(type='bandpass',freqmin=0.2, freqmax= 10)
                    
                    st_n.detrend(type='linear')
                    st_n.detrend(type='demean')
                    st_n.filter(type='bandpass',freqmin=0.2, freqmax= 10)
                    
                    #ensure all combined components are same length
                    lens = [len(st_z[0].data),len(st_n[0].data),len(st_e[0].data)]
                    min_len = np.min(lens)
                    s_point = int(5*sr)
                    e_point = int(min_len - 5*sr)
                    
                    #combine Z,N,E and calibrate
                    st_cal= cal_s*np.sqrt(st_z[0].data[s_point:e_point]**2 + st_n[0].data[s_point:e_point]**2 + st_e[0].data[s_point:e_point]**2 )
                    
                                   
                    # calculate Seismic energy
                    
                    dl=len(st_cal)
                    p = np.linspace(0, dl/sr, num=dl) 
                            
                    y2= np.square(st_cal)
                                   
                    y_int2 = integrate.cumtrapz(y2, p, initial=0) # total energy
                    EI = y_int2[-1]
    
                    E_seis= consts*(r*r)*EI                               
                
            # save to matrix (in J)
                    Energies_it[e,2] = E_seis

                except:
                    Energies[e,1] = 0
                    Energies[e,2] = 0
                    print(" Cant find event at: ", t1)


np.savetxt("/Users/william/Documents/Fuego_catalogue/analysis/Fuego_inf_trem_corrected.csv", Energies_it,delimiter=",",header="Time,Inf,Seis")


















