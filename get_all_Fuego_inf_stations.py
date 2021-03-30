#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 08:15:32 2020

@author: root
"""

def get_all_Fuego_inf_stations(timestamp):
    
    from obspy.clients.earthworm import Client
    from obspy import UTCDateTime
    from obspy import Stream
    import numpy as np

    
    station_activity = np.zeros(shape=(0,3))
    
    num = 0
    
    t1 = UTCDateTime(timestamp) #the format is year:day_of_the_year:month
    t2 = t1 + 24*60*60 + 1200
#%%
    sta = 'VF01' # STATION VF01
    cha = 'HDF' # CHANNEL - inf
    net = 'XZ'  # Fuego volcano
    loc = ''    # location
    
    client = Client('138.253.112.23', 16022)
    
    t1_2s = UTCDateTime(1526774400)
    t2_2s = t1_2s + 2
    
    st_blank = Stream()
    st_blank = client.get_waveforms(net, sta, loc, cha, t1_2s , t2_2s)
    st_blank.detrend(type='linear')
    st_blank.detrend(type='demean')
    st_blank.filter(type='bandpass',freqmin=0.5, freqmax=5)
    
    
#%% VF01    
    try:  
        ID = 1
        sta = 'VF01' # STATION VF01
        cha = 'HDF' # CHANNEL - inf
        net = 'XZ'  # Fuego volcano
        loc = ''    # location, 
            
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        sr = st1[0].stats.sampling_rate
        
        if sum(abs(st1[0].data)) > 10 and st1[0].stats.npts > 7200*sr:
            station_activity = np.lib.pad(station_activity, ((0,1),(0,0)), 'constant', constant_values=(0))
    
            station_activity[num][0] = ID
            station_activity[num][1] = st1[0].stats.starttime
            station_activity[num][2] = st1[0].stats.endtime
            num+=1
        
        else:
            st1 = st_blank
            
    except:
        st1 = st_blank
            
            
        
#%% VF02   
    try:  
        ID = 2
        sta = 'VF02' # STATION VF01
        cha = 'HDF' # CHANNEL - inf
        net = 'XZ'  # Fuego volcano
        loc = ''    # location, 
    
        st2 = Stream()
        st2 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        st2.detrend(type='linear')
        st2.detrend(type='demean')
        
        break_test=st2
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        sr = st2[0].stats.sampling_rate
        
        if sum(abs(st2[0].data)) > 10 and st2[0].stats.npts > 7200*sr:
            station_activity = np.lib.pad(station_activity, ((0,1),(0,0)), 'constant', constant_values=(0))
            station_activity[num][0] = ID
            station_activity[num][1] = st2[0].stats.starttime
            station_activity[num][2] = st2[0].stats.endtime
            num+=1
            
        
        else:
            st2 = st_blank
            
    except:
        st2 = st_blank        
        
        
#%% VF03   
    try:  
        ID = 3
        sta = 'VF03' # STATION VF01
        cha = 'HDF' # CHANNEL - inf
        net = 'XZ'  # Fuego volcano
        loc = ''    # location, 
        
        st3 = Stream()
        st3 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        st3.detrend(type='linear')
        st3.detrend(type='demean')
        
        break_test=st3
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        sr = st3[0].stats.sampling_rate
        
        if sum(abs(st3[0].data)) > 10 and st3[0].stats.npts > 7200*sr:
            station_activity = np.lib.pad(station_activity, ((0,1),(0,0)), 'constant', constant_values=(0))
            station_activity[num][0] = ID
            station_activity[num][1] = st3[0].stats.starttime
            station_activity[num][2] = st3[0].stats.endtime
            num+=1
            
        
        else:
            st3 = st_blank
            
    except:
        st3 = st_blank        
        
        
        
 #%% VF04  
    try:  
        ID = 4
        sta = 'VF04' # STATION VF01
        cha = 'HDF' # CHANNEL - inf
        net = 'XZ'  # Fuego volcano
        loc = ''    # location, 

        st4 = Stream()
        st4 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        st4.detrend(type='linear')
        st4.detrend(type='demean')
        
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        sr = st4[0].stats.sampling_rate
        
        if sum(abs(st4[0].data)) > 10 and st4[0].stats.npts > 7200*sr:
            station_activity = np.lib.pad(station_activity, ((0,1),(0,0)), 'constant', constant_values=(0))
            station_activity[num][0] = ID
            station_activity[num][1] = st4[0].stats.starttime
            station_activity[num][2] = st4[0].stats.endtime
            num+=1
            
        
        else:
            st4 = st_blank
            
    except:
        st4 = st_blank        
 
               
 #%% VF05   
    try:  
        ID = 5
        sta = 'VF05' # STATION VF01
        cha = 'HDF' # CHANNEL - inf
        net = 'XZ'  # Fuego volcano
        loc = ''    # location, 

        st5 = Stream()
        st5 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        st5.detrend(type='linear')
        st5.detrend(type='demean')
        
        break_test=st5
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        sr = st5[0].stats.sampling_rate
        
        
        if sum(abs(st5[0].data)) > 10 and st5[0].stats.npts > 7200*sr:
            station_activity = np.lib.pad(station_activity, ((0,1),(0,0)), 'constant', constant_values=(0))
            station_activity[num][0] = ID
            station_activity[num][1] = st5[0].stats.starttime
            station_activity[num][2] = st5[0].stats.endtime
            num+=1
            
        
        else:
            st5 = st_blank
            
    except:
        st5 = st_blank        
        
        
 #%% VF06  
    try:  
        ID = 6
        sta = 'VF06' # STATION VF01
        cha = 'HDF' # CHANNEL - inf
        net = 'XZ'  # Fuego volcano
        loc = ''    # location, 

        st6 = Stream()
        st6 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        st6.detrend(type='linear')
        st6.detrend(type='demean')
        
        break_test=st6
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        sr = st6[0].stats.sampling_rate
        
        
        if sum(abs(st6[0].data)) > 10 and st6[0].stats.npts > 7200*sr:
            station_activity = np.lib.pad(station_activity, ((0,1),(0,0)), 'constant', constant_values=(0))
            station_activity[num][0] = ID
            station_activity[num][1] = st6[0].stats.starttime
            station_activity[num][2] = st6[0].stats.endtime
            num+=1
            
        
        else:
            st6 = st_blank
            
    except:
        st6 = st_blank        
        





#%% FG8

    try:  
        ID = 7
        sta = 'FG8' # STATION VF01
        cha = 'BDF' # CHANNEL - inf
        net = 'GI'  # Fuego volcano
        loc = '01'    # location, 

        st7 = Stream()
        st7 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        st7.detrend(type='linear')
        st7.detrend(type='demean')
        
        break_test=st7
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        sr = st7[0].stats.sampling_rate
        
        
        if sum(abs(st7[0].data)) > 10 and st7[0].stats.npts > 7200*sr:
            station_activity = np.lib.pad(station_activity, ((0,1),(0,0)), 'constant', constant_values=(0))
            station_activity[num][0] = ID
            station_activity[num][1] = st7[0].stats.starttime
            station_activity[num][2] = st7[0].stats.endtime
            num+=1
            
        
        else:
            st7 = st_blank
            
    except:
        st7 = st_blank   


    try:  
        ID = 8
        sta = 'FG8' # STATION VF01
        cha = 'BDF' # CHANNEL - inf
        net = 'GI'  # Fuego volcano
        loc = '02'    # location, 

        st8 = Stream()
        st8 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        st8.detrend(type='linear')
        st8.detrend(type='demean')
        
        break_test=st8
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        sr = st8[0].stats.sampling_rate
        
        
        if sum(abs(st8[0].data)) > 10 and st8[0].stats.npts > 7200*sr:
            station_activity = np.lib.pad(station_activity, ((0,1),(0,0)), 'constant', constant_values=(0))
            station_activity[num][0] = ID
            station_activity[num][1] = st8[0].stats.starttime
            station_activity[num][2] = st8[0].stats.endtime
            num+=1
            
        
        else:
            st8 = st_blank
            
    except:
        st8 = st_blank   
        
    
    
    try:  
        ID = 9
        sta = 'FG8' # STATION VF01
        cha = 'BDF' # CHANNEL - inf
        net = 'GI'  # Fuego volcano
        loc = '03'    # location, 

        st9 = Stream()
        st9 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        st9.detrend(type='linear')
        st9.detrend(type='demean')
        
        break_test=st9
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        sr = st9[0].stats.sampling_rate
        
        
        if sum(abs(st9[0].data)) > 10 and st9[0].stats.npts > 7200*sr:
            station_activity = np.lib.pad(station_activity, ((0,1),(0,0)), 'constant', constant_values=(0))
            station_activity[num][0] = ID
            station_activity[num][1] = st9[0].stats.starttime
            station_activity[num][2] = st9[0].stats.endtime
            num+=1
            
        
        else:
            st9 = st_blank
            
    except:
        st9 = st_blank   

      
#%% FG10

    try:  
        ID = 10
        sta = 'FG10' # STATION VF01
        cha = 'BDF' # CHANNEL - inf
        net = 'GI'  # Fuego volcano
        loc = '01'    # location, 

        st10 = Stream()
        st10 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        st10.detrend(type='linear')
        st10.detrend(type='demean')
        
        break_test=st10
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        sr = st10[0].stats.sampling_rate
        
        
        if sum(abs(st10[0].data)) > 10 and st10[0].stats.npts > 7200*sr:
            station_activity = np.lib.pad(station_activity, ((0,1),(0,0)), 'constant', constant_values=(0))
            station_activity[num][0] = ID
            station_activity[num][1] = st10[0].stats.starttime
            station_activity[num][2] = st10[0].stats.endtime
            num+=1
            
        
        else:
            st10 = st_blank
            
    except:
        st10 = st_blank   


    try:  
        ID = 11
        sta = 'FG10' # STATION VF01
        cha = 'BDF' # CHANNEL - inf
        net = 'GI'  # Fuego volcano
        loc = '02'    # location, 

        st11 = Stream()
        st11 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        st11.detrend(type='linear')
        st11.detrend(type='demean')
        
        break_test=st11
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        sr = st11[0].stats.sampling_rate
        
        
        if sum(abs(st11[0].data)) > 10 and st11[0].stats.npts > 7200*sr:
            station_activity = np.lib.pad(station_activity, ((0,1),(0,0)), 'constant', constant_values=(0))
            station_activity[num][0] = ID
            station_activity[num][1] = st11[0].stats.starttime
            station_activity[num][2] = st11[0].stats.endtime
            num+=1
            
        
        else:
            st11 = st_blank
            
    except:
        st11 = st_blank   
        
    
    
    try:  
        ID = 12
        sta = 'FG10' # STATION VF01
        cha = 'BDF' # CHANNEL - inf
        net = 'GI'  # Fuego volcano
        loc = '03'    # location, 

        st12 = Stream()
        st12 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        st12.detrend(type='linear')
        st12.detrend(type='demean')
        
        break_test=st12
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        sr = st12[0].stats.sampling_rate
        
        
        if sum(abs(st12[0].data)) > 10 and st12[0].stats.npts > 7200*sr:
            station_activity = np.lib.pad(station_activity, ((0,1),(0,0)), 'constant', constant_values=(0))
            station_activity[num][0] = ID
            station_activity[num][1] = st12[0].stats.starttime
            station_activity[num][2] = st12[0].stats.endtime
            num+=1
            
        
        else:
            st12 = st_blank
            
    except:
        st12 = st_blank   





#%% FG11

    try:  
        ID = 13
        sta = 'FG11' # STATION VF01
        cha = 'HDF' # CHANNEL - inf
        net = 'GI'  # Fuego volcano
        loc = '01'    # location, 

        st13 = Stream()
        st13 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        st13.detrend(type='linear')
        st13.detrend(type='demean')
        
        break_test=st13
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        sr = st13[0].stats.sampling_rate
        
        
        if sum(abs(st13[0].data)) > 10 and st13[0].stats.npts > 7200*sr:
            station_activity = np.lib.pad(station_activity, ((0,1),(0,0)), 'constant', constant_values=(0))
            station_activity[num][0] = ID
            station_activity[num][1] = st13[0].stats.starttime
            station_activity[num][2] = st13[0].stats.endtime
            num+=1
            
        
        else:
            st13 = st_blank
            
    except:
        st13 = st_blank   


    try:  
        ID = 14
        sta = 'FG11' # STATION VF01
        cha = 'HDF' # CHANNEL - inf
        net = 'GI'  # Fuego volcano
        loc = '02'    # location, 

        st14 = Stream()
        st14 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        st14.detrend(type='linear')
        st14.detrend(type='demean')
        
        break_test=st14
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        sr = st14[0].stats.sampling_rate
        
        
        if sum(abs(st14[0].data)) > 10 and st14[0].stats.npts > 7200*sr:
            station_activity = np.lib.pad(station_activity, ((0,1),(0,0)), 'constant', constant_values=(0))
            station_activity[num][0] = ID
            station_activity[num][1] = st14[0].stats.starttime
            station_activity[num][2] = st14[0].stats.endtime
            num+=1
            
        
        else:
            st14 = st_blank
            
    except:
        st14 = st_blank   
        
    
    
    try:  
        ID = 15
        sta = 'FG11' # STATION VF01
        cha = 'HDF' # CHANNEL - inf
        net = 'GI'  # Fuego volcano
        loc = '04'    # location, 

        st15 = Stream()
        st15 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        st15.detrend(type='linear')
        st15.detrend(type='demean')
        
        break_test=st15
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        sr = st15[0].stats.sampling_rate
        
        
        if sum(abs(st15[0].data)) > 10 and st15[0].stats.npts > 7200*sr:
            station_activity = np.lib.pad(station_activity, ((0,1),(0,0)), 'constant', constant_values=(0))
            station_activity[num][0] = ID
            station_activity[num][1] = st15[0].stats.starttime
            station_activity[num][2] = st15[0].stats.endtime
            num+=1
            
        
        else:
            st15 = st_blank
            
    except:
        st15 = st_blank   
      
    try:  
        ID = 16
        sta = 'FG11' # STATION VF01
        cha = 'HDF' # CHANNEL - inf
        net = 'GI'  # Fuego volcano
        loc = '05'    # location, 

        st16 = Stream()
        st16 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        st16.detrend(type='linear')
        st16.detrend(type='demean')
        
        break_test=st16
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        sr = st16[0].stats.sampling_rate
        
        
        if sum(abs(st16[0].data)) > 10 and st16[0].stats.npts > 7200*sr:
            station_activity = np.lib.pad(station_activity, ((0,1),(0,0)), 'constant', constant_values=(0))
            station_activity[num][0] = ID
            station_activity[num][1] = st16[0].stats.starttime
            station_activity[num][2] = st16[0].stats.endtime
            num+=1
            
        
        else:
            st16 = st_blank
            
    except:
        st16 = st_blank 



#%% FG12

    try:  
        ID = 17
        sta = 'FG12' # STATION VF01
        cha = 'BDF' # CHANNEL - inf
        net = 'GI'  # Fuego volcano
        loc = '01'    # location, 

        st17 = Stream()
        st17 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        st17.detrend(type='linear')
        st17.detrend(type='demean')
        
        break_test=st17
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        sr = st17[0].stats.sampling_rate
        
        
        if sum(abs(st17[0].data)) > 10 and st17[0].stats.npts > 7200*sr:
            station_activity = np.lib.pad(station_activity, ((0,1),(0,0)), 'constant', constant_values=(0))
            station_activity[num][0] = ID
            station_activity[num][1] = st17[0].stats.starttime
            station_activity[num][2] = st17[0].stats.endtime
            num+=1
            
        
        else:
            st17 = st_blank
            
    except:
        st17 = st_blank   


    try:  
        ID = 18
        sta = 'FG12' # STATION VF01
        cha = 'BDF' # CHANNEL - inf
        net = 'GI'  # Fuego volcano
        loc = '02'    # location, 

        st18 = Stream()
        st18 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        st18.detrend(type='linear')
        st18.detrend(type='demean')
        
        break_test=st18
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        sr = st18[0].stats.sampling_rate
        
        
        if sum(abs(st18[0].data)) > 10 and st18[0].stats.npts > 7200*sr:
            station_activity = np.lib.pad(station_activity, ((0,1),(0,0)), 'constant', constant_values=(0))
            station_activity[num][0] = ID
            station_activity[num][1] = st18[0].stats.starttime
            station_activity[num][2] = st18[0].stats.endtime
            num+=1
            
        
        else:
            st18 = st_blank
            
    except:
        st18 = st_blank   
        
    
    
    try:  
        ID = 19
        sta = 'FG12' # STATION VF01
        cha = 'BDF' # CHANNEL - inf
        net = 'GI'  # Fuego volcano
        loc = '03'    # location, 

        st19 = Stream()
        st19 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        st19.detrend(type='linear')
        st19.detrend(type='demean')
        
        break_test=st19
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        sr = st19[0].stats.sampling_rate
        
        
        if sum(abs(st19[0].data)) > 10 and st19[0].stats.npts > 7200*sr:
            station_activity = np.lib.pad(station_activity, ((0,1),(0,0)), 'constant', constant_values=(0))
            station_activity[num][0] = ID
            station_activity[num][1] = st19[0].stats.starttime
            station_activity[num][2] = st19[0].stats.endtime
            num+=1
            
        
        else:
            st19 = st_blank
            
    except:
        st19 = st_blank   


#%% FG13

    try:  
        ID = 20
        sta = 'FG13' # STATION VF01
        cha = 'BDF' # CHANNEL - inf
        net = 'GI'  # Fuego volcano
        loc = '01'    # location, 

        st20 = Stream()
        st20 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        st20.detrend(type='linear')
        st20.detrend(type='demean')
        
        break_test=st20
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        sr = st20[0].stats.sampling_rate
        
        
        if sum(abs(st20[0].data)) > 10 and st20[0].stats.npts > 7200*sr:
            station_activity = np.lib.pad(station_activity, ((0,1),(0,0)), 'constant', constant_values=(0))
            station_activity[num][0] = ID
            station_activity[num][1] = st20[0].stats.starttime
            station_activity[num][2] = st20[0].stats.endtime
            num+=1
            
        
        else:
            st20 = st_blank
            
    except:
        st20 = st_blank   


    try:  
        ID = 21
        sta = 'FG13' # STATION VF01
        cha = 'BDF' # CHANNEL - inf
        net = 'GI'  # Fuego volcano
        loc = '02'    # location, 

        st21 = Stream()
        st21 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        st21.detrend(type='linear')
        st21.detrend(type='demean')
        
        break_test=st21
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        sr = st21[0].stats.sampling_rate
        
        
        if sum(abs(st21[0].data)) > 10 and st21[0].stats.npts > 7200*sr:
            station_activity = np.lib.pad(station_activity, ((0,1),(0,0)), 'constant', constant_values=(0))
            station_activity[num][0] = ID
            station_activity[num][1] = st21[0].stats.starttime
            station_activity[num][2] = st21[0].stats.endtime
            num+=1
            
        
        else:
            st21 = st_blank
            
    except:
        st21 = st_blank   
        
    
    
    try:  
        ID = 22
        sta = 'FG13' # STATION VF01
        cha = 'BDF' # CHANNEL - inf
        net = 'GI'  # Fuego volcano
        loc = '03'    # location, 

        st22 = Stream()
        st22 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        st22.detrend(type='linear')
        st22.detrend(type='demean')
        
        break_test=st22
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        sr = st22[0].stats.sampling_rate
        
        
        if sum(abs(st22[0].data)) > 10 and st22[0].stats.npts > 7200*sr:
            station_activity = np.lib.pad(station_activity, ((0,1),(0,0)), 'constant', constant_values=(0))
            station_activity[num][0] = ID
            station_activity[num][1] = st22[0].stats.starttime
            station_activity[num][2] = st22[0].stats.endtime
            num+=1
            
        
        else:
            st22 = st_blank
            
    except:
        st22 = st_blank   



#%% FG15

    try:  
        ID = 23
        sta = 'FG15' # STATION VF01
        cha = 'BDF' # CHANNEL - inf
        net = 'GI'  # Fuego volcano
        loc = '01'    # location, 

        st23 = Stream()
        st23 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        st23.detrend(type='linear')
        st23.detrend(type='demean')
        
        break_test=st23
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        sr = st23[0].stats.sampling_rate
        
        
        if sum(abs(st23[0].data)) > 10 and st23[0].stats.npts > 7200*sr:
            station_activity = np.lib.pad(station_activity, ((0,1),(0,0)), 'constant', constant_values=(0))
            station_activity[num][0] = ID
            station_activity[num][1] = st23[0].stats.starttime
            station_activity[num][2] = st23[0].stats.endtime
            num+=1
            
        
        else:
            st23 = st_blank
            
    except:
        st23 = st_blank   


    try:  
        ID = 24
        sta = 'FG15' # STATION VF01
        cha = 'BDF' # CHANNEL - inf
        net = 'GI'  # Fuego volcano
        loc = '02'    # location, 

        st24 = Stream()
        st24 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        st24.detrend(type='linear')
        st24.detrend(type='demean')
        
        break_test=st24
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        sr = st24[0].stats.sampling_rate
        
        
        if sum(abs(st24[0].data)) > 10 and st24[0].stats.npts > 7200*sr:
            station_activity = np.lib.pad(station_activity, ((0,1),(0,0)), 'constant', constant_values=(0))
            station_activity[num][0] = ID
            station_activity[num][1] = st24[0].stats.starttime
            station_activity[num][2] = st24[0].stats.endtime
            num+=1
            
        
        else:
            st24 = st_blank
            
    except:
        st24 = st_blank   
        
    
    
    try:  
        ID = 25
        sta = 'FG15' # STATION VF01
        cha = 'BDF' # CHANNEL - inf
        net = 'GI'  # Fuego volcano
        loc = '03'    # location, 

        st25 = Stream()
        st25 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        st25.detrend(type='linear')
        st25.detrend(type='demean')
        
        break_test=st25
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        sr = st25[0].stats.sampling_rate
        
        
        if sum(abs(st25[0].data)) > 10 and st25[0].stats.npts > 7200*sr:
            station_activity = np.lib.pad(station_activity, ((0,1),(0,0)), 'constant', constant_values=(0))
            station_activity[num][0] = ID
            station_activity[num][1] = st25[0].stats.starttime
            station_activity[num][2] = st25[0].stats.endtime
            num+=1
            
        
        else:
            st25 = st_blank
            
    except:
        st25 = st_blank   


    try:  
        ID = 26
        sta = 'FG15' # STATION VF01
        cha = 'BDF' # CHANNEL - inf
        net = 'GI'  # Fuego volcano
        loc = '04'    # location, 

        st26 = Stream()
        st26 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        st26.detrend(type='linear')
        st26.detrend(type='demean')
        
        break_test=st26
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        sr = st26[0].stats.sampling_rate
        
        
        if sum(abs(st26[0].data)) > 10 and st26[0].stats.npts > 7200*sr:
            station_activity = np.lib.pad(station_activity, ((0,1),(0,0)), 'constant', constant_values=(0))
            station_activity[num][0] = ID
            station_activity[num][1] = st26[0].stats.starttime
            station_activity[num][2] = st26[0].stats.endtime
            num+=1
            
        
        else:
            st26 = st_blank
            
    except:
        st26 = st_blank   


    try:  
        ID = 27
        sta = 'FG15' # STATION VF01
        cha = 'BDF' # CHANNEL - inf
        net = 'GI'  # Fuego volcano
        loc = '05'    # location, 

        st27 = Stream()
        st27 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        st27.detrend(type='linear')
        st27.detrend(type='demean')
        
        break_test=st27
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        sr = st27[0].stats.sampling_rate
        
        
        if sum(abs(st27[0].data)) > 10 and st27[0].stats.npts > 7200*sr:
            station_activity = np.lib.pad(station_activity, ((0,1),(0,0)), 'constant', constant_values=(0))
            station_activity[num][0] = ID
            station_activity[num][1] = st27[0].stats.starttime
            station_activity[num][2] = st27[0].stats.endtime
            num+=1
            
        
        else:
            st27 = st_blank
            
    except:
        st27 = st_blank   
        
    
    
    try:  
        ID = 28
        sta = 'FG15' # STATION VF01
        cha = 'BDF' # CHANNEL - inf
        net = 'GI'  # Fuego volcano
        loc = '06'    # location, 

        st28 = Stream()
        st28 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        st28.detrend(type='linear')
        st28.detrend(type='demean')
        
        break_test=st28
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        sr = st28[0].stats.sampling_rate
        
        
        if sum(abs(st28[0].data)) > 10 and st28[0].stats.npts > 7200*sr:
            station_activity = np.lib.pad(station_activity, ((0,1),(0,0)), 'constant', constant_values=(0))
            station_activity[num][0] = ID
            station_activity[num][1] = st28[0].stats.starttime
            station_activity[num][2] = st28[0].stats.endtime
            num+=1
            
        
        else:
            st28 = st_blank
            
    except:
        st28 = st_blank   




#%% FV01    
    try:  
        ID = 29
        sta = 'FV01' # STATION VF01
        cha = 'HDF' # CHANNEL - inf
        net = 'GI'  # Fuego volcano
        loc = ''    # location, 
            
        st29 = Stream()
        st29 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        st29.detrend(type='linear')
        st29.detrend(type='demean')
        
        break_test=st29
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        sr = st29[0].stats.sampling_rate
        
        if sum(abs(st29[0].data)) > 10 and st29[0].stats.npts > 7200*sr:
            station_activity = np.lib.pad(station_activity, ((0,1),(0,0)), 'constant', constant_values=(0))
    
            station_activity[num][0] = ID
            station_activity[num][1] = st29[0].stats.starttime
            station_activity[num][2] = st29[0].stats.endtime
            num+=1
        
        else:
            st29 = st_blank
            
    except:
        st29 = st_blank
            

#%% FV02   
    try:  
        ID = 30
        sta = 'FV02' # STATION VF01
        cha = 'HDF' # CHANNEL - inf
        net = 'GI'  # Fuego volcano
        loc = ''    # location, 
            
        st30 = Stream()
        st30 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        st30.detrend(type='linear')
        st30.detrend(type='demean')
        
        break_test=st30
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        sr = st30[0].stats.sampling_rate
        
        if sum(abs(st30[0].data)) > 10 and st30[0].stats.npts > 7200*sr:
            station_activity = np.lib.pad(station_activity, ((0,1),(0,0)), 'constant', constant_values=(0))
    
            station_activity[num][0] = ID
            station_activity[num][1] = st30[0].stats.starttime
            station_activity[num][2] = st30[0].stats.endtime
            num+=1
        
        else:
            st30 = st_blank
            
    except:
        st30 = st_blank
            
            
#%% FV03   
    try:  
        ID = 31
        sta = 'FV03' # STATION VF01
        cha = 'HDF' # CHANNEL - inf
        net = 'GI'  # Fuego volcano
        loc = ''    # location, 
            
        st31 = Stream()
        st31 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        st31.detrend(type='linear')
        st31.detrend(type='demean')
        
        break_test=st31
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        sr = st31[0].stats.sampling_rate
        
        if sum(abs(st31[0].data)) > 10 and st31[0].stats.npts > 7200*sr:
            station_activity = np.lib.pad(station_activity, ((0,1),(0,0)), 'constant', constant_values=(0))
    
            station_activity[num][0] = ID
            station_activity[num][1] = st31[0].stats.starttime
            station_activity[num][2] = st31[0].stats.endtime
            num+=1
        
        else:
            st31 = st_blank
            
    except:
        st31 = st_blank


#%% FV04   
    try:  
        ID = 32
        sta = 'FV04' # STATION VF01
        cha = 'HDF' # CHANNEL - inf
        net = 'GI'  # Fuego volcano
        loc = ''    # location, 
            
        st32 = Stream()
        st32 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        st32.detrend(type='linear')
        st32.detrend(type='demean')
        
        break_test=st32
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        sr = st32[0].stats.sampling_rate
        
        if sum(abs(st32[0].data)) > 10 and st32[0].stats.npts > 7200*sr:
            station_activity = np.lib.pad(station_activity, ((0,1),(0,0)), 'constant', constant_values=(0))
    
            station_activity[num][0] = ID
            station_activity[num][1] = st32[0].stats.starttime
            station_activity[num][2] = st32[0].stats.endtime
            num+=1
        
        else:
            st32 = st_blank
        
    except:
        st32 = st_blank
        

#%% FV08   
    try:  
        ID = 33
        sta = 'FV08' # STATION VF01
        cha = 'HDF' # CHANNEL - inf
        net = 'GI'  # Fuego volcano
        loc = ''    # location, 
        
        st33 = Stream()
        st33 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        st33.detrend(type='linear')
        st33.detrend(type='demean')
        
        break_test=st33
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        sr = st33[0].stats.sampling_rate
        
        if sum(abs(st33[0].data)) > 10 and st33[0].stats.npts > 7200*sr:
            station_activity = np.lib.pad(station_activity, ((0,1),(0,0)), 'constant', constant_values=(0))
    
            station_activity[num][0] = ID
            station_activity[num][1] = st33[0].stats.starttime
            station_activity[num][2] = st33[0].stats.endtime
            num+=1
        
        else:
            st33 = st_blank
            
    except:
        st33 = st_blank





        
  #%%      give all traces and list of active stations
        
        
    return(st1,st2,st3,st4,st5,st6,st7,st8,st9,st10,st11,st12,st13,st14,st15,st16,st17,st18,st19,st20,st21,st22,st23,st24,st25,st26,st27,st28,st29,st30,st31,st32,st33,station_activity)    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        