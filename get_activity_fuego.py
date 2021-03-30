#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 15:16:55 2019

@author: root
"""

def get_activity_fuego(day):

    from obspy.clients.earthworm import Client
    from obspy import UTCDateTime
    from obspy import Stream
    import numpy as np
    
    year1=2018
    month1=3
    day1=10
    hour1=0
    minute1=0
    second1=0
    
    num=0
    nums=0
    numa=0
    
    t0 = UTCDateTime(year1, month1, day1, hour1, minute1, second1) #the format is year:day_of_the_year:month
    t1 = t0 + day*24*60*60
    t2 = t1 + 23*60*60 + 59*60 +59.999 #UTCDateTime(year2, month2, day2, hour2, minute2, second2) # notice we have here 10 minutes, but we can select our times. 
        
    print(UTCDateTime(t1))
    
    
    port = '138.253.113.19' # ip, port - ip's 138.253.113.19 or 138.253.112.23
#%% FG10 seismic    
    try:  
        
        sta = 'FG10' # STATION LB01
        cha = 'BHZ' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        
        if np.mean(abs(st1[0].data)) > 2 :
            
            FG10s=1
            num += 1
            nums += 1
        else:
            FG10s=0
                  
    except: 
        
        FG10s=0
        
#%% FG11 Seismic
    try:  
        
        sta = 'FG11' # STATION LB01
        cha = 'BHZ' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
       
        if np.mean(abs(st1[0].data)) > 2 :
            
            FG11s=1
            num += 1
            nums += 1
        else:
            FG11s=0

            
    except:  

        FG11s=0

#%% FG12 Seismic
    try:  
        
        sta = 'FG12' # STATION LB01
        cha = 'BHZ' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = '00'    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        if np.mean(abs(st1[0].data)) > 2 :
            
            FG12s=1
            num += 1
            nums += 1
        else:
            FG12s=0
            
    except: 
        
        FG12s=0
        
#%% FG13 Seismic
    try:  
        
        sta = 'FG13' # STATION LB01
        cha = 'BHZ' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = '00'    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        if np.mean(abs(st1[0].data)) > 2 :
            
            FG13s=1
            num += 1
            nums += 1
        else:
            FG13s=0
            
    except:
        

        FG13s=0
        
#%% FG14 Seismic
    try:  
        
        sta = 'FG14' # STATION LB01
        cha = 'BHZ' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = '01'    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        if np.mean(abs(st1[0].data)) > 2 :
            
            FG14s=1
            num += 1
            nums += 1
        else:
            FG14s=0
            
    except: 
        

        FG14s=0
        

#%% FG16 Seismic
    try:  
        
        sta = 'FG16' # STATION LB01
        cha = 'BHZ' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = '00'    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        if np.mean(abs(st1[0].data)) > 2 :
            
            FG16s=1
            num += 1
            nums += 1
        else:
            FG16s=0
          
    except: 
        

        FG16s=0


#%% FG8 Seismic
    try:  
        
        sta = 'FG8' # STATION LB01
        cha = 'BHZ' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = '00'    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        if np.mean(abs(st1[0].data)) > 2 :
            
            FG8s=1
            num += 1
            nums += 1
        else:
            FG8s=0
                   
    except:     

        FG8s=0
        
#%% FG3 Seismic
    try:  
        
        sta = 'FG3' # STATION LB01
        cha = 'SHZ' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = '01'    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        if np.mean(abs(st1[0].data)) > 2 :
            
            FG3s=1
            num += 1
            nums += 1
        else:
            FG3s=0

    except: 

        FG3s=0



#%%
        




#%% FG10 inf  
    try:  
        
        sta = 'FG10' # STATION LB01
        cha = 'BDF' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = '01'    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        if np.mean(abs(st1[0].data)) > 2 :
            
            FG10_1=1
            num += 1
            numa += 1
        else:
            FG10_1=0
                  
    except: 
        
        FG10_1=0
        
    try:  
        
        sta = 'FG10' # STATION LB01
        cha = 'BDF' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = '02'    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        if np.mean(abs(st1[0].data)) > 2 :
            
            FG10_2=1
            num += 1
            numa += 1
        else:
            FG10_2=0
                  
    except: 
        
        FG10_2=0
        
        
    try:  
        
        sta = 'FG10' # STATION LB01
        cha = 'BDF' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = '03'    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        if np.mean(abs(st1[0].data)) > 2 :
            
            FG10_3=1
            num += 1
            numa += 1
        else:
            FG10_3=0
                  
    except: 
        
        FG10_3=0
        
#%% FG11 inf
    try:  
        
        sta = 'FG11' # STATION LB01
        cha = 'BDF' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = '01'    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
       
        if np.mean(abs(st1[0].data)) > 2 :
            
            FG11_1=1
            num += 1
            numa += 1
        else:
            FG11_1=0
        
            
    except:  

        FG11_1=0
        
    try:  
        
        sta = 'FG11' # STATION LB01
        cha = 'BDF' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = '02'    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
       
        if np.mean(abs(st1[0].data)) > 2 :
            
            FG11_2=1
            num += 1
            numa += 1
        else:
            FG11_2=0
            
    except:  

        FG11_2=0
        
    try:  
        
        sta = 'FG11' # STATION LB01
        cha = 'BDF' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = '04'    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
       
        if np.mean(abs(st1[0].data)) > 2 :
            
            FG11_4=1
            num += 1
            numa += 1
        else:
            FG11_4=0
            
    except:  

        FG11_4=0
    
    try:  
        
        sta = 'FG11' # STATION LB01
        cha = 'BDF' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = '05'    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
       
        if np.mean(abs(st1[0].data)) > 2 :
            
            FG11_5=1
            num += 1
            numa += 1
        else:
            FG11_5=0
            
    except:  

        FG11_5=0

#%% FG12 inf
    try:  
        
        sta = 'FG12' # STATION LB01
        cha = 'BDF' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = '01'    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        if np.mean(abs(st1[0].data)) > 2 :
            
            FG12_1=1
            num += 1
            numa += 1
        else:
            FG12_1=0
            
    except: 
        
        FG12_1=0
    
    try:  
        
        sta = 'FG12' # STATION LB01
        cha = 'BDF' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = '02'    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        if np.mean(abs(st1[0].data)) > 2 :
            
            FG12_2=1
            num += 1
            numa += 1
        else:
            FG12_2=0
            
    except: 
        
        FG12_2=0
        
    try:  
        
        sta = 'FG12' # STATION LB01
        cha = 'BDF' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = '03'    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        if np.mean(abs(st1[0].data)) > 2 :
            
            FG12_3=1
            num += 1
            numa += 1
        else:
            FG12_3=0
            
    except: 
        
        FG12_3=0
        
#%% FG13 inf
    try:  
        
        sta = 'FG13' # STATION LB01
        cha = 'BDF' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = '01'    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        if np.mean(abs(st1[0].data)) > 2 :
            
            FG13_1=1
            num += 1
            numa += 1
        else:
            FG13_1=0
            
    except:
        
        FG13_1=0
        
    try:  
        
        sta = 'FG13' # STATION LB01
        cha = 'BDF' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = '02'    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        if np.mean(abs(st1[0].data)) > 2 :
            
            FG13_2=1
            num += 1
            numa += 1
        else:
            FG13_2=0
            
    except:
        
        FG13_2=0
        
    try:  
        
        sta = 'FG13' # STATION LB01
        cha = 'BDF' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = '03'    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        if np.mean(abs(st1[0].data)) > 2 :
            
            FG13_3=1
            num += 1
            numa += 1
        else:
            FG13_3=0
            
    except:
        
        FG13_3=0
        

#%% FG15 inf
    try:  
        
        sta = 'FG15' # STATION LB01
        cha = 'BDF' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = '01'    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        if np.mean(abs(st1[0].data)) > 2 :
            
            FG15_1=1
            num += 1
            numa += 1
        else:
            FG15_1=0
          
    except: 
        
        FG15_1=0
        
    try:  
        
        sta = 'FG15' # STATION LB01
        cha = 'BDF' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = '02'    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        if np.mean(abs(st1[0].data)) > 2 :
            
            FG15_2=1
            num += 1
            numa += 1
        else:
            FG15_2=0
          
    except: 
        
        FG15_2=0
        
    try:  
        
        sta = 'FG15' # STATION LB01
        cha = 'BDF' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = '03'    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        if np.mean(abs(st1[0].data)) > 2 :
            
            FG15_3=1
            num += 1
            numa += 1
        else:
            FG15_3=0
          
    except: 
        
        FG15_3=0
        
    try:  
        
        sta = 'FG15' # STATION LB01
        cha = 'BDF' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = '04'    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        if np.mean(abs(st1[0].data)) > 2 :
            
            FG15_4=1
            num += 1
            numa += 1
        else:
            FG15_4=0
          
    except: 
        
        FG15_4=0

    try:  
        
        sta = 'FG15' # STATION LB01
        cha = 'BDF' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = '05'    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        if np.mean(abs(st1[0].data)) > 2 :
            
            FG15_5=1
            num += 1
            numa += 1
        else:
            FG15_5=0
          
    except: 
        
        FG15_5=0
        
    try:  
        
        sta = 'FG15' # STATION LB01
        cha = 'BDF' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = '06'    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        if np.mean(abs(st1[0].data)) > 2 :
            
            FG15_6=1
            num += 1
            numa += 1
        else:
            FG15_6=0
          
    except: 
        
        FG15_6=0


#%% FG8 inf
    try:  
        
        sta = 'FG8' # STATION LB01
        cha = 'BDF' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = '01'    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        if np.mean(abs(st1[0].data)) > 2 :
            
            FG8_1=1
            num += 1
            numa += 1
        else:
            FG8_1=0
          
    except: 
        
        FG8_1=0

    try:  
        
        sta = 'FG8' # STATION LB01
        cha = 'BDF' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = '02'    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        if np.mean(abs(st1[0].data)) > 2 :
            
            FG8_2=1
            num += 1
            numa += 1
        else:
            FG8_2=0
          
    except: 
        
        FG8_2=0
        
    try:  
        
        sta = 'FG8' # STATION LB01
        cha = 'BDF' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = '03'    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        if np.mean(abs(st1[0].data)) > 2 :
            
            FG8_3=1
            num += 1
            numa += 1
        else:
            FG8_3=0
          
    except: 
        
        FG8_3=0

#%% FV01 inf
    try:  
        
        sta = 'FV01' # STATION LB01
        cha = 'HDF' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        if np.mean(abs(st1[0].data)) > 2 :
            
            FV01a=1
            num += 1
            numa += 1
        else:
            FV01a=0
          
    except: 
        
        FV01a=0

#%% FV02 inf
    try:  
        
        sta = 'FV02' # STATION LB01
        cha = 'HDF' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        if np.mean(abs(st1[0].data)) > 2 :
            
            FV02a=1
            num += 1
            numa += 1
        else:
            FV02a=0
          
    except: 
        
        FV02a=0
        
#%% FV03 inf
    try:  
        
        sta = 'FV03' # STATION LB01
        cha = 'HDF' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        if np.mean(abs(st1[0].data)) > 2 :
            
            FV03a=1
            num += 1
            numa += 1
        else:
            FV03a=0
          
    except: 
        
        FV03a=0

#%% FV04 inf
    try:  
        
        sta = 'FV04' # STATION LB01
        cha = 'HDF' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        if np.mean(abs(st1[0].data)) > 2 :
            
            FV04a=1
            num += 1
            numa += 1
        else:
            FV04a=0
          
    except: 
        
        FV04a=0



#%% VF01 inf
    try:  
        
        sta = 'VF01' # STATION LB01
        cha = 'HDF' # CHANNEL - Vertical
        net = 'XZ'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        if np.mean(abs(st1[0].data)) > 2 :
            
            VF01a=1
            num += 1
            numa += 1
        else:
            VF01a=0
          
    except: 
        
        VF01a=0


#%% VF02 inf
    try:  
        
        sta = 'VF02' # STATION LB01
        cha = 'HDF' # CHANNEL - Vertical
        net = 'XZ'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        if np.mean(abs(st1[0].data)) > 2 :
            
            VF02a=1
            num += 1
            numa += 1
        else:
            VF02a=0
          
    except: 
        
        VF02a=0
        
        
#%% VF03 inf
    try:  
        
        sta = 'VF03' # STATION LB01
        cha = 'HDF' # CHANNEL - Vertical
        net = 'XZ'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        if np.mean(abs(st1[0].data)) > 2 :
            
            VF03a=1
            num += 1
            numa += 1
        else:
            VF03a=0
          
    except: 
        
        VF03a=0

#%% VF04 inf
    try:  
        
        sta = 'VF04' # STATION LB01
        cha = 'HDF' # CHANNEL - Vertical
        net = 'XZ'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        if np.mean(abs(st1[0].data)) > 2 :
            
            VF04a=1
            num += 1
            numa += 1
        else:
            VF04a=0
          
    except: 
        
        VF04a=0

#%% VF05 inf
    try:  
        
        sta = 'VF05' # STATION LB01
        cha = 'HDF' # CHANNEL - Vertical
        net = 'XZ'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        if np.mean(abs(st1[0].data)) > 2 :
            
            VF05a=1
            num += 1
            numa += 1
        else:
            VF05a=0
          
    except: 
        
        VF05a=0

#%% VF06 inf
    try:  
        
        sta = 'VF06' # STATION LB01
        cha = 'HDF' # CHANNEL - Vertical
        net = 'XZ'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        if np.mean(abs(st1[0].data)) > 2 :
            
            VF06a=1
            num += 1
            numa += 1
        else:
            VF06a=0
          
    except: 
        
        VF06a=0



    #%% return all stations
    
    return(num,nums,numa,FG3s,FG8s,FG10s,FG11s,FG12s,FG13s,FG14s,FG16s,FG8_1,FG8_2,FG8_3,FG10_1,FG10_2,FG10_3,FG11_1,FG11_2,FG11_4,FG11_5,FG12_1,FG12_2,FG12_3,FG13_1,FG13_2,FG13_3,FG15_1,FG15_2,FG15_3,FG15_4,FG15_5,FG15_6,FV01a,FV02a,FV03a,FV04a,VF01a,VF02a,VF03a,VF04a,VF05a,VF06a)
    
    
    
    
    
    
    
    
    
    
    
    