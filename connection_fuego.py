#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 13 11:19:55 2018

@author: william
"""

import os
from collections import OrderedDict
import numpy as np
import obspy
import obspy.signal
import scipy.signal as sgn
import matplotlib.pyplot as plt 
import matplotlib.mlab as mlab
from scipy import integrate
from obspy.clients.earthworm import Client
from obspy import UTCDateTime
from obspy.signal.trigger import plot_trigger, trigger_onset
from obspy.signal.trigger import classic_sta_lta, recursive_sta_lta
from obspy import Stream
from numpy import argmax

# inf calibs
c1 = 1.22E-04
c2 = 6.10E-05
c3 = 2.78E-04
c4 = 6.25E-06
c5 = 2.50E-06

#for seismics
c6 = 3.30E-09
c7 = 1.33E-09


sta = 'FG12' # STATION 
cha = 'BHZ' # CHANNEL
net = 'GI'  # 
loc = '00'    # location, it depends mostly of which network you are in. 
client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23

t1 = UTCDateTime(2019,1,24,21,36,0) #the format is year:day_of_the_year:month
#t1 =UTCDateTime(1560601280)  
t2 = t1 +  5*60 

#t3 = t2 + 3*24*60*60

st = Stream()
st = client.get_waveforms(net, sta, loc, cha, t1 , t2 + 2)
#st[0].plot(color='r',starttime=t1, endtime=t2)
tr = st[0] 

tr.data = tr.data * c6
#
#st2 = Stream()
#st2 = client.get_waveforms(net, sta, loc, cha, t2 - 2 , t3)

#tr2= st2[0]

#%%
#tr = tr1 + tr2

sr = st[0].stats.sampling_rate
tr.detrend(type='linear')
tr.detrend(type='demean')
tr.filter(type='bandpass',freqmin=0.1, freqmax=6) #((sr/2)-(sr/20)))
#tr.plot(type='dayplot',starttime=t1, endtime=t2)
tr.plot(color='r',starttime=t1, endtime=t2)

#tr.spectrogram(log=False, dbscale=True)

#ti = np.linspace(0,len(tr.data)/sr,len(tr.data))
#f, t, Sxx = sgn.spectrogram(tr.data,fs=50)
#plt.figure(figsize=(12,3))
#plt.pcolormesh(t, f, Sxx, shading='gouraud')
#plt.ylabel('Frequency [Hz]')
#plt.xlabel('Time [sec]')
##plt.ylim([-10,10])
#plt.plot(ti,tr.data)
#plt.show()

### STATION, CHANNEL (DDF --> 400 Hz), NETWWORK AND LOCATION CODES 
#sta = 'FG12' # STATION 
#cha = 'BHZ' # CHANNEL
#net = 'GI'  # 
#loc = '00'    # location, it depends mostly of which network you are in. 
#
#
#sts = Stream()
#sts = client.get_waveforms(net, sta, loc, cha, t1 , t2)
##print(st)
#
##2014-11-25 20:20:16 (UTC)
#trs = sts[0]
#trs.detrend(type='linear')
#trs.detrend(type='demean')
##trs.filter(type='bandpass',freqmin=10, freqmax=20)
#trs.plot(color='k',starttime=t1, endtime=t2)

#tr.spectrogram(log=False, dbscale=True)


#%%  fft



start = t1
end=t2
sf = tr.stats.sampling_rate
#
dom,cf, bwid50 = freq_info_NqV(tr,start,end,sf)

mean_amp = np.mean(abs(tr.data))
med_amp = np.median(abs(tr.data))

print('')
print('central f =', cf, 'Hz')
print('dominant f =', dom, 'Hz')
print('bandwith 50% =', bwid50, 'Hz')

print('mean amp =', mean_amp)
print('median amp =', med_amp)

#
#
#
#
#%% get spectrum 

#frequency info
sf=st[0].stats.sampling_rate
nq=sf/2

end = t2
start = t1
 
window=end-start
tr_data=tr.data
m=np.mean(tr_data)
tr_data = tr_data-m
famp = abs(np.fft.fft(tr_data))

 # bandwidth 50%
values = famp
searchval = max(famp)/2
ii = np.where(values >= searchval)[0]

above=[]
for p in range(0,len(ii)):
    if ii[p]/window < nq:
        above.append(ii[p]/window)        
bwid50=max(above)-min(above)  

#%%

#

Fs = tr.stats.sampling_rate;  # sampling rate
Ts = 1/Fs; # sampling interval
t = np.arange(0,(t2-t1),Ts) # time vector

y = tr.data

n = len(y) # length of the signal
k = np.arange(n)
T = n/Fs
frq = k/T # two sides frequency range
#frq = frq[range(n/2)] # one side frequency range

Y = famp/(window*Fs) # fft computing and normalization
#Y = Y[range(n/2)]


#
#x_cf = [cf,cf]
#yf = [0,100]
#

maxf=max(abs(Y))
hf = max(abs(Y))/2
cf_a = Y[int(cf/frq[1])]

ypf=[0,max(abs(Y))]
ycf=[0,cf_a]
xpeak=[dom,dom]
cff=[cf,cf]
y50=[hf,hf]
bwl=[min(above),max(above)]


#plt.figure(99)
#plt.plot(frq,abs(Y),'r')
#plt.xlabel('Freq (Hz)')
#plt.ylabel('normalised amplitude')
#plt.xlim([0, 5])

plt.figure(100)
plt.plot(frq,abs(Y),'r')
plt.plot(cff,ycf,'k')
plt.plot(xpeak,ypf,'b')
plt.plot(bwl,y50,'g')
plt.plot(min(above),hf,'go')
plt.plot(max(above),hf,'go')
plt.plot(dom,maxf,'bo')
plt.plot(cf,cf_a,'ko')
plt.xlabel('Frequency (Hz)')
plt.ylabel('|Amplitude|')
plt.title('fft')
plt.legend(('Frequency spectra','Central frequency','Dominant frequency','50% bandwidth'))
plt.xlim([0, 10])


##  range ratio

fps=int(len(famp)/(sf/0.5))
fpe=int(len(famp)/(sf/1.5))

mps=int(len(famp)/(sf/1.5))
mpe=int(len(famp)/(sf/3))

lps=int(len(famp)/(sf/3))
lpe=int(len(famp)/(sf/4.5))
#
f_peak_m = np.mean(famp[fps:fpe])/(window*sf)
m_peak_m = np.mean(famp[mps:mpe])/(window*sf)
l_peak_m = np.mean(famp[lps:lpe])/(window*sf)
peak_r = m_peak_m/((f_peak_m+l_peak_m)/2)

peak2_r = l_peak_m/m_peak_m

peaks_r = f_peak_m/l_peak_m

print("")
print("range ratio = {:.3}".format(peak_r) )
print("peaks ratio = {:.3}".format(peaks_r) )
print("2nd peak ratio = {:.3}".format(peak2_r) )
print("")





#plt.figure(101)
#plt.plot(frq,abs(Y),'r')
##plt.plot(cff,ycf,'k')
##plt.plot(xpeak,ypf,'b')
##plt.plot(bwl,y50,'g')
##plt.plot(min(above),hf,'go')
##plt.plot(max(above),hf,'go')
##plt.plot(dom,maxf,'bo')
##plt.plot(cf,cf_a,'ko')
#plt.xlabel('Frequency (Hz)')
#plt.ylabel('|Amplitude|')
#plt.title('fft')
##plt.legend(('Frequency spectra','Central frequency','Dominant frequency','50% bandwidth'))
#plt.xlim([0, 10])
#


#
##%% Infrasound fft
##
#tr = sts[0]
#start = t1
#end=t2
#sf= tr.stats.sampling_rate
##
#dom,cf, bwid50 = freq_info_NqV(tr,start,end,sf)
#
#print('')
#print('Inf central f =', cf, 'Hz')
#print('Inf dominant f =', dom, 'Hz')
#print('Inf bandwith 50% =', bwid50, 'Hz')
#
#
#
#
#
##%% get spectrum
#
#sf=st[0].stats.sampling_rate
#nq=sf/2
#
##frequency info
#  
#
#end = t2
#start = t1
# 
#window=end-start
#tr_data=tr.data
#m=np.mean(tr_data)
#tr_data = tr_data-m
#famp = abs(np.fft.fft(tr_data))
#
#
## dominant f
#peak= argmax(abs(famp))/window
#if peak > nq:
#    peak = sf-peak        
#
## center f 
#hal =sum(famp)/4 #half of the first half of the complete fft (ie. less than nq)
#num=0
#misf_i=100000000000000
#for t in range(0,len(famp),1):
#    n=famp[t]                  
#    num += n
#    misf = abs(num-hal)
#    if misf < misf_i and t != 0:
#        misf_i = misf
#        cf=t/window
#
## bandwidth 50%
#values = famp
#searchval = max(famp)/2
#ii = np.where(values >= searchval)[0]
#
#above=[]
#for p in range(0,len(ii)):
#    if ii[p]/window < nq:
#        above.append(ii[p]/window)        
#bwid50=max(above)-min(above)    
##
## double peak ratio
#fps=int(len(famp)/(sf/0.5))
#midp=int(len(famp)/(sf/2))
#spe=int(len(famp)/(sf/4.5))
#
#f_peak_m = np.mean(famp[fps:midp])/(window*sf)
#s_peak_m = np.mean(famp[midp:spe])/(window*sf)
#peak_r = f_peak_m/s_peak_m
#
##print('')
##print('first peak mean = ', f_peak_m)
##print('second peak mean = ', s_peak_m)
##print('ratio between peaks = ', peak_r)
##print('')
#
## 2Hz drop ratio
#
#fpbs=int(len(famp)/(sf/0.5))
#fpbe=int(len(famp)/(sf/1.8))
#
#spbs=int(len(famp)/(sf/2.2))
#spbe=int(len(famp)/(sf/3.5))
#
#peaks_ma = np.mean(famp[fpbs:fpbe])/(window*sf)
#peaks_mb = np.mean(famp[spbs:spbe])/(window*sf)
#peaks_m = (peaks_ma + peaks_mb)/2
#
#trough_m = np.mean(famp[fpbe:spbs])/(window*sf)
#
#trough_r = peaks_m/trough_m
#
#print('peaks mean = ', peaks_m)
#print('trough mean = ', trough_m)
#print('ratio between peaks and 2Hz dip = ', trough_r)
#print('')
#
#
##%%
#
#
#Fs = tr.stats.sampling_rate;  # sampling rate
#Ts = 1/Fs; # sampling interval
#t = np.arange(0,(t2-t1),Ts) # time vector
#
#y = tr.data
#
#n = len(y) # length of the signal
#k = np.arange(n)
#T = n/Fs
#frq = k/T # two sides frequency range
##frq = frq[range(n/2)] # one side frequency range
#
#Y = famp/(window*Fs) # fft computing and normalization
##Y = Y[range(n/2)]
#
#
##
##x_cf = [cf,cf]
##yf = [0,100]
##
#
#maxf=max(abs(Y))
#hf = max(abs(Y))/2
#cf_a = Y[int(cf/frq[1])]
#
#ypf=[0,max(abs(Y))]
#ycf=[0,cf_a]
#xpeak=[peak,peak]
#cff=[cf,cf]
#y50=[hf,hf]
#bwl=[min(above),max(above)]
#
#
##plt.figure(99)
##plt.plot(frq,abs(Y),'r')
##plt.xlabel('Freq (Hz)')
##plt.ylabel('normalised amplitude')
##plt.xlim([0, 5])
#
#plt.figure(101)
#plt.plot(frq,abs(Y),'r')
##plt.plot(cff,ycf,'k')
##plt.plot(xpeak,ypf,'b')
##plt.plot(bwl,y50,'g')
##plt.plot(min(above),hf,'go')
##plt.plot(max(above),hf,'go')
##plt.plot(peak,maxf,'bo')
##plt.plot(cf,cf_a,'ko')
#plt.xlabel('Frequency (Hz)')
#plt.ylabel('|Amplitude|')
#plt.title('Infrasound fft')
##plt.legend(('Frequency spectra','Central frequency','Dominant frequency','50% bandwidth'))
#plt.xlim([0, 25])
#



#%%
##
#from obspy.signal.trigger import classic_sta_lta, recursive_sta_lta
#from obspy.signal.trigger import plot_trigger, trigger_onset
#
#
#
#trace=st[0]
#
##window endpoints
#start= t1 + 2 #time window start 
#end= t2
##end=sample[0].stats.endtime 
#trs = trace.slice(starttime = start  , endtime= end) #cut out sample waveform with same window length as chosen event
##trs.filter("bandpass", freqmin=fmin,freqmax=fmax)
##trs_e = obspy.signal.filter.envelope(trs.data)
##print('reference waveform')
#
#
#trs.plot(type='relative',color='b')#, starttime=start , endtime=end)
#trace = st[0]
#
#sr = trace.stats.sampling_rate
#nsta=int(10*sr)                                      #2
#nlta=int(120*sr)                                     #20
#stream1=trace.data
#cft=recursive_sta_lta(stream1, nsta, nlta)
#trig_on=10                                         #8
#trig_off=0.01                                        #0.2
#plot_trigger(trace, cft, trig_on, trig_off) 
#
#on_off = trigger_onset(cft,trig_on,trig_off)
#
#
#
#for x in range(0,len(on_off)):
#    if ((on_off[x,1] - on_off[x,0])/sr) > 120:
#        
#        t_start = t1 + on_off[x,0]/sr  
#        t_end = t1 + on_off[x,1]/sr 
#        event1 = tr.slice(starttime = t_start , endtime= t_end )
#        
#        
#        tr1 = event1
#        sf = tr1.stats.sampling_rate
#        dom,cf, bwid50 = freq_info_NqV(tr1,t_start,t_end,sf)
#        med_dat = np.median(abs(event1.data))
#        
#        if dom > 1 and cf > 2:
#            
##            print('')
##            print('central f =', cf, 'Hz')
##            print('dominant f =', dom, 'Hz')
##            print('bandwith 50% =', bwid50, 'Hz')
##            print('median amp =', med_dat)
#            print(event1.stats.starttime,'to', event1.stats.endtime )
#            event1.plot(color='r')
#            
#
#
#
#
#
#
#
#

































    