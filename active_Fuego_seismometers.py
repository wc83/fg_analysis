#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 15:47:20 2019

@author: root
"""
#kill code

import numpy as np
import matplotlib.pyplot as plt

num_active=np.zeros(shape=(0,44))
first = 0
last = 820
for x in range(first,last): 
    print('day',x+1,'of', last)
    num,nums,numa,FG3s,FG8s,FG10s,FG11s,FG12s,FG13s,FG14s,FG16s,FG8_1,FG8_2,FG8_3,FG10_1,FG10_2,FG10_3,FG11_1,FG11_2,FG11_4,FG11_5,FG12_1,FG12_2,FG12_3,FG13_1,FG13_2,FG13_3,FG15_1,FG15_2,FG15_3,FG15_4,FG15_5,FG15_6,FV01a,FV02a,FV03a,FV04a,VF01a,VF02a,VF03a,VF04a,VF05a,VF06a =get_activity_fuego(x)
    
    num_active=np.lib.pad(num_active, ((0,1),(0,0)), 'constant', constant_values=(0))
    
    num_active[x][0]=x+1
    num_active[x][1]=num
    num_active[x][2]=nums
    num_active[x][3]=FG3s
    num_active[x][4]=FG8s
    num_active[x][5]=FG10s
    num_active[x][6]=FG11s
    num_active[x][7]=FG12s
    num_active[x][8]=FG13s
    num_active[x][9]=FG14s
    num_active[x][10]=FG16s
    
    num_active[x][11]=numa
    
    num_active[x][12]=FG8_1
    num_active[x][13]=FG8_2
    num_active[x][14]=FG8_3
    
    num_active[x][15]=FG10_1
    num_active[x][16]=FG10_2
    num_active[x][17]=FG10_3
    
    num_active[x][18]=FG11_1
    num_active[x][19]=FG11_2
    num_active[x][20]=FG11_4
    num_active[x][21]=FG11_5
    
    num_active[x][22]=FG12_1
    num_active[x][23]=FG12_2
    num_active[x][24]=FG12_3
    
    num_active[x][25]=FG13_1
    num_active[x][26]=FG13_2
    num_active[x][27]=FG13_3
    
    num_active[x][28]=FG15_1
    num_active[x][29]=FG15_2
    num_active[x][30]=FG15_3
    num_active[x][31]=FG15_4
    num_active[x][32]=FG15_5
    num_active[x][33]=FG15_6
    
    num_active[x][34]=FV01a
    num_active[x][35]=FV02a
    num_active[x][36]=FV03a
    num_active[x][37]=FV04a
    num_active[x][38]=VF01a
    num_active[x][39]=VF02a
    num_active[x][40]=VF03a
    num_active[x][41]=VF04a
    num_active[x][42]=VF05a
    num_active[x][43]=VF06a
    

#%% plot data
    
#plt.figure()
#plt.plot(num_active[:,0],num_active[:,1])
#plt.xlabel('Day')
#plt.ylabel('Number of active stations')
#plt.title('All Station Activity')
#plt.ylim([0,10])


#plt.figure()
#plt.plot(num_active[:,0],num_active[:,4])
#plt.xlabel('Day')
##plt.ylabel('Number of active stations')
#plt.title('FG3 Station Activity')
#plt.ylim([0,2])
#
#plt.figure()
#plt.plot(num_active[:,0],num_active[:,5])
#plt.xlabel('Day')
##plt.ylabel('Number of active stations')
#plt.title('FG8 Station Activity')
#plt.ylim([0,2])
#
#plt.figure()
#plt.plot(num_active[:,0],num_active[:,6])
#plt.xlabel('Day')
##plt.ylabel('Number of active stations')
#plt.title('FG10 Station Activity')
#plt.ylim([0,2])
#
#plt.figure()
#plt.plot(num_active[:,0],num_active[:,7])
#plt.xlabel('Day')
##plt.ylabel('Number of active stations')
#plt.title('FG11 Station Activity')
#plt.ylim([0,2])
#
#plt.figure()
#plt.plot(num_active[:,0],num_active[:,8])
#plt.xlabel('Day')
##plt.ylabel('Number of active stations')
#plt.title('FG12 Station Activity')
#plt.ylim([0,2])
#
#plt.figure()
#plt.plot(num_active[:,0],num_active[:,9])
#plt.xlabel('Day')
##plt.ylabel('Number of active stations')
#plt.title('FG13 Station Activity')
#plt.ylim([0,2])
#
#plt.figure()
#plt.plot(num_active[:,0],num_active[:,10])
#plt.xlabel('Day')
##plt.ylabel('Number of active stations')
#plt.title('FG14 Station Activity')
#plt.ylim([0,2])
#
#plt.figure()
#plt.plot(num_active[:,0],num_active[:,11])
#plt.xlabel('Day')
##plt.ylabel('Number of active stations')
#plt.title('FG16 Station Activity')
#plt.ylim([0,2])


#%%

plt.figure()
plt.plot(num_active[:,0],num_active[:,1])
plt.xlabel('Day')
plt.ylabel('Number of active stations')
plt.title('Station Activity')
plt.ylim([0,30])

plt.figure()
plt.plot(num_active[:,0],num_active[:,11])
plt.xlabel('Day')
plt.ylabel('Number of active Acoustic stations')
plt.title('Acoustic Station Activity')
plt.ylim([0,20])


plt.figure()
plt.plot(num_active[:,0],num_active[:,2])
plt.xlabel('Day')
plt.ylabel('Number of active seismic stations')
plt.title('Seismic Station Activity')
plt.ylim([0,10])

#%% save data
     
#np.savetxt("/Users/william/Documents/scanner/all_stations/active_stations_fuego_all_instruments_v2.csv", num_active,delimiter=",",header="day,num,nums,numa,FG3s,FG8s,FG10s,FG11s,FG12s,FG13s,FG14s,FG16s,FG8a,FG10a,FG11a,FG12a,FG13a,FG15a,FV01a,FV02a,FV03a,FV04a,VF01a,VF02a,VF03a,VF04a,VF05a,VF06a")




