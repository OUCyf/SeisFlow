#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 30 23:17:54 2021

@author: yf
"""
import obspy
import numpy as np
import matplotlib.pyplot as plt


def fft(signal,sampling_rate):
    
    npts=len(signal) # 采样点长度
    freqs = np.linspace(0, sampling_rate, npts)
    fft_y=np.fft.fft(signal) 
    abs_y=np.abs(fft_y)                # 取复数的绝对值，即复数的模(双边频谱)
    angle_y=np.angle(fft_y)            #取复数的角度
    normalization_y=abs_y/max(abs_y) 

    return freqs, normalization_y




# path = "/Users/yf/3.Project/2.Velocity_Change/2013106/*SAC"
path = "/Users/yf/share2/cluster/datapool2/yinf/2.Project/1.yunnan_vel_change/1.DATA/SAC_DATA/G1/2013/G1.53266/2013106/*SAC"
path_mseed = "/Users/yf/share2/datapool/yinf/4.Data_yf/Binchuan.40T.data/2011.2016/data/RATE0100/G1/2013/G1.53266/2013106/*.mseed"

st = obspy.read(path)
st.filter('bandpass', freqmin=0.5, freqmax=10, corners=4, zerophase=True)
st_mseed = obspy.read(path_mseed)

##
s_sac = st[2].data
rate_sac = 20 #采样频率为Hz

s_mseed = st_mseed[2].data
rate_mseed = 100 #采样频率为Hz

fre_sac,y_sac = fft(s_sac,rate_sac)
fre_mseed,y_mseed = fft(s_mseed,rate_mseed)





# plot sac
plt.style.use('ggplot')
fig = plt.figure(figsize=(10, 5))
ax = fig.add_axes([0, 0, 0.5, 0.5])
ax.plot(fre_sac,y_sac)
ax.set_ylim(0, 1)
# ax.set_xlim(0, sampling_rate/2)
ax.set_xlim(0.0, 10)
ax.set_xlabel("Frequency(Hz)",fontsize=20)
ax.set_ylabel("Amp",fontsize=20)
ax.set_title('Spectrum analysis',fontsize=20,color='k')


# plot mseed
ax = fig.add_axes([0.6, 0, 0.5, 0.5])
ax.plot(fre_mseed,y_mseed)
ax.set_ylim(0, 0.02)
# ax.set_xlim(0, sampling_rate/2)
ax.set_xlim(0.000, 10)
ax.set_xlabel("Frequency(Hz)",fontsize=20)
ax.set_ylabel("Amp",fontsize=20)
ax.set_title('Spectrum analysis',fontsize=20,color='k')
