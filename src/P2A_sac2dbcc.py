#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 20:07:31 2021

Input
    1. Generate tmp file: find $(pwd) -name \*.SAC > ../tmp.sac

    2. Run this script

    3. KILL: ps -aux|grep python|grep -v grep|gawk '{print $2}' |xargs kill -9

Output
    1.DB_output format: [ network_name/year/network_name.station_name/yearday/sac ]
    
    
@author: yf
"""

import os
import sys
import glob
import obspy
import numpy as np
import time as time_pac
from obspy.core import UTCDateTime
from mpi4py import MPI

os.environ["OMP_NUM_THREADS"] = "1" # export OMP_NUM_THREADS=1

def CPU_split(Source_num,MPI_n):
    '''
    The number of all earthquakes is evenly divided into N calculation cores 
    (MPI_N), and the remaining earthquake numbers are divided into the first 
    several cores
    '''

    cpu_int = Source_num // MPI_n
    cpu_rem = Source_num % MPI_n
    cpu_each_num = np.zeros(shape=(MPI_n,1) )
    
    for i in range(0,MPI_n,1):
        if i < cpu_rem:
            cpu_each_num[i]=cpu_int+1
        else:
            cpu_each_num[i]=cpu_int
    CPU_chunk = np.cumsum(cpu_each_num)
    
    CPU_splits=[]
    for i in range(0,MPI_n,1):
        if i==0:
            index_start = 1
        else:
            index_start = round( CPU_chunk[i-1]+1 )
            
        index_end = round( CPU_chunk[i] )
            
        CPU_splits.append([index_start,index_end])


    return CPU_splits





#%% 1. input parameters

# tmp.mseed fiel path
TMP_PATH = '/cluster/datapool2/yinf/2.Project/1.yunnan_vel_change/1.DATA/SAC_DATA/tmp_server.sac'

# mseed_ID
BEGIN_sac_ID = '/cluster/datapool2/yinf/2.Project/1.yunnan_vel_change/1.DATA/SAC_DATA/G1/2016/G1.53034/2016104/G1.53034.01.SHN.D.2016.104.00.00.00.SAC'

# ouput
OUT_PATH = '/cluster/datapool2/yinf/2.Project/1.yunnan_vel_change/1.DATA/SeisNoise_DATA'

# file of sta and lon
STA_LAT_LON_PATH = "/cluster/datapool2/yinf/2.Project/1.yunnan_vel_change/1.DATA/sta.all"

# filter
freqmin = 0.5
freqmax = 8                                                                    # hz
sampling_rate = 20                                                               # targeted sampling rate (hz)

# MPI_n
MPI_n = 10








#%% 2. MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# 2.1 read and split
if rank == 0:
    # Read TMP_FILE file
    print("rank: ",str(rank),flush=True)
    print(f"Now begin to calculate with {MPI_n} cores...",flush=True)
    TMP_FILE = []
    with open(TMP_PATH,'r') as f:
        for line in f:
            TMP_FILE.append(line.strip("\n"))

    BEGIN_ID = TMP_FILE.index(BEGIN_sac_ID)
    END_ID   = len(TMP_FILE)
    sacfile_num = END_ID-BEGIN_ID
    
    CPU_splits = CPU_split(sacfile_num, MPI_n)
    
    
    STA_INFO={}                                                             # store a station corresponding to the source
    with open(STA_LAT_LON_PATH,'r') as f:
        for line in f:
            line = line.strip()
            line = line.strip('\t')
            line = line.strip('\n')
            line = line.split()
            sta_name = line[0]
            info = line[1:]
            STA_INFO.update({ sta_name: info})


else:
    TMP_FILE,CPU_splits,STA_INFO = [None for _ in range(3)]

# 2.2 broadcast the variables
TMP_FILE = comm.bcast(TMP_FILE,root=0)
CPU_splits = comm.bcast(CPU_splits,root=0)
STA_INFO = comm.bcast(STA_INFO ,root=0)



#%% 3. loop
begin_ID = CPU_splits[rank][0]-1
end_ID = CPU_splits[rank][1]

tt0=time_pac.time()
for i in range(begin_ID,end_ID,1):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    
    try:
        # 3.1 read mseed
        t0=time_pac.time()
        tr = obspy.read(TMP_FILE[i])
        
        NET  =tr[0].stats.network
        STA = tr[0].stats.station
        LOC = tr[0].stats.location
        CHN = tr[0].stats.channel
        tr[0].stats.sac["stlo"] = STA_INFO[STA][0]
        tr[0].stats.sac["stla"] = STA_INFO[STA][1]
        
        # tr[0].stats.sac["stel"] = STA_INFO[STA][2]
        
        
        # 3.2 滤波降采样
        tr.filter('bandpass', freqmin=freqmin, freqmax=freqmax, corners=4, zerophase=True)  # bandpass
        # tr.resample(sampling_rate,window='hanning',no_filter=True, strict_length=True)         # resample
        
        # 3.3 write sac
        time = obspy.UTCDateTime(tr[0].stats.starttime)
        year = "{:04d}".format(time.year)
        julday = "{:03d}".format(time.julday)
        hour = "{:02d}".format(time.hour) 
        minute = "{:02d}".format(time.minute) 
        second = "{:02d}".format(time.second)
        
        UTCtime_begin = UTCDateTime(year=time.year, julday=time.julday)
        UTCtime_end = UTCtime_begin + 60*60*24
        tr.trim(UTCtime_begin, UTCtime_end, pad=True,fill_value=0)
        
        sac_path = os.path.join(OUT_PATH,
                                year,
                                year+julday)
        if not os.path.exists(sac_path):
            try:
                os.makedirs(sac_path)
            except OSError as reason:
                print("Rank = "+str(rank)+": now ID "+str(i)+ " || "+ str(end_ID)+" error~")
                print("    mkdir sac_path error with file:"+TMP_FILE[i])
                print('    reason: ' + str(reason))
                sys.exit()

        sac_name = os.path.join(sac_path, 
                                NET+'.'+
                                STA+'.'+
                                LOC+'.'+
                                CHN+'.'+
                                'D'+'.'+  # what's the D?
                                year+'.'+
                                julday+'.'+
                                hour+'.'+
                                minute+'.'+
                                second+'.'+
                                'SAC')
        tr.write(sac_name, format='SAC')
        
        
        tt1=time_pac.time()
        print("Rank = "+str(rank)+": now ID "+str(i)+ " || "+ str(end_ID) +" || "+'time takes '+str(tt1-tt0)+' s'+" ok~")
        
    except Exception:
        tt1=time_pac.time()
        print("Rank = "+str(rank)+": now ID "+str(i)+ " || "+ str(end_ID) +" || "+'time takes '+str(tt1-tt0)+' s'+" skip~")
        print("     skipping with file:"+TMP_FILE[i]);continue



#%% 4.end
comm.barrier()
if rank == 0:
    print("\n\n*****************************************************************",flush=True)
    print("Successful !\nAll of the sac files have been computed.\n\n",flush=True)
    sys.exit()


