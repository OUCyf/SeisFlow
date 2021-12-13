#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 20:07:31 2021

Input
    1. Generate tmp file: find $(pwd) -name \*mseed > ../tmp.mseed

    2. Run this script

    3. KILL: ps -aux|grep python|grep -v grep|gawk '{print $2}' |xargs kill -9

Output
    1.DB_output format: [ network_name/year/network_name.station_name/yearday/sac ]
    
    
@author: yf
"""

import os
import sys
import obspy
import numpy as np
import time as time_pac
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
TMP_PATH = './tmp_2011_2015.mseed'

# mseed_ID
BEGIN_mseed_ID = '/datapool/yinf/4.Data_yf/Binchuan.40T.data/2011.2016/data/RATE0100/G1/2013/G1.53266/2013106/G1.53266.01.SHN.D.2013.106.00.00.00.mseed'

# ouput
OUT_PATH = './SAC_DATA1'

# Resp file path
RESP_PATH = "./resp"

# some station with 200hz rate
STA_200hz = [ "CKT0"]

# remove_response
f1 = 0.001
f2 = 0.003
f3 = 45
f4 = 50

# resample
freqmin = 0.003
freqmax = 40                                                                    # hz
sampling_rate = 100                                                               # targeted sampling rate (hz)

# MPI_n
MPI_n = 30



#%% 2. MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# 2.1 read and split
if rank == 0:
    # Read TMP_FILE file
    print("rank: ",str(rank),flush=True)
    print(f"Now begin to calculate Green's function with {MPI_n} cores...",flush=True)
    TMP_FILE = []
    with open(TMP_PATH,'r') as f:
        for line in f:
            TMP_FILE.append(line.strip("\n"))

    BEGIN_ID = TMP_FILE.index(BEGIN_mseed_ID)
    END_ID   = len(TMP_FILE)
    mseedfile_num = END_ID-BEGIN_ID
    
    CPU_splits = CPU_split(mseedfile_num, MPI_n)

else:
    TMP_FILE,CPU_splits = [None for _ in range(2)]

# 2.2 broadcast the variables
TMP_FILE = comm.bcast(TMP_FILE,root=0)
CPU_splits = comm.bcast(CPU_splits,root=0)




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
        
        # NOTE: if reference station in "RATE0100", then skip
        if (TMP_FILE[i].split('/')[-6] == "RATE0100") & (tr[0].stats.station in STA_200hz):
            continue
        
        
        # 3.2 read resp
        t1=time_pac.time();print("t1:"+str(round(t1-t0,4)))
        NET  =tr[0].stats.network
        STA = tr[0].stats.station
        LOC = tr[0].stats.location
        CHN = tr[0].stats.channel
        RESP_FILE_NAME = "RESP" + '.' + NET + '.' + STA + '.'+LOC +'.'+CHN
        RESP_FILE_PATH = os.path.join(RESP_PATH, RESP_FILE_NAME) 
        inv = obspy.read_inventory(RESP_FILE_PATH)
        
        # 3.3 merge the data and fill the gap with zeros
        t2=time_pac.time();print("t2:"+str(round(t2-t1,4)))
        tr.merge(method=1, fill_value=0)
        
        # 3.4 remean taper
        t3=time_pac.time();print("t3:"+str(round(t3-t2,4)))
        tr.detrend(type='demean')
        tr.detrend(type='simple')
        tr.taper(max_percentage=0.05)
        
        #3.5 remove response
        t4=time_pac.time();print("t4:"+str(round(t4-t3,4)))
        tr.remove_response(inventory=inv, pre_filt=[f1,f2,f3,f4], output="VEL", plot=False)
        
        #3.6 滤波降采样
        t5=time_pac.time();print("t5:"+str(round(t5-t4,4)))
        tr.filter('bandpass', freqmin=freqmin, freqmax=freqmax, corners=4, zerophase=True)  # bandpass
        tr.resample(sampling_rate,window='hanning',no_filter=True, strict_length=True)         # resample

        # 3.7 write sac
        t6=time_pac.time();print("t6:"+str(round(t6-t5,4)))
        time = obspy.UTCDateTime(tr[0].stats.starttime)
        year = "{:04d}".format(time.year)
        julday = "{:03d}".format(time.julday)
        hour = "{:02d}".format(time.hour) 
        minute = "{:02d}".format(time.minute) 
        second = "{:02d}".format(time.second)
        
        sac_path = os.path.join(OUT_PATH,
                                NET,
                                year,
                                NET+'.'+STA,
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
        t7=time_pac.time();print("t7:"+str(round(t7-t6,4)))
        # print(str(round(t1-t0),4),
        #       str(round(t2-t1),4),
        #       str(round(t3-t2),4),
        #       str(round(t4-t3),4),
        #       str(round(t5-t4),4),
        #       str(round(t6-t5),4),
        #       str(round(t7-t6),4),flush=True)
        
        
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
    print("Successful !\nAll of the mseed files have been computed.\n\n",flush=True)
    sys.exit()


