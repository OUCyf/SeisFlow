#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 30 23:01:58 2021

Input
    1.Data format: reftek data with folder of one network 
                    [ year/yearday/unit(速采)/0 and 1 file ]
    2.component order: BHZ,BHE,BHN
    2.(optional) sample_rate, network, location

    # 3.(optional) info.list which includes 会先检查在此时间段内是否包含此台站，然后核对经纬度信息。
    #                 Good: 台站名称对应，且经纬度信息对应
    #                 Check: 台站名称和经纬度有一个对应不上 （会转换，）
    #                 Bad: 台站名称和经纬度都对应不上 （不会转换）
    #                 [ network_name/station_name/location/stla/stlo/altitude/time_lapses ]



Output
    1.DB_output format: [ network_name/year/network_name.station_name/yearday/sac ]


@author: yf
"""

import os
import sys
import glob
import json5
import obspy
from obspy import Stream


#%% 1. input parameters
# path = '/Users/yf/share1/DATA/yinfu/1.Data/bichun_db/2011'
# out_path = '/Users/yf/share1/DATA/yinfu/1.Data/bichun_db/out/'
path = '/DATA/yinfu/1.Data/bichun_db/2011'
out_path = '/DATA/yinfu/1.Data/bichun_db/out/'

component_order = ['BHZ','BHE','BHN']

start_mode = False
year_start = '2011'
yearday_start = '2011020'
unit_start = '9A8F'

end_mode = False
year_end = '2011'
yearday_end = '2011025'
unit_end = '9A8F'

network = 'G1'              # optional
location = '01'             # optional
# sample_rate = 50           # optional targeted sampling rate (hz)  desample





#%% 2. prepare
# 2.1 weather does optional var exist?
network_exists = 'network' in locals() or 'network' in globals()
location_exists = 'location' in locals() or 'location' in globals()
sample_rate_exists = 'sample_rate' in locals() or 'sample_rate' in globals()

# 2.2 generate list info
info_list = []
all_year = sorted(glob.glob(path))

for i in range(0,len(all_year),1):                                              # 循环年
    all_yearday = sorted(glob.glob(all_year[i]+'/*'))
    year_st = all_year[i].split('/')[-1]
    
    for j in range(0,len(all_yearday),1):                                       # 循环天
        all_unit = sorted(glob.glob(all_yearday[j]+'/*'))
        yearday_st = all_yearday[j].split('/')[-1]
        print("Making info list: [year] ",year_st," [yearday] ",yearday_st)
        
        for k in range(0,len(all_unit),1):
            all_reftek = sorted(glob.glob(all_unit[k]+'/1/*'))
            if all_reftek!=[]:                                                  # 如果 1 文件里存在文件
                info_list += [ all_unit[k]+'/1' ]
            else:
                continue
print("\n\n*****************************")
print("totall num of info list: ", len(info_list),"/station's day (unit)")
print("*****************************")


# 2.3 start and end path
if start_mode==True:
    path_start = os.path.join(os.path.dirname(path),year_start,yearday_start,unit_start,'1')
    start_index = info_list.index(path_start)
else:
    path_start = info_list[0]
    start_index = 0

if end_mode==True:
    path_end = os.path.join(os.path.dirname(path),year_end,yearday_end,unit_end,'1')
    end_index = info_list.index(path_end)
else:
    path_end = info_list[-1]
    end_index = len(info_list)

print("\n\n*****************************")
print("path_start: ",path_start,'\n')
print("path_end: ",path_end,'\n')
print("totall num of transform: ",end_index-start_index ,"/station's day (unit)")
print("*****************************\n\n")



#%% 3. transform
stats_sac = obspy.core.util.attribdict.AttribDict()                             # generate sac AttribDict in obspy
sta_info={}

for i in range(start_index,end_index+1,1):
    print('Totall be trans days:',i-start_index,'  ||  ','now path:',info_list[i],'\n')
    st_init = Stream()                                                          # inintal empty stream
    all_reftek = sorted(glob.glob(info_list[i]+'/*'))                           # files of one day of one station

    # a. 读取一天的数据 and set channels
    for j in range(0,len(all_reftek),1):
        
        try:
            if network_exists and location_exists:
                st = obspy.read(all_reftek[j], network=network, location=location)
            else:
                st = obspy.read(all_reftek[j])
        except:
            print('OMG！Reading reftek files went wrong！')
            continue
        
        for n in range(0,3,1):
            st[n].stats.channel = component_order[n]
            
        st_init += st

    # b. 合并一天的数据 and get position
    # merge the data and fill the gap with zeros
    try:
        st_init.merge(method=1, fill_value=0)
    except:
        print('OMG！Merge reftek files went wrong！')
        continue
    
    if sample_rate_exists:
        st_init.resample(sample_rate,no_filter=False, strict_length=False)

    position = st_init[0].stats.reftek130['position']
    # print(position)
    if not bool(position.strip()):
        position = st_init[1].stats.reftek130['position']
    if not bool(position.strip()):
        position = st_init[2].stats.reftek130['position']
    if not bool(position.strip()):
        print("this day's position went wrong with GPS! NO TRANS!!\n")
        continue

    # c. write sac data
    stats_sac.stla = float(position[1:10])*0.01     # N
    stats_sac.stlo = float(position[11:20])*0.01    # E
    for qq in range(0,3,1):
        st_init[qq].stats.sac = stats_sac
    
    network_name = st_init[0].stats.network
    station_name = st_init[0].stats.station
    location_name = st_init[0].stats.location
    channel_name = st_init[0].stats.channel
    time = obspy.UTCDateTime(st_init[0].stats.starttime)
    year = "{:04d}".format(time.year)
    julday = "{:03d}".format(time.julday)
    hour = "{:02d}".format(time.hour) 
    minute = "{:02d}".format(time.minute) 
    second = "{:02d}".format(time.second)
    
    sac_path = os.path.join(out_path,
                            network_name,
                            year,
                            network_name+'.'+station_name,
                            julday)
    
    if not os.path.exists(sac_path):
        try:
            os.makedirs(sac_path)
        except OSError as reason:
            print('OMG！Something went wrong！' + str(reason))
            sys.exit()
            
    sac_name = os.path.join(sac_path, 
                            network_name+'.'+
                            station_name+'.'+
                            location_name+'.'+
                            channel_name+'.'+
                            'D'+'.'+                                            # what's the D?
                            year+'.'+
                            julday+'.'+
                            hour+'.'+
                            minute+'.'+
                            second+'.'+
                            'sac')
    
    st_init.write(sac_name, format='SAC')

    # d. mark info
    one_station_info = [year,julday,stats_sac.stla,stats_sac.stlo,info_list[i]]
    # print(station_name,one_station_info)
    if station_name not in list(sta_info.keys()):
        one_list = []
        one_list.append(one_station_info)
        sta_info.update( {station_name: one_list} )
    else:
        sta_info[station_name].append(one_station_info)




#%% 4. output
parameter_path = os.path.join(out_path,'parameter.txt')
fout = open(parameter_path, 'w') 
fout.write( str("totall num of info list: "+str(len(info_list))+"/station's day (unit)"+"\n") )
fout.write( str("path_start: "+path_start+"\n") )
fout.write( str("path_end: "+path_end+"\n") )
fout.write( str("totall num of transform: "+str(end_index-start_index)+"/station's day (unit)"+"\n") )
fout.close()


sta_info_path = os.path.join(out_path,'sta_info.json')  #[year]\t[julday]\t[stla]\t[stlo]\t[reftek_file_path]
with open(sta_info_path,'w') as f1:
    json5.dump(sta_info, f1, indent=2)
f1.close()


info_list_path = os.path.join(out_path,'info_list.txt')
fout = open(info_list_path, 'w') 
for i in range(0,len(info_list),1):
    fout.write( str(info_list[i]) )
    fout.write('\n')
fout.close()




