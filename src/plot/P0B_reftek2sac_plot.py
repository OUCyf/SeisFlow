#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 30 23:01:58 2021

Input

Output

@author: yf
"""

# obspy-scan /bay_mobil/mobil/20090622/1081019/*_1.*

import datetime as dt
import matplotlib.dates as mdates
import os
import json5
import matplotlib.pyplot as plt
import numpy as np
from obspy.core import UTCDateTime


#%% 1. input parameters
# sta_info_path  = '/Users/yf/share1/DATA/yinfu/1.Data/bichun_db/out/sta_info.json'
# plot_output_path = '/Users/yf/share1/DATA/yinfu/1.Data/bichun_db/out/'

sta_info_path  = '/Users/yf/3.Project/2.Velocity_Change/1.Work/sta_info.json'
plot_output_path = '/Users/yf/3.Project/2.Velocity_Change/1.Work'


fig_format = 'pdf'



#%% 2. read file
with open(sta_info_path, 'r',encoding='utf-8') as f1:
    sta_info_r = json5.load(f1)
sta_info = sta_info_r.copy()


#%% 3. plot
fig1, ax = plt.subplots(1, 1, dpi=800,figsize=(12, 6))

labels = list(sta_info.keys())
data = list(sta_info.values())

category_colors = plt.get_cmap('RdYlGn')(
    np.linspace(0.15, 0.85, len(labels)))

for i in range(0,len(labels),1):
    color = category_colors[i]
    for j in range(0,len(data[i])-1,1):
        # a. begin date
        yearday = data[i][j][0]+data[i][j][1]
        UTCtime = UTCDateTime(yearday, iso8601=True)
        dttime = dt.datetime(UTCtime.year, UTCtime.month, UTCtime.day)
        bdate = mdates.date2num(dttime)
        
        # b. end date
        yearday = data[i][j+1][0]+data[i][j+1][1]
        UTCtime = UTCDateTime(yearday, iso8601=True)
        dttime = dt.datetime(UTCtime.year, UTCtime.month, UTCtime.day)
        edate = mdates.date2num(dttime)
        
        # c. barh plot
        ax.barh(labels[i], edate - bdate, left=bdate, height=0.8, align='center',color=color)
        ax.axis('tight')
        
        
# Major ticks every 6 months.
# fmt_half_year = mdates.MonthLocator(interval=1)
fmt_half_year = mdates.DayLocator(interval=2)
ax.xaxis.set_major_locator(fmt_half_year)

# Minor ticks every month.
# fmt_month = mdates.MonthLocator(interval=1)
fmt_month = mdates.DayLocator(interval=1)
ax.xaxis.set_minor_locator(fmt_month)

# Text in the x axis will be displayed in 'YYYY-mm' format.
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))

ax.set_ylabel("Station",fontsize=17)
ax.set_xlabel("Time",fontsize=17)
ax.grid(True)
ax.xaxis_date()
fig1.autofmt_xdate()



#%% 4. save fig
figurename1 = os.path.join(plot_output_path,'continuity.'+fig_format)
fig1.savefig(figurename1, dpi=800, format=fig_format)



#%% 5. plot stla stlo dev


# fig2, ax = plt.subplots(1, 1, dpi=800,figsize=(12, 6))


# for i in range(0,len(labels),1):
#     color = category_colors[i]
    
#     for j in range(0,len(data[i])-1,1):
#         # a. 
#         stla = data[i][j][2]
#         stlo = data[i][j][3]
#         UTCtime = UTCDateTime(yearday, iso8601=True)
#         dttime = dt.datetime(UTCtime.year, UTCtime.month, UTCtime.day)
#         bdate = mdates.date2num(dttime)
        
#         # b. end date
#         yearday = data[i][j+1][0]+data[i][j+1][1]
#         UTCtime = UTCDateTime(yearday, iso8601=True)
#         dttime = dt.datetime(UTCtime.year, UTCtime.month, UTCtime.day)
#         edate = mdates.date2num(dttime)
        
#         # c. barh plot
#         ax.barh(labels[i], edate - bdate, left=bdate, height=0.8, align='center',color=color)
#         ax.axis('tight')





