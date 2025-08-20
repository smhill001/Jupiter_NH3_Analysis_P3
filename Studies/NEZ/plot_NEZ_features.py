# -*- coding: utf-8 -*-
"""
Created on Mon Aug 18 11:58:23 2025

@author: smhil
"""
import csv
from numpy import genfromtxt
import matplotlib.pyplot as pl
import numpy as np
from datetime import datetime
import matplotlib.dates as mdates

pth="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Studies/NEZ/"
infile="2024-2025 Mean Sys1 0-360 15N-15S blobs fNH3 Tracker.csv"

fig,axs=pl.subplots(2,3,figsize=(12,6), dpi=150, facecolor="white",sharex=True)
fig.suptitle("fNH3 Time Series")
data_array = np.genfromtxt(pth+infile, delimiter=',',dtype=None)
titles=np.array([['fNH3 Long. Centroid','fNH3 Maximum','Mean Cloud Pressure'],
                 ['fNH3 Lat. Centroid','fNH3 Area','Minimum Cloud Pressure']])
ylabels=np.array([['Longitude (deg)','fNH3 (ppm)','Pressure (mb)'],
                  ['PG Latitude (deg)','Area (sq. deg.)','Pressure (mb)']])
columns=np.array([[7,12,13],[6,4,14]])

date_form = mdates.DateFormatter("%Y-%m-%d")

pathmapplots="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Studies/NEZ/"

for pl_col in range(0,2):
    for pl_row in range(0,3):
        col=columns[pl_col,pl_row]
        title=titles[pl_col,pl_row]
        for i in range(1,11):
            filtered_data = data_array[data_array[:,0] == str(i).encode('ascii')]
            t=[]
            for j in range(0,len(filtered_data)):
                t.append(datetime.strptime(filtered_data[j,3].decode("ascii")[:19], "%Y-%m-%dT%H:%M:%S"))
            axs[pl_col,pl_row].scatter(t,
                       np.array(filtered_data[:,col],dtype=float),s=5)
            
        axs[pl_col,pl_row].set_title(title)
        axs[pl_col,pl_row].set_ylabel(ylabels[pl_col,pl_row],fontsize=10)
        axs[pl_col,pl_row].tick_params('x', labelsize=8)
        axs[pl_col,pl_row].xaxis.set_major_formatter(date_form)
        axs[pl_col,pl_row].xaxis.set_major_locator(mdates.MonthLocator(interval=1))
        print(ylabels[pl_col,pl_row])
        if "mb" in ylabels[pl_col,pl_row]:
            axs[pl_col,pl_row].invert_yaxis()


fig.subplots_adjust(left=0.05, bottom=0.05, right=0.98, top=0.90,
            wspace=0.25, hspace=0.15)     
fig.savefig(pathmapplots+"2024-25 fNH3 time series.png",dpi=300)
