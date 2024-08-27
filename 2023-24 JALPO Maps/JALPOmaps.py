"""
Created on Mon Jul  1 13:32:12 2024

@author: smhil
"""

from matplotlib.pyplot import imread
import pylab as pl
import numpy as np
from os import listdir

pathRGB="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/2023-24 Maps/"

RGBfiles=listdir(pathRGB)

fig1,axs1=pl.subplots(4,4,figsize=(12.0,6.0), dpi=150, facecolor="white",
                      sharex=True,sharey=True)
subplot=0
for RGBfile in RGBfiles:
    print("## ",RGBfile)
    if '.jpg' in RGBfile:
        print()
        print(RGBfile)
        print()
        RGB=imread(pathRGB+RGBfile)
        axs1[np.mod(subplot,4),int(subplot/4)].imshow(RGB[2790:2874,975:1185,:],
                   origin='upper',extent=[225,175,-5,15],aspect="equal")
        axs1[np.mod(subplot,4),int(subplot/4)].set_title(RGBfile,fontsize=10)
        axs1[np.mod(subplot,4),int(subplot/4)].set_xlim(225,175)
        axs1[np.mod(subplot,4),int(subplot/4)].set_xticks(np.linspace(225,175,11),
                                                          minor=False)
        axs1[np.mod(subplot,4),int(subplot/4)].set_ylim(-5,15)
        axs1[np.mod(subplot,4),int(subplot/4)].set_yticks(np.linspace(-5,15,5),
                                                          minor=False)
        axs1[np.mod(subplot,4),int(subplot/4)].tick_params(axis='both', 
                                                            which='major', 
                                                            labelsize=8)
        axs1[np.mod(subplot,4),int(subplot/4)].grid(linewidth=0.2)

        
        subplot=subplot+1

for i in range(0,4):
    axs1[3,i].set_xlabel("Sys. 2 Longitude (deg)",fontsize=9)
    axs1[i,0].set_ylabel("PG Latitude (deg)")



pathmapplots="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Studies/NEZ/"

fig1.subplots_adjust(left=0.05, bottom=0.05, right=0.97, top=0.93,
            wspace=0.10, hspace=0.03)     
fig1.savefig(pathmapplots+"JALPOmaps.png",dpi=300)
