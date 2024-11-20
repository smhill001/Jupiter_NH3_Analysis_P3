# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 16:51:26 2024

@author: smhil
"""

from matplotlib.pyplot import imread
import pylab as pl
import numpy as np
from os import listdir
import copy
import scipy.ndimage as ndimage
from imageio import imwrite
import get_WINJupos_ephem as WJ_ephem


pathRGB="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/5micron/"

RGBfiles=listdir(pathRGB)
RGBfiles=["j22mapsL3_08.jpg","j22mapsL3_09a.jpg","j22mapsL3_09b.jpg","j22mapsL3_11.jpg","j22mapsL3_13.jpg"]

yranges={"j22mapsL3_08.jpg":[1350,1580],
         "j22mapsL3_09a.jpg":[2360,2590],
         "j22mapsL3_09b.jpg":[1603,1833],
         "j22mapsL3_11.jpg":[1099,1329],
         "j22mapsL3_13.jpg":[ 340, 570]}

outnames={"j22mapsL3_08.jpg":"2022-08-02-0000_0-Jupiter_5micron_RGB_CM2_L360_MAP-BARE",
         "j22mapsL3_09a.jpg":"2022-08-16-1200_0-Jupiter_5micron_RGB_CM2_L360_MAP-BARE",
         "j22mapsL3_09b.jpg":"2022-08-18-0000_0-Jupiter_5micron_RGB_CM2_L360_MAP-BARE",
         "j22mapsL3_11.jpg":"2022-09-03-1200_0-Jupiter_5micron_RGB_CM2_L360_MAP-BARE",
         "j22mapsL3_13.jpg":"2022-09-29-1200_0-Jupiter_5micron_RGB_CM2_L360_MAP-BARE"}

fig1,axs1=pl.subplots(5,1,figsize=(6.0,6.0), dpi=150, facecolor="white",
                      sharex=True,sharey=True)
fig2,axs2=pl.subplots(5,1,figsize=(4.0,6.0), dpi=150, facecolor="white",
                      sharex=True,sharey=True)
#lon_back_plane=np.zeros(230,939,dtype=float)
colat=np.arange(30,150,1)
rows=np.array((1.0-(np.cos(colat*np.pi/180.)/np.cos(30*np.pi/180.)))*115,dtype=int)
mymap=np.zeros((180,360,3),dtype=float)
#print(colat)
#print()
#print(rows)
subplot=0
for RGBfile in RGBfiles:
    print("## ",RGBfile)
    if '.jpg' in RGBfile:
        print()
        print(RGBfile)
        print()
        RGB=imread(pathRGB+RGBfile)
        print(RGB.shape)
        
        outname=outnames[RGBfile]
        RGBtime=(outname[0:10]+"_"+outname[11:13]+":"+outname[13:15]+":00")
        eph=WJ_ephem.get_WINJupos_ephem(RGBtime,planet="Jupiter")
        #time.sleep(5)
        RGB_CM1=float(eph[0].strip())
        RGB_CM2=float(eph[1].strip())
        RGB_CM3=float(eph[2].strip())

        y0=yranges[RGBfile][0]
        y1=yranges[RGBfile][1]
        RGBpatch=copy.deepcopy(RGB[y0:y1,294:1231,:]) #-> row, col, depth
        axs1[subplot].imshow(RGBpatch,origin='upper')

        NH3LonLims=[360.,0.]
        axs1[subplot].xlim=[360-NH3LonLims[0],360-NH3LonLims[1]]
        axs1[subplot].set_xticks(np.linspace(450,0,31), minor=False)
        axs1[4].xticklabels=np.array(np.mod(np.linspace(450,0,31),360))


        longscaledmap=ndimage.zoom(RGBpatch, (1, 0.3842, 1))
        print(longscaledmap.shape)
        colatindx=30
        for row in rows:
            mymap[colatindx,:,:]=longscaledmap[row,:,:]
            #print(colatindx,row,np.max(longscaledmap[row,:,:]),np.max(mymap[colatindx,:,:]))

            colatindx=colatindx+1
        
        mymap=np.roll(mymap,-120,axis=1)
        RGBroll=RGB_CM2-RGB_CM3
        print("########## RGBroll=",RGBroll)
        mymap=np.roll(mymap,int(-RGBroll),axis=1)

        axs2[subplot].imshow(np.array(mymap,dtype=int),origin='upper')
        #axs2[subplot].xlim=[360-NH3LonLims[0],360-NH3LonLims[1]]
        axs2[subplot].set_xticks(np.linspace(360,0,25), minor=False)
        axs2[4].xticklabels=np.array(np.mod(np.linspace(360,0,25),360))
        axs2[subplot].tick_params(axis='both',which='major',labelsize=7)


        imwrite(pathRGB+outnames[RGBfile]+'.png', mymap)#.astype(np.uint16))

        #+60 deg= row 0
        #0 deg = row 114
        #-60 deg = row 229
        #co-latitude runs from 30 to 150
        #co-latitude is 90 - latitude
        #120 long deg = 1232 and 293 (scale= 360/(1232-293) = 0.3834 deg/pix)
        
            #RGBroll=RGB_CM2-RGB_CM1
            #RGB=np.roll(RGB,int(RGBroll),axis=1)


        #axs1[subplot].set_title(RGBfile,fontsize=10)
        #axs1[subplot].set_xlim(225,175)
        #axs1[subplot].set_xticks(np.linspace(225,175,11), minor=False)
        #axs1[subplot].set_ylim(-5,15)
        #axs1[subplot].set_yticks(np.linspace(-5,15,5),minor=False)
        #axs1[subplot].tick_params(axis='both',which='major',labelsize=8)
        #axs1[subplot].grid(linewidth=0.2)

        
        subplot=subplot+1

#for i in range(0,4):
#    axs1[3,i].set_xlabel("Sys. 2 Longitude (deg)",fontsize=9)
#    axs1[i,0].set_ylabel("PG Latitude (deg)")



#pathmapplots="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Studies/NEZ/"
pathmapplots=pathRGB

fig1.subplots_adjust(left=0.05, bottom=0.05, right=0.97, top=0.93,
            wspace=0.10, hspace=0.03)     
fig1.savefig(pathmapplots+"5 micron maps.png",dpi=300)

fig2.subplots_adjust(left=0.05, bottom=0.05, right=0.97, top=0.93,
            wspace=0.10, hspace=0.03)     
fig2.savefig(pathmapplots+"5 micron maps remapped.png",dpi=300)
