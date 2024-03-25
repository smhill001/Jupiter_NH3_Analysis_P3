# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 07:06:42 2024

@author: smhil
"""


import numpy as np
import matplotlib.pyplot as pl
from astropy.io import fits

path='C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/'
hist=np.loadtxt(path+'comparehistogram_nohdr.txt',usecols=range(5))

pathFITS='C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Analysis Data/L3 FITS/'
#CH4file="2022-10-09-0401_5-Jupiter_Img_L3OCH4_S0.fits"
#CH4file="2022-10-13-0345_5-Jupiter_Img_L3OCH4_S0.fits"
#CH4file="2022-10-20-0440_4-Jupiter_Img_L3OCH4_S0.fits"
CH4file="2022-09-19-0352_3-Jupiter_Img_L3OCH4_S0.fits" #VLT

CH4hdulist=fits.open(pathFITS+CH4file)
CH4hdulist.info()
CH4hdr=CH4hdulist[0].header
CH4data=CH4hdulist[0].data
CH4hdulist.close()

b=np.arange(-0.1,0.4,0.001)
histCH4, bin_edges = np.histogram(CH4data.flatten(),bins=b)

#NH3file="2022-10-09-0401_5-Jupiter_Img_L3ONH3_S0.fits"
#NH3file="2022-10-13-0345_5-Jupiter_Img_L3ONH3_S0.fits"
#NH3file="2022-10-20-0440_4-Jupiter_Img_L3ONH3_S0.fits"
NH3file="2022-09-19-0352_3-Jupiter_Img_L3ONH3_S0.fits" #VLT

NH3hdulist=fits.open(pathFITS+NH3file)
NH3hdulist.info()
NH3hdr=NH3hdulist[0].header
NH3data=NH3hdulist[0].data
NH3hdulist.close()

b=np.arange(-0.1,0.4,0.001)
histNH3, bin_edges = np.histogram(NH3data.flatten(),bins=b)

fig,axs=pl.subplots(2,2,figsize=(6.0,6.0), dpi=150, facecolor="white",
                    sharex=True,sharey=True)

axs[0,0].plot(hist[:,0],hist[:,1],label="VLT Irwin CH4")
axs[0,0].plot(hist[:,0],hist[:,2],label="SCT Irwin CH4")
axs[0,0].plot(b[0:len(b)-1]+0.0005,histCH4,label="SCT Hill CH4")
axs[0,0].legend(fontsize=8)
axs[0,1].plot(hist[:,0],hist[:,3])
axs[0,1].plot(hist[:,0],hist[:,4])
axs[0,1].plot(b[0:len(b)-1]+0.0005,histNH3)


axs[1,0].plot(hist[:,0],hist[:,1])
axs[1,0].plot(hist[:,0]+0.04,hist[:,2])
axs[1,0].plot(b[0:len(b)-1]+0.0005+0.04,histCH4)
axs[1,1].plot(hist[:,0],hist[:,3])
axs[1,1].plot(hist[:,0]+0.012,hist[:,4])
axs[1,1].plot(b[0:len(b)-1]+0.0005+0.012,histNH3)

for i in range(0,2):
    for j in range(0,2):
        axs[i,j].set_xlim(0.0,0.3)
        axs[i,j].set_ylim(0,2500)
