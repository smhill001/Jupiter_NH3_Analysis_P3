# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 19:53:33 2023

@author: smhil
"""

import sys
drive='C:'
sys.path.append(drive+'/Astronomy/Python Play')
sys.path.append(drive+'/Astronomy/Python Play/Util_P3')
sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy_P3')
#sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy')
sys.path.append(drive+'/Astronomy/Python Play/SPLibraries_P3')
#import scipy.ndimage as nd
from matplotlib.pyplot import imread
import pylab as pl
import numpy as np
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
from datetime import datetime
import ephem
import EWLibV006_P3 as EWL
import plot_TEXES_Groups_P3 as PTG
import NH3_Filter_Library_P3 as NH3
from astropy.convolution import convolve, Box1DKernel

smooth=True
      
belt={"SSTB":[-39.6,-36.2],
      "STB":[-32.4,-27.1],
      "SEB":[-19.7,-7.2],
      "NEB":[6.9,17.4],
      "NTB":[24.2,31.4],
      "NNTB":[35.4,39.6]}

zone={"STZ":[-36.2,-32.4],
      "STrZ":[-27.1,-19.7],
      "EZ":[-7.2,6.9],
      "NTrZ":[17.4,24.2],
      "NTZ":[31.4,35.4]}

path="/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/"
   
figavgprof,axsavgprof=pl.subplots(1,1,figsize=(6.0,4.0), dpi=150, facecolor="white")


PTG.plot_Teifel(axsavgprof,clr='0.5',width=3.)

CMOS2020EW=np.loadtxt(path+
                          '2020 CMOS_NH3_Meridian_EW.csv',
                          usecols=range(3),delimiter=",")
if smooth:
    axsavgprof.plot(CMOS2020EW[:,0],convolve(CMOS2020EW[:,1],Box1DKernel(3),boundary='extend'),
                    color="r",linewidth=1.0,label="2020 Celestron 11")
else:
    axsavgprof.plot(CMOS2020EW[:,0],CMOS2020EW[:,1],color="r",linewidth=1.0,label="2020 Celestron 11")
       

CMOS2021EW=np.loadtxt(path+
                          '2021 CMOS_NH3_Meridian_EW.csv',
                          usecols=range(3),delimiter=",")
if smooth:
    axsavgprof.plot(CMOS2021EW[:,0],convolve(CMOS2021EW[:,1],Box1DKernel(3),boundary='extend'),
                    color="b",linewidth=1.0,label="2021 Celestron 11")
else:
    axsavgprof.plot(CMOS2021EW[:,0],CMOS2021EW[:,1],color="b",linewidth=1.0,label="2021 Celestron 11")

PTG.plot_VLTMUSEandC11_EW_profiles(axsavgprof,"C11 2022",LonRng=1.,CalModel='SCT-Obs-Final',
                                   clr='k',width=1.5,smooth=smooth)
PTG.plot_VLTMUSEandC11_EW_profiles(axsavgprof,"VLTMUSE 2022",LonRng=1.,CalModel='VLT-Obs-Final',
                                   clr='k',width=1.5,style='dashed',smooth=smooth)
#PTG.plot_VLTMUSEandC11_EW_profiles(axsavgprof,"VLTMUSE 2022",clr='k',width=1.5,style='dashed',smooth=True)


for zb in belt:
    print(zb,belt[zb])
    axsavgprof.fill_between([belt[zb][0],belt[zb][1]],np.array([0.,0.]),np.array([1000.,1000.]),
                            color="0.5",alpha=0.2)
    axsavgprof.annotate(zb,xy=[np.mean(belt[zb]),0.01],ha="center")
for zb in zone:
    axsavgprof.annotate(zb,xy=[np.mean(zone[zb]),0.05],ha="center")

# Plot layout details and labeling
axsavgprof.set_title("Ammonia Absorption Profiles")
axsavgprof.grid(linewidth=0.2)
axsavgprof.set_xlim(-45.,45.)
axsavgprof.set_ylim(0.0,1.0)
axsavgprof.set_xticks(np.linspace(-45.,45.,7), minor=False)
axsavgprof.set_yticks(np.linspace(0.0,1.0,6), minor=False)
axsavgprof.tick_params(axis='both', which='major', labelsize=8)
axsavgprof.set_xlabel("Latitude (deg)",fontsize=10)
axsavgprof.set_ylabel("Equivalent Width (nm)",fontsize=10)
axsavgprof.legend(fontsize=8,loc=2)
figavgprof.subplots_adjust(left=0.09, right=0.97, top=0.93, bottom=0.11)
  
figavgprof.savefig(drive+path+"Jupiter-NH3_MeridProfile_Paper.png",dpi=300)