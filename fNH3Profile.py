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
import extract_profile as EP

smooth=False
      
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
   
figavgprof,axsavgprof=pl.subplots(1,1,figsize=(6.0,6.0), dpi=150, facecolor="white")


plevel=0.752910
PTG.plot_TEXES_Groups(axsavgprof,clr='C2',prs=plevel,mult=1000000.)
plevel=0.657540
PTG.plot_TEXES_Groups(axsavgprof,clr='C0',prs=plevel,mult=1000000.)
plevel=0.574240
PTG.plot_TEXES_Groups(axsavgprof,clr='C4',prs=plevel,mult=1000000.)

PTG.plot_VLTMUSEandC11_fNH3_profiles(axsavgprof,"SCT 2022",LonRng=1.,CalModel='SCT-Obs-Final',
                                   clr='k',width=1.5,smooth=smooth)
PTG.plot_VLTMUSEandC11_fNH3_profiles(axsavgprof,"VLTMUSE 2022",LonRng=1.,CalModel='VLT-Obs-Final',
                                   clr='k',width=1.0,style='dashed',smooth=smooth)
PTG.plot_VLTMUSEandC11_fNH3_profiles(axsavgprof,"SCT 2023",LonRng=1.,CalModel='SCT-Obs-Final',
                                   clr='r',width=1.5,smooth=smooth)
#PTG.plot_VLTMUSEandC11_EW_profiles(axsavgprof,"VLTMUSE 2022",clr='k',width=1.5,style='dashed',smooth=True)


for zb in belt:
    print(zb,belt[zb])
    axsavgprof.fill_between([belt[zb][0],belt[zb][1]],np.array([0.,0.]),np.array([1000.,1000.]),
                            color="0.5",alpha=0.2)
    axsavgprof.annotate(zb,xy=[np.mean(belt[zb]),0.01],ha="center")
for zb in zone:
    axsavgprof.annotate(zb,xy=[np.mean(zone[zb]),0.05],ha="center")

# Plot layout details and labeling
axsavgprof.set_xlim(-30.,30.)
axsavgprof.set_xlabel("Planetographic Latitude (deg)",fontsize=10)
axsavgprof.set_ylim(0.,300.)
axsavgprof.set_ylabel("Ammonia Abundance (ppm)",fontsize=10)
axsavgprof.legend(loc=1,ncol=2, borderaxespad=0.,prop={'size':8})
axsavgprof.grid(linewidth=0.2)
axsavgprof.tick_params(axis='both', which='major', labelsize=8)

figavgprof.subplots_adjust(left=0.10, bottom=0.12, right=0.98, top=0.92)  
figavgprof.savefig(path+"GRS_NH3_AbundanceProfile_AVG_2022.png",dpi=300)

"""pth='C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Analysis Data/fNH3/'
fn='2023-08-15-1112_8-Jupiter_fNH3.fits'
Lats,AvgMerid,StdMerid=EP.extract_profile(pth,fn,LonRng=10.)
axsavgprof.plot(Lats,AvgMerid*1e6,color='C5')

fn='2023-08-16-1130_1-Jupiter_fNH3.fits'
Lats,AvgMerid,StdMerid=EP.extract_profile(pth,fn,LonRng=10.)
axsavgprof.plot(Lats,AvgMerid*1e6,color='C6')"""