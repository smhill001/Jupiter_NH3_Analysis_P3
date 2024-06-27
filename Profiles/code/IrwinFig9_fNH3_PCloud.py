# -*- coding: utf-8 -*-
"""
Created on Sun Jun 23 11:44:27 2024

@author: smhil
"""

"""
Created on Tue Jun  4 12:39:15 2024

@author: smhil
"""
import sys
sys.path.append('c:/Astronomy/Python Play')
import pylab as pl
import numpy as np
from astropy.convolution import convolve
import plot_TEXES_Groups_P3 as PTG
from astropy.convolution import Box1DKernel
sys.path.append('./Maps')
import Profile_L3 as PFL3

prof="Zonal"

SCT22NH3,VLT22NH3,SCT23NH3=PFL3.Profile_L3(param="fNH3",profile=prof,
                                                      ProfileHalfWidth=1, 
                                                      LatPlotLims=[30,150],
                                                      ZonePlotHalfWidth=60,
                                                      smooth=False)
SCT22CH4,VLT22CH4,SCT23CH4=PFL3.Profile_L3(param="PCld",profile=prof,
                                                      ProfileHalfWidth=1,
                                                      LatPlotLims=[30,150],
                                                      ZonePlotHalfWidth=60,
                                                      smooth=False)

figamfVLT22,axsamfVLT22=pl.subplots(1,2,figsize=(8.0,4.0), dpi=150, facecolor="white")
figamfVLT22.suptitle("VLT22 f(NH3) and PCloud")
axsamfVLT22[0].scatter(VLT22NH3['Amf'],VLT22NH3['Pro'],label="fNH3")
axsamfVLT22[0].scatter(VLT22CH4['Amf'],VLT22CH4['Pro'],label="PCloud")
axsamfVLT22[0].scatter(VLT22NH3['Amf'],VLT22NH3['Pro']*10.,label="Scaled fNH3")
axsamfVLT22[0].scatter(VLT22NH3['Amf'],VLT22NH3['Pro']*VLT22NH3['Amf']**1.0,s=5,label="Corr. fNH3")
axsamfVLT22[0].scatter(VLT22CH4['Amf'],VLT22CH4['Pro']*VLT22CH4['Amf']**0.3,s=5,label="Corr. PCloud")

axsamfVLT22[0].tick_params(axis='both', which='major', labelsize=8)
axsamfVLT22[0].set_xlabel("One-Way Airmass Factor",fontsize=10)
axsamfVLT22[0].set_ylabel("f(NH3) (ppm) and PCloud (mb)",fontsize=10)
axsamfVLT22[0].legend(fontsize=6,loc="lower center")
axsamfVLT22[0].set_xlim(1,3)
axsamfVLT22[0].set_ylim(0.0,1200.0)
axsamfVLT22[0].set_yticks(np.linspace(0.0,1200.0,7), minor=False)
axsamfVLT22[0].grid(linewidth=0.2)
axsamfVLT22[0].legend(loc="best",fontsize=8,ncol=2)


axsamfVLT22[1].scatter(VLT22NH3['Amf'],VLT22CH4['Pro']/VLT22NH3['Pro'],label="PCloud/fNH3")

axsamfVLT22[1].tick_params(axis='both', which='major', labelsize=8)
axsamfVLT22[1].set_ylabel("Ratio (PCloud/f(NH3)",fontsize=10)
axsamfVLT22[1].set_xlabel("One-Way Airmass Factor",fontsize=10)
axsamfVLT22[1].legend(fontsize=6,loc="lower center")
axsamfVLT22[1].set_xlim(1,3)
axsamfVLT22[1].set_xticks(np.linspace(1,3,5), minor=False)
axsamfVLT22[1].set_ylim(0.0,20.0)
axsamfVLT22[1].set_yticks(np.linspace(0.0,20.0,5), minor=False)
axsamfVLT22[1].grid(linewidth=0.2)
axsamfVLT22[0].legend(loc="best",fontsize=8,ncol=2)

figamfVLT22.subplots_adjust(left=0.09, right=0.97, top=0.93, bottom=0.11)

path="c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/"
figamfVLT22.savefig(path+"Profiles/output/IrwinFig9VLT22fNH3.png",dpi=300)


figamfSCT22,axsamfSCT22=pl.subplots(1,2,figsize=(8.0,4.0), dpi=150, facecolor="white")
figamfSCT22.suptitle("SCT22 f(NH3) and PCloud")

axsamfSCT22[0].scatter(SCT22NH3['Amf'],SCT22NH3['Pro'],label="fNH3")
axsamfSCT22[0].scatter(SCT22CH4['Amf'],SCT22CH4['Pro'],label="PCloud")
axsamfSCT22[0].scatter(SCT22NH3['Amf'],SCT22NH3['Pro']*10,label="Scaled fNH3")
axsamfSCT22[0].scatter(SCT22NH3['Amf'],SCT22NH3['Pro']*SCT22NH3['Amf']**1.0,s=5,label="Corr. fNH3")
axsamfSCT22[0].scatter(SCT22CH4['Amf'],SCT22CH4['Pro']*SCT22CH4['Amf']**0.3,s=5,label="Corr. PCloud")

axsamfSCT22[0].tick_params(axis='both', which='major', labelsize=8)
axsamfSCT22[0].set_xlabel("One-Way Airmass Factor",fontsize=10)
axsamfSCT22[0].set_ylabel("f(NH3) (ppm) and PCloud (mb)",fontsize=10)
axsamfSCT22[0].legend(fontsize=6,loc="lower center")
axsamfSCT22[0].set_xlim(1,3)
axsamfSCT22[0].set_ylim(0.0,1200.0)
axsamfSCT22[0].set_yticks(np.linspace(0.0,1200.0,7), minor=False)
axsamfSCT22[0].grid(linewidth=0.2)
axsamfSCT22[0].legend(loc="best",fontsize=8,ncol=2)


axsamfSCT22[1].scatter(SCT22NH3['Amf'],SCT22CH4['Pro']/SCT22NH3['Pro'],label="PCloud/fNH3")

axsamfSCT22[1].tick_params(axis='both', which='major', labelsize=8)
axsamfSCT22[1].set_ylabel("Ratio (PCloud/f(NH3)",fontsize=10)
axsamfSCT22[1].set_xlabel("One-Way Airmass Factor",fontsize=10)
axsamfSCT22[1].legend(fontsize=6,loc="lower center")
axsamfSCT22[1].set_xlim(1,3)
axsamfSCT22[1].set_xticks(np.linspace(1,3,5), minor=False)
axsamfSCT22[1].set_ylim(0.0,20.0)
axsamfSCT22[1].set_yticks(np.linspace(0.0,20.0,5), minor=False)
axsamfSCT22[1].grid(linewidth=0.2)

figamfSCT22.subplots_adjust(left=0.09, right=0.97, top=0.93, bottom=0.11)

path="c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/"
figamfSCT22.savefig(path+"Profiles/output/IrwinFig9SCT22fNH3.png",dpi=300)

figamfSCT23,axsamfSCT23=pl.subplots(1,2,figsize=(8.0,4.0), dpi=150, facecolor="white")
figamfSCT23.suptitle("SCT22 f(NH3) and PCloud")

axsamfSCT23[0].scatter(SCT23NH3['Amf'],SCT23NH3['Pro'],label="fNH3")
axsamfSCT23[0].scatter(SCT23CH4['Amf'],SCT23CH4['Pro'],label="PCloud")
axsamfSCT23[0].scatter(SCT23NH3['Amf'],SCT23NH3['Pro']*10,label="Scaled NH3")
axsamfSCT23[0].scatter(SCT23NH3['Amf'],SCT23NH3['Pro']*SCT23NH3['Amf']**1.0,s=5,label="Corr. fNH3")
axsamfSCT23[0].scatter(SCT23CH4['Amf'],SCT23CH4['Pro']*SCT23CH4['Amf']**0.3,s=5,label="Corr. PCloud")

axsamfSCT23[0].tick_params(axis='both', which='major', labelsize=8)
axsamfSCT23[0].set_xlabel("One-Way Airmass Factor",fontsize=10)
axsamfSCT23[0].set_ylabel("f(NH3) (ppm) and PCloud (mb)",fontsize=10)
axsamfSCT23[0].legend(fontsize=6,loc="lower center")
axsamfSCT23[0].set_xlim(1,3)
axsamfSCT23[0].set_ylim(0.0,1200.0)
axsamfSCT23[0].set_yticks(np.linspace(0.0,1200.0,7), minor=False)
axsamfSCT23[0].grid(linewidth=0.2)
axsamfSCT22[0].legend(loc="best",fontsize=8,ncol=2)


axsamfSCT23[1].scatter(SCT23NH3['Amf'],SCT23CH4['Pro']/SCT23NH3['Pro'],label="PCloud/fNH3")

axsamfSCT23[1].tick_params(axis='both', which='major', labelsize=8)
axsamfSCT23[1].set_ylabel("Ratio (PCloud/f(NH3)",fontsize=10)
axsamfSCT23[1].set_xlabel("One-Way Airmass Factor",fontsize=10)
axsamfSCT23[1].legend(fontsize=6,loc="lower center")
axsamfSCT23[1].set_xlim(1,3)
axsamfSCT23[1].set_xticks(np.linspace(1,3,5), minor=False)
axsamfSCT23[1].set_ylim(0.0,20.0)
axsamfSCT23[1].set_yticks(np.linspace(0.0,20.0,5), minor=False)
axsamfSCT23[1].grid(linewidth=0.2)

figamfSCT23.subplots_adjust(left=0.09, right=0.97, top=0.93, bottom=0.11)

path="c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/"
figamfSCT23.savefig(path+"Profiles/output/IrwinFig9SCT23fNH3.png",dpi=300)


