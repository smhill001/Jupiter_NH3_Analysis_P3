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
import Profile_L2 as PFL2

prof="Zonal"

SCT20NH3,SCT21NH3,SCT22NH3,VLT22NH3,SCT23NH3=PFL2.Profile_L2(band="NH3",profile=prof,
                                                      ProfileHalfWidth=1, 
                                                      LatPlotLims=[30,150],
                                                      ZonePlotHalfWidth=60,
                                                      smooth=False)
SCT22CH4,VLT22CH4,SCT23CH4=PFL2.Profile_L2(band="CH4",profile=prof,
                                                      ProfileHalfWidth=1,
                                                      LatPlotLims=[30,150],
                                                      ZonePlotHalfWidth=60,
                                                      smooth=False)

figamfVLT22,axsamfVLT22=pl.subplots(1,2,figsize=(8.0,4.0), dpi=150, facecolor="white")
axsamfVLT22[0].scatter(VLT22NH3['Amf'],VLT22NH3['Pro'],label="NH3")
axsamfVLT22[0].scatter(VLT22CH4['Amf'],VLT22CH4['Pro'],label="CH4")
axsamfVLT22[0].scatter(VLT22NH3['Amf'],VLT22NH3['Pro']*2.6,label="Scaled NH3")
axsamfVLT22[0].scatter(VLT22NH3['Amf'],VLT22NH3['Pro']*VLT22NH3['Amf']**1.0,s=5,label="Corr. NH3")
axsamfVLT22[0].scatter(VLT22CH4['Amf'],VLT22CH4['Pro']*VLT22CH4['Amf']**0.3,s=5,label="Corr. NH3")

axsamfVLT22[0].tick_params(axis='both', which='major', labelsize=8)
axsamfVLT22[0].set_xlabel("One-Way Airmass Factor",fontsize=10)
axsamfVLT22[0].set_ylabel("Equivalent Width (nm)",fontsize=10)
axsamfVLT22[0].set_xlim(1,3)
axsamfVLT22[0].set_xticks(np.linspace(1,3,5), minor=False)

axsamfVLT22[0].set_ylim(0.0,2.0)
axsamfVLT22[0].set_yticks(np.linspace(0.0,2.0,5), minor=False)
axsamfVLT22[0].grid(linewidth=0.2)
axsamfVLT22[0].legend(loc="best",fontsize=8,ncol=2)

axsamfVLT22[1].scatter(VLT22NH3['Amf'],VLT22CH4['Pro']/VLT22NH3['Pro'],label="CH4/NH3")

axsamfVLT22[1].tick_params(axis='both', which='major', labelsize=8)
axsamfVLT22[1].set_ylabel("CH4/NH3 Ratio",fontsize=10)
axsamfVLT22[1].set_xlabel("One-Way Airmass Factor",fontsize=10)
axsamfVLT22[1].set_xlim(1,3)
axsamfVLT22[1].set_xticks(np.linspace(1,3,5), minor=False)
axsamfVLT22[1].set_ylim(0.0,6.0)
axsamfVLT22[1].set_yticks(np.linspace(0.0,6.0,7), minor=False)
axsamfVLT22[1].grid(linewidth=0.2)

figamfVLT22.suptitle("VLT22 Absorption")
figamfVLT22.subplots_adjust(left=0.09, right=0.97, top=0.93, bottom=0.11)

path="c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/"
figamfVLT22.savefig(path+"Profiles/output/IrwinFig9VLT22.png",dpi=300)

figamfSCT22,axsamfSCT22=pl.subplots(1,2,figsize=(8.0,4.0), dpi=150, facecolor="white")
axsamfSCT22[0].scatter(SCT22NH3['Amf'],SCT22NH3['Pro'],label="NH3")
axsamfSCT22[0].scatter(SCT22CH4['Amf'],SCT22CH4['Pro'],label="CH4")
axsamfSCT22[0].scatter(SCT22NH3['Amf'],SCT22NH3['Pro']*2.6,label="Scaled NH3")
axsamfSCT22[0].scatter(SCT22NH3['Amf'],SCT22NH3['Pro']*SCT22NH3['Amf']**1.0,s=5,label="Corr. NH3")
axsamfSCT22[0].scatter(SCT22CH4['Amf'],SCT22CH4['Pro']*SCT22CH4['Amf']**0.3,s=5,label="Corr. CH4")

axsamfSCT22[0].tick_params(axis='both', which='major', labelsize=8)
axsamfSCT22[0].set_xlabel("One-Way Airmass Factor",fontsize=10)
axsamfSCT22[0].set_ylabel("Equivalent Width (nm)",fontsize=10)
axsamfSCT22[0].set_xlim(1,3)
axsamfVLT22[0].set_xticks(np.linspace(1,3,5), minor=False)
axsamfSCT22[0].set_ylim(0.0,2.0)
axsamfSCT22[0].set_yticks(np.linspace(0.0,2.0,5), minor=False)
axsamfSCT22[0].grid(linewidth=0.2)
axsamfSCT22[0].legend(loc="best",fontsize=8,ncol=2)


axsamfSCT22[1].scatter(SCT22NH3['Amf'],SCT22CH4['Pro']/SCT22NH3['Pro'],label="CH4/NH3")

axsamfSCT22[1].tick_params(axis='both', which='major', labelsize=8)
axsamfSCT22[1].set_ylabel("CH4/NH3 Ratio",fontsize=10)
axsamfSCT22[1].set_xlabel("One-Way Airmass Factor",fontsize=10)
axsamfSCT22[1].set_xlim(1,3)
axsamfSCT22[1].set_xticks(np.linspace(1,3,5), minor=False)
axsamfSCT22[1].set_ylim(0.0,6.0)
axsamfSCT22[1].set_yticks(np.linspace(0.0,6.0,7), minor=False)
axsamfSCT22[1].grid(linewidth=0.2)

figamfSCT22.suptitle("SCT22 Absorption")
figamfSCT22.subplots_adjust(left=0.09, right=0.97, top=0.93, bottom=0.11)

path="c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/"
figamfSCT22.savefig(path+"Profiles/output/IrwinFig9SCT22.png",dpi=300)

figamfSCT23,axsamfSCT23=pl.subplots(1,2,figsize=(8.0,4.0), dpi=150, facecolor="white")
axsamfSCT23[0].scatter(SCT23NH3['Amf'],SCT23NH3['Pro'],label="NH3")
axsamfSCT23[0].scatter(SCT23CH4['Amf'],SCT23CH4['Pro'],label="CH4")
axsamfSCT23[0].scatter(SCT23NH3['Amf'],SCT23NH3['Pro']*2.6,label="Scaled NH3")
axsamfSCT23[0].scatter(SCT23NH3['Amf'],SCT23NH3['Pro']*SCT23NH3['Amf']**1.0,s=5,label="Corr. NH3")
axsamfSCT23[0].scatter(SCT23CH4['Amf'],SCT23CH4['Pro']*SCT23CH4['Amf']**0.3,s=5,label="Corr. CH4")

axsamfSCT23[0].tick_params(axis='both', which='major', labelsize=8)
axsamfSCT23[0].set_xlabel("One-Way Airmass Factor",fontsize=10)
axsamfSCT23[0].set_ylabel("Equivalent Width (nm)",fontsize=10)
axsamfSCT23[0].set_xlim(1,3)
axsamfSCT23[0].set_xticks(np.linspace(1,3,5), minor=False)
axsamfSCT23[0].set_ylim(0.0,2.0)
axsamfSCT23[0].set_yticks(np.linspace(0.0,2.0,5), minor=False)
axsamfSCT23[0].grid(linewidth=0.2)
axsamfSCT23[0].legend(loc="best",fontsize=8,ncol=2)


axsamfSCT23[1].scatter(SCT23NH3['Amf'],SCT23CH4['Pro']/SCT23NH3['Pro'],label="CH4/NH3")

axsamfSCT23[1].tick_params(axis='both', which='major', labelsize=8)
axsamfSCT23[1].set_ylabel("CH4/NH3 Ratio",fontsize=10)
axsamfSCT23[1].set_xlabel("One-Way Airmass Factor",fontsize=10)
axsamfSCT23[1].set_xlim(1,3)
axsamfSCT23[1].set_xticks(np.linspace(1,3,5), minor=False)
axsamfSCT23[1].set_ylim(0.0,6.0)
axsamfSCT23[1].set_yticks(np.linspace(0.0,6.0,7), minor=False)
axsamfSCT23[1].grid(linewidth=0.2)

figamfSCT23.suptitle("SCT23 Absorption")
figamfSCT23.subplots_adjust(left=0.09, right=0.97, top=0.93, bottom=0.11)

path="c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/"
figamfSCT23.savefig(path+"Profiles/output/IrwinFig9SCT23.png",dpi=300)


