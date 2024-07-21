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
from scipy.optimize import curve_fit

def sorted_reduced(x,y,xmin,xmax):
    """
    Parameters
    ----------
    x : TYPE
        Array of x-values to be ordered by increasing value
    y : TYPE
        Array of y-values to be ordered according to x-values by increasing value
    xmin : TYPE
        Minimum x-value to include.
    xmax : TYPE
        Maximum x-value to include.

    Returns
    -------
    xreduced : TYPE
        Array of x-values order by increasing value and bounded by xmin and xmax.
    yreduced : TYPE
        Array of y-values corresponding to sorted, bounded x-values

    """
    ###VLT22 special treatment
    idx_sort = x.argsort()
    xsort=x[idx_sort]
    ysort=y[idx_sort]

    reducedindices=[i for i, e in enumerate(np.array(xsort)) if (xmin<=e<=xmax)]
    xreduced=xsort[reducedindices]
    yreduced=ysort[reducedindices]

    return(xreduced,yreduced)

def lin(x, a, b):
    return a + (b * x)
def exp(x, a, b):
    return a * np.exp(b * x)
def power(x, a, b):
    return a * x**b


prof="Zonal"
path="c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/"

###############################################################################
# COMPUTE NH3 AND CH4 ABSORPTION PROFILES ACROSS ALL DATA
###############################################################################
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

###############################################################################
# VLT22: SET UP PLOTS
###############################################################################
figamfVLT22,axsamfVLT22=pl.subplots(1,2,figsize=(8.0,4.0), dpi=150, facecolor="white")
figamfVLT22.suptitle("VLT22 Absorption")

###############################################################################
# VLT22: PLOT NH3 and CH4 ABSORPTION DATA and FITS vs JOVIAN AIRMASS
###############################################################################
NH3V22Amf,NH3V22Pro=sorted_reduced(VLT22NH3['Amf'],VLT22NH3['Pro'],1.0,2.5)
VLT22NH3param, VLT22NH3param_cov = curve_fit(power,NH3V22Amf,NH3V22Pro)
fity=power(NH3V22Amf,VLT22NH3param[0],VLT22NH3param[1])
axsamfVLT22[0].plot(NH3V22Amf,fity,label="Fit")
axsamfVLT22[0].scatter(VLT22NH3['Amf'],VLT22NH3['Pro'],label="NH3")

CH4V22Amf,CH4V22Pro=sorted_reduced(VLT22CH4['Amf'],VLT22CH4['Pro'],1.0,2.5)
VLT22CH4param, VLT22CH4param_cov = curve_fit(power,CH4V22Amf,CH4V22Pro)
fity=power(CH4V22Amf,VLT22CH4param[0],VLT22CH4param[1])
axsamfVLT22[0].plot(CH4V22Amf,fity,label="Fit")
axsamfVLT22[0].scatter(VLT22CH4['Amf'],VLT22CH4['Pro'],label="CH4")

###############################################################################
# VLT22: ADJUSTED SCATTER PLOTS OF NH3 and CH4 VERSUS JOVIAN AIRMASS
###############################################################################
axsamfVLT22[0].scatter(VLT22NH3['Amf'],VLT22NH3['Pro']*2.6,label="Scaled NH3")
axsamfVLT22[0].scatter(VLT22NH3['Amf'],VLT22NH3['Pro']*VLT22NH3['Amf']**0.93,s=5,label="Corr. NH3")
axsamfVLT22[0].scatter(VLT22CH4['Amf'],VLT22CH4['Pro']*VLT22CH4['Amf']**0.43,s=5,label="Corr. NH3")

axsamfVLT22[0].tick_params(axis='both', which='major', labelsize=8)
axsamfVLT22[0].set_xlabel("One-Way Airmass Factor",fontsize=10)
axsamfVLT22[0].set_ylabel("Equivalent Width (nm)",fontsize=10)
axsamfVLT22[0].set_xlim(1,3)
axsamfVLT22[0].set_xticks(np.linspace(1,3,5), minor=False)
axsamfVLT22[0].set_ylim(0.0,2.0)
axsamfVLT22[0].set_yticks(np.linspace(0.0,2.0,5), minor=False)
axsamfVLT22[0].grid(linewidth=0.2)
axsamfVLT22[0].legend(loc="best",fontsize=8,ncol=2)

###############################################################################
# VLT22: SCATTER PLOT OF CH4/NH3 RATIO VERSUS JOVIAN AIRMASS
###############################################################################
RV22Amf,RV22Pro=NH3V22Amf,CH4V22Pro/NH3V22Pro
VLT22Rparam, VLT22Rparam_cov = curve_fit(power,RV22Amf,RV22Pro)
fity=lin(RV22Amf,VLT22Rparam[0],VLT22Rparam[1])
axsamfVLT22[1].plot(NH3V22Amf,fity,label="Fit")
axsamfVLT22[1].scatter(RV22Amf,RV22Pro,label="CH4/NH3")

axsamfVLT22[1].plot(NH3V22Amf,power(CH4V22Amf,VLT22CH4param[0],VLT22CH4param[1]) \
                    /power(NH3V22Amf,VLT22NH3param[0],VLT22NH3param[1]),label="Fit2")

axsamfVLT22[1].tick_params(axis='both', which='major', labelsize=8)
axsamfVLT22[1].set_ylabel("CH4/NH3 Ratio",fontsize=10)
axsamfVLT22[1].set_xlabel("One-Way Airmass Factor",fontsize=10)
axsamfVLT22[1].set_xlim(1,3)
axsamfVLT22[1].set_xticks(np.linspace(1,3,5), minor=False)
axsamfVLT22[1].set_ylim(0.0,6.0)
axsamfVLT22[1].set_yticks(np.linspace(0.0,6.0,7), minor=False)
axsamfVLT22[1].grid(linewidth=0.2)

###############################################################################
# VLT22: ADJUST SUBPLOTS and SAVE FIGURE
###############################################################################
figamfVLT22.subplots_adjust(left=0.09, right=0.97, top=0.93, bottom=0.11)
figamfVLT22.savefig(path+"Profiles/output/IrwinFig9VLT22.png",dpi=300)

###############################################################################
###############################################################################
# SCT22: SET UP PLOTS
###############################################################################
figamfSCT22,axsamfSCT22=pl.subplots(1,2,figsize=(8.0,4.0), dpi=150, facecolor="white")
figamfSCT22.suptitle("SCT22 Absorption")

###############################################################################
# SCT22: PLOT NH3 and CH4 ABSORPTION DATA and FITS vs JOVIAN AIRMASS
###############################################################################
NH3S22Amf,NH3S22Pro=sorted_reduced(SCT22NH3['Amf'],SCT22NH3['Pro'],1.0,3.0)
SCT22NH3param, SCT22NH3param_cov = curve_fit(power,NH3S22Amf,NH3S22Pro)
fity=power(NH3S22Amf,SCT22NH3param[0],SCT22NH3param[1])
axsamfSCT22[0].plot(NH3S22Amf,fity,label="Fit")
axsamfSCT22[0].scatter(SCT22NH3['Amf'],SCT22NH3['Pro'],label="NH3")

CH4S22Amf,CH4S22Pro=sorted_reduced(SCT22CH4['Amf'],SCT22CH4['Pro'],1.0,3.0)
SCT22CH4param, SCT22CH4param_cov = curve_fit(power,CH4S22Amf,CH4S22Pro)
fity=power(CH4S22Amf,SCT22CH4param[0],SCT22CH4param[1])
axsamfSCT22[0].plot(CH4S22Amf,fity,label="Fit")
axsamfSCT22[0].scatter(SCT22CH4['Amf'],SCT22CH4['Pro'],label="NH3")

###############################################################################
# SCT22: ADJUSTED SCATTER PLOTS OF NH3 and CH4 VERSUS JOVIAN AIRMASS
###############################################################################
axsamfSCT22[0].scatter(SCT22NH3['Amf'],SCT22NH3['Pro']*2.6,label="Scaled NH3")
axsamfSCT22[0].scatter(SCT22NH3['Amf'],SCT22NH3['Pro']*SCT22NH3['Amf']**0.93,s=5,label="Corr. NH3")
axsamfSCT22[0].scatter(SCT22CH4['Amf'],SCT22CH4['Pro']*SCT22CH4['Amf']**0.43,s=5,label="Corr. CH4")

axsamfSCT22[0].tick_params(axis='both', which='major', labelsize=8)
axsamfSCT22[0].set_xlabel("One-Way Airmass Factor",fontsize=10)
axsamfSCT22[0].set_ylabel("Equivalent Width (nm)",fontsize=10)
axsamfSCT22[0].set_xlim(1,3)
axsamfVLT22[0].set_xticks(np.linspace(1,3,5), minor=False)
axsamfSCT22[0].set_ylim(0.0,2.0)
axsamfSCT22[0].set_yticks(np.linspace(0.0,2.0,5), minor=False)
axsamfSCT22[0].grid(linewidth=0.2)
axsamfSCT22[0].legend(loc="best",fontsize=8,ncol=2)

###############################################################################
# SCT22: SCATTER PLOT OF CH4/NH3 RATIO VERSUS JOVIAN AIRMASS
###############################################################################
axsamfSCT22[1].scatter(SCT22NH3['Amf'],SCT22CH4['Pro']/SCT22NH3['Pro'],label="CH4/NH3")
axsamfSCT22[1].tick_params(axis='both', which='major', labelsize=8)
axsamfSCT22[1].set_ylabel("CH4/NH3 Ratio",fontsize=10)
axsamfSCT22[1].set_xlabel("One-Way Airmass Factor",fontsize=10)
axsamfSCT22[1].set_xlim(1,3)
axsamfSCT22[1].set_xticks(np.linspace(1,3,5), minor=False)
axsamfSCT22[1].set_ylim(0.0,6.0)
axsamfSCT22[1].set_yticks(np.linspace(0.0,6.0,7), minor=False)
axsamfSCT22[1].grid(linewidth=0.2)

###############################################################################
# SCT22: ADJUST SUBPLOTS and SAVE FIGURE
###############################################################################
figamfSCT22.subplots_adjust(left=0.09, right=0.97, top=0.93, bottom=0.11)
figamfSCT22.savefig(path+"Profiles/output/IrwinFig9SCT22.png",dpi=300)

###############################################################################
###############################################################################
# VCT23: SET UP PLOTS
###############################################################################
figamfSCT23,axsamfSCT23=pl.subplots(1,2,figsize=(8.0,4.0), dpi=150, facecolor="white")
figamfSCT23.suptitle("SCT23 Absorption")

###############################################################################
# SCT23: PLOT NH3 and CH4 ABSORPTION DATA and FITS vs JOVIAN AIRMASS
###############################################################################
NH3S23Amf,NH3S23Pro=sorted_reduced(SCT23NH3['Amf'],SCT23NH3['Pro'],1.0,3.0)
SCT23NH3param, SCT23NH3param_cov = curve_fit(power,NH3S23Amf,NH3S23Pro)
fity=power(NH3S23Amf,SCT23NH3param[0],SCT23NH3param[1])
axsamfSCT23[0].plot(NH3S23Amf,fity,label="Fit")
axsamfSCT23[0].scatter(SCT23NH3['Amf'],SCT23NH3['Pro'],label="NH3")

CH4S23Amf,CH4S23Pro=sorted_reduced(SCT23CH4['Amf'],SCT23CH4['Pro'],1.0,3.0)
SCT23CH4param, SCT23CH4param_cov = curve_fit(power,CH4S23Amf,CH4S23Pro)
fity=power(CH4S23Amf,SCT23CH4param[0],SCT23CH4param[1])
axsamfSCT23[0].plot(CH4S23Amf,fity,label="Fit")
axsamfSCT23[0].scatter(SCT23CH4['Amf'],SCT23CH4['Pro'],label="NH3")

###############################################################################
# SCT23: ADJUSTED SCATTER PLOTS OF NH3 and CH4 VERSUS JOVIAN AIRMASS
###############################################################################
axsamfSCT23[0].scatter(SCT23NH3['Amf'],SCT23NH3['Pro']*2.6,label="Scaled NH3")
axsamfSCT23[0].scatter(SCT23NH3['Amf'],SCT23NH3['Pro']*SCT23NH3['Amf']**0.93,s=5,label="Corr. NH3")
axsamfSCT23[0].scatter(SCT23CH4['Amf'],SCT23CH4['Pro']*SCT23CH4['Amf']**0.43,s=5,label="Corr. CH4")

axsamfSCT23[0].tick_params(axis='both', which='major', labelsize=8)
axsamfSCT23[0].set_xlabel("One-Way Airmass Factor",fontsize=10)
axsamfSCT23[0].set_ylabel("Equivalent Width (nm)",fontsize=10)
axsamfSCT23[0].set_xlim(1,3)
axsamfSCT23[0].set_xticks(np.linspace(1,3,5), minor=False)
axsamfSCT23[0].set_ylim(0.0,2.0)
axsamfSCT23[0].set_yticks(np.linspace(0.0,2.0,5), minor=False)
axsamfSCT23[0].grid(linewidth=0.2)
axsamfSCT23[0].legend(loc="best",fontsize=8,ncol=2)

###############################################################################
# SCT23: SCATTER PLOT OF CH4/NH3 RATIO VERSUS JOVIAN AIRMASS
###############################################################################
axsamfSCT23[1].scatter(SCT23NH3['Amf'],SCT23CH4['Pro']/SCT23NH3['Pro'],label="CH4/NH3")
axsamfSCT23[1].tick_params(axis='both', which='major', labelsize=8)
axsamfSCT23[1].set_ylabel("CH4/NH3 Ratio",fontsize=10)
axsamfSCT23[1].set_xlabel("One-Way Airmass Factor",fontsize=10)
axsamfSCT23[1].set_xlim(1,3)
axsamfSCT23[1].set_xticks(np.linspace(1,3,5), minor=False)
axsamfSCT23[1].set_ylim(0.0,6.0)
axsamfSCT23[1].set_yticks(np.linspace(0.0,6.0,7), minor=False)
axsamfSCT23[1].grid(linewidth=0.2)

###############################################################################
# SCT23: ADJUST SUBPLOTS and SAVE FIGURE
###############################################################################
figamfSCT23.subplots_adjust(left=0.09, right=0.97, top=0.93, bottom=0.11)
figamfSCT23.savefig(path+"Profiles/output/IrwinFig9SCT23.png",dpi=300)

print()
print("######################################################################")
print("NH3 ##################################################################")
print(VLT22NH3param, VLT22NH3param_cov)
print(SCT22NH3param, SCT22NH3param_cov)
print(SCT23NH3param, SCT23NH3param_cov)
print("CH4 ##################################################################")
print(VLT22CH4param, VLT22CH4param_cov)
print(SCT22CH4param, SCT22CH4param_cov)
print(SCT23CH4param, SCT23CH4param_cov)
print("######################################################################")
print()