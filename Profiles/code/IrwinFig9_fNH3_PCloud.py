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
import IrwinLib as IL

prof="Zonal"
path="c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/"

###############################################################################
# COMPUTE fNH3 AND PCloud (USUALLY ZONAL) PROFILES ACROSS ALL DATA
###############################################################################
SCT22NH3,VLT22NH3,SCT23NH3,SCT24NH3=PFL3.Profile_L3(param="fNH3",profile=prof,
                                                      ProfileHalfWidth=1, 
                                                      LatPlotLims=[30,150],
                                                      ZonePlotHalfWidth=60,
                                                      smooth=False)
SCT22CH4,VLT22CH4,SCT23CH4,SCT24CH4=PFL3.Profile_L3(param="PCld",profile=prof,
                                                      ProfileHalfWidth=1,
                                                      LatPlotLims=[30,150],
                                                      ZonePlotHalfWidth=60,
                                                      smooth=False)

###############################################################################
# VLT22: SET UP PLOTS
###############################################################################
figamfVLT22,axsamfVLT22=pl.subplots(1,2,figsize=(8.0,4.0), dpi=150, facecolor="white")
figamfVLT22.suptitle("VLT22 PCloud and fNH3")

VLT22NH3param, VLT22NH3param_cov,VLT22NH3R2,VLT22CH4param, VLT22CH4param_cov,VLT22CH4R2,VLT22RatioR2,VLT22fig=\
    IL.make_fit_plot(VLT22NH3['Amf'],VLT22NH3['Pro'],VLT22CH4['Amf'],
                     VLT22CH4['Pro'],"L3",figamfVLT22,axsamfVLT22,path)



###############################################################################
###############################################################################
# SCT22: SET UP PLOTS
###############################################################################
figamfSCT22,axsamfSCT22=pl.subplots(1,2,figsize=(8.0,4.0), dpi=150, facecolor="white")
figamfSCT22.suptitle("SCT22 PCloud and fNH3")

SCT22NH3param, SCT22NH3param_cov,SCT22NH3R2,SCT22CH4param, SCT22CH4param_cov,SCT22CH4R2,SCT22RatioR2,SCT22fig=\
    IL.make_fit_plot(SCT22NH3['Amf'],SCT22NH3['Pro'],SCT22CH4['Amf'],
                     SCT22CH4['Pro'],"L3",figamfSCT22,axsamfSCT22,path)



###############################################################################
###############################################################################
# SCT23: SET UP PLOTS
###############################################################################
figamfSCT23,axsamfSCT23=pl.subplots(1,2,figsize=(8.0,4.0), dpi=150, facecolor="white")
figamfSCT23.suptitle("SCT23 PCloud and fNH3")

SCT23NH3param, SCT23NH3param_cov,SCT23NH3R2,SCT23CH4param, SCT23CH4param_cov,SCT23CH4R2,SCT23RatioR2,SCT23fig=\
    IL.make_fit_plot(SCT23NH3['Amf'],SCT23NH3['Pro'],SCT23CH4['Amf'],
                     SCT23CH4['Pro'],"L3",figamfSCT23,axsamfSCT23,path)

###############################################################################
###############################################################################
# SCT24: SET UP PLOTS
###############################################################################
figamfSCT24,axsamfSCT24=pl.subplots(1,2,figsize=(8.0,4.0), dpi=150, facecolor="white")
figamfSCT24.suptitle("SCT24 PCloud and fNH3")

SCT24NH3param, SCT24NH3param_cov,SCT24NH3R2,SCT24CH4param, SCT24CH4param_cov,SCT24CH4R2,SCT24RatioR2,SCT24fig=\
    IL.make_fit_plot(SCT24NH3['Amf'],SCT24NH3['Pro'],SCT24CH4['Amf'],
                     SCT24CH4['Pro'],"L3",figamfSCT24,axsamfSCT24,path)

print()
print("######################################################################")
print("NH3 ##################################################################")
print(VLT22NH3param, VLT22NH3param_cov)
print(VLT22NH3R2)
print(SCT22NH3param, SCT22NH3param_cov)
print(SCT22NH3R2)
print(SCT23NH3param, SCT23NH3param_cov)
print(SCT23NH3R2)
print("CH4 ##################################################################")
print(VLT22CH4param, VLT22CH4param_cov)
print(VLT22CH4R2)
print(SCT22CH4param, SCT22CH4param_cov)
print(SCT22CH4R2)
print(SCT23CH4param, SCT23CH4param_cov)
print(SCT23CH4R2)
print("######################################################################")
print()

L3fitdata = open(path+'/Profiles/output/L3fitdata.csv', 'w')
tmp="Observatory,Year,Gas,Constant,ExponentCH4,R2\n"
L3fitdata.write(tmp)

tmp="VLT,2022,NH3,"+"{:.4f}".format(VLT22NH3param[0])+","+"{:.3f}".format(VLT22NH3param[1])+","+"{:.3f}".format(VLT22NH3R2)+"\n"
L3fitdata.write(tmp)
tmp="SCT,2022,NH3,"+"{:.4f}".format(SCT22NH3param[0])+","+"{:.3f}".format(SCT22NH3param[1])+","+"{:.3f}".format(SCT22NH3R2)+"\n"
L3fitdata.write(tmp)
tmp="SCT,2023,NH3,"+"{:.4f}".format(SCT23NH3param[0])+","+"{:.3f}".format(SCT23NH3param[1])+","+"{:.3f}".format(SCT23NH3R2)+"\n"
L3fitdata.write(tmp)
tmp="VLT,2022,CH4,"+"{:.4f}".format(VLT22CH4param[0])+","+"{:.3f}".format(VLT22CH4param[1])+","+"{:.3f}".format(VLT22CH4R2)+"\n"
L3fitdata.write(tmp)
tmp="SCT,2022,CH4,"+"{:.4f}".format(SCT22CH4param[0])+","+"{:.3f}".format(SCT22CH4param[1])+","+"{:.3f}".format(SCT22CH4R2)+"\n"
L3fitdata.write(tmp)
tmp="SCT,2023,CH4,"+"{:.4f}".format(SCT23CH4param[0])+","+"{:.3f}".format(SCT23CH4param[1])+","+"{:.3f}".format(SCT23CH4R2)+"\n"
L3fitdata.write(tmp)

tmp="VLT,2022,CH4/NH3,"+"{:.4f}".format(VLT22CH4param[0]/VLT22NH3param[0])+","+\
    "{:.3f}".format(VLT22NH3param[1]-VLT22CH4param[1])+","+"{:.3f}".format(VLT22RatioR2)+"\n"
L3fitdata.write(tmp)
tmp="SCT,2022,CH4/NH3,"+"{:.4f}".format(SCT22CH4param[0]/SCT22NH3param[0])+","+\
    "{:.3f}".format(SCT22NH3param[1]-SCT22CH4param[1])+","+"{:.3f}".format(SCT22RatioR2)+"\n"
L3fitdata.write(tmp)
tmp="SCT,2023,CH4/NH3,"+"{:.4f}".format(SCT23CH4param[0]/SCT23NH3param[0])+","+\
    "{:.3f}".format(SCT23NH3param[1]-SCT23CH4param[1])+","+"{:.3f}".format(SCT23RatioR2)+"\n"
L3fitdata.write(tmp)


L3fitdata.close()
