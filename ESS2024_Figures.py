"""
Created on Fri Apr 26 10:17:21 2024

@author: smhil
"""

import sys
import os
drive='c:'
sys.path.append(drive+'/Astronomy/Python Play')
sys.path.append(drive+'/Astronomy/Python Play/Util_P3')
sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy')
sys.path.append('./Molecular Absorption/code')
sys.path.append('./Services')
sys.path.append('./Maps')
sys.path.append('./Profiles/code')
sys.path.append('./Studies/GRS')

import numpy as np
from astropy.io import fits
import Map_Jup_Atm_P3 as MJA
import get_albedo_continua_crossections as GACC
import Profile_L2 as PL2
import Profile_L3 as PL3
import image_array_new as IA
import grs_zoom_maps as GZM
import get_obs_list as GOL
import shutil

pathout='C:/Astronomy/Projects/SAS 2021 Ammonia/Publications/2024-ESS-Data/'
sourcefiles=GOL.get_obs_list()
pathFITS='C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Analysis Data/L3 FITS/'
pathFITS2='C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Analysis Data/L2 FITS/'

CH4suffix="-Jupiter_Map_L3PCld_S0"
NH3suffix="-Jupiter_Map_L3fNH3_S0"

CH4suffixI2="-Jupiter_Img_L2TCH4"
NH3suffixI2="-Jupiter_Img_L2TNH3"
CH4suffixI3="-Jupiter_Img_L3PCld_S0"
NH3suffixI3="-Jupiter_Img_L3fNH3_S0"

###############################################################################
# FIGURE 1
###############################################################################
Albedo,Continua,CH4,NH3,NH3_Lutz_Owen_1980,fig_molecules,ax_molecules= \
    GACC.get_albedo_continua_crossections(600,680,9,0.0,0.7,8,LutzPlot=False,
                                          Crossect=True, ModList=[])

fig_molecules.savefig(pathout+'Fig1/Fig1.png',dpi=300)

###############################################################################
# FIGURE 2
###############################################################################
fig2=IA.image_array_new(obsdate="20221009UTa",target="Jupiter",imagetype='Img',
                        contour=False)

fig2.savefig(pathout+'Fig2/Fig2.png',dpi=300)

try:
    TCH4file=sourcefiles["20221009UTa"]['CH4file'][0:17]+CH4suffixI2+\
            sourcefiles["20221009UTa"]['Metadata']['Variation']+".fits"
    variation=sourcefiles["20221009UTa"]['Metadata']['Variation']
except:
    TCH4file=sourcefiles["20221009UTa"]['CH4file'][0:17]+CH4suffixI2+".fits"
    variation=""
try:
    TNH3file=sourcefiles["20221009UTa"]['NH3file'][0:17]+NH3suffixI2+\
            sourcefiles["20221009UTa"]['Metadata']['Variation']+".fits"
    variation=sourcefiles["20221009UTa"]['Metadata']['Variation']
except:
    TNH3file=sourcefiles["20221009UTa"]['NH3file'][0:17]+NH3suffixI2+".fits"
    variation=""

shutil.copyfile(pathFITS2+TCH4file,pathout+"FITS/"+TCH4file)
shutil.copyfile(pathFITS2+TNH3file,pathout+"FITS/"+TNH3file)

try:
    PCloudfileI=sourcefiles["20221009UTa"]['CH4file'][0:17]+CH4suffixI3+\
            sourcefiles["20221009UTa"]['Metadata']['Variation']+".fits"
    variation=sourcefiles["20221009UTa"]['Metadata']['Variation']
except:
    PCloudfileI=sourcefiles["20221009UTa"]['CH4file'][0:17]+CH4suffixI3+".fits"
    variation=""
try:
    fNH3fileI=sourcefiles["20221009UTa"]['NH3file'][0:17]+NH3suffixI3+\
            sourcefiles["20221009UTa"]['Metadata']['Variation']+".fits"
    variation=sourcefiles["20221009UTa"]['Metadata']['Variation']
except:
    fNH3fileI=sourcefiles["20221009UTa"]['NH3file'][0:17]+NH3suffixI3+".fits"
    variation=""
    
shutil.copyfile(pathFITS+PCloudfileI,pathout+"FITS/"+PCloudfileI)
shutil.copyfile(pathFITS+fNH3fileI,pathout+"FITS/"+fNH3fileI)
###############################################################################
# FIGURE 3
###############################################################################
fig3,axsfig3=PL2.Profile_L2(band="NH3",profile="Meridional",
                                              ProfileHalfWidth=1,LatPlotLims=[45,135],
                                              ZonePlotHalfWidth=45,smooth=False)

fig3.savefig(pathout+'Fig3/Fig3.png',dpi=300)

###############################################################################
# FIGURE 4
###############################################################################
fig4,axsfig4=PL3.Profile_L3(param="fNH3",profile="Meridional",ProfileHalfWidth=1,
               LatPlotLims=[60,120],ZonePlotHalfWidth=45,smooth=False,
               inset=True)

fig4.savefig(pathout+'Fig4/Fig4.png',dpi=300)

###############################################################################
# FIGURE 5
###############################################################################
fig1,axs1,fig2,axs2,fig3,axs3=MJA.Map_Jup_Atm_P3(obskey="20221009UTa",
                                                 imagetype='Map', 
                                                 Smoothing=False,
                                                 LatLims=[60,120],LonRng=30, 
                                                 CMpref='subobs',LonSys='2',
                                                 showbands=False,coef=[0.,0.],
                                                 subproj='paper')
"""
if Level=="L2":
    CH4suffix="-Jupiter_"+imagetype+"_L2TCH4"
    NH3suffix="-Jupiter_"+imagetype+"_L2TNH3"
    pathFITS='C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Analysis Data/L2 FITS/'
elif Level=="L3":
    CH4suffix="-Jupiter_Map_L3PCld_S0"
    NH3suffix="-Jupiter_Map_L3fNH3_S0"
"""                  
         
                
dx,dy=0,-8
BRs={155:6,137:10,133:7,120:8}
N=1
for BR in BRs:
    axs1[0].arrow(BR-dx,BRs[BR]-dy,dx,dy,color='w',width=0.8,length_includes_head=True,
                  head_starts_at_zero=True,head_width=3,head_length=4,zorder=5)
    axs1[0].text(BR-dx,BRs[BR]-dy,N,color='w',fontsize=12,weight='bold',
                 horizontalalignment='center',zorder=5)
    axs1[1].arrow(BR-dx,BRs[BR]-dy,dx,dy,color='w',width=0.8,length_includes_head=True,
                  head_starts_at_zero=True,head_width=3,head_length=4,zorder=5)
    axs1[1].text(BR-dx,BRs[BR]-dy,N,color='w',fontsize=12,weight='bold',
                 horizontalalignment='center',zorder=5)
    axs2[0].arrow(BR-dx,BRs[BR]-dy,dx,dy,color='w',width=0.8,length_includes_head=True,
                  head_starts_at_zero=True,head_width=3,head_length=4,zorder=5)
    axs2[0].text(BR-dx,BRs[BR]-dy,N,color='w',fontsize=12,weight='bold',
                 horizontalalignment='center',zorder=5)
    axs2[1].arrow(BR-dx,BRs[BR]-dy,dx,dy,color='w',width=0.8,length_includes_head=True,
                  head_starts_at_zero=True,head_width=3,head_length=4,zorder=5)
    axs2[1].text(BR-dx,BRs[BR]-dy,N,color='w',fontsize=12,weight='bold',
                 horizontalalignment='center',zorder=5)
    axs3[0].arrow(BR-dx,BRs[BR]-dy,dx,dy,color='w',width=0.8,length_includes_head=True,
                  head_starts_at_zero=True,head_width=3,head_length=4,zorder=5)
    axs3[0].text(BR-dx,BRs[BR]-dy,N,color='w',fontsize=12,weight='bold',
                 horizontalalignment='center',zorder=5)
    N=N+1
    
fig1.savefig(pathout+'Fig5/Fig5a.png',dpi=300)
fig2.savefig(pathout+'Fig5/Fig5b.png',dpi=300)
fig3.savefig(pathout+'Fig5/Fig5c.png',dpi=300)

try:
    PCloudfile=sourcefiles["20221009UTa"]['CH4file'][0:17]+CH4suffix+\
            sourcefiles["20221009UTa"]['Metadata']['Variation']+".fits"
    variation=sourcefiles["20221009UTa"]['Metadata']['Variation']
except:
    PCloudfile=sourcefiles["20221009UTa"]['CH4file'][0:17]+CH4suffix+".fits"
    variation=""
try:
    fNH3file=sourcefiles["20221009UTa"]['NH3file'][0:17]+NH3suffix+\
            sourcefiles["20221009UTa"]['Metadata']['Variation']+".fits"
    variation=sourcefiles["20221009UTa"]['Metadata']['Variation']
except:
    fNH3file=sourcefiles["20221009UTa"]['NH3file'][0:17]+NH3suffix+".fits"
    variation=""

shutil.copyfile(pathFITS+fNH3file,pathout+"FITS/"+fNH3file)
shutil.copyfile(pathFITS+PCloudfile,pathout+"FITS/"+PCloudfile)
###############################################################################
# FIGURE 6
###############################################################################
fig1,axs1,fig2,axs2,fig3,axs3=MJA.Map_Jup_Atm_P3(obskey="20220919UTa",
                                                 imagetype='Map', 
                                                 Smoothing=False,
                                                 LatLims=[90,130],LonRng=20, 
                                                 CMpref=25,LonSys='2',
                                                 showbands=False,coef=[0.,0.],
                                                 subproj='paper')
fig1.savefig(pathout+'Fig6/Fig6b.png',dpi=300)

try:
    fNH3file=sourcefiles["20220919UTa"]['NH3file'][0:17]+NH3suffix+\
            sourcefiles["20220919UTa"]['Metadata']['Variation']+".fits"
    variation=sourcefiles["20220919UTa"]['Metadata']['Variation']
except:
    fNH3file=sourcefiles["20220919UTa"]['NH3file'][0:17]+NH3suffix+".fits"
    variation=""

shutil.copyfile(pathFITS+fNH3file,pathout+"FITS/"+fNH3file)
    
fig1,axs1,fig2,axs2,fig3,axs3=MJA.Map_Jup_Atm_P3(obskey="20220919UTb",
                                                 imagetype='Map', 
                                                 Smoothing=False,
                                                 LatLims=[90,130],LonRng=20, 
                                                 CMpref=25,LonSys='2',
                                                 showbands=False,coef=[0.,0.],
                                                 subproj='paper')
fig1.savefig(pathout+'Fig6/Fig6c.png',dpi=300)

try:
    fNH3file=sourcefiles["20220919UTb"]['NH3file'][0:17]+NH3suffix+\
            sourcefiles["20220919UTb"]['Metadata']['Variation']+".fits"
    variation=sourcefiles["20220919UTb"]['Metadata']['Variation']
except:
    fNH3file=sourcefiles["20220919UTb"]['NH3file'][0:17]+NH3suffix+".fits"
    variation=""
    
shutil.copyfile(pathFITS+fNH3file,pathout+"FITS/"+fNH3file)

###############################################################################
# FIGURE 7
###############################################################################
fig1,axs1,fig2,axs2,fig3,axs3=MJA.Map_Jup_Atm_P3(obskey="20221013UTa",
                                                 imagetype='Map', 
                                                 Smoothing=False,
                                                 LatLims=[90,130],LonRng=20, 
                                                 CMpref=25,LonSys='2',
                                                 showbands=False,coef=[0.,0.],
                                                 subproj='paper')
fig1.savefig(pathout+'Fig7/Fig7a.png',dpi=300)

try:
    fNH3file=sourcefiles["20221013UTa"]['NH3file'][0:17]+NH3suffix+\
            sourcefiles["20221013UTa"]['Metadata']['Variation']+".fits"
    variation=sourcefiles["20221013UTa"]['Metadata']['Variation']
except:
    fNH3file=sourcefiles["20221013UTa"]['NH3file'][0:17]+NH3suffix+".fits"
    variation=""
    
shutil.copyfile(pathFITS+fNH3file,pathout+"FITS/"+fNH3file)

fig1,axs1,fig2,axs2,fig3,axs3=MJA.Map_Jup_Atm_P3(obskey="20221020UTa",
                                                 imagetype='Map', 
                                                 Smoothing=False,
                                                 LatLims=[90,130],LonRng=20, 
                                                 CMpref=25,LonSys='2',
                                                 showbands=False,coef=[0.,0.],
                                                 subproj='paper')
fig1.savefig(pathout+'Fig7/Fig7b.png',dpi=300)

try:
    fNH3file=sourcefiles["20221020UTa"]['NH3file'][0:17]+NH3suffix+\
            sourcefiles["20221020UTa"]['Metadata']['Variation']+".fits"
    variation=sourcefiles["20221020UTa"]['Metadata']['Variation']
except:
    fNH3file=sourcefiles["20221020UTa"]['NH3file'][0:17]+NH3suffix+".fits"
    variation=""
    
shutil.copyfile(pathFITS+fNH3file,pathout+"FITS/"+fNH3file)
###############################################################################
# FIGURE 8
###############################################################################
fig1,MapAvg=GZM.grs_zoom_maps(year=2022)
fig1.savefig(pathout+'Fig8/Fig8a.png',dpi=300)
hdu = fits.PrimaryHDU(MapAvg.astype(np.float32))

hdul = fits.HDUList([hdu])

hdul[0].header['BITPIX']=-32
hdul[0].header['DATE-OBS']='2022 Avg 5'#.replace('_','T')
hdul[0].header['AUTHOR']='Hill, S. M.'
hdul[0].header['FILENAME']='2022_fNH3_Avg5_C0_Sys2_N0-S40_Lon005-045'

hdul[0].header['OBJECT']='Jupiter'

hdul[0].header['TELESCOP']='SCT'
hdul[0].header['INSTRUME']='ASI120MM'
hdul[0].header['BUNIT']='ppm'
hdul[0].header['CALIBRA']='VLT-Filter'
hdul[0].header['VERSION']=('TBD','TBD')
hdul[0].header['CTYPE1']=('Sys. 3 Longitude','deg')
hdul[0].header['CTYPE2']=('PG Latitude','deg')

try:
    os.remove(pathout+'FITS/'+hdul[0].header['FILENAME']+'.fits')
except: 
    print("file doesn't exist")
hdul.writeto(pathout+'FITS/'+hdul[0].header['FILENAME']+'.fits')
hdul.close()


fig1,MapAvg=GZM.grs_zoom_maps(year=2023)
fig1.savefig(pathout+'Fig8/Fig8b.png',dpi=300)
hdu = fits.PrimaryHDU(MapAvg.astype(np.float32))

hdul = fits.HDUList([hdu])

hdul[0].header['BITPIX']=-32
hdul[0].header['DATE-OBS']='2023 Avg 10'#.replace('_','T')
hdul[0].header['AUTHOR']='Hill, S. M.'
hdul[0].header['FILENAME']='2023_fNH3_Avg10_C0_Sys2_N0-S40_Lon025-065'

hdul[0].header['OBJECT']='Jupiter'

hdul[0].header['TELESCOP']='SCT'
hdul[0].header['INSTRUME']='ASI120MM'
hdul[0].header['BUNIT']='ppm'
hdul[0].header['CALIBRA']='VLT-Filter'
hdul[0].header['VERSION']=('TBD','TBD')
hdul[0].header['CTYPE1']=('Sys. 3 Longitude','deg')
hdul[0].header['CTYPE2']=('PG Latitude','deg')

try:
    os.remove(pathout+'FITS/'+hdul[0].header['FILENAME']+'.fits')
except: 
    print("file doesn't exist")
hdul.writeto(pathout+'FITS/'+hdul[0].header['FILENAME']+'.fits')
hdul.close()