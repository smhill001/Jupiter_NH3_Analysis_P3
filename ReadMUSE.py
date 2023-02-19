# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 14:21:05 2023

@author: smhil
"""

import sys
drive='c:'
sys.path.append(drive+'/Astronomy/Python Play')
sys.path.append(drive+'/Astronomy/Python Play/Util_P3')
sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy')

import os
from matplotlib.pyplot import imread
import pylab as pl
import numpy as np
from imageio import imwrite
from numpy import inf
from re import search
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
from astropy.io import fits
import RetrievalLibrary as RL


path='c:/Astronomy/Projects/SAS 2021 Ammonia/VLT MUSE/'
MUSEfile='2022-09-19_obs18_proj.fits'
MUSEhdulist=fits.open(path+MUSEfile)
MUSEhdulist.info()
MUSEhdr=MUSEhdulist[0].header
MUSEdata=MUSEhdulist[1].data
MUSEhdulist.close()

dateobs=MUSEhdr['DATE-OBS']
fn=dateobs[0:10]+'-'+dateobs[11:13]+dateobs[14:16]+'_'+str(float(dateobs[17:19])/60.)[2]+'-Jupiter-'

print(dateobs)
print(fn)
#print(MUSEhdr)
wavelength=np.linspace(470.,935.13,3722)

F620_wvs=[615.,625.]
F620_idx=[np.argmin(abs(F620_wvs[0]-wavelength)),np.argmin(abs(F620_wvs[1]-wavelength))]
MUSE620=np.mean(MUSEdata[F620_idx[0]:F620_idx[1],:,:],axis=0)
MUSE620scaled=np.nan_to_num(65535.*MUSE620*100.)
MUSE620scaled[MUSE620scaled<=0.]=0.0
MUSE620abs16bit = np.flipud(MUSE620scaled.astype(np.uint16))
imwrite(path+fn+'MUSE620.png', MUSE620abs16bit)

F632_wvs=[627.,637.]
F632_idx=[np.argmin(abs(F632_wvs[0]-wavelength)),np.argmin(abs(F632_wvs[1]-wavelength))]
MUSE632=np.mean(MUSEdata[F632_idx[0]:F632_idx[1],:,:],axis=0)
MUSE632scaled=np.nan_to_num(65535.*MUSE632*100.)
MUSE632scaled[MUSE632scaled<=0.]=0.0
MUSE632abs16bit = np.flipud(MUSE632scaled.astype(np.uint16))
imwrite(path+fn+'MUSE632.png', MUSE632abs16bit)

F647_wvs=[642.,652.]
F647_idx=[np.argmin(abs(F647_wvs[0]-wavelength)),np.argmin(abs(F647_wvs[1]-wavelength))]
MUSE647=np.mean(MUSEdata[F647_idx[0]:F647_idx[1],:,:],axis=0)
MUSE647scaled=np.nan_to_num(65535.*MUSE647*100.)
MUSE647scaled[MUSE647scaled<=0.]=0.0
MUSE647abs16bit = np.flipud(MUSE647scaled.astype(np.uint16))
imwrite(path+fn+'MUSE647.png', MUSE647abs16bit)

F656_wvs=[651.,661.]
F656_idx=[np.argmin(abs(F656_wvs[0]-wavelength)),np.argmin(abs(F656_wvs[1]-wavelength))]
MUSE656=np.mean(MUSEdata[F656_idx[0]:F656_idx[1],:,:],axis=0)
MUSE656scaled=np.nan_to_num(65535.*MUSE656*100.)
MUSE656scaled[MUSE656scaled<=0.]=0.0
MUSE656abs16bit = np.flipud(MUSE656scaled.astype(np.uint16))
imwrite(path+fn+'MUSE656.png', MUSE656abs16bit)

fig,axs=pl.subplots(2,2,figsize=(7.0,4.5), dpi=150, facecolor="white",
                    sharey=True,sharex=True)      
#fig.suptitle(NH3time.replace("_"," ")+", CM2="+str(int(CM2))+", Calibration = "+CalModel+", Data "+smthtitle,x=0.5,ha='center',color='k')


axs[0,0].imshow(MUSE620,"gray",origin='lower')
axs[0,1].imshow(MUSE632,"gray",origin='lower')
axs[1,0].imshow(MUSE647,"gray",origin='lower')
axs[1,1].imshow(MUSE656,"gray",origin='lower')

F450_wvs=[471.,480.]
F450_idx=[np.argmin(abs(F450_wvs[0]-wavelength)),np.argmin(abs(F450_wvs[1]-wavelength))]
MUSE450=np.mean(MUSEdata[F450_idx[0]:F450_idx[1],:,:],axis=0)
MUSE450scaled=np.nan_to_num(65535.*MUSE450*100.)
MUSE450scaled[MUSE450scaled<=0.]=0.0
MUSE450abs16bit = np.flipud(MUSE450scaled.astype(np.uint16))
imwrite(path+fn+'MUSE450.png', MUSE450abs16bit)

F550_wvs=[530.,580.]
F550_idx=[np.argmin(abs(F550_wvs[0]-wavelength)),np.argmin(abs(F550_wvs[1]-wavelength))]
MUSE550=np.mean(MUSEdata[F550_idx[0]:F550_idx[1],:,:],axis=0)
MUSE550scaled=np.nan_to_num(65535.*MUSE550*100.)
MUSE550scaled[MUSE550scaled<=0.]=0.0
MUSE550abs16bit = np.flipud(MUSE550scaled.astype(np.uint16))
imwrite(path+fn+'MUSE550.png', MUSE550abs16bit)

F650_wvs=[630.,680.]
F650_idx=[np.argmin(abs(F650_wvs[0]-wavelength)),np.argmin(abs(F650_wvs[1]-wavelength))]
MUSE650=np.mean(MUSEdata[F650_idx[0]:F650_idx[1],:,:],axis=0)
MUSE650scaled=np.nan_to_num(65535.*MUSE650*100.)
MUSE650scaled[MUSE450scaled<=0.]=0.0
MUSE650abs16bit = np.flipud(MUSE650scaled.astype(np.uint16))
imwrite(path+fn+'MUSE650.png', MUSE650abs16bit)

fig1,axs1=pl.subplots(2,2,figsize=(7.0,4.5), dpi=150, facecolor="white",
                    sharey=True,sharex=True)      
#fig.suptitle(NH3time.replace("_"," ")+", CM2="+str(int(CM2))+", Calibration = "+CalModel+", Data "+smthtitle,x=0.5,ha='center',color='k')


axs1[0,0].imshow(MUSE450,"gray",origin='lower')
axs1[0,1].imshow(MUSE550,"gray",origin='lower')
axs1[1,0].imshow(MUSE650,"gray",origin='lower')
