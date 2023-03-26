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
import GeneralSpecUtils_P3 as GSU
import scipy


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

filterwavelength=['620','632','647','656','658','672','730','889','940']
filterdata={'620':{'transfile':'620CH4/620CH4_Transmission.txt',
                   'filtname':'620CH4','filtwdth':10.},
             '632':{'transfile':'632OI/632OI_Transmission.txt',
                    'filtname':'632OI','filtwdth':10.},
             '647':{'transfile':'647CNT/647CNT_Transmission.txt',
                    'filtname':'647NH3','filtwdth':10.},
             '656':{'transfile':'656HIA/656HIA_Transmission.txt',
                    'filtname':'656HIA','filtwdth':10.},
             '658':{'transfile':'658NII/658NII_Transmission.txt',
                    'filtname':'658NII','filtwdth':5.},
             '672':{'transfile':'672SII/672SII_Transmission.txt',
                    'filtname':'672SII','filtwdth':10.},
             '730':{'transfile':'730OII/730OII_Transmission.txt',
                    'filtname':'730OII','filtwdth':10.},
             '889':{'transfile':'889CH4/889CH4_Transmission.txt',
                    'filtname':'889CH4','filtwdth':10.},
             '940':{'transfile':'940NIR/940NIR_Transmission.txt',
                    'filtname':'940NIR','filtwdth':10.}}

filtpath='C:/Astronomy/Projects/Techniques/InstrumentPerformance-P3/Filters/'
filterdata['620']['FiltTrans']=np.loadtxt(filtpath+filterdata['620']['transfile'],usecols=range(2))

#F620_wvs=[617.5,622.5]
F620_wvs=[615.,625.]
F620_idx=[np.argmin(abs(F620_wvs[0]-wavelength)),np.argmin(abs(F620_wvs[1]-wavelength))]
MUSE620=np.mean(MUSEdata[F620_idx[0]:F620_idx[1],:,:],axis=0)
MUSE620scaled=np.nan_to_num(65535.*MUSE620*100.)
MUSE620scaled[MUSE620scaled<=0.]=0.0
MUSE620abs16bit = np.flipud(MUSE620scaled.astype(np.uint16))
imwrite(path+fn+'MUSE620-10nm.png', MUSE620abs16bit)

#F632_wvs=[629.5,634.5]
F632_wvs=[627.,637.]
F632_idx=[np.argmin(abs(F632_wvs[0]-wavelength)),np.argmin(abs(F632_wvs[1]-wavelength))]
MUSE632=np.mean(MUSEdata[F632_idx[0]:F632_idx[1],:,:],axis=0)
MUSE632scaled=np.nan_to_num(65535.*MUSE632*100.)
MUSE632scaled[MUSE632scaled<=0.]=0.0
MUSE632abs16bit = np.flipud(MUSE632scaled.astype(np.uint16))
imwrite(path+fn+'MUSE632-10nm.png', MUSE632abs16bit)

#F647_wvs=[644.5,649.5]
F647_wvs=[642.,652.]
F647_idx=[np.argmin(abs(F647_wvs[0]-wavelength)),np.argmin(abs(F647_wvs[1]-wavelength))]
MUSE647=np.mean(MUSEdata[F647_idx[0]:F647_idx[1],:,:],axis=0)
MUSE647scaled=np.nan_to_num(65535.*MUSE647*100.)
MUSE647scaled[MUSE647scaled<=0.]=0.0
MUSE647abs16bit = np.flipud(MUSE647scaled.astype(np.uint16))
imwrite(path+fn+'MUSE647-10nm.png', MUSE647abs16bit)

F656_wvsBLU=[651.0,654.0]
#F656_wvs=[653.5,658.5]
#F656_wvs=[651.0,661.0]
F656_idxBLU=[np.argmin(abs(F656_wvsBLU[0]-wavelength)),np.argmin(abs(F656_wvsBLU[1]-wavelength))]
MUSE656BLU=np.mean(MUSEdata[F656_idxBLU[0]:F656_idxBLU[1],:,:],axis=0)
F656_wvsRED=[658.0,661.0]
F656_idxRED=[np.argmin(abs(F656_wvsRED[0]-wavelength)),np.argmin(abs(F656_wvsRED[1]-wavelength))]
MUSE656RED=np.mean(MUSEdata[F656_idxRED[0]:F656_idxRED[1],:,:],axis=0)
MUSE656=(MUSE656BLU+MUSE656RED)/2.

MUSE656scaled=np.nan_to_num(65535.*MUSE656*100.)
MUSE656scaled[MUSE656scaled<=0.]=0.0
MUSE656abs16bit = np.flipud(MUSE656scaled.astype(np.uint16))
imwrite(path+fn+'MUSE656-6nmSplit.png', MUSE656abs16bit)

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

MUSESpec=np.nanmean(MUSEdata[:,:,:],axis=(1,2))
print("MUSESpec.shape=",MUSESpec.shape,MUSESpec)
WaveGrid,SignalonGrid=GSU.uniform_wave_grid(wavelength,MUSESpec,Extend=False,Fine=False)
MuseSpecGrid=np.column_stack((WaveGrid,SignalonGrid))
print("MuseSpecGrid.shape=",MuseSpecGrid.shape)
RefPath="c:/Astronomy/Python Play/SPLibraries/SpectralReferenceFiles/ReferenceLibrary/"
G2V=np.loadtxt(RefPath+"g2v.dat", dtype=float, usecols=(0,1))
print("G2V.shape=",G2V.shape)
AlbedoSpec=GSU.SpectrumMath(MuseSpecGrid,G2V,"Divide")

fig2,axs2=pl.subplots(1,1,figsize=(7.0,4.5), dpi=150, facecolor="white",
                    sharey=True,sharex=True)      
axs2.plot(wavelength,MUSESpec)
axs2.plot(WaveGrid,SignalonGrid)
axs2.plot(AlbedoSpec[:,0],AlbedoSpec[:,1])

F620_idx=[np.argmin(abs(F620_wvs[0]-WaveGrid)),np.argmin(abs(F620_wvs[1]-WaveGrid))]
MUSE620Albedo=np.mean(AlbedoSpec[F620_idx[0]:F620_idx[1],1],axis=0)
print(MUSE620Albedo)
F632_idx=[np.argmin(abs(F632_wvs[0]-WaveGrid)),np.argmin(abs(F632_wvs[1]-WaveGrid))]
MUSE632Albedo=np.mean(AlbedoSpec[F632_idx[0]:F632_idx[1],1],axis=0)
print(MUSE632Albedo)
F647_idx=[np.argmin(abs(F647_wvs[0]-WaveGrid)),np.argmin(abs(F647_wvs[1]-WaveGrid))]
MUSE647Albedo=np.mean(AlbedoSpec[F647_idx[0]:F647_idx[1],1],axis=0)
print(MUSE647Albedo)

F656_idxRED=[np.argmin(abs(F656_wvsRED[0]-WaveGrid)),np.argmin(abs(F656_wvsRED[1]-WaveGrid))]
MUSE656AlbedoRED=np.mean(AlbedoSpec[F656_idxRED[0]:F656_idxRED[1],1],axis=0)
F656_idxBLU=[np.argmin(abs(F656_wvsBLU[0]-WaveGrid)),np.argmin(abs(F656_wvsBLU[1]-WaveGrid))]
MUSE656AlbedoBLU=np.mean(AlbedoSpec[F656_idxBLU[0]:F656_idxBLU[1],1],axis=0)
MUSE656Albedo=(MUSE656AlbedoRED+MUSE656AlbedoBLU)/2.
print(MUSE656Albedo)
