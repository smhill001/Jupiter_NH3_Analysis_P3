# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 22:24:59 2023

@author: smhil
"""
from astropy.io import fits
import numpy as np
import pylab as pl

path="C:/Astronomy/Projects/Planets/Jupiter/Imaging Data/20220919UT/"
CH4hdulist=fits.open(path+'2022-09-19-0352_3-Jupiter-620CH4AbsMap.fits')
CH4hdulist.info()
CH4hdr=CH4hdulist[0].header
CH4data=CH4hdulist[0].data
CH4hdulist.close()

NH3hdulist=fits.open(path+'2022-09-19-0352_3-Jupiter-647NH3AbsMap.fits')
NH3hdulist.info()
NH3hdr=NH3hdulist[0].header
NH3data=NH3hdulist[0].data
NH3hdulist.close()

CH4_tau=-np.log(CH4data*0.893)
NH3_tau=-np.log(NH3data*0.962)

CH4_Ncol=1000*CH4_tau/0.454
NH3_Ncol=1000*NH3_tau/3.129

fNH3=(1.81e-3)*NH3_Ncol/CH4_Ncol

pl.imshow(fNH3)
