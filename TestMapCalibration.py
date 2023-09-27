# -*- coding: utf-8 -*-
"""
Created on Sun Sep 17 20:54:46 2023

@author: smhil
"""
import os
import numpy as np
from astropy.io import fits


pathin="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Analysis Data/L2 FITS/"
fn="2023-09-06-1203_2-Jupiter_620CH4AbsMap.fits"
CH4hdulist=fits.open(pathin+fn)
CH4hdulist.info()
CH4hdr=CH4hdulist[0].header
CH4data=CH4hdulist[0].data
CH4hdulist.close()

tau=-np.log(CH4data)


hdu = fits.PrimaryHDU(tau.astype(np.float32))
hdul = fits.HDUList([hdu])
try:
    os.remove(pathin+'CH4tau.fits')
except: 
    print("file doesn't exist")
hdul.writeto(pathin+'CH4tau.fits')
hdul.close()

