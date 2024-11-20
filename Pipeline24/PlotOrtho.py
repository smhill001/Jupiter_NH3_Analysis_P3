# -*- coding: utf-8 -*-
"""
Created on Sat Sep 21 22:13:39 2024

@author: smhil
"""

import sys
drive='c:'
sys.path.append(drive+'/Astronomy/Python Play')
sys.path.append(drive+'/Astronomy/Python Play/Util_P3')
sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy')
sys.path.append('./Services')

import pylab as pl
from astropy.io import fits
sys.path.append('./Services')
import get_obs_list as getlist

import get_WINJupos_ephem as WJ_ephem
import numpy as np
import planetmapper


pth="C:/Astronomy/Projects/Planets/Saturn/Imaging Data/20240919UT/"
fn="2024-09-19-0327_3-Saturn_685NIR-FlatStack600-WV2x20-Aligned_map.fits"

test=fits.open(pth+fn)
test.info()
testhdr=test[0].header
testdata=test[0].data
test.close()

#pl.imshow(testdata[0,:,:])

planetmapper.BodyXY.plot_map(testdata[0,:,:],projection="orthographic")
