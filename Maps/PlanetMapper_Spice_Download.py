# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 18:00:39 2023

@author: smhil
"""
import sys
#sys.path.append('C:/Astronomy/Projects/SAS 2021 Ammonia/Visualization-and-Analysis/')

import socket
import planetmapper as pm
from planetmapper.kernel_downloader import download_urls
hostname = socket.gethostname()
print(hostname)
from config_VA import spice_path

pm.set_kernel_path(spice_path[hostname])
#pm.set_kernel_path("/mnt/data/spice_kernels/")
download_urls('https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/')
download_urls('https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/')
#download_urls('https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/earth_200101_990825_predict.bpc')
download_urls('https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430.bsp')
download_urls('https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/jup365.bsp')

# Need to make sure that I have the right time intervals for JUNO kernals
download_urls('https://naif.jpl.nasa.gov/pub/naif/JUNO/kernels/spk/spk_ref_251119_281001_260129.bsp')
download_urls('https://naif.jpl.nasa.gov/pub/naif/JUNO/kernels/spk/spk_ref_250724_281001_250724.bsp')
download_urls('https://naif.jpl.nasa.gov/pub/naif/JUNO/kernels/spk/spk_ref_250526_281001_250505.bsp')
download_urls('https://naif.jpl.nasa.gov/pub/naif/JUNO/kernels/spk/spk_ref_231110_251016_231110.bsp')
#download_urls('https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/sat459.bsp')
