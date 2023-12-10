# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 18:00:39 2023

@author: smhil
"""

from planetmapper.kernel_downloader import download_urls
download_urls('https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430.bsp')
download_urls('https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/jup365.bsp')
