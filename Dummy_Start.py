#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 09:27:47 2026

@author: shill
"""


#def dummy():
from scipy.signal import convolve2d
import sys
sys.path.append("/mnt/data/git_repos/Jupiter_NH3_Analysis_P3/Services/")
sys.path.append("/mnt/data/git_repos/Jupiter_NH3_Analysis_P3/Maps/")
sys.path.append("/mnt/data/git_repos/Jupiter_NH3_Analysis_P3/HST/")
sys.path.append("/mnt/data/git_repos/Jupiter_NH3_Analysis_P3/Winds/")
sys.path.append("/mnt/data/git_repos/Visualization-and-Analysis/")
sys.path.append("/mnt/data/git_repos/Data-Management-and-Access/processes/")
import get_spice_ephem
import convert_system3_to_I_II_spice
import PlanetMapper_Spice_Furnish
import read_HST_GO
import make_HST_fits_script
