# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 13:55:35 2023

    First break-off code from AmmoniaTest.py
        
        COMPUTE and PLOT FILTER TRANSMISSIONS CONVOLVED WITH DISK-INTEGRATED 
            ALBEDO AND CONTINUUM
        COMPUTE and PLOT K_eff, l_eff, AND WEIGHTING FUNCTION CALCULATIONS
        COMPUTE and PLOT **RAYLEIGH CANCELING** AND GAS ABSORPTION ONLY
            WEIGHTING FUNCTIONS
        COMPUTE ESTIMATED JUPITER ABSORPTION USING GALILEAN MOONS


@author: Steven Hill
"""
import sys
sys.path.append('c:/Astronomy/Python Play')
sys.path.append('c:/Astronomy/Python Play/Util_P3')
sys.path.append('c:/Astronomy/Python Play/SPLibraries_P3')
sys.path.append('c:/Astronomy/Python Play/SpectroPhotometry/Spectroscopy_P3')
import matplotlib.pyplot as pl
import numpy as np
from copy import deepcopy
#from scipy import interpolate
import GeneralSpecUtils_P3 as GSU
#from numpy import genfromtxt
import NH3_Filter_Library_P3 as NFL
#import AmmoniaTest_P3 as AT

###############################################################################
# LOAD JOVIAN DISK-INTEGRATEDALBEDO DATA FROM KARKOSCHKA, 1994 (DATA FROM 1993)
###############################################################################
#Plot Layout Configuration
ContinuumModel='All'
x0,x1,xtks=600.,680.,9
y0,y1,ytks=0.0,0.7,8
Albedo,Continua,CH4,NH3=NFL.Get_Albedo_and_Absorption(x0,x1,xtks,y0,y1,ytks,
                                                              Crossect=True)
########## END OF FIRST FUNCTION AND PLOT ##############

###############################################################################
# RETRIEVE FILTER TRANSMISSIONS FROM MASTER FILTER LIBRARY
#   !!!! NEED TO ADD TELESCOPE THROUGHPUT CALIBRATION HERE
###############################################################################
#filterwavelength=['620','632','647','656','658','672','730','889','940']
filterwavelength=['620','632','647','656']

###############################################################################
# PLOT FILTER TRANSMISSIONS CONVOLVED WITH DISK-INTEGRATED ALBEDO AND CONTINUUM
###############################################################################
Telescopes=["SCT","VLT"]
Models=['Spline1','Spline2','Piecewise1','Piecewise2','Linear 2pt.']
pth='c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/'
filtereffectivedata = open(pth+'JupiterFilterPerformance.csv', 'w')
tmp="Wavelength (nm),Filter Name,Telescope, Model,Type,Transmission,Tau\n"
filtereffectivedata.write(tmp)

for Tele in Telescopes:
    for Model in Models:
    
        Continuum_Albedo=np.zeros((Continua[Model]['WaveGrid'].size,2))
        Continuum_Albedo[:,0]=Continua[Model]['WaveGrid']
        Continuum_Albedo[:,1]=Continua[Model]['Albedo']
        
        filterdata,axsFilt=NFL.compute_filter_Jupiter_transmissions(x0,x1,xtks,y0,y1,ytks,
                                                             filterwavelength,
                                                             Albedo,Continuum_Albedo,
                                                             Model,
                                                             Telescope=Tele)
        
        
        for filtr in filterwavelength:
            print(filtr)
            tmp=filtr+","+filterdata[filtr]['filtname']+","+Tele+','+Model+','\
                +str(filterdata[filtr]['TransInt'])+","\
                    +str(filterdata[filtr]['Tau_Albedo'])+"\n"
            filtereffectivedata.write(tmp)
filtereffectivedata.close()
          
