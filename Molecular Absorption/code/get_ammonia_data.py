# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 15:53:50 2023

@author: smhil
"""
def get_ammonia_data(Source="Irwin"):
    import sys
    sys.path.append('c:/Astronomy/Python Play')
    sys.path.append('c:/Astronomy/Python Play/Util_P3')
    sys.path.append('c:/Astronomy/Python Play/SPLibraries_P3')
    sys.path.append('c:/Astronomy/Python Play/SpectroPhotometry/Spectroscopy_P3')
    import numpy as np
    import GeneralSpecUtils_P3 as GSU
    from numpy import genfromtxt
    projpath="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/"

    fn='Lutz&Owen1980_Figure5_AmmoniaCrossSection.csv'
    pathin="Molecular Absorption/Cross Sections/NH3/"
    NH3_Lutz_Owen_1980 = np.array(genfromtxt(projpath+pathin+fn, delimiter=','))
    
    fn='Exomol_NH3.csv'
    NH3_Exomol = np.array(genfromtxt(projpath+pathin+fn, delimiter=','))
    #print("***",NH3_Exomol.shape)
    wv,csec=GSU.uniform_wave_grid(NH3_Exomol[:,0]*1000.,NH3_Exomol[:,1],Extend=False,Fine=False)
    NH3_Exomol_regrid=np.transpose(np.array([wv,csec]))
    
    if Source=="Irwin":
        NH3=NH3_Exomol_regrid
    elif Source=="Lutz":
        NH3=NH3_Lutz_Owen_1980
        
    return(NH3)
