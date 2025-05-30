# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 15:34:15 2023

@author: smhil
"""
def get_MUSE_albedos():
    import sys
    sys.path.append('c:/Astronomy/Python Play')
    sys.path.append('c:/Astronomy/Python Play/Util_P3')
    sys.path.append('c:/Astronomy/Python Play/SPLibraries_P3')
    sys.path.append('c:/Astronomy/Python Play/SpectroPhotometry/Spectroscopy_P3')
    import numpy as np
    import GeneralSpecUtils_P3 as GSU
    projpath="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/"
    MUSEpath="MUSE/data output/"
    ###############################################################################
    # LOAD JOVIAN DISK-INTEGRATEDALBEDO DATA FROM MUSE OBSERVATIONS)
    ###############################################################################
    
    temp_20220730 = np.fromfile(file=projpath+MUSEpath+"2022-07-30-0729_8_MUSE_Albedo.txt",
                                         dtype=float, count=-1, sep=" ")
    nrows=int(temp_20220730.size/2)
    MUSE_20220730=np.reshape(temp_20220730,[nrows,2])

    temp_20220919 = np.fromfile(file=projpath+MUSEpath+"2022-09-19-0352_3_MUSE_Albedo.txt",
                                         dtype=float, count=-1, sep=" ")
    nrows=int(temp_20220919.size/2)
    MUSE_20220919=np.reshape(temp_20220919,[nrows,2])
        
    temp_20230812 = np.fromfile(file=projpath+MUSEpath+"2023-08-12-0948_2_MUSE_Albedo.txt",
                                         dtype=float, count=-1, sep=" ")
    nrows=int(temp_20230812.size/2)
    MUSE_20230812=np.reshape(temp_20230812,[nrows,2])

    temp_20230815 = np.fromfile(file=projpath+MUSEpath+"2023-08-15-0734_0_MUSE_Albedo.txt",
                                         dtype=float, count=-1, sep=" ")
    nrows=int(temp_20230815.size/2)
    MUSE_20230815=np.reshape(temp_20230815,[nrows,2])

    temp_20220816 = np.fromfile(file=projpath+MUSEpath+"2023-08-16-0743_5_MUSE_Albedo.txt",
                                         dtype=float, count=-1, sep=" ")
    nrows=int(temp_20220816.size/2)
    MUSE_20230816=np.reshape(temp_20220816,[nrows,2])

    temp_20230923 = np.fromfile(file=projpath+MUSEpath+"2023-09-23-0559_2_MUSE_Albedo.txt",
                                         dtype=float, count=-1, sep=" ")
    nrows=int(temp_20230923.size/2)
    MUSE_20230923=np.reshape(temp_20230923,[nrows,2])

    return(MUSE_20220730,MUSE_20220919,MUSE_20230812,MUSE_20230815,
           MUSE_20230816,MUSE_20230923)

