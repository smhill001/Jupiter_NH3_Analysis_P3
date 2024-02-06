# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 15:34:15 2023

@author: smhil
"""
def get_karkoschka_data(Type='Jupiter'):
    import sys
    sys.path.append('c:/Astronomy/Python Play')
    sys.path.append('c:/Astronomy/Python Play/Util_P3')
    sys.path.append('c:/Astronomy/Python Play/SPLibraries_P3')
    sys.path.append('c:/Astronomy/Python Play/SpectroPhotometry/Spectroscopy_P3')
    import numpy as np
    import GeneralSpecUtils_P3 as GSU
    projpath="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/"
    
    ###############################################################################
    # LOAD JOVIAN DISK-INTEGRATEDALBEDO DATA FROM KARKOSCHKA, 1994 (DATA FROM 1993)
    ###############################################################################
    col={"Jupiter":3,"CH4":2}
    CH4path="Molecular Absorption/Cross Sections/CH4/"
    print('hi')
    print(projpath)
    print(projpath+CH4path+"1993.tab.txt")
    Jupiter_Karkoschka1993 = np.fromfile(file=projpath+CH4path+"1993.tab.txt",
                                         dtype=float, count=-1, sep=" ")
    kark1993nrows=int(Jupiter_Karkoschka1993.size/8)
    Jupiter_Karkoschka1993=np.reshape(Jupiter_Karkoschka1993,[kark1993nrows,8])
    
    KarkRef1993=np.zeros((kark1993nrows,2))
    KarkRef1993[:,0]=Jupiter_Karkoschka1993[:,0]
    KarkRef1993[:,1]=Jupiter_Karkoschka1993[:,col[Type]] #Albedo
    WaveGrid,SignalonGrid=GSU.uniform_wave_grid(KarkRef1993[:,0],
                                                KarkRef1993[:,1],
                                                Extend=False,Fine=False)
    
    KarkRef1993_regrid=np.zeros((WaveGrid.size,2))
    KarkRef1993_regrid[:,0]=WaveGrid
    KarkRef1993_regrid[:,1]=SignalonGrid
    
    return(KarkRef1993_regrid)

