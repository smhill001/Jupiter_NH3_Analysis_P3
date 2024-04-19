# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 11:17:50 2023

@author: smhil
"""

def get_albedo_continua_crossections(x0,x1,xtks,y0,y1,ytks,
                                     LutzPlot=True,Crossect=True,plotMUSE=True,
                                     ModList=['Spline1','Spline2','Piecewise1',
                                              'Piecewise2','Linear 2pt.']):
    """
    PURPOSE: 
        1) Retrieves a reference albedo for Jupiter (Karkoschka, 1994)
        2) Computes the continuum models provided in the model list
        3) Retrieves gas absorption cross sections for CH4 and NH3
        4) Creates a composite plot including the reference albedo, the
           continuum models, and the gas absorption cross sections
    
    Calls:
        ->get_karkoschka_data
        ->get_ammonia_data
        ->get_continuum_model
    Called by:
        JupiterFilterPerformance->
        AmmoniaTest_P3.py->
        NH3_Filter_Library_P3.py/SpectralModeling->
    Example:
        get_albedo_continua_crossections(600,680,9,0,0.7,8,LutzPlot=True,Crossect=True, ModList=[])
    """
    import sys
    sys.path.append('c:/Astronomy/Python Play')
    sys.path.append('c:/Astronomy/Python Play/Util_P3')
    sys.path.append('c:/Astronomy/Python Play/SPLibraries_P3')
    sys.path.append('c:/Astronomy/Python Play/SpectroPhotometry/Spectroscopy_P3')
    import matplotlib.pyplot as pl
    import numpy as np
    import sys
    sys.path.append('./Services')
    import get_karkoschka_data
    import get_ammonia_data
    import get_continuum_model
    import get_MUSE_albedos as GMA
    import get_filter_base_dict as GFBD

    projpath="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/"

    # LOAD JOVIAN DISK-INTEGRATEDALBEDO DATA FROM KARKOSCHKA, 1994 (DATA FROM 1993)
    Albedo=get_karkoschka_data.get_karkoschka_data(Type='Jupiter')
    if plotMUSE==True:
        MUSE_20220730,MUSE_20220919=GMA.get_MUSE_albedos()
    CH4=get_karkoschka_data.get_karkoschka_data(Type='CH4')
    WaveGrid=CH4[:,0]
    NH3=get_ammonia_data.get_ammonia_data(Source="Irwin")
    NH3_Lutz_Owen_1980=get_ammonia_data.get_ammonia_data(Source="Lutz")
    
    # COMPUTE AND LOAD CONTINUUM MODELS
    Continua=get_continuum_model.get_continuum_model(Albedo,
                                                     ModelList=ModList)
    
    #SET UP PLOT
    fig_molecules,ax_molecules=pl.subplots(3,1,figsize=(5.0,6.0), dpi=150, facecolor="white",
                          sharex=True)
    #fig_molecules.suptitle("Albedo, Continuum Models, and Molecular Crossections")
    ax_molecules[0].set_xlim(x0,x1)
    ax_molecules[0].set_xticks(np.linspace(x0,x1,xtks, endpoint=True))
    ax_molecules[0].set_ylim(y0,y1)
    ax_molecules[0].set_yticks(np.linspace(y0,y1,ytks, endpoint=True))
    ax_molecules[0].grid(linewidth=0.2)
    ax_molecules[0].tick_params(axis='both', which='major', labelsize=8)
    ax_molecules[0].set_ylabel("Albedo",fontsize=8)
    
    # PLOT ALBEDO AND CONTINUUM MODELS IN THE LIST
    ax_molecules[0].plot(Albedo[:,0],Albedo[:,1],label='Karkoschka, 1994',linewidth=1.0,color='C0')
    if plotMUSE==True:
        ax_molecules[0].plot(MUSE_20220730[:,0],MUSE_20220730[:,1]*0.70,label='MUSE 2022-07-30',linewidth=1.0,color='C1')
        ax_molecules[0].plot(MUSE_20220919[:,0],MUSE_20220919[:,1]*0.63,label='MUSE 2022-09-19',linewidth=1.0,color='C2')

    for ModelName in ModList:
        ax_molecules.plot(Continua[ModelName]['WaveGrid'],
                          Continua[ModelName]['Albedo'],label=ModelName,
                          linewidth=0.5)
    ax_molecules[0].legend(fontsize=6, loc="upper right")


    path,filterdict=GFBD.get_filter_base_dict("SCT")
    T620=filterdict['620']['FiltTrans']
    T632=filterdict['632']['FiltTrans']
    T647=filterdict['647']['FiltTrans']
    T656=filterdict['656']['FiltTrans']

    ax_molecules[1].set_ylim(0,1)
    ax_molecules[1].tick_params(axis='y', which='major', labelsize=8)
    ax_molecules[1].plot(T620[:,0],T620[:,1]/np.nanmax(T620[:,1]),label='620 nm',linewidth=1.0)
    ax_molecules[1].plot(T632[:,0],T632[:,1]/np.nanmax(T632[:,1]),label='632 nm',linewidth=1.0)
    ax_molecules[1].plot(T647[:,0],T647[:,1]/np.nanmax(T647[:,1]),label='647 nm',linewidth=1.0)
    ax_molecules[1].plot(T656[:,0],T656[:,1]/np.nanmax(T656[:,1]),label='656 nm',linewidth=1.0)
    ax_molecules[1].set_ylim(0,1.)
    ax_molecules[1].set_ylabel("Normalized Transmission",fontsize=8)
    ax_molecules[1].legend(fontsize=6,loc="upper right")
    ax_molecules[1].grid(linewidth=0.2)



    #OVERPLOT MOLECULAR ABSORPTION CROSSECTIONS   
    #if Crossect:
    ax_molecules[2].ticklabel_format(axis='y')
    ax_molecules[2].tick_params(axis='both', which='major', labelsize=8)
    ax_molecules[2].set_yscale('log')
    ax_molecules[2].set_ylim(1e-3,1e2)
    ax_molecules[2].set_ylabel("Absorption Coefficient 1/(km-atm)", fontsize=8)#,color="green")
    
    ax_molecules[2].plot(CH4[:,0],CH4[:,1],label='CH4 (Karkoschka, 1994) ',linewidth=1.0,color='C0')
    ax_molecules[2].plot(NH3[:,0],NH3[:,1],label='NH3 (ExoMol) ',linewidth=1.0,color='C1')
    
    if LutzPlot:
        ax_molecules[2].plot(NH3_Lutz_Owen_1980[:,0],NH3_Lutz_Owen_1980[:,1],
                             label='NH3 (Lutz & Owen, 1980) ',linewidth=0.5,color='C3')
    ax_molecules[2].legend(fontsize=6, loc="upper right")

    ax_molecules[2].grid(linewidth=0.2)
    ax_molecules[2].set_xlabel("Wavelength (nm)",fontsize=8)
    
    
    fig_molecules.subplots_adjust(left=0.12, right=0.95, top=0.98, bottom=0.06,
                                  wspace=0.1,hspace=0.1)
          

    outpath="Molecular Absorption/data output/"

    fig_molecules.savefig(projpath+outpath+\
                          'Albedo_Continua_Crossections.png',dpi=320)
    
    return(Albedo,Continua,CH4,NH3,NH3_Lutz_Owen_1980)