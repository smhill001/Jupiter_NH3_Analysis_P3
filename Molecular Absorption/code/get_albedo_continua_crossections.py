# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 11:17:50 2023

@author: smhil
"""

def get_albedo_continua_crossections(x0,x1,xtks,y0,y1,ytks,
                                     LutzPlot=True,Crossect=True,
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
    projpath="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/"

    # LOAD JOVIAN DISK-INTEGRATEDALBEDO DATA FROM KARKOSCHKA, 1994 (DATA FROM 1993)
    Albedo=get_karkoschka_data.get_karkoschka_data(Type='Jupiter')
    CH4=get_karkoschka_data.get_karkoschka_data(Type='CH4')
    WaveGrid=CH4[:,0]
    NH3=get_ammonia_data.get_ammonia_data(Source="Irwin")
    NH3_Lutz_Owen_1980=get_ammonia_data.get_ammonia_data(Source="Lutz")
    
    # COMPUTE AND LOAD CONTINUUM MODELS
    Continua=get_continuum_model.get_continuum_model(Albedo,
                                                     ModelList=ModList)
    
    #SET UP PLOT
    fig_molecules,ax_molecules=pl.subplots(1,1,figsize=(6.0,3.5), dpi=150, facecolor="white",
                          sharex=True)
    ax_molecules.set_title("Albedo, Continuum Models, and Molecular Crossections")
    ax_molecules.set_xlim(x0,x1)
    ax_molecules.set_xticks(np.linspace(x0,x1,xtks, endpoint=True))
    ax_molecules.set_ylim(y0,y1)
    ax_molecules.set_yticks(np.linspace(y0,y1,ytks, endpoint=True))
    ax_molecules.grid(linewidth=0.2)
    ax_molecules.tick_params(axis='both', which='major', labelsize=8)
    ax_molecules.set_ylabel("Albedo",color="black")
    
    # PLOT ALBEDO AND CONTINUUM MODELS IN THE LIST
    ax_molecules.plot(Albedo[:,0],Albedo[:,1],label='Jupiter Albedo (Karkoschka, 1994)',linewidth=1.0,color='C0')
    for ModelName in ModList:
        ax_molecules.plot(Continua[ModelName]['WaveGrid'],
                          Continua[ModelName]['Albedo'],label=ModelName,
                          linewidth=0.5)
        
    #OVERPLOT MOLECULAR ABSORPTION CROSSECTIONS   
    if Crossect:
        axs1b = ax_molecules.twinx()  # instantiate a second axes that shares the same x-axis
        axs1b.ticklabel_format(axis='y')
        axs1b.tick_params(axis='y', which='major', labelsize=8)
        axs1b.set_yscale('log')
        axs1b.set_ylim(1e-4,1e3)
        axs1b.set_ylabel("Absorption Coefficient 1/(km-atm)")#,color="green")
        
        axs1b.plot(CH4[:,0],CH4[:,1],label='CH4 Abs. Coef. (Karkoschka, 1994) ',linewidth=1.0,color='C2')
        axs1b.plot(NH3[:,0],NH3[:,1],label='NH3 Abs. Coef. (ExoMol) ',linewidth=1.0,color='C3')
        if LutzPlot:
            axs1b.plot(NH3_Lutz_Owen_1980[:,0],NH3_Lutz_Owen_1980[:,1],label='NH3 Abs. Coef. (Lutz & Owen, 1980) ',linewidth=0.5,color='C3')
        axs1b.legend(fontsize=7, loc=1)

    ax_molecules.set_xlabel("Wavelength (nm)")
    ax_molecules.legend(fontsize=7, loc=2)
    fig_molecules.subplots_adjust(left=0.10, right=0.90, top=0.9, bottom=0.14)
    
    outpath="Molecular Absorption/data output/"

    fig_molecules.savefig(projpath+outpath+\
                          'Albedo_Continua_Crossections.png',dpi=320)
    
    return(Albedo,Continua,CH4,NH3,NH3_Lutz_Owen_1980)