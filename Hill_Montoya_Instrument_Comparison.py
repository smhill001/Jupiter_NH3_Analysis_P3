# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 10:08:45 2022

@author: smhil
"""

def Hill_Montoya_Instrument_Comparison():
    import sys
    sys.path.append('c:/Astronomy/Python Play')
    sys.path.append('c:/Astronomy/Python Play/Util_P3')
    sys.path.append('c:/Astronomy/Python Play/SPLibraries_P3')
    sys.path.append('c:/Astronomy/Python Play/SpectroPhotometry/Spectroscopy_P3')
    import matplotlib.pyplot as pl
    import numpy as np

    path_Montoyoa="c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Montoya/"
    Montoya_NUV=path_Montoyoa+"Filter-U.txt"
    Montoya_BLU=path_Montoyoa+"Filter-B.txt"
    Montoya_GRN=path_Montoyoa+"Filter-G.txt"
    Montoya_RED=path_Montoyoa+"Filter-R.txt"
    Montoya_CH4=path_Montoyoa+"Filter-CH4.txt"
    
    path_Hill="c:/Astronomy/Projects/Techniques/InstrumentPerformance-P3/Filters/"
    Hill_NUV=path_Hill+"380NUV/380NUV_Transmission.txt"
    Hill_BLU=path_Hill+"450BLU/450BLU_Transmission.txt"
    Hill_GRN=path_Hill+"550GRN/550GRN_Transmission.txt"
    Hill_RED=path_Hill+"650RED/650RED_Transmission.txt"
    Hill_CH4=path_Hill+"889CH4/889CH4_Transmission.txt"


    ###### Retrieve filter transmissions and convovle with disk integrated albedoes
    M_NUV = np.fromfile(file=Montoya_NUV, dtype=float, count=-1, sep='\t')    
    M_NUV=np.reshape(M_NUV,[int(M_NUV.size/2),2])
    M_BLU = np.fromfile(file=Montoya_BLU, dtype=float, count=-1, sep='\t')    
    M_BLU=np.reshape(M_BLU,[int(M_BLU.size/2),2])
    M_GRN = np.fromfile(file=Montoya_GRN, dtype=float, count=-1, sep='\t')    
    M_GRN=np.reshape(M_GRN,[int(M_GRN.size/2),2])
    M_RED = np.fromfile(file=Montoya_RED, dtype=float, count=-1, sep='\t')    
    M_RED=np.reshape(M_RED,[int(M_RED.size/2),2])
    M_CH4 = np.fromfile(file=Montoya_CH4, dtype=float, count=-1, sep='\t')    
    M_CH4=np.reshape(M_CH4,[int(M_CH4.size/2),2])

    H_NUV = np.fromfile(file=Hill_NUV, dtype=float, count=-1, sep='\t')    
    H_NUV=np.reshape(H_NUV,[int(H_NUV.size/4),4])
    H_BLU = np.fromfile(file=Hill_BLU, dtype=float, count=-1, sep='\t')    
    H_BLU=np.reshape(H_BLU,[int(H_BLU.size/4),4])
    H_GRN = np.fromfile(file=Hill_GRN, dtype=float, count=-1, sep='\t')    
    H_GRN=np.reshape(H_GRN,[int(H_GRN.size/4),4])
    H_RED = np.fromfile(file=Hill_RED, dtype=float, count=-1, sep='\t')    
    H_RED=np.reshape(H_RED,[int(H_RED.size/4),4])
    H_CH4 = np.fromfile(file=Hill_CH4, dtype=float, count=-1, sep='\t')    
    H_CH4=np.reshape(H_CH4,[int(H_CH4.size/2),2])
    
    fig,axs=pl.subplots(1,1,figsize=(6.5, 2.5), dpi=150, facecolor="white")

    x0=300.
    x1=1000.
    xtks=15
    y0=0.0
    y1=1.2
    ytks=7

    # Set x limits
    axs.set_xlim(x0,x1)
    # Set x ticks
    axs.set_xticks(np.linspace(x0,x1,xtks, endpoint=True))
    # Set y limits
    axs.set_ylim(y0,y1)
    axs.set_yticks(np.linspace(y0,y1,ytks, endpoint=True))
    # Set y ticks
    axs.grid(linewidth=0.2)
    axs.tick_params(axis='both', which='major', labelsize=8)
    axs.set_ylabel("Transmission",fontsize=8,color="black")
    axs.set_xlabel("Wavelength (nm)",fontsize=8)
    axs.plot(M_NUV[:,0],M_NUV[:,1]/100.,'--',label='Montoya NUV',linewidth=1,color='m')
    axs.plot(H_NUV[:,0],H_NUV[:,1]/100.,label='Hill NUV',linewidth=1,color='m')
    
    axs.plot(M_BLU[:,0],M_BLU[:,1]/100.,'--',label='Montoya BLU',linewidth=1,color='b')
    axs.plot(H_BLU[:,0],H_BLU[:,1]/1.2,label='Hill BLU',linewidth=1,color='b')

    axs.plot(M_GRN[:,0],M_GRN[:,1]/100.,'--',label='Montoya GRN',linewidth=1,color='g')
    axs.plot(H_GRN[:,0],H_GRN[:,1]/1.1,label='Hill GRN',linewidth=1,color='g')

    axs.plot(M_RED[:,0],M_RED[:,1]/100.,'--',label='Montoya RED',linewidth=1,color='r')
    axs.plot(H_RED[:,0],H_RED[:,1]/1.5,label='Hill RED',linewidth=1,color='r')

    axs.plot(M_CH4[:,0],M_CH4[:,1]/100.,'--',label='Montoya CH4',linewidth=1,color='k')
    axs.plot(H_CH4[:,0],H_CH4[:,1],label='Hill CH4',linewidth=1,color='k')

    axs.legend(loc=0,ncol=5, borderaxespad=0.,prop={'size':6})
    pl.subplots_adjust(left=0.08, bottom=0.15, right=0.98, top=0.92,
            wspace=None, hspace=None)

    pl.savefig('c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Filter Comparison.png',dpi=320)
