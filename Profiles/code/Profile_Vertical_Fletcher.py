# -*- coding: utf-8 -*-

def Profile_Vertical_Fletcher():
    """
    Created on Fri Nov 24 11:55:45 2023

    @author: smhil
    """
    import sys
    drive='C:'
    sys.path.append(drive+'/Astronomy/Python Play')
    sys.path.append(drive+'/Astronomy/Python Play/Util_P3')
    sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy_P3')
    sys.path.append(drive+'/Astronomy/Python Play/SPLibraries_P3')
    import pylab as pl
    import scipy
    import numpy as np
    import plot_TEXES_Groups_P3 as PTG
    import matplotlib.ticker as ticker



    pth="c:/Astronomy/Projects/SAS 2021 Ammonia/GeminiTEXES2017/ZonalResults/"
    pressure = np.fromfile(file=pth+"zmean_g1_retnh3_pressure.txt", dtype=float, count=-1, sep=" ")
    pressure=np.array(pressure)
    size=pressure.size
    dataavg=np.zeros(size)
    datastd=np.zeros(size)
    CH4=np.ones(size)*1810.
    figavgprof,axsavgprof=pl.subplots(1,1,figsize=(6.0,6.0), dpi=150, facecolor="white")

    #plevel=0.752910
    for i in np.arange(0,size):
        plevel=pressure[i]
        scaled_data_mean,scaled_data_std=PTG.plot_TEXES_Groups(axsavgprof,
                                                               clr='C2',
                                                               prs=plevel,
                                                               mult=1000000.)
        dataavg[i]=np.nanmean(scaled_data_mean)
        #datastd[i]=scaled_data_std
        
    axsavgprof.set_yscale('log')
    axsavgprof.set_ylim([0.001,100.])
    axsavgprof.invert_yaxis()
    
    figvertprof,axsvertprof=pl.subplots(1,1,figsize=(6.0,6.0), dpi=150, facecolor="white")
    axsvertprof.plot(dataavg,pressure,label='Ammonia - Fletcher etal. [2020]')
    axsvertprof.plot(CH4,pressure,label='Methane - citation')
    print(dataavg)
    #Upper Haze
    axsvertprof.fill_between([1,10000],[0.6,0.6],[0.2,0.2],color='0.9',
                             linewidth=[0,0],label='Haze')
    #Sheet Cloud
    axsvertprof.fill_between([1,10000],[0.65,0.65],[0.67,0.67],color='0.5',
                             linewidth=[0,0],label='Sheet Cloud')

    axsvertprof.set_xscale('log')
    axsvertprof.set_xlim([1,10000.])
    axsvertprof.set_yscale('log')
    axsvertprof.set_ylim([0.1,1.])
    axsvertprof.yaxis.set_major_formatter(ticker.FormatStrFormatter('%1.1f'))
    axsvertprof.yaxis.set_minor_formatter(ticker.FormatStrFormatter('%1.1f'))
    axsvertprof.invert_yaxis()
    
    axsvertprof.set_title('Abundances and Clouds')
    axsvertprof.set_xlabel('Abundance (ppm)')
    axsvertprof.set_ylabel('Pressure (bar)')
    axsvertprof.legend()
    path="c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Analysis Data/Profiles/"

    figvertprof.savefig(path+"Profile Vertical Fletcher.png",dpi=300)

