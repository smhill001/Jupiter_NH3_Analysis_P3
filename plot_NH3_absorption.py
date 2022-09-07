# -*- coding: utf-8 -*-
"""
Created on Sun Aug  7 13:00:33 2022

@author: smhil
"""

def plot_NH3_absorption(ax,clr='C0'):
    """
    PURPOSE:    This code reads and plots the NH3 laboratory measured
                absorption of the 645nm Lutz & Owen, 1980.
    """
    import numpy as np
    from numpy import genfromtxt
    import GeneralSpecUtils_P3 as GSU
    fn='Lutz&Owen1980_Figure5_AmmoniaCrossSection.txt'
    pth="c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/"
    Lutz = np.array(genfromtxt(pth+fn, delimiter=','))
    startwvnum=1.0e7/660.0
    wavelength=1.0e7/(Lutz[:,0]+startwvnum)
    temp=np.flip(np.array([wavelength,1000.*Lutz[:,1]]).transpose(),0)
    WaveGrid,SignalonGrid=GSU.uniform_wave_grid(temp[:,0],temp[:,1],Extend=False,Fine=False)
    print(np.nanmax(SignalonGrid))
    NH3_K_Lutz=np.array([WaveGrid,SignalonGrid])
    print(np.nanmax(NH3_K_Lutz[1,:]))

    ax.scatter(NH3_K_Lutz[0,:],NH3_K_Lutz[1,:],label='Lutz & Owen, 1980',color=clr)
    np.savetxt(pth+"foo.csv", NH3_K_Lutz.transpose(), delimiter=",")    