# -*- coding: utf-8 -*-
def Profile_L2_Variations(band="CH4",profile="Meridional",LonRng=1,
                          Variation="VLT 2022 Var"):
    """
    Created on Thu Sep  7 11:17:51 2023
    
    Plots absorption annual profile averages from multiple years for
    both methane and ammonia.
    
    @author: smhil
    """
    import sys
    drive='C:'
    sys.path.append(drive+'/Astronomy/Python Play')
    import pylab as pl
    import numpy as np
    from astropy.convolution import convolve
    import plot_TEXES_Groups_P3 as PTG
    from astropy.convolution import Box1DKernel
    
    smooth=False
          
    belt={"SSTB":[-39.6,-36.2],
          "STB":[-32.4,-27.1],
          "SEB":[-19.7,-7.2],
          "NEB":[6.9,17.4],
          "NTB":[24.2,31.4],
          "NNTB":[35.4,39.6]}
    
    zone={"STZ":[-36.2,-32.4],
          "STrZ":[-27.1,-19.7],
          "EZ":[-7.2,6.9],
          "NTrZ":[17.4,24.2],
          "NTZ":[31.4,35.4]}
    
    path="/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/"
       
    figavgprof,axsavgprof=pl.subplots(1,1,figsize=(6.0,4.0), dpi=150, facecolor="white")
    
    if band=="NH3" and profile=="Meridional":
        PTG.plot_Teifel(axsavgprof,clr='0.5',width=3.)
   
    PTG.plot_profile_L2(axsavgprof,Variation,LonRng=LonRng,
                        profile=profile,clr='k',width=1.5,band=band,style='dashed',smooth=smooth)
    
    if profile=="Meridional":
        for zb in belt:
            #print(zb,belt[zb])
            axsavgprof.fill_between([belt[zb][0],belt[zb][1]],np.array([0.,0.]),np.array([1000.,1000.]),
                                    color="0.5",alpha=0.2)
            axsavgprof.annotate(zb,xy=[np.mean(belt[zb]),0.01],ha="center")
        for zb in zone:
            axsavgprof.annotate(zb,xy=[np.mean(zone[zb]),0.05],ha="center")
    
    # Plot layout details and labeling
    axsavgprof.grid(linewidth=0.2)
    axsavgprof.set_xlim(-45.,45.)
    axsavgprof.set_ylim(0.0,1.0)
    axsavgprof.set_xticks(np.linspace(-45.,45.,7), minor=False)
    if band=="CH4":
        axsavgprof.set_title("Methane 619 nm Absorption Profiles")
        axsavgprof.set_yticks(np.linspace(0.0,2.0,5), minor=False)
    elif band=="NH3":
        axsavgprof.set_title("Ammonia 646 nm Absorption Profiles")
        axsavgprof.set_yticks(np.linspace(0.0,1.0,5), minor=False)
    if profile=="Meridional":
        xlabel="Planetographic Latitude (deg)"
    elif profile=="Zonal":
        xlabel="Longitude from CM (deg)"
    axsavgprof.tick_params(axis='both', which='major', labelsize=8)
    axsavgprof.set_xlabel(xlabel,fontsize=10)
    axsavgprof.set_ylabel("Equivalent Width (nm)",fontsize=10)
    axsavgprof.legend(fontsize=8,loc=2)
    figavgprof.subplots_adjust(left=0.09, right=0.97, top=0.93, bottom=0.11)
      
    figavgprof.savefig(drive+path+"Analysis Data/Profiles/Profile_"+band+"_"+profile+"_Absorption.png",dpi=300)