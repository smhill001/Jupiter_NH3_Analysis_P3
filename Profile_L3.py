# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 10:26:54 2023

@author: smhil
"""
def Profile_L3(param="PCloud",profile="Meridional",LonRng=1):
    import sys
    drive='C:'
    sys.path.append(drive+'/Astronomy/Python Play')
    sys.path.append(drive+'/Astronomy/Python Play/Util_P3')
    sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy_P3')
    sys.path.append(drive+'/Astronomy/Python Play/SPLibraries_P3')
    import pylab as pl
    import numpy as np
    import plot_TEXES_Groups_P3 as PTG
    
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
    
    path="/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Analysis Data/Profiles/"
       
    figavgprof,axsavgprof=pl.subplots(1,1,figsize=(6.0,6.0), dpi=150, facecolor="white")
    
    if param=="fNH3" and profile=="Meridional":
        plevel=0.752910
        PTG.plot_TEXES_Groups(axsavgprof,clr='C2',prs=plevel,mult=1000000.)
        plevel=0.657540
        PTG.plot_TEXES_Groups(axsavgprof,clr='C0',prs=plevel,mult=1000000.)
        plevel=0.574240
        PTG.plot_TEXES_Groups(axsavgprof,clr='C4',prs=plevel,mult=1000000.)

    
    PTG.plot_profiles_L3(axsavgprof,"SCT 2022",LonRng=LonRng,profile=profile,
                             clr='k',width=1.5,smooth=smooth,param=param)
    PTG.plot_profiles_L3(axsavgprof,"VLTMUSE 2022",LonRng=LonRng,profile=profile,
                             clr='k',width=1.0,smooth=smooth,param=param,
                             style='dashed')
    PTG.plot_profiles_L3(axsavgprof,"SCT 2023",LonRng=LonRng,profile=profile,
                             clr='r',width=1.5,smooth=smooth,param=param)
    
    if profile=="Meridional":
        for zb in belt:
            print(zb,belt[zb])
            axsavgprof.fill_between([belt[zb][0],belt[zb][1]],np.array([0.,0.]),
                                    np.array([1000.,1000.]),color="0.5",alpha=0.2)
            axsavgprof.annotate(zb,xy=[np.mean(belt[zb]),0.01],ha="center")
        for zb in zone:
            axsavgprof.annotate(zb,xy=[np.mean(zone[zb]),0.05],ha="center")
    
    # Plot layout details and labeling
    axsavgprof.set_xlim(-45.,45.)
    axsavgprof.set_xlabel("Planetographic Latitude (deg)",fontsize=10)
    if param=="PCloud":
        axsavgprof.set_title("Effective Cloud-Top Pressure Profiles")
        axsavgprof.set_ylim(400.,1000.)
        axsavgprof.set_ylabel("Effective Cloud-Top Pressure (mb)",fontsize=10)
        axsavgprof.set_yticks(np.linspace(400,1000,7), minor=False)
        axsavgprof.invert_yaxis()
    elif param=="fNH3":
        axsavgprof.set_title("Column-Integrated Ammonia Abundance Profiles")
        axsavgprof.set_ylim(0.,200.)
        axsavgprof.set_ylabel("Column-Integrated Ammonia Abundance (ppm)",fontsize=10)
        axsavgprof.set_yticks(np.linspace(0,300,7), minor=False)

    axsavgprof.legend(loc=1,ncol=2, borderaxespad=0.,prop={'size':8})
    axsavgprof.grid(linewidth=0.2)
    axsavgprof.tick_params(axis='both', which='major', labelsize=8)
    
    figavgprof.subplots_adjust(left=0.12, bottom=0.12, right=0.98, top=0.92)  
    figavgprof.savefig(path+"Profile_"+param+"_"+profile+".png",dpi=300)