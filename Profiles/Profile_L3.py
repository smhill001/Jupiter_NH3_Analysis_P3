# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 10:26:54 2023

@author: smhil
"""
def Profile_L3(param="fNH3",profile="Meridional",ProfileHalfWidth=45,
               LatPlotLims=[60,120],ZonePlotHalfWidth=45,smooth=False):

    #(param="PCloud",profile="Meridional",LonRng=1):
    import sys
    drive='C:'
    sys.path.append(drive+'/Astronomy/Python Play')
    sys.path.append(drive+'/Astronomy/Python Play/Util_P3')
    sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy_P3')
    sys.path.append(drive+'/Astronomy/Python Play/SPLibraries_P3')
    import pylab as pl
    import numpy as np
    import plot_TEXES_Groups_P3 as PTG
             
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
       
    figavgprof,axsavgprof=pl.subplots(1,1,figsize=(6.0,4.0), dpi=150, facecolor="white")
    
    #ax=axsavgprof
    ax='None'
    if param=="fNH3" and profile=="Meridional":
        plevel=0.752910
        PTG.plot_TEXES_Groups(ax,clr='C2',prs=plevel,mult=1000000.)
        plevel=0.657540
        PTG.plot_TEXES_Groups(ax,clr='r',prs=plevel,mult=1000000.)
        plevel=0.574240
        PTG.plot_TEXES_Groups(ax,clr='C4',prs=plevel,mult=1000000.)
    
    
    LatsSCT22,OutProSCT22,OutStdSCT22=PTG.plot_profile_L3(axsavgprof,"SCT 2022 L3",ProfileHalfWidth=ProfileHalfWidth,
                        LatPlotLims=LatPlotLims,ZonePlotHalfWidth=ZonePlotHalfWidth,
                        profile=profile,clr='k',width=1.5,param=param,smooth=smooth)
    LatsVLT22,OutProVLT22,OutStdVLT22=PTG.plot_profile_L3(axsavgprof,"VLTMUSE 2022",ProfileHalfWidth=ProfileHalfWidth,
                        LatPlotLims=LatPlotLims,ZonePlotHalfWidth=ZonePlotHalfWidth,
                        profile=profile,clr='k',width=1.5,param=param,smooth=smooth,
                        style='dashed')
    LatsSCT23,OutProSCT23,OutStdSCT23=PTG.plot_profile_L3(axsavgprof,"SCT 2023",ProfileHalfWidth=ProfileHalfWidth,
                        LatPlotLims=LatPlotLims,ZonePlotHalfWidth=ZonePlotHalfWidth,
                        profile=profile,clr='r',width=1.5,param=param,smooth=smooth)

    if profile=="Meridional":
        for zb in belt:
            print(zb,belt[zb])
            axsavgprof.fill_between([belt[zb][0],belt[zb][1]],np.array([0.,0.]),
                                    np.array([2000.,2000.]),color="0.5",alpha=0.2)
            axsavgprof.annotate(zb,xy=[np.mean(belt[zb]),0.01],ha="center")
        for zb in zone:
            axsavgprof.annotate(zb,xy=[np.mean(zone[zb]),0.05],ha="center")
    
    # Plot layout details and labeling
    if profile=="Meridional":
        axsavgprof.set_xlim(90-LatPlotLims[1],90-LatPlotLims[0])
        axsavgprof.set_xlabel("Planetographic Latitude (deg)",fontsize=10)
    if profile=="Zonal":
        axsavgprof.set_xlim(-ZonePlotHalfWidth,ZonePlotHalfWidth)
        axsavgprof.set_xlabel("Longitude from Sys. II CM (deg)",fontsize=10)

    if param=="PCloud":
        axsavgprof.set_title("Effective Cloud-Top Pressure Profiles")
        axsavgprof.set_ylim(0.,1100.)
        axsavgprof.set_ylabel("Effective Cloud-Top Pressure (mb)",fontsize=10)
        #axsavgprof.set_yticks(np.linspace(400,1100,8), minor=False)
        axsavgprof.invert_yaxis()
    elif param=="fNH3":
        axsavgprof.set_title("Column-Integrated Ammonia Abundance Profiles")
        axsavgprof.set_ylim(0.,200.)
        axsavgprof.set_ylabel("Column-Integrated Ammonia Abundance (ppm)",fontsize=10)
        axsavgprof.set_yticks(np.linspace(0,200,9), minor=False)

    axsavgprof.legend(loc=1,ncol=2, borderaxespad=0.,prop={'size':8})
    axsavgprof.grid(linewidth=0.2)
    axsavgprof.tick_params(axis='both', which='major', labelsize=8)
    
    axsavgprof.annotate("ProfileHalfWidth="+str(ProfileHalfWidth),(0.01,0.01),
                        xycoords='subfigure fraction',size=8)
    axsavgprof.annotate("Smoothing="+str(smooth),(0.01,0.04),
                        xycoords='subfigure fraction',size=8)
    
    figavgprof.subplots_adjust(left=0.12, bottom=0.12, right=0.98, top=0.92)  
    figavgprof.savefig(path+"Profile_"+param+"_"+profile+".png",dpi=300)

    ###########################################################################
    # Plot change in profile relative to 2022 profile
    ###########################################################################
    
    figresid,axsresid=pl.subplots(1,1,figsize=(6.0,4.0), dpi=150, facecolor="white")
    proresid=OutProSCT23-OutProSCT22
    axsresid.plot(LatsVLT22,proresid,label='2023 minus 2022')
    stdresid=np.sqrt(OutStdSCT22**2+OutStdSCT23**2)
    axsresid.fill_between(LatsVLT22,proresid-stdresid,proresid+stdresid,
                          color='C0',alpha=.1)

    if profile=="Meridional":
        axsresid.set_xlim(90-LatPlotLims[1],90-LatPlotLims[0])
        axsresid.set_xlabel("Planetographic Latitude (deg)",fontsize=10)
    if profile=="Zonal":
        axsresid.set_xlim(-ZonePlotHalfWidth,ZonePlotHalfWidth)
        axsresid.set_xlabel("Longitude from Sys. II CM (deg)",fontsize=10)

    axsresid.set_title("Change in Ammonia Abundance Profiles")
    axsresid.set_ylim(-40.,40.)
    axsresid.set_ylabel("Column-Integrated Ammonia Abundance (ppm)",fontsize=10)
    axsresid.set_yticks(np.linspace(-40,40,9), minor=False)
    axsresid.grid(linewidth=0.2)
    axsresid.tick_params(axis='both', which='major', labelsize=8)
    axsresid.legend(fontsize=8)

    if profile=="Meridional":
        for zb in belt:
            print(zb,belt[zb])
            axsresid.fill_between([belt[zb][0],belt[zb][1]],np.array([-2000.,-2000.]),
                                    np.array([2000.,2000.]),color="0.5",alpha=0.2)
            axsresid.annotate(zb,xy=[np.mean(belt[zb]),-14.9],ha="center")
        for zb in zone:
            axsresid.annotate(zb,xy=[np.mean(zone[zb]),-14.0],ha="center")
            
    figresid.subplots_adjust(left=0.12, bottom=0.12, right=0.98, top=0.92)  
    figresid.savefig(path+"Residuals_"+param+"_"+profile+".png",dpi=300)
