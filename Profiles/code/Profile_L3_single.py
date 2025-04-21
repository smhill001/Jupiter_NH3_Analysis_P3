def Profile_L3_single(param="fNH3",profile="Meridional",ProfileHalfWidth=45,
               LatPlotLims=[60,120],ZonePlotHalfWidth=45,smooth=False,
               inset=True,collection="2022 CMOS",colat=90):
    """
    PURPOSE:    To compute and plot averaged meridional and zonal profiles over
                a single set of observing sessions
                
    CALLS:      plot_profile_L3_granular
                
    Parameters
    ----------
    param : TYPE, optional
        DESCRIPTION. The default is "fNH3".
    profile : TYPE, optional
        DESCRIPTION. The default is "Meridional".
    ProfileHalfWidth : TYPE, optional
        DESCRIPTION. The default is 45.
    LatPlotLims : TYPE, optional
        DESCRIPTION. The default is [60,120].
    ZonePlotHalfWidth : TYPE, optional
        DESCRIPTION. The default is 45.
    smooth : TYPE, optional
        DESCRIPTION. The default is False.
    inset : TYPE, optional
        DESCRIPTION. The default is True.
    Batch0 : TYPE, optional
        DESCRIPTION. The default is "2022 CMOS".
    Batch1 : TYPE, optional
        DESCRIPTION. The default is "2023 CMOS".

    Returns
    -------
    None.

    """
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
    import plot_profile_L3_granular as PPL3G
    
    ###########################################################################
    #! Set up belt and zone boundaries (need to make this a common service!)             
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
    
    pathout="/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Profiles/output/"
    
    ###########################################################################
    # Set up figures for average profiles and spaghetti plots versus longitude
    # (or latitude if meridional?), and scatter plots versus air mass and air
    # mass spaghetti plots.
    
    figavgprof,axsavgprof=pl.subplots(1,1,figsize=(6.0,4.0), dpi=150, facecolor="white")
    figspg,axsspg=pl.subplots(1,1,figsize=(8.0,6.0), dpi=150,
                                    sharex=True,sharey=True,facecolor="white")

    figamf,axsamf=pl.subplots(1,1,figsize=(6.0,6.0), dpi=150, facecolor="white")
    figamfspg,axsamfspg=pl.subplots(1,1,figsize=(8.0,6.0), dpi=150,
                                    sharex=True,sharey=True,facecolor="white")   
   
    ###########################################################################
    # Compute profiles for the SCT for 2022, 2023, 2024 and for the VLT for 
    # 2022. Requires four separate calls to plot_profile_L3.
    LatsSCT22,OutProSCT22,OutStdSCT22,OutamfSCT22,NumS22=PPL3G.plot_profile_L3_granular(axsspg,
                        axsamfspg,collection,ProfileHalfWidth=ProfileHalfWidth,
                        LatPlotLims=LatPlotLims,ZonePlotHalfWidth=ZonePlotHalfWidth,
                        profile=profile,clr='k',width=1.,param=param,smooth=smooth,
                        colat=colat)
            
    axsavgprof.plot(LatsSCT22,OutProSCT22,color='k',linewidth=1.,linestyle='solid',
            #label=reference+' (Avg. '+str(Num)+')')  
            label='SCT 2022 (Avg. ')#+str(Num)+')')  

    axsavgprof.fill_between(LatsSCT22, OutProSCT22-OutStdSCT22, OutProSCT22+OutStdSCT22,
                    color='k',alpha=.05)
    
    axsamf.scatter(OutamfSCT22,OutProSCT22,s=5,label='2022 SCT (Avg. '+str(NumS22)+')')

    ###########################################################################
    # Plot layout details and labeling
    #if profile=="Meridional":
    if param=="PCld":
        axsavgprof.set_title("Cloud Pressure Profiles")
        axsavgprof.set_ylim(1200.,2200.)
        axsavgprof.set_ylabel("Cloud Pressure (mb)",fontsize=10)
        axsamf.set_title("Cloud Pressure Profiles")
        axsamf.set_ylim(1200.,2200.)
        axsamf.set_ylabel("Cloud Pressure (mb)",fontsize=10)
        #axsavgprof.set_yticks(np.linspace(400,1100,8), minor=False)
        yb=2180 #Belt and Zone label locations
        yz=2140
        axsavgprof.invert_yaxis()
        axsamf.invert_yaxis()
    elif param=="fNH3":
        axsavgprof.set_title("Ammonia Abundance Profiles")
        axsavgprof.set_ylim(0.,200.)
        axsavgprof.set_ylabel("Ammonia Abundance (ppm)",fontsize=10)
        axsavgprof.set_yticks(np.linspace(0,200.,9), minor=False)
        axsamf.set_title("Ammonia Abundance Profiles")
        axsamf.set_ylim(0.,200.)
        axsamf.set_ylabel("Ammonia Abundance (ppm)",fontsize=10)
        axsamf.set_yticks(np.linspace(0,200,9), minor=False)
        yb=2  #Belt and Zone label locations
        yz=15
        
    if profile=="Zonal":
        axsavgprof.set_xlim(-ZonePlotHalfWidth,ZonePlotHalfWidth)
        axsavgprof.set_xlabel("Longitude from Sys. II CM (deg)",fontsize=10)
        axsamf.set_xlim(1,3)
        axsamf.set_xlabel("One-Way Air Mass",fontsize=10)

    if profile=="Meridional":
            axsavgprof.set_xlim(90-LatPlotLims[1],90-LatPlotLims[0])
            axsavgprof.set_xlabel("Planetographic Latitude (deg)",fontsize=10)
            axsamf.set_xlim(1,3)
            axsamf.set_xlabel("One-Way Air Mass",fontsize=10)
            for zb in belt:
                print(zb,belt[zb])
                axsavgprof.fill_between([belt[zb][0],belt[zb][1]],np.array([0.,0.]),
                                        np.array([5000.,5000.]),color="0.5",alpha=0.2)
                axsavgprof.annotate(zb,xy=[np.mean(belt[zb]),yb],ha="center")
            for zb in zone:
                axsavgprof.annotate(zb,xy=[np.mean(zone[zb]),yz],ha="center")
        

    axsavgprof.legend(loc=1,ncol=2, borderaxespad=0.,prop={'size':8})
    axsavgprof.grid(linewidth=0.2)
    axsavgprof.tick_params(axis='both', which='major', labelsize=8)
    
    axsavgprof.annotate("ProfileHalfWidth="+str(ProfileHalfWidth),(0.01,0.01),
                        xycoords='subfigure fraction',size=8)
    axsavgprof.annotate("Smoothing="+str(smooth),(0.01,0.04),
                        xycoords='subfigure fraction',size=8)
    
    figavgprof.subplots_adjust(left=0.12, bottom=0.12, right=0.98, top=0.94)  
    figavgprof.savefig(pathout+"Profile_"+param+"_"+profile+".png",dpi=300)

    axsamf.legend(loc=1,ncol=2, borderaxespad=0.,prop={'size':8})

    SCT22={'Lats':LatsSCT22,'Pro':OutProSCT22,'Std':OutStdSCT22,'Amf':OutamfSCT22}
   
    return(SCT22)#,figavgprof,axsavgprof)
