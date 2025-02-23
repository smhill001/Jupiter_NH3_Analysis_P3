def Profile_L3(param="fNH3",profile="Meridional",ProfileHalfWidth=45,
               LatPlotLims=[60,120],ZonePlotHalfWidth=45,smooth=False,
               inset=True,Batch0="2022 CMOS",Batch1="2023 CMOS"):

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
    
    path="/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Profiles/output/"
       
    figavgprof,axsavgprof=pl.subplots(1,1,figsize=(6.0,4.0), dpi=150, facecolor="white")
    figspg,axsspg=pl.subplots(2,2,figsize=(8.0,6.0), dpi=150,
                                    sharex=True,sharey=True,facecolor="white")

    figamf,axsamf=pl.subplots(1,1,figsize=(6.0,6.0), dpi=150, facecolor="white")
    figamfspg,axsamfspg=pl.subplots(2,2,figsize=(8.0,6.0), dpi=150,
                                    sharex=True,sharey=True,facecolor="white")   
   
    LatsSCT22,OutProSCT22,OutStdSCT22,OutamfSCT22=PTG.plot_profile_L3(axsavgprof,axsamf,axsspg[0,0],axsamfspg[0,0],Batch0,ProfileHalfWidth=ProfileHalfWidth,
                        LatPlotLims=LatPlotLims,ZonePlotHalfWidth=ZonePlotHalfWidth,
                        profile=profile,clr='k',width=1.,param=param,smooth=smooth)
    LatsVLT22,OutProVLT22,OutStdVLT22,OutamfVLT22=PTG.plot_profile_L3(axsavgprof,axsamf,axsspg[0,1],axsamfspg[0,1],"2022 VLTMUSE",ProfileHalfWidth=ProfileHalfWidth,
                        LatPlotLims=LatPlotLims,ZonePlotHalfWidth=ZonePlotHalfWidth,
                        profile=profile,clr='k',width=0.5,param=param,smooth=smooth,
                        style='dashed')
    LatsSCT23,OutProSCT23,OutStdSCT23,OutamfSCT23=PTG.plot_profile_L3(axsavgprof,axsamf,axsspg[1,0],axsamfspg[1,0],Batch1,ProfileHalfWidth=ProfileHalfWidth,
                        LatPlotLims=LatPlotLims,ZonePlotHalfWidth=ZonePlotHalfWidth,
                        profile=profile,clr='C0',width=2.5,param=param,smooth=smooth)
    LatsSCT24,OutProSCT24,OutStdSCT24,OutamfSCT24=PTG.plot_profile_L3(axsavgprof,axsamf,axsspg[1,1],axsamfspg[1,1],"2024 CMOS",ProfileHalfWidth=ProfileHalfWidth,
                        LatPlotLims=LatPlotLims,ZonePlotHalfWidth=ZonePlotHalfWidth,
                        profile=profile,clr='C1',width=2.5,param=param,smooth=smooth)
    if inset==True and param=="fNH3" and profile=="Meridional":
        plevel=0.752910
        PTG.plot_TEXES_Groups(axsavgprof,clr='C2',prs=plevel,mult=1.0e6)
        plevel=0.657540
        PTG.plot_TEXES_Groups(axsavgprof,clr='C0',prs=plevel,mult=1.0e6)
        plevel=0.574240
        PTG.plot_TEXES_Groups(axsavgprof,clr='C4',prs=plevel,mult=1.0e6)
        if inset==True:
            x1,x2,y1,y2=-30.,10.,110,135.
            axsins=axsavgprof.inset_axes([0.07,0.5,0.45,0.20],xlim=(x1,x2),ylim=(y1,y2))
            axsavgprof.indicate_inset_zoom(axsins,edgecolor="black")
            PTG.plot_TEXES_Groups(axsins,clr='C0',prs=plevel,mult=1.0e6)
            axsins.tick_params(axis='both', which='major', labelsize=8)
            axsins.set_yticks(np.linspace(110,135,6),)   

        LatsSCT22,OutProSCT22,OutStdSCT22,OutamfSCT22=PTG.plot_profile_L3(axsins,axsamf,axsspg[0,0],axsamfspg[0,0],Batch0,ProfileHalfWidth=ProfileHalfWidth,
                            LatPlotLims=LatPlotLims,ZonePlotHalfWidth=ZonePlotHalfWidth,
                            profile=profile,clr='k',width=1.5,param=param,smooth=smooth)
        LatsVLT22,OutProVLT22,OutStdVLT22,OutamfVLT22=PTG.plot_profile_L3(axsins,axsamf,axsspg[0,1],axsamfspg[0,1],"2022 VLTMUSE",ProfileHalfWidth=ProfileHalfWidth,
                            LatPlotLims=LatPlotLims,ZonePlotHalfWidth=ZonePlotHalfWidth,
                            profile=profile,clr='k',width=1.5,param=param,smooth=smooth,
                            style='dashed')
        LatsSCT23,OutProSCT23,OutStdSCT23,OutamfSCT23=PTG.plot_profile_L3(axsins,axsamf,axsspg[1,0],axsamfspg[1,0],Batch1,ProfileHalfWidth=ProfileHalfWidth,
                            LatPlotLims=LatPlotLims,ZonePlotHalfWidth=ZonePlotHalfWidth,
                            profile=profile,clr='r',width=1.5,param=param,smooth=smooth)
        LatsSCT24,OutProSCT24,OutStdSCT24,OutamfSCT24=PTG.plot_profile_L3(axsins,axsamf,axsspg[1,1],axsamfspg[1,1],"2024 CMOS",ProfileHalfWidth=ProfileHalfWidth,
                            LatPlotLims=LatPlotLims,ZonePlotHalfWidth=ZonePlotHalfWidth,
                            profile=profile,clr='C3',width=1.5,param=param,smooth=smooth)


    # Plot layout details and labeling
    if profile=="Meridional":
        axsavgprof.set_xlim(90-LatPlotLims[1],90-LatPlotLims[0])
        axsavgprof.set_xlabel("Planetographic Latitude (deg)",fontsize=10)
        axsamf.set_xlim(1,3)
        axsamf.set_xlabel("One-Way Air Mass",fontsize=10)
    if profile=="Zonal":
        axsavgprof.set_xlim(-ZonePlotHalfWidth,ZonePlotHalfWidth)
        axsavgprof.set_xlabel("Longitude from Sys. II CM (deg)",fontsize=10)
        axsamf.set_xlim(1,3)
        axsamf.set_xlabel("One-Way Air Mass",fontsize=10)

    if param=="PCld":
        axsavgprof.set_title("Effective Cloud-Top Pressure Profiles")
        axsavgprof.set_ylim(1200.,2200.)
        axsavgprof.set_ylabel("Effective Cloud-Top Pressure (mb)",fontsize=10)
        axsamf.set_title("Effective Cloud-Top Pressure Profiles")
        axsamf.set_ylim(1200.,2200.)
        axsamf.set_ylabel("Effective Cloud-Top Pressure (mb)",fontsize=10)
        #axsavgprof.set_yticks(np.linspace(400,1100,8), minor=False)
        yb=2180
        yz=2140
        axsavgprof.invert_yaxis()
        axsamf.invert_yaxis()
    elif param=="fNH3":
        axsavgprof.set_title("Column-Integrated Ammonia Abundance Profiles")
        axsavgprof.set_ylim(0.,200.)
        axsavgprof.set_ylabel("Column-Integrated Ammonia Abundance (ppm)",fontsize=10)
        axsavgprof.set_yticks(np.linspace(0,200.,9), minor=False)
        axsamf.set_title("Column-Integrated Ammonia Abundance Profiles")
        axsamf.set_ylim(0.,200.)
        axsamf.set_ylabel("Column-Integrated Ammonia Abundance (ppm)",fontsize=10)
        axsamf.set_yticks(np.linspace(0,200,9), minor=False)
        yb=2
        yz=15
        
    if profile=="Meridional":
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
    figavgprof.savefig(path+"Profile_"+param+"_"+profile+".png",dpi=300)

    ###########################################################################
    # Plot change in profile relative to 2022 profile
    ###########################################################################
    
    figresid,axsresid=pl.subplots(1,1,figsize=(6.0,4.0), dpi=150, facecolor="white")
    proresid=OutProSCT23-OutProSCT22
    axsresid.plot(LatsVLT22,proresid,label="2023 SCT minus 2022 SCT",linewidth=2.5)
    stdresid=np.sqrt(OutStdSCT22**2+OutStdSCT23**2)
    axsresid.fill_between(LatsVLT22,proresid-stdresid,proresid+stdresid,
                          color='C0',alpha=.1)
    
    proresid=OutProSCT24-OutProSCT23
    axsresid.plot(LatsSCT23,proresid,label="2024 SCT minus 2023 SCT",linewidth=2.5)
    stdresid=np.sqrt(OutStdSCT23**2+OutStdSCT24**2)
    axsresid.fill_between(LatsSCT23,proresid-stdresid,proresid+stdresid,
                          color='C1',alpha=.1)

    if profile=="Meridional":
        axsresid.set_xlim(90-LatPlotLims[1],90-LatPlotLims[0])
        axsresid.set_xlabel("Planetographic Latitude (deg)",fontsize=10)
    if profile=="Zonal":
        axsresid.set_xlim(-ZonePlotHalfWidth,ZonePlotHalfWidth)
        axsresid.set_xlabel("Longitude from Sys. II CM (deg)",fontsize=10)

    axsresid.grid(linewidth=0.2)
    axsresid.tick_params(axis='both', which='major', labelsize=8)
    axsresid.legend(fontsize=8)

    if param=="PCld":
        axsresid.set_title("Change in Cloud Top Pressure")
        axsresid.set_ylim(-250.,250.)
        axsresid.set_ylabel("Effective Cloud Top Pressure (mb)",fontsize=10)
        axsresid.set_yticks(np.linspace(-250,250,11), minor=False)
        axsresid.invert_yaxis()
        yb=245
        yz=225
        annotate_sign=1
    elif param=="fNH3":
        axsresid.set_title("Change in Ammonia Abundance")
        axsresid.set_ylim(-40.,40.)
        axsresid.set_ylabel("Column-Integrated Ammonia Abundance (ppm)",fontsize=10)
        axsresid.set_yticks(np.linspace(-40,40,9), minor=False)
        yb=-39
        yz=-34
        
    if profile=="Meridional":
        for zb in belt:
            print(zb,belt[zb])
            axsresid.fill_between([belt[zb][0],belt[zb][1]],np.array([-5000.,-5000.]),
                                    np.array([5000.,5000.]),color="0.5",alpha=0.2)
            axsresid.annotate(zb,xy=[np.mean(belt[zb]),yb],ha="center")
        for zb in zone:
            axsresid.annotate(zb,xy=[np.mean(zone[zb]),yz],ha="center")
    
    figresid.subplots_adjust(left=0.12, bottom=0.12, right=0.98, top=0.92)  
    figresid.savefig(path+"Residuals_"+param+"_"+profile+".png",dpi=300)

    SCT22={'Lats':LatsSCT22,'Pro':OutProSCT22,'Std':OutStdSCT22,'Amf':OutamfSCT22}
    VLT22={'Lats':LatsVLT22,'Pro':OutProVLT22,'Std':OutStdVLT22,'Amf':OutamfVLT22}
    SCT23={'Lats':LatsSCT23,'Pro':OutProSCT23,'Std':OutStdSCT23,'Amf':OutamfSCT23}
    SCT24={'Lats':LatsSCT24,'Pro':OutProSCT24,'Std':OutStdSCT24,'Amf':OutamfSCT24}

    #return(figavgprof,axsavgprof)
   
    return(SCT22,VLT22,SCT23,SCT24)#,figavgprof,axsavgprof)
