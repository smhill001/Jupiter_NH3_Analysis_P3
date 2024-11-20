def Profile_L2(band="NH3",profile="Meridional",ProfileHalfWidth=45,
               LatPlotLims=[30,150],ZonePlotHalfWidth=45,smooth=False):

    
    
    import sys
    sys.path.append('c:/Astronomy/Python Play')
    import pylab as pl
    import numpy as np
    from astropy.convolution import convolve
    import plot_TEXES_Groups_P3 as PTG
    from astropy.convolution import Box1DKernel
              
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
    
    path="c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/"
       
    figavgprof,axsavgprof=pl.subplots(1,1,figsize=(6.0,4.0), dpi=150, facecolor="white")
    figspg,axsspg=pl.subplots(2,3,figsize=(12.0,6.0), dpi=150,
                                    sharex=True,sharey=True,facecolor="white")
    #figavgprof.suptitle=("Test")
    figamf,axsamf=pl.subplots(1,1,figsize=(6.0,4.0), dpi=150, facecolor="white")
    figamfspg,axsamfspg=pl.subplots(2,3,figsize=(12.0,6.0), dpi=150,
                                    sharex=True,sharey=True,facecolor="white")

    if band=="NH3" and profile=="Meridional":
        PTG.plot_Teifel(axsavgprof,clr='0.5',width=4.)
        
    if band=="NH3":
        LatsSCT20,OutProSCT20,OutStdSCT20,OutamfSCT20=PTG.plot_profile_L2(axsavgprof,axsamf,axsspg[0,0],axsamfspg[0,0],"2020 CMOS",ProfileHalfWidth=ProfileHalfWidth,
                            LatPlotLims=LatPlotLims,ZonePlotHalfWidth=ZonePlotHalfWidth,
                            profile=profile,clr='C0',width=1.0,band=band,smooth=smooth)
        LatsSCT21,OutProSCT21,OutStdSCT21,OutamfSCT21=PTG.plot_profile_L2(axsavgprof,axsamf,axsspg[0,1],axsamfspg[0,1],"2021 CMOS",ProfileHalfWidth=ProfileHalfWidth,
                            LatPlotLims=LatPlotLims,ZonePlotHalfWidth=ZonePlotHalfWidth,
                            profile=profile,clr='C1',width=1.0,band=band,smooth=smooth)
        
    LatsSCT22,OutProSCT22,OutStdSCT22,OutamfSCT22=PTG.plot_profile_L2(axsavgprof,axsamf,axsspg[0,2],axsamfspg[0,2],"2022 CMOS",ProfileHalfWidth=ProfileHalfWidth,
                        LatPlotLims=LatPlotLims,ZonePlotHalfWidth=ZonePlotHalfWidth,
                        profile=profile,clr='k',width=1.5,band=band,smooth=smooth)
    LatsVLT22,OutProVLT22,OutStdVLT22,OutamfVLT22=PTG.plot_profile_L2(axsavgprof,axsamf,axsspg[1,0],axsamfspg[1,0],"2022 VLTMUSE",ProfileHalfWidth=ProfileHalfWidth,
                        LatPlotLims=LatPlotLims,ZonePlotHalfWidth=ZonePlotHalfWidth,
                        profile=profile,clr='k',width=1.5,band=band,style='dashed',smooth=smooth)
    LatsSCT23,OutProSCT23,OutStdSCT23,OutamfSCT23=PTG.plot_profile_L2(axsavgprof,axsamf,axsspg[1,1],axsamfspg[1,1],"2023 CMOS",ProfileHalfWidth=ProfileHalfWidth,
                        LatPlotLims=LatPlotLims,ZonePlotHalfWidth=ZonePlotHalfWidth,
                        profile=profile,clr='C3',width=1.0,band=band,smooth=smooth)
    LatsVLT23,OutProVLT23,OutStdVLT23,OutamfVLT23=PTG.plot_profile_L2(axsavgprof,axsamf,axsspg[1,0],axsamfspg[1,0],"2023 VLTMUSE",ProfileHalfWidth=ProfileHalfWidth,
                        LatPlotLims=LatPlotLims,ZonePlotHalfWidth=ZonePlotHalfWidth,
                        profile=profile,clr='C4',width=1.0,band=band,smooth=smooth)
    LatsSCT24,OutProSCT24,OutStdSCT24,OutamfSCT24=PTG.plot_profile_L2(axsavgprof,axsamf,axsspg[1,0],axsamfspg[1,0],"2024 CMOS",ProfileHalfWidth=ProfileHalfWidth,
                        LatPlotLims=LatPlotLims,ZonePlotHalfWidth=ZonePlotHalfWidth,
                        profile=profile,clr='C5',width=1.0,band=band,smooth=smooth)
    
    if profile=="Meridional":
        for zb in belt:
            #print(zb,belt[zb])
            axsavgprof.fill_between([belt[zb][0],belt[zb][1]],np.array([0.,0.]),np.array([1000.,1000.]),
                                    color="0.5",alpha=0.2)
            axsavgprof.annotate(zb,xy=[np.mean(belt[zb]),0.01],ha="center")
        for zb in zone:
            axsavgprof.annotate(zb,xy=[np.mean(zone[zb]),0.05],ha="center")
    
    # Plot layout details and labeling
    
    if band=="CH4":
        y0,y1,ny=0.0,2.5,6
    elif band=="NH3":
        y0,y1,ny=0.0,1.0,6
        
    if profile=="Meridional":
        xlabel="Planetographic Latitude (deg)"
        xlims=[90-LatPlotLims[1],90-LatPlotLims[0]]
        #axsavgprof.set_xlim(90-LatPlotLims[1],90-LatPlotLims[0])
    elif profile=="Zonal":
        xlabel="Longitude from CM (deg)"
        xlims=[-ZonePlotHalfWidth,ZonePlotHalfWidth]
        #axsavgprof.set_xlim(-ZonePlotHalfWidth,ZonePlotHalfWidth)
        

    # FORMAT SUMMARY PROFILE PLOT
    axsavgprof.set_title(band+" "+profile+" Absorption Profiles")
    axsavgprof.set_xlim(xlims[0],xlims[1])
    axsavgprof.grid(linewidth=0.2)
    axsavgprof.set_ylim(y0,y1)
    axsavgprof.set_yticks(np.linspace(y0,y1,ny), minor=False)
    axsavgprof.tick_params(axis='both', which='major', labelsize=8)
    axsavgprof.set_xlabel(xlabel,fontsize=10)
    axsavgprof.set_ylabel("Equivalent Width (nm)",fontsize=10)
    axsavgprof.legend(fontsize=8,loc="upper right")
    
    figavgprof.subplots_adjust(left=0.09, right=0.97, top=0.93, bottom=0.11)
    figavgprof.savefig(path+"Profiles/output/Profile_"+band+"_"+profile+"_Absorption.png",dpi=300)
    figavgprof.savefig(path+"Profiles/output/Profile_"+band+"_"+profile+"_Absorption.pdf",dpi=300)

    #FORMAT SUMMARY AIRMASS PLOT
    axsamf.set_title(band+" "+profile+" Absorption Profiles")
    axsamf.set_xlim(1,3)
    axsamf.set_xticks(np.linspace(1,3,5), minor=False)
    axsamf.grid(linewidth=0.2)
    axsamf.set_ylim(y0,y1)
    axsamf.set_yticks(np.linspace(y0,y1,ny), minor=False)   
    axsamf.tick_params(axis='both', which='major', labelsize=8)
    axsamf.set_xlabel("One-Way Airmass Factor",fontsize=10)
    axsamf.set_ylabel("Equivalent Width (nm)",fontsize=10)
    axsamf.legend(fontsize=6,loc="upper right")
    
    figamf.subplots_adjust(left=0.09, right=0.97, top=0.93, bottom=0.11)
    figamf.savefig(path+"Profiles/output/AirMassProfile_"+band+"_"+profile+"_Absorption.png",dpi=300)

    #FORMAT ARRAY/ANNUAL PLOTS
    figspg.suptitle(band+" "+profile+" Absorption Profiles")
    
    axsspg[0,0].set_xlim(-ZonePlotHalfWidth,ZonePlotHalfWidth)
    axsspg[0,0].set_ylim(y0,y1)
    axsspg[0,0].set_yticks(np.linspace(y0,y1,ny), minor=False)
    
    figspg.delaxes(axsspg[1,2])
    figspg.subplots_adjust(left=0.05, bottom=0.08, right=0.98, top=0.90)  
    figspg.savefig(path+"Profiles/output/Spaghetti_"+band+"_"+profile+".png",dpi=300)

    
    figamfspg.suptitle(band+" "+profile+" Absorption Profiles")
    
    axsamfspg[0,0].set_xlim(1,3)
    axsamfspg[0,0].set_xticks(np.linspace(1,3,5), minor=False)
    axsamfspg[0,0].set_ylim(y0,y1)
    axsamfspg[0,0].set_yticks(np.linspace(y0,y1,ny), minor=False)
   
    for i in range(0,2):
        axsspg[i,0].set_ylabel("Equivalent Width (nm)",fontsize=8)
        axsamfspg[i,0].set_ylabel("Equivalent Width (nm)",fontsize=8)

    for i in range(0,3):
        axsspg[1,i].set_xlabel(xlabel,fontsize=8)
        axsamfspg[1,i].set_xlabel("One-Way Airmass Factor",fontsize=8)
    
    figamfspg.delaxes(axsamfspg[1,2])
    figamfspg.subplots_adjust(left=0.05, bottom=0.08, right=0.98, top=0.90)  
    figamfspg.savefig(path+"Profiles/output/AMF_scatter_"+band+"_"+profile+".png",dpi=300)

    ###########################################################################
    # Plot change in profile relative to 2022 profile
    ###########################################################################
    
    figresid,axsresid=pl.subplots(1,1,figsize=(6.0,4.0), dpi=150, facecolor="white")
    proresid=OutProSCT23-OutProSCT22
    axsresid.plot(LatsVLT22,proresid,label='2023 SCT minus 2022 SCT')
    stdresid=np.sqrt(OutStdSCT22**2+OutStdSCT23**2)
    axsresid.fill_between(LatsSCT22,proresid-stdresid,proresid+stdresid,
                          color='C0',alpha=.1)

    proresid=OutProSCT24-OutProSCT23
    axsresid.plot(LatsSCT22,proresid,label='2024 SCT minus 2023 SCT')
    stdresid=np.sqrt(OutStdSCT23**2+OutStdSCT24**2)
    axsresid.fill_between(LatsSCT22,proresid-stdresid,proresid+stdresid,
                          color='C1',alpha=.1)

    #axsresid.plot(LatsVLT22,CMOS2021EW[:,1]-OutProSCT22)
    #axsresid.plot(LatsVLT22,CMOS2020EW[:,1]-OutProSCT22)
    if profile=="Meridional":
        axsresid.set_xlim(90-LatPlotLims[1],90-LatPlotLims[0])
    if profile=="Zonal":
        axsresid.set_xlim(-ZonePlotHalfWidth,ZonePlotHalfWidth)
        
    axsresid.set_xlabel(xlabel,fontsize=10)

    axsresid.set_title("Change in "+band+" Absorption Profiles")
    axsresid.set_ylim(-0.2,0.2)
    axsresid.set_ylabel("Equivalent Width (nm)",fontsize=10)
    axsresid.set_yticks(np.linspace(-0.2,0.2,9), minor=False)
    axsresid.grid(linewidth=0.2)
    axsresid.tick_params(axis='both', which='major', labelsize=8)
    axsresid.legend(fontsize=8)

    if profile=="Meridional":
        for zb in belt:
            print(zb,belt[zb])
            axsresid.fill_between([belt[zb][0],belt[zb][1]],np.array([-2000.,-2000.]),
                                    np.array([2000.,2000.]),color="0.5",alpha=0.2)
            axsresid.annotate(zb,xy=[np.mean(belt[zb]),-0.19],ha="center")
        for zb in zone:
            axsresid.annotate(zb,xy=[np.mean(zone[zb]),-0.17],ha="center")
            
    figresid.subplots_adjust(left=0.12, bottom=0.12, right=0.98, top=0.92)  
    figresid.savefig(path+"Profiles/output/Residuals_"+band+"_"+profile+".png",dpi=300)

    if band=="NH3":
        SCT20={'Lats':LatsSCT20,'Pro':OutProSCT20,'Std':OutStdSCT20,'Amf':OutamfSCT20}
        SCT21={'Lats':LatsSCT21,'Pro':OutProSCT21,'Std':OutStdSCT21,'Amf':OutamfSCT21}
    SCT22={'Lats':LatsSCT22,'Pro':OutProSCT22,'Std':OutStdSCT22,'Amf':OutamfSCT22}
    VLT22={'Lats':LatsVLT22,'Pro':OutProVLT22,'Std':OutStdVLT22,'Amf':OutamfVLT22}
    SCT23={'Lats':LatsSCT23,'Pro':OutProSCT23,'Std':OutStdSCT23,'Amf':OutamfSCT23}
    VLT23={'Lats':LatsVLT23,'Pro':OutProVLT23,'Std':OutStdVLT23,'Amf':OutamfVLT23}

    #return(LatsVLT22,OutProVLT22,OutStdVLT22,OutamfVLT22)#,figavgprof,axsavgprof)
    
    if band=="NH3":
        return(SCT20,SCT21,SCT22,VLT22,SCT23,VLT23)#,figavgprof,axsavgprof)
    else:
        return(SCT22,VLT22,SCT23,VLT23)#,figavgprof,axsavgprof)
    
    #return(figavgprof,axsavgprof)