def continousmapscript(collections=['AGU2025'],LonSys='1',lats=[75,105],LonLims=[0,360],
                       localmax=False,segment=False,variance=False,ctbls=['terrain_r','Blues'],
                       proj='AGU2025'):
    """
    Created on Sun Aug 25 15:04:22 2024
    
    @author: smhil
    """
    import pylab as pl
    import time
    import get_map_collection as gmc
    import L4_Jup_Map_Plot as L4MP
    import make_L4_cont_maps as M4CM


    ts=time.time()

    maps2022=["20220810-20220812","20220828-20220901","20220904-20220905",
              "20220912-20220913","20220919-20220919","20221009-20221013",
              "20221019-20221021"]
    
    maps2023=["20230815-20230818","20230827-20230831",
              "20230905-20230906","20230922-20230929","20231005-20231006",
              "20231015-20231019","20231022-20231026","20231103-20231107",
              "20231110-20231110","20231112-20231113","20231115-20231115",
              "20231128-20231129","20231206-20231207","20231217-20231218",
              "20231229-20231229","20240129-20240202","20240229-20240301"]
    
    maps2024=["20241006-20241010","20241027-20241027","20241105-20241105",
              "20241115-20241115","20241118-20241118","20241128-20241129",
              "20241202-20241203","20241205-20241205",#"20241231-20241231",
              "20250106-20250106","20250116-20250117","20250127-20250127",
              "20250128-20250128","20250129-20250129","20250302-20250302"]
    
    maps2024GRS=["20241006-20241010","20241027-20241027","20241105-20241105",
              "20241115-20241115","20241118-20241118","20241128-20241129",
              "20241202-20241203","20241205-20241205",#"20241231-20241231",
              "20250106-20250106","20250117-20250117","20250127-20250127",
              #"20250128-20250128",
              "20250129-20250129","20250302-20250302"]
    
    maps2024SEBOutbreak=["20241105-20241105",
              "20241115-20241115","20241118-20241118","20241128-20241129",
              "20241202-20241203","20241205-20241205",#"20241231-20241231",
              "20250106-20250106","20250117-20250117",#"20250127-20250127",
              "20250128-20250128","20250129-20250129"]#,"20250302-20250302"]
    
    maps2024NTBOutbreak=[
              "20250106-20250106","20250117-20250117","20250127-20250127",
              "20250128-20250128","20250129-20250129","20250302-20250302"]
    
    maps2024NEZ=["20241105-20241105",
              "20241115-20241115","20241118-20241118","20241128-20241129",
              "20241202-20241203","20241205-20241205",
              "20250106-20250106","20250116-20250117","20250127-20250128",
              "20250128-20250129","20250302-20250302"]

    maps2024NEZG=["20241006-20241010","20241027-20241027","20241105-20241105",
              "20241115-20241115","20241118-20241118","20241128-20241128","20241129-20241129",
              "20241202-20241202","20241203-20241203","20241205-20241205",
              "20250106-20250106","20250116-20250116","20250117-20250117","20250127-20250127",
              "20250128-20250128","20250129-20250129","20250302-20250302"]
    
    maps2024NEZNovDec=["20241128-20241129","20241202-20241203","20241205-20241205"]

    maps2025=["20250919-20250919","20250925-20250926","20250930-20251001","20251016-20251017"]
    
    
    AGU2025=["20220730-20220730","20220818-20220818","20220904-20220905","20220925-20220925",
             "20231015-20231019","20231112-20231113","20231217-20231218","20240129-20240202",
             "20240919-20240919","20240925-20240929",
             "20241022-20241023","20241105-20241105","20241202-20241203",
             "20250106-20250106","20250128-20250129","20250309-20250309",
             "20250919-20250919","20250925-20250926","20250930-20251001","20251016-20251017"]
    AGU2025IRTF={"20220730-20220730":"20220725-20220726",
                 "20220818-20220818":"20220817-20220818",
                 "20220904-20220905":False,
                 "20220925-20220925":"20220929-20220929",
                 "20231015-20231019":"20231014-20231015",
                 "20231112-20231113":False,
                 "20231217-20231218":False,
                 "20240129-20240202":"20240205-20240205",
                 "20240919-20240919":"20240920-20240920",
                 "20240925-20240929":"20240920-20240920",
                 "20241022-20241023":"20241022-20241022",
                 "20241105-20241105":False,
                 "20241202-20241203":False,
                 "20250106-20250106":False,
                 "20250128-20250129":False,
                 "20250309-20250309":"20250308-20250308",
                 "20250919-20250919":"20250913-20250913",
                 "20250925-20250926":False,
                 "20250930-20251001":False,
                 "20251016-20251017":False}

    AGU2025CH4889={"20220730-20220730":"20220802-20220803",
                 "20220818-20220818":"20220816-20220816",
                 "20220904-20220905":False,
                 "20220925-20220925":"20220929-20220929",
                 "20231015-20231019":"20231017-20231017",
                 "20231112-20231113":False,
                 "20231217-20231218":False,
                 "20240129-20240202":"20240131-20240201",
                 "20240919-20240919":False,
                 "20240925-20240929":False,
                 "20241022-20241023":False,
                 "20241105-20241105":False,
                 "20241202-20241203":False,
                 "20250106-20250106":False,
                 "20250128-20250129":False,
                 "20250309-20250309":"20250310-20250310",
                 "20250919-20250919":"20250919-20250920",
                 "20250925-20250926":False,
                 "20250930-20251001":False,
                 "20251016-20251017":False}

    pathmapplots="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Studies/"+proj+"/"
    
    if int(lats[0])<90:
        latstr=str(90-lats[0])+"N"
    if int(lats[0])==90:
        latstr=str(90-lats[0])
    if int(lats[0])>90:
        latstr=str(lats[0]-90)+"S"
        
    if int(lats[1])<90:
        latstr=latstr+"-"+str(90-lats[1])+"N"
    if int(lats[1])==90:
        latstr=latstr+"-"+str(90-lats[1])
    if int(lats[1])>90:
        latstr=latstr+"-"+str(lats[1]-90)+"S"
        
    lonstr=str(LonLims[0])+"-"+str(LonLims[1])
    ###############################################################################
    # 2023
    ############################################################################### 
    for collection in collections:
        if collection=='2022':
            maps=maps2022
        if collection=='2023':
            maps=maps2023
        if collection=='2024':
            maps=maps2024
        if collection=='2024 SEB Outbreak':
            maps=maps2024SEBOutbreak
        if collection=='2024 NTB Outbreak':
            maps=maps2024NTBOutbreak
        if collection=='2024 GRS':
            maps=maps2024GRS
        if collection=='2024 NEZ':
            maps=maps2024NEZ
        if collection=='2024 NEZ Granular':
            maps=maps2024NEZG
        if collection=='2024 NEZ Nov-Dec':
            maps=maps2024NEZNovDec
        if collection=='2025':
            maps=maps2025
        if collection=='AGU2025':
            maps=AGU2025
            mapsIRTF=AGU2025IRTF
            mapsCH4889=AGU2025CH4889
        
        xfig=6.0
        yfig=6.0
        
        aspectratio=(LonLims[1]-LonLims[0])/(lats[1]-lats[0])
        aspect_ratio_map = {
                            1:     [3.0, 6.0],
                            4/3.:  [3.5, 6.0],
                            2:     [4.5, 6.0],
                            3:     [3.0], #e.g. 120x360
                            4:     [7.05, 6.0],
                            6:     [8.5, 5.0],
                            9:     [10, 4.5],
                            12:    [1.0] #e.g. 30x360
                            }
        yfactor = aspect_ratio_map[aspectratio]
      
        
        nrows=len(maps)
        fig23NH3,axs23NH3=pl.subplots(nrows,1,figsize=(xfig,yfig*yfactor[0]), dpi=150, facecolor="white",
                              sharex=True,sharey=True)   
        fig23NH3.subplots_adjust(left=0.10, bottom=0.08, right=0.98, top=0.94,
                    wspace=0.25, hspace=0.08)     
        fig23NH3.suptitle("Ammonia Abundance (ppm)")
        axs23NH3[nrows-1].set_xlabel("System "+LonSys+" Longitude (deg)",fontsize=8)        
        
        
        fig23CH4,axs23CH4=pl.subplots(nrows,1,figsize=(xfig,yfig*yfactor[0]), dpi=150, facecolor="white",
                              sharex=True,sharey=True)   
        fig23CH4.subplots_adjust(left=0.10, bottom=0.08, right=0.98, top=0.94,
                    wspace=0.25, hspace=0.08)     
        fig23CH4.suptitle("Effective Cloud-Top Pressure (mb)")
        axs23CH4[nrows-1].set_xlabel("System "+LonSys+" Longitude (deg)",fontsize=8)
        
        
        fig23RGB,axs23RGB=pl.subplots(nrows,1,figsize=(xfig,yfig*yfactor[0]), dpi=150, facecolor="white",
                              sharex=True,sharey=True)   
        fig23RGB.subplots_adjust(left=0.10, bottom=0.08, right=0.98, top=0.94,
                    wspace=0.25, hspace=0.08)     
        fig23RGB.suptitle("Visual Context")
        axs23RGB[nrows-1].set_xlabel("System "+LonSys+" Longitude (deg)",fontsize=8)
        
        
        fig23IRTF,axs23IRTF=pl.subplots(nrows,1,figsize=(xfig,yfig*yfactor[0]), dpi=150, facecolor="white",
                              sharex=True,sharey=True)   
        fig23IRTF.subplots_adjust(left=0.10, bottom=0.08, right=0.98, top=0.94,
                    wspace=0.25, hspace=0.08)     
        fig23IRTF.suptitle("5-um Radiance (Log10, uncalibrated)")
        axs23IRTF[nrows-1].set_xlabel("System "+LonSys+" Longitude (deg)",fontsize=8)
        
        fig23CH4889,axs23CH4889=pl.subplots(nrows,1,figsize=(xfig,yfig*yfactor[0]), dpi=150, facecolor="white",
                              sharex=True,sharey=True)   
        fig23CH4889.subplots_adjust(left=0.10, bottom=0.08, right=0.98, top=0.94,
                    wspace=0.25, hspace=0.08)     
        fig23CH4889.suptitle("Methane 889-nm Radiance (unscaled)")
        axs23CH4889[nrows-1].set_xlabel("System "+LonSys+" Longitude (deg)",fontsize=8)
        
        
        
        #######################################################################
        # DO SCATTER PLOTS
        #######################################################################
        figNH3vCloud,axNH3vCloud=pl.subplots(1,1,figsize=(6,6), dpi=150, facecolor="white")
        figNH3vCloud.suptitle("Cloud-Top Pressure vs Ammonia Abundance: 2022-2025",
                      fontsize=14,x=0.5,ha='center',color='k')

        axNH3vCloud.tick_params(axis='both', which='major', labelsize=10)
        axNH3vCloud.set_title("Cloud-Top Pressure vs Ammonia Abundance: 2022-2025",fontsize=12)

        figNH3vIRTF,axNH3vIRTF=pl.subplots(1,1,figsize=(6,6), dpi=150, facecolor="white")
        figNH3vIRTF.suptitle("IRTF Radiance vs Ammonia Abundance: 2022-2025",
                      fontsize=14,x=0.5,ha='center',color='k')
        axNH3vIRTF.tick_params(axis='both', which='major', labelsize=10)
        axNH3vIRTF.set_title("5-um Radiance vs Ammonia Abundance: 2022-2025",fontsize=12)

        figIRTFvPCld,axIRTFvPCld=pl.subplots(1,1,figsize=(6,6), dpi=150, facecolor="white")
        figIRTFvPCld.suptitle("Cloud-Top Pressure vs 5-um Radiance: 2022-2025",
                      fontsize=14,x=0.5,ha='center',color='k')
        axIRTFvPCld.tick_params(axis='both', which='major', labelsize=10)
        axIRTFvPCld.set_title("Cloud-Top Pressure vs 5-um Radiance: 2022-2025",fontsize=12)
    
        fig889vPCld,ax889vPCld=pl.subplots(1,1,figsize=(6,6), dpi=150, facecolor="white")
        fig889vPCld.suptitle("Cloud-Top Pressure vs 889-nm Radiance: 2022-2025",
                      fontsize=14,x=0.5,ha='center',color='k')
        ax889vPCld.tick_params(axis='both', which='major', labelsize=10)
        ax889vPCld.set_title("Cloud-Top Pressure vs 5-um Radiance: 2022-2025",fontsize=12)
    
        counter=0
        countmax=len(maps)
        cb=False
        for mp in maps:
            #if len(collection)>4:
            #    ROI,obslist,CM=NRC.NEDF_ROI_collections(collection=mp)
            #else:
            #    obslist,dummy=gmc.get_map_collection(mp)
            obslist,dummy=gmc.get_map_collection(mp)

            M4CM.make_L4_cont_maps(collection=mp,obskeys=False,
                                  FiveMicron=False,Five_obskey='',IRTFcollection=False,
                                  lats=[0,180],LonLims=[0,360],variance=True,proj='../../Data/L4 FITS (cont maps)/',
                                  bare_maps=False,cb=False,LimbCorrection=True,
                                  lonhalfwidth=45,boxcar=9)

            print("------------ collection,mp=",collection,mp)
            L4MP.L4_Jup_Map_Plot(collection=mp,IRTFcollection=mapsIRTF[mp],
                                CH4889collection=mapsCH4889[mp],LonSys=LonSys,
                                  lats=lats,LonLims=LonLims,figsz=[6.0,6.0],ROI=False,
                                  variance=False,localmax=False,segment=False,
                                  proj='AGU2025',ctbls=['terrain_r','Blues'],
                                  cont=False,bare_maps=False,cb=False,
                                  axNH3=axs23NH3[counter],
                                  axCH4=axs23CH4[counter],
                                  axRGB=axs23RGB[counter],
                                  axIRTF=axs23IRTF[counter],
                                  axCH4889=axs23CH4889[counter],
                                  axNH3vCloud=axNH3vCloud,
                                  axNH3vIRTF=axNH3vIRTF,
                                  axIRTFvPCld=axIRTFvPCld,
                                  ax889vPCld=ax889vPCld,
                                  counter=counter,countmax=countmax)
            #pl.show()
                
            counter=counter+1
            cb=False
        
    fig23NH3.savefig(pathmapplots+collection+" NH3 Stack Sys"+LonSys+" "+lonstr+" "+latstr+" map.png",dpi=300)
    fig23CH4.savefig(pathmapplots+collection+" CH4 Stack Sys"+LonSys+" "+lonstr+" "+latstr+" map.png",dpi=300)
    fig23RGB.savefig(pathmapplots+collection+" RGB Stack Sys"+LonSys+" "+lonstr+" "+latstr+" map.png",dpi=300)
    fig23IRTF.savefig(pathmapplots+collection+" IRTF Stack Sys"+LonSys+" "+lonstr+" "+latstr+" map.png",dpi=300)
    fig23CH4889.savefig(pathmapplots+collection+" 889CH4 Stack Sys"+LonSys+" "+lonstr+" "+latstr+" map.png",dpi=300)
    figNH3vCloud.savefig(pathmapplots+collection+" ContMap PCldvfNH3 Scatter.png",dpi=300)
    figNH3vIRTF.savefig(pathmapplots+collection+" ContMap IRTFvfNH3 Scatter.png",dpi=300)
    figIRTFvPCld.savefig(pathmapplots+collection+" ContMap IRTFvPCld Scatter.png",dpi=300)
    fig889vPCld.savefig(pathmapplots+collection+" ContMap 889vPCld Scatter.png",dpi=300)
    
    elapsed=time.time()-ts
    print("elapsed time=",elapsed)