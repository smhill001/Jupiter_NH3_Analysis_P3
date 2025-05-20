def continousmapscript(collections=['2024'],LonSys='1',lats=[75,105],LonLims=[0,360],
                       localmax=False,variance=False):
    """
    Created on Sun Aug 25 15:04:22 2024
    
    @author: smhil
    """
    import pylab as pl
    import MakeContiguousMap as MCM
    import time
    import NEDF_ROI_collections as NRC
    import get_map_collection as gmc


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
              "20250106-20250106","20250117-20250117","20250127-20250127",
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
              "20241115-20241115","20241118-20241118","20241128-20241128","20241129-20241129",
              "20241202-20241202","20241203-20241203","20241205-20241205",
              "20250106-20250106","20250117-20250117","20250127-20250127",
              "20250128-20250128","20250129-20250129","20250302-20250302"]

    pathmapplots="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Studies/maps/"
    
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
            
        nrows=len(maps)
        fig23NH3,axs23NH3=pl.subplots(nrows,1,figsize=(6.0,6.0), dpi=150, facecolor="white",
                              sharex=True,sharey=True)   
        fig23NH3.subplots_adjust(left=0.10, bottom=0.08, right=0.98, top=0.94,
                    wspace=0.25, hspace=0.08)     
        fig23NH3.suptitle("Ammonia Abundance (ppm)")
        axs23NH3[nrows-1].set_xlabel("System "+LonSys+" Longitude (deg)",fontsize=8)        
        
        fig23CH4,axs23CH4=pl.subplots(nrows,1,figsize=(6.0,6.0), dpi=150, facecolor="white",
                              sharex=True,sharey=True)   
        fig23CH4.subplots_adjust(left=0.10, bottom=0.08, right=0.98, top=0.94,
                    wspace=0.25, hspace=0.08)     
        fig23CH4.suptitle("Effective Cloud-Top Pressure (mb)")
        axs23CH4[nrows-1].set_xlabel("System "+LonSys+" Longitude (deg)",fontsize=8)
        
        fig23RGB,axs23RGB=pl.subplots(nrows,1,figsize=(6.0,6.0), dpi=150, facecolor="white",
                              sharex=True,sharey=True)   
        fig23RGB.subplots_adjust(left=0.10, bottom=0.08, right=0.98, top=0.94,
                    wspace=0.25, hspace=0.08)     
        fig23RGB.suptitle("Visual Context")
        axs23RGB[nrows-1].set_xlabel("System "+LonSys+" Longitude (deg)",fontsize=8)
    
    
        counter=0
        for mp in maps:
            #if len(collection)>4:
            #    ROI,obslist,CM=NRC.NEDF_ROI_collections(collection=mp)
            #else:
            #    obslist,dummy=gmc.get_map_collection(mp)
            obslist,dummy=gmc.get_map_collection(mp)


            MCM.MakeContiguousMap(axs23NH3[counter],axs23CH4[counter],axs23RGB[counter],
                                                                       obslist,
                                                                       collection=mp,
                                                                       LonSys=LonSys,
                                                                       lats=lats,
                                                                       LonLims=LonLims,
                                                                       localmax=localmax,
                                                                       variance=variance)
            #pl.show()
            counter=counter+1
    
        fig23NH3.savefig(pathmapplots+collection+" NH3 Stack Sys"+LonSys+" "+lonstr+" "+latstr+"_map.png",dpi=300)
        fig23CH4.savefig(pathmapplots+collection+" CH4 Stack Sys"+LonSys+" "+lonstr+" "+latstr+"_map.png",dpi=300)
        fig23RGB.savefig(pathmapplots+collection+" RGB Stack Sys"+LonSys+" "+lonstr+" "+latstr+"_map.png",dpi=300)
    
    elapsed=time.time()-ts
    print("elapsed time=",elapsed)