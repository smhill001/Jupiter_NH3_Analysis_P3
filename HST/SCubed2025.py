def SCubed2025(obskeyHST,LonSys,ROI_ID=False,makefits=False,HST=True,SCT=False):
    
    import sys
    sys.path.append('C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/HST/')
    sys.path.append('C:/Astronomy/Projects/SAS 2021 Ammonia/Visualization-and-Analysis/')
    import L3_Jup_Map_Plot_V2 as L3MP
    import read_HST_GO as HGO

    ROI=False
    segment=False

    ###############################################################################
    ###############################################################################
    #
    # PJ 77 - HST data 2025-10-16
    #
    ###############################################################################
    ###############################################################################
    # 20251015UTa HST
    ###############################################################################
    if obskeyHST=='20251016UTa':
        collection="20251016-20251016"
        if LonSys=='3':        
            plotoptions=["contours"]
            CoLatLims=[45,135]
            LonRng=45
            CMpref=255
            #HGO.HSTGO_process_and_plot("20251016UTa",[45,135],[210,300],LonSys='3')
        if LonSys=='1':        
            #plotoptions=["contours","scatter","wave","resid","correl"]
            plotoptions=["scatter","resid","correl"]
            CoLatLims=[75,105]
            LonRng=45
            CMpref=280
            if ROI_ID=="NEZ Hot Spot A1":
                ROI={"Hot Spot":[82,84,267.0,2.0],
                     "Gyre":[85,87,270.0,3.0],
                     "Cloud Plume":[82,84,262.0,2.0],
                     "Reference":[75,77,270.0,5.0]}
                CoLatLims=[75,90]
                LonRng=15
                CMpref=270
            if ROI_ID=="NEZ Hot Spot A2":
                ROI={"Hot Spot":[81,82,297.0,2.0],
                     "Gyre":[83,86,299.0,4.0],
                     "Cloud Plume":[81,83,293.0,2.0],
                     "Reference":[75,77,295.0,5.0]}
                CoLatLims=[75,90]
                LonRng=15
                CMpref=295
            #HGO.HSTGO_process_and_plot("20251016UTa",[75,105],[235,325],LonSys='1')
    if obskeyHST=='20251016UTc':
        collection="20251017-20251017"
        if LonSys=='3':        
            plotoptions=["contours"]
            CoLatLims=[45,135]
            LonRng=30
            CMpref=30
            #HGO.HSTGO_process_and_plot("20251016UTc",[45,135],[0,60],LonSys='3')
        if LonSys=='1':        
            plotoptions=["scatter","resid","correl"]#,"contours"]
            CoLatLims=[75,105]
            LonRng=45
            CMpref=45
            
            #!!!! ROI is custom to SCubed proprosal and 20251016UTc
            #!!!! NEZ Hot Spot
            if ROI_ID=="NEZ Hot Spot":
                ROI={"Hot Spot":[82,83,14.0,2.0],
                     "Gyre":[84,86,15.0,3.0],
                     "Cloud Plume":[82,84,5.0,3.0],
                     "Reference":[76,78,15,4.0]}
                CoLatLims=[75,90]
                LonRng=15
                CMpref=15
                #segment=True
            elif ROI_ID=="SED":
            #!!!! SEZ South Equatorial Disturbance
                ROI={"Hot Spot":[98,100,61.0,2.0],
                     "Gyre":[95,97,63.0,4.0],
                     "Cloud Plume":[96,98,53.0,2.0],
                     "Reference":[104,105,60,8.0]}
                CoLatLims=[90,105]
                LonRng=15
                CMpref=60
            
    if obskeyHST=='20251016UTf':
        collection="20251016-20251017"
        if LonSys=='3':        
            plotoptions=["contours"]
            CoLatLims=[45,135]
            LonRng=45
            CMpref=90
            #HGO.HSTGO_process_and_plot("20251016UTc",[45,135],[0,60],LonSys='3')
        if LonSys=='1':        
            plotoptions=["contours","scatter","wave","resid","correl"]
            CoLatLims=[75,105]
            LonRng=45
            CMpref=105
    ###########################################################################        
    if obskeyHST=='20251120UTa':
        collection="20251116-20251116"
        if LonSys=='3':        
            plotoptions=["contours"]
            CoLatLims=[45,135]
            LonRng=30
            CMpref=135
            #HGO.HSTGO_process_and_plot("20251016UTc",[45,135],[0,60],LonSys='3')
        if LonSys=='1':        
            plotoptions=["contours","scatter","wave","resid","correl"]
            CoLatLims=[75,105]
            LonRng=30
            CMpref=65
    if obskeyHST=='20251120UTb':
        collection="20251119-20251119"
        if LonSys=='3':        
            plotoptions=["contours"]
            CoLatLims=[45,135]
            LonRng=30
            CMpref=300
            #HGO.HSTGO_process_and_plot("20251016UTc",[45,135],[0,60],LonSys='3')
        if LonSys=='1':        
            plotoptions=["contours","scatter","wave","resid","correl"]
            CoLatLims=[75,105]
            LonRng=30
            CMpref=230
    if obskeyHST=='20251120UTc':
        collection="20251119-20251119"
        if LonSys=='3':        
            plotoptions=["contours"]
            CoLatLims=[45,135]
            LonRng=30
            CMpref=50
            #HGO.HSTGO_process_and_plot("20251016UTc",[45,135],[0,60],LonSys='3')
        if LonSys=='1':        
            plotoptions=["contours","scatter","wave","resid","correl"]
            CoLatLims=[75,105]
            LonRng=30
            CMpref=330
            
    ###############################################################################
    ###############################################################################

    if makefits:
        LonLimsWest=[CMpref-LonRng,CMpref+LonRng]
        HGO.HSTGO_process_and_plot(obskeyHST,CoLatLims,LonLimsWest,LonSys=LonSys)
        
    if HST:
        L3MP.L3_Jup_Map_Plot_V2(obskey=obskeyHST, 
                           CoLatLims=CoLatLims,LonRng=LonRng,CMpref=CMpref,LonSys=LonSys, 
                           subproj='SCubed 2025/'+obskeyHST,plotoptions=plotoptions,
                           dataversion='H',smoothcont=5,segment=segment,
                           ROI=ROI)

    if SCT:
        L3MP.L3_Jup_Map_Plot_V2(obskey=collection, 
                       CoLatLims=CoLatLims,LonRng=LonRng,CMpref=CMpref,LonSys=LonSys, 
                       subproj='SCubed 2025/'+obskeyHST,plotoptions=plotoptions,
                       dataversion=2,smoothcont=0,ROI=ROI)


