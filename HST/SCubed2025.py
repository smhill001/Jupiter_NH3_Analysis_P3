def SCubed2025(obskeyHST,LonSys,makefits=False,HST=True,SCT=False):
    
    import L3_Jup_Map_Plot_V2 as L3MP
    import sys
    sys.path.append('C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/HST/')
    import read_HST_GO as HGO

    #!!!! ROI is custom to SCubed proprosal and 20251016UTc
    #ROI={"Hot Spot":[82,84,15.0,3.0],
    #     "Gyre":[85,87,15.0,3.0],
    #     "Cloud Plume":[82,84,5.0,3.0],
    #     "NEB Reference":[76,78,15,4.0]}
    
    ROI=False

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
            plotoptions=["contours","scatter","wave","resid","correl"]
            CoLatLims=[75,105]
            LonRng=45
            CMpref=280
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
            plotoptions=["scatter","wave","resid","correl"]
            #CoLatLims=[75,105]
            #LonRng=45
            #CMpref=45
            #FOR SCUBED HOT-SPOT PLOT:
            CoLatLims=[75,90]
            LonRng=15
            CMpref=15
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
        HGO.HSTGO_process_and_plot(obskeyHST,CoLatLims,LonLimsWest,LonSys='1')
        
    if HST:
        L3MP.L3_Jup_Map_Plot_V2(obskey=obskeyHST, 
                           CoLatLims=CoLatLims,LonRng=LonRng,CMpref=CMpref,LonSys=LonSys, 
                           subproj='SCubed 2025/'+obskeyHST,plotoptions=plotoptions,
                           dataversion='H',smoothcont=5,
                           ROI=ROI)

    if SCT:
        L3MP.L3_Jup_Map_Plot_V2(obskey=collection, 
                       CoLatLims=CoLatLims,LonRng=LonRng,CMpref=CMpref,LonSys=LonSys, 
                       subproj='SCubed 2025/'+obskeyHST,plotoptions=plotoptions,
                       dataversion=2,smoothcont=0,ROI=ROI)


