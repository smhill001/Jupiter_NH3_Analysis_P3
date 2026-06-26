# -*- coding: utf-8 -*-

def SCubed2025(obskeyHST,LonSys,ROI_ID=False,makefits=False,HST=True,SCT=False):
    
    import sys
    import socket
    hostname = socket.gethostname()
    from config_VA import Host_path
    #sys.path.append('C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/HST/')
    #sys.path.append('C:/Astronomy/Projects/SAS 2021 Ammonia/Visualization-and-Analysis/')
    #sys.path.append()
    import L3_Jup_Map_Plot_V2 as L3MP
    #import read_HST_GO as HGO

    ROI=False
    segment=False
    compare=False

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
            if ROI_ID=="-NEZ-A1":
                ROI={"Hot Spot":[82,84,267.0,2.0],
                     "Gyre":[85,87,270.0,3.0],
                     "Cloud Plume":[82,84,262.0,2.0],
                     "Reference":[75,77,270.0,5.0]}
                CoLatLims=[75,90]
                LonRng=15
                CMpref=270
            if ROI_ID=="-NEZ-A2":
                ROI={"Hot Spot":[81,82,297.0,2.0],
                     "Gyre":[83,86,299.0,4.0],
                     "Cloud Plume":[81,83,293.0,2.0],
                     "Reference":[75,77,295.0,5.0]}
                CoLatLims=[75,90]
                LonRng=15
                CMpref=295
            #HGO.HSTGO_process_and_plot("20251016UTa",[75,105],[235,325],LonSys='1')
    ###############################################################################
    # 20251015UTc HST
    ###############################################################################
    if obskeyHST=='20251016UTc':
        collection="20251017-20251017"
        if LonSys=='3':        
            plotoptions=[""]#"contours"]
            CoLatLims=[45,135]
            LonRng=30
            CMpref=30
            #HGO.HSTGO_process_and_plot("20251016UTc",[45,135],[0,60],LonSys='3')
        if LonSys=='1':        
            plotoptions=["scatter","resid","correl"]#,"contours","wave"]
            CoLatLims=[75,105]
            LonRng=45
            CMpref=45
            segment=False
            #!!!! ROI is custom to SCubed proprosal and 20251016UTc
            #!!!! NEZ Hot Spot
            if ROI_ID=="-NEZ-East":
                ROI={"Hot Spot":[82,83,14.0,2.0],
                     "Gyre":[84,86,15.0,3.0],
                     "Cloud Plume":[82,84,5.0,3.0],
                     "Reference":[76,78,15,4.0]}
                CoLatLims=[75,90]
                LonRng=15
                CMpref=15
                segment=False
                compare=False
            if ROI_ID=="-NEZ-West":
                ROI={"Hot Spot":[83,84,53.0,2.0],
                     "Gyre":[84,85,53.0,3.0],
                     "Cloud Plume":[84,86,34.0,3.0],
                     "Reference":[76,78,45,4.0]}
                CoLatLims=[75,90]
                LonRng=15
                CMpref=45
                segment=False
            elif ROI_ID=="-SED":
            #!!!! SEZ South Equatorial Disturbance
                ROI={"Hot Spot":[98,100,61.0,2.0],
                     "Gyre":[95,97,63.0,4.0],
                     "Cloud Plume":[96,98,53.0,2.0],
                     "Reference":[104,105,60,8.0]}
                CoLatLims=[90,105]
                LonRng=15
                CMpref=60
    ###############################################################################
    # 20251015UTf HST
    ###############################################################################

    if obskeyHST=='20251016UTf':
        collection="20251016-20251017"
        if LonSys=='3':        
            plotoptions=["surface"]#"contours"]
            CoLatLims=[45,135]
            LonRng=45
            CMpref=90
            #HGO.HSTGO_process_and_plot("20251016UTc",[45,135],[0,60],LonSys='3')
        if LonSys=='1':        
            plotoptions=["contours","scatter","wave","resid","correl"]
            CoLatLims=[75,105]
            LonRng=45
            CMpref=105
            if ROI_ID=="-NEZ-East":
                plotoptions=["scatter","resid","correl"]
                ROI={"Hot Spot":[83,85,96.0,3.0],
                     "Gyre":[85,87,98.0,4.0],
                     "Cloud Plume":[82,84,79.0,3.0],
                     "Reference":[76,78,90.0,4.0]}
                CoLatLims=[75,90]
                LonRng=15
                CMpref=90
                segment=True
            if ROI_ID=="SED": #Doesn't work because of limited eastern data boundary, but that may be fixable
                plotoptions=["scatter","resid","correl"]
                ROI={"Hot Spot":[98,100,61.0,2.0],
                     "Gyre":[95,97,63.0,4.0],
                     "Cloud Plume":[96,98,53.0,2.0],
                     "Reference":[104,105,60,8.0]}
                CoLatLims=[90,105]
                LonRng=15
                CMpref=60
                segment=True
            if ROI_ID=="-NEZ-West":
                plotoptions=["scatter","resid","correl"]
                ROI={"Hot Spot":[82,84,125.0,3.0],
                     "Gyre":[84,86,126.0,3.0],
                     "Cloud Plume":[82,84,113.0,3.0],
                     "Reference":[76,78,120.0,4.0]}
                CoLatLims=[75,90]
                LonRng=15
                CMpref=120
                segment=True

    ###############################################################################
    ###############################################################################
    #
    # PJ 78 - HST data 2025-11-20
    #
    ###############################################################################
    ###############################################################################
    # 20251120UTa HST
    ###############################################################################  
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
    ###############################################################################
    # 20251120UTb HST
    ###############################################################################  
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
    ###############################################################################
    # 20251120UTc HST
    ###############################################################################  
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
    #
    # PJ 79 - HST data 2025-12-21 through 2025-12-23
    #
    ###############################################################################
    ###############################################################################
    # 20251221UTa HST
    ###############################################################################  
    if obskeyHST=='20251221UTa':
        collection="20251222-20251222"
        if LonSys=='3':        
            plotoptions=["contours"]
            CoLatLims=[45,135]
            LonRng=30
            CMpref=95
            #HGO.HSTGO_process_and_plot("20251016UTc",[45,135],[0,60],LonSys='3')
        if LonSys=='1':        
            plotoptions=["contours","scatter","wave","resid","correl"]
            CoLatLims=[75,105]
            LonRng=30
            CMpref=250

    ###############################################################################
    # 20251221UTb HST
    ###############################################################################  
    if obskeyHST=='20251221UTb':
        collection="20251222-20251222"
        if LonSys=='3':        
            plotoptions=["contours"]
            CoLatLims=[45,135]
            LonRng=30
            CMpref=120
            #HGO.HSTGO_process_and_plot("20251016UTc",[45,135],[0,60],LonSys='3')
        if LonSys=='1':        
            plotoptions=["contours","scatter","wave","resid","correl"]
            CoLatLims=[75,105]
            LonRng=30
            CMpref=285

    ###############################################################################
    # 20251221UTc HST
    ###############################################################################  
    if obskeyHST=='20251221UTc':
        collection="20251222-20251222"
        if LonSys=='3':        
            plotoptions=["contours"]
            CoLatLims=[45,135]
            LonRng=30
            CMpref=310
            #HGO.HSTGO_process_and_plot("20251016UTc",[45,135],[0,60],LonSys='3')
        if LonSys=='1':        
            plotoptions=["contours","scatter","wave","resid","correl"]
            CoLatLims=[75,105]
            LonRng=30
            CMpref=125

    ###############################################################################
    ###############################################################################
    ###############################################################################

    if makefits:
        LonLimsWest=[CMpref-LonRng,CMpref+LonRng]
        HGO.HSTGO_process_and_plot(obskeyHST,CoLatLims,LonLimsWest,LonSys=LonSys)
        
    if HST:
        L3MP.L3_Jup_Map_Plot_V2(obskey=obskeyHST, 
                                CoLatLims=CoLatLims,LonRng=LonRng,CMpref=CMpref,LonSys=LonSys, 
                                subproj='SCubed 2025/'+obskeyHST,plotoptions=plotoptions,
                                dataversion='H',smoothcont=5,segment=segment,ROI_ID=ROI_ID,
                                ROI=ROI,compare=compare)

    if SCT:
        L3MP.L3_Jup_Map_Plot_V2(obskey=collection, 
                                CoLatLims=CoLatLims,LonRng=LonRng,CMpref=CMpref,LonSys=LonSys, 
                                subproj='SCubed 2025/'+obskeyHST,plotoptions=plotoptions,
                                dataversion=2,smoothcont=0,segment=segment,ROI_ID=ROI_ID,
                                ROI=ROI,compare=compare)
