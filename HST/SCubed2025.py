def SCubed2025(obskeyHST):
    
    import L4_Jup_Map_Plot as L4MP
    import read_HST_GO as HGO
    
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
        #Sys 3
        L4MP.L4_Jup_Map_Plot(collection="20251016-20251016",IRTFcollection=False, 
                        CH4889collection=False,LonSys='3',lats=[45,135],
                        LonLims=[210,300],proj='SCubed 2025/20251015UTa')
        HGO.HSTGO_process_and_plot("20251015UTa",[45,135],[210,300],LonSys='3')
        
        #Sys 1
        L4MP.L4_Jup_Map_Plot(collection="20251016-20251016",IRTFcollection=False, 
                        CH4889collection=False,LonSys='1',lats=[75,105],
                        LonLims=[235,325],proj='SCubed 2025/20251015UTa')
        HGO.HSTGO_process_and_plot("20251015UTa",[75,105],[235,325],LonSys='1')
    
    ###############################################################################
    # 20251016UTb HST
    ###############################################################################
    if obskeyHST=='20251016UTb':
        #Sys 3
        L4MP.L4_Jup_Map_Plot(collection="20251016-20251016",IRTFcollection=False, 
                        CH4889collection=False,LonSys='3',lats=[45,135],
                        LonLims=[0,60],proj='SCubed 2025/20251016UTb')
        HGO.HSTGO_process_and_plot("20251016UTb",[45,135],[0,60],LonSys='3')
        
        #Sys 1
        L4MP.L4_Jup_Map_Plot(collection="20251016-20251017",IRTFcollection=False, 
                        CH4889collection=False,LonSys='1',lats=[75,105],
                        LonLims=[0,90],proj='SCubed 2025/20251016UTb')
        HGO.HSTGO_process_and_plot("20251016UTb",[75,105],[0,90],LonSys='1')
    
    ###############################################################################
    ###############################################################################
    #
    # PJ 78 - HST data 2025-11-20
    #
    ###############################################################################
    ###############################################################################
    # 20251120UTb HST
    ###############################################################################
    if obskeyHST=='20251120UTb':
        #Sys 3
        L4MP.L4_Jup_Map_Plot(collection="20251119-20251119",IRTFcollection=False,
                        CH4889collection=False,LonSys='3',lats=[45,135],
                        LonLims=[270,330],proj='SCubed 2025/20251120UTb')
        HGO.HSTGO_process_and_plot("20251120UTb",[45,135],[270,330],LonSys='3')
        
        
        #Sys 1
        L4MP.L4_Jup_Map_Plot(collection="20251119-20251119",IRTFcollection=False,
                        CH4889collection=False,LonSys='1',lats=[75,105],
                        LonLims=[200,260],proj='SCubed 2025/20251120UTb')
        HGO.HSTGO_process_and_plot("20251120UTb",[75,105],[200,260],LonSys='1')
