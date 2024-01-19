def grs_zoom_maps(year):
    """
    Created on Fri Jan  5 07:32:06 2024
    
    @author: smhil
    """
    
    import Map_Jup_Atm_2022_P3 as MJP
    
    GRS2022=['20220730UTa','20220810UTa','20220818UTa','20220828UTa','20220904UTa',
              '20220919UTa','20220919UTb','20221013UTa','20221020UTa','20230113UTa']
    
    GRS2023=["20230831UTa","20230831UTb","20230905UTa","20230924UTa","20231005UTa",
             "20231015UTa","20231017UTa","20231017UTb","20231022UTa","20231103UTa",
             "20231110UTb","20231110UTc","20231113UTa","20231113UTb","20231115UTa",
             "20231115UTb"]
    
    if year==2022:
        GRS=GRS2022
        CM=25
    elif year==2023:
        GRS=GRS2023
        CM=45
        
    for obskey in GRS:
        
        MJP.Map_Jup_Atm_2022_P3(obskey=obskey,imagetype='Map', Smoothing=False,
                            LatLims=[90,135],LonRng=30,CMpref=CM,LonSys='2',
                            showbands=False,coef=[0.,0.],subproj='GRS')
