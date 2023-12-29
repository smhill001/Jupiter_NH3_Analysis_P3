def AmmoniaMapsScript_P3(reference,Level='L3'):
    """
    Created on Sat Dec 04 22:23:20 2021
    
    reference must begin with a year, e.g., "2022 CMOS"
    @author: Steven Hill
    """
    import sys
    import numpy as np
    import AmmoniaMaps_P3 as AMP3
    sys.path.append('./Services')
    import make_L2_abs_data as make_L2
    import make_L3_env_data as make_L3
    import get_batch_lists as GBL
    import Map_Jup_Atm_2022_P3 as MapJup

    ###########################################################################
    # Get the data set to batch process
    ###########################################################################         
    DataSets=GBL.get_batch_lists()

    if int(reference[0:4])>2021:       
        #######################################################################
        # If the year is 2022 or later, all processing options exist
        # Select calibration based on telescope
        #######################################################################         
        if "SCT" or "CMOS" in reference:
            Cal="SCT-Obs-Final"
        elif "VLT" in reference:
            Cal="VLT-Obs-Final"
            
        #######################################################################
        # Loop over processing options
        #######################################################################         
        Frst=True
        for date in DataSets[reference]:
            if len(date)==11:
                version=date[10]
                dataset=date[0:10].replace('-','')+'UT'+version
            else:
                dataset=date.replace('-','')+'UT'
                
            ###################################################################
            # Execute Processing from L1 to L2
            ###################################################################
            if Level=='L2':
                make_L2.make_l2_abs_data(obsdate=dataset,
                                         target="Jupiter",
                                         imagetype='Map')
                
            ###################################################################
            # Execute Processing from L2 to L3
            ###################################################################
            elif Level=='L3':
                make_L3.make_L3_env_data(obsdate=dataset,target="Jupiter",
                                     imagetype='Map',CalModel=Cal,
                                     Smoothing=False,First=Frst)
                
            ###################################################################
            # Execute Processing from L3 to map of NH3, PCloud, and context
            ###################################################################
            elif Level=='plots':
                MapJup.Map_Jup_Atm_2022_P3(obskey=dataset,imagetype='Map',
                                           Smoothing=False,LatLims=[45,135],LonRng=45,
                                           CMpref='subobs',LonSys='1',showbands=False,
                                           coef=[0.0,0.0],subproj='Dark Features')
            Frst=False
                
    else:
        #######################################################################
        # If the year is 2021 or earlier, use old code to process absorption
        # data only, e.g., L2.
        #######################################################################         
        Lats,AvgMeridEW,StdZoneEW=AMP3.AmmoniaMaps_P3(DateSelection=DataSets[reference],
                                                      cont=False)   
        pathout="/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/"
        np.savetxt(pathout+reference+'_NH3_Meridian_EW.csv', 
                   np.transpose([Lats,AvgMeridEW,StdZoneEW]), fmt='%3.5f', 
                   delimiter=',')
