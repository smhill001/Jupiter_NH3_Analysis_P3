def Study_Change_2022():
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
            # Execute Processing from L3 to map of NH3, PCloud, and context
            ###################################################################
            MapJup.Map_Jup_Atm_2022_P3(obskey=dataset,imagetype='Map',
                                       Smoothing=False,LatLims=[45,135],LonRng=45,
                                       delta_CM2=0,LonSys='2',showbands=False,
                                       coef=[0.65,0.25])
            Frst=False
