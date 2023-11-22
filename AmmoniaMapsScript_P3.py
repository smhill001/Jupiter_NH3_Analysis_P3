# -*- coding: utf-8 -*-
def AmmoniaMapsScript_P3(reference,Level='L3'):
    """
    Created on Sat Dec 04 22:23:20 2021
    
    reference must begin with a year, e.g., "2022 CMOS"
    @author: Steven Hill
    """
    import sys
    import numpy as np
    import AmmoniaMaps_P3 as AMP3
    import Retrieve_Jup_Atm_2022_P3 as RJA
    sys.path.append('./Services')
    import make_L2_abs_data as make_L2
    import make_L3_env_data as make_L3
    import get_batch_lists as GBL

    DataSets=GBL.get_batch_lists()

    print("int(reference[0:4])=",int(reference[0:4]))
    if int(reference[0:4])>2021:
        print("int(reference[0:4])=",int(reference[0:4]))
        
        if "SCT" or "CMOS" in reference:
            Cal="SCT-Obs-Final"
        elif "VLT" in reference:
            Cal="VLT-Obs-Final"
            
        for date in DataSets[reference]:
            if len(date)==11:
                version=date[10]
                dataset=date[0:10].replace('-','')+'UT'+version
            else:
                dataset=date.replace('-','')+'UT'
            if Level=='L3':
                #RJA.Retrieve_Jup_Atm_2022_P3(obsdate=dataset,target="Jupiter",
                #                         imagetype='Map',CalModel=Cal,
                #                         Smoothing=True,LatLims=[45,135],LonRng=45,
                #                         delta_CM2=0,showbands=False,LonSys='2')
                make_L3.make_L3_env_data(obsdate=dataset,target="Jupiter",
                                     imagetype='Map',CalModel=Cal,
                                     Smoothing=True,LonSys='2')
            elif Level=='L2':
                make_L2.make_l2_abs_data(obsdate=dataset,
                                         target="Jupiter",
                                         imagetype='Map')
                
    else:
        Lats,AvgMeridEW,StdZoneEW=AMP3.AmmoniaMaps_P3(DateSelection=DataSets[reference],
                                                      cont=False)   
        pathout="/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/"
        np.savetxt(pathout+reference+'_NH3_Meridian_EW.csv', 
                   np.transpose([Lats,AvgMeridEW,StdZoneEW]), fmt='%3.5f', 
                   delimiter=',')
