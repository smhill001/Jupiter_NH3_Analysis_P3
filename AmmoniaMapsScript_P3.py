# -*- coding: utf-8 -*-
"""
Created on Sat Dec 04 22:23:20 2021

@author: Steven Hill
"""
def AmmoniaMapsScript_P3(reference):
    import numpy as np
    import AmmoniaMaps_P3 as AMP3
    import Retrieve_Jup_Atm_2022_P3 as RJA
    DataSets={"2020-21 CMOS":["2020-07-20","2020-07-29","2020-07-30","2020-07-31",
               "2021-06-22","2021-07-08","2021-07-19","2021-07-20",
               "2021-09-10","2021-09-13","2021-09-15","2021-09-19",
               "2021-09-20","2021-09-23","2021-09-26","2021-09-27",
               "2021-10-17","2021-10-19","2021-10-22","2021-11-30",
               "2021-12-02","2021-12-03"],
    "2020 CMOS":["2020-07-20","2020-07-29","2020-07-30","2020-07-31"],
    "2020 CCD":["2020-09-02","2020-09-03","2020-09-04","2020-09-13",
                  "2020-09-15"],
    "2021 CMOS":["2021-06-22","2021-07-08","2021-07-19","2021-07-20",
                   "2021-09-10","2021-09-13","2021-09-15","2021-09-19",
                   "2021-09-20","2021-09-23","2021-09-26","2021-09-27",
                   "2021-10-17","2021-10-19","2021-10-22","2021-11-30",
                   "2021-12-02","2021-12-03"],
    "2022 CMOS":["2022-08-10","2022-08-12","2022-08-18","2022-08-28",
                 "2022-08-30","2022-09-01","2022-09-04","2022-09-05",
                 "2022-09-12","2022-09-13","2022-09-19a","2022-09-19b",
                 "2022-09-25","2022-10-09a","2022-10-09b","2022-10-13",
                 "2022-10-19","2022-10-20","2022-10-21","2023-01-13"],
    "2023 CMOS":["2023-08-15","2023-08-16","2023-08-17","2023-08-18",
                 "2023-08-27a","2023-08-27b","2023-08-30","2023-08-31a",
                 "2023-08-31b","2023-09-05a","2023-09-05b","2023-09-06a",
                 "2023-09-06b"]}

    print("int(reference[0:4])=",int(reference[0:4]))
    if int(reference[0:4])>2021:
        print("int(reference[0:4])=",int(reference[0:4]))

        for date in DataSets[reference]:
            if len(date)==11:
                version=date[10]
                dataset=date[0:10].replace('-','')+'UT'+version
            else:
                dataset=date.replace('-','')+'UT'
            RJA.Retrieve_Jup_Atm_2022_P3(obsdate=dataset,target="Jupiter",
                                     imagetype='Map',CalModel='SCT-Obs-Final',
                                     Smoothing=True,LatLims=[45,135],LonRng=45,
                                     delta_CM2=0,showbands=False,LonSys='2')
    else:
        Lats,AvgMeridEW,StdZoneEW=AMP3.AmmoniaMaps_P3(DateSelection=DataSets[reference],
                                                      cont=False)   
        pathout="/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/"
        np.savetxt(pathout+reference+'_NH3_Meridian_EW.csv', 
                   np.transpose([Lats,AvgMeridEW,StdZoneEW]), fmt='%3.5f', 
                   delimiter=',')
