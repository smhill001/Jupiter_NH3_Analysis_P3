# -*- coding: utf-8 -*-
"""
Created on Sat Dec 04 22:23:20 2021

@author: Steven Hill
"""
def AmmoniaMapsScript_P3(reference):
    import numpy as np
    import AmmoniaMaps_P3 as AMP3
    
    DataSets={"2020-21 CMOS":["2020-07-20","2020-07-29","2020-07-30","2020-07-31",
               "2021-06-22","2021-07-08","2021-07-19","2021-07-20",
               "2021-09-10","2021-09-13","2021-09-15","2021-09-19",
               "2021-09-20","2021-09-23","2021-09-26","2021-09-27",
               "2021-10-17","2021-10-19","2021-10-22","2021-11-30",
               "2021-12-02","2021-12-03"],
    "2020 CMOS":["2020-07-20","2020-07-29","2020-07-30","2020-07-31"],
    "2021 CMOS":["2021-06-22","2021-07-08","2021-07-19","2021-07-20",
                   "2021-09-10","2021-09-13","2021-09-15","2021-09-19",
                   "2021-09-20","2021-09-23","2021-09-26","2021-09-27",
                   "2021-10-17","2021-10-19","2021-10-22","2021-11-30",
                   "2021-12-02","2021-12-03"],
    "2020 CCD":["2020-09-02","2020-09-03","2020-09-04","2020-09-13",
                  "2020-09-15"]}

    Lats,AvgMeridEW,StdZoneEW=AMP3.AmmoniaMaps_P3(DateSelection=DataSets[reference],cont=False)
    
    pathout="/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/"
    np.savetxt(pathout+reference+'_NH3_Meridian_EW.csv', np.transpose([Lats,AvgMeridEW,StdZoneEW]), fmt='%3.5f', delimiter=',')
