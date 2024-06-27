# -*- coding: utf-8 -*-
"""
Created on Fri Jan  5 07:32:06 2024

@author: smhil
"""

import Map_Jup_Atm_P3 as MJP

batch={'1':{'LonRng':25,'CMpref':200}}

data={'1':["2023-08-18a","2023-08-27a","2023-08-27b","2023-09-05a","2023-09-05b",
         "2023-09-22a","2023-09-29a","2023-10-06a","2023-10-06b","2023-10-15b",
         "2023-11-28a""2023-12-07a","2023-12-07b","2023-12-07c","2023-12-18f",
         "2024-01-29a","2024-01-31a","2024-01-31b","2024-02-02b"]}

for i in data['1']:
    obskey=i[0:10].replace('-','')+'UT'+i[10]
    print(obskey)
    MJP.Map_Jup_Atm_P3(obskey=obskey,imagetype='Map',Smoothing=False,
                   LatLims=[75,95],LonRng=batch['1']['LonRng'],
                   CMpref=batch['1']['CMpref'],LonSys='1',
                   showbands=False,coef=[0.,0.],subproj='NEZ',figxy=[8.,2.5])
"""
MJP.Map_Jup_Atm_2022_P3(obskey="20230830UTa",imagetype='Map', Smoothing=False,
                    LatLims=[75,95],LonRng=20,CMpref=5,LonSys='1',
                    showbands=False,coef=[0.,0.],subproj='NEZ')

MJP.Map_Jup_Atm_2022_P3(obskey="20230922UTa",imagetype='Map', Smoothing=False,
                    LatLims=[75,95],LonRng=20,CMpref=205,LonSys='1',
                    showbands=False,coef=[0.,0.],subproj='NEZ')

MJP.Map_Jup_Atm_2022_P3(obskey="20230922UTa",imagetype='Map', Smoothing=False,
                    LatLims=[75,95],LonRng=20,CMpref=170,LonSys='1',
                    showbands=False,coef=[0.,0.],subproj='NEZ')

MJP.Map_Jup_Atm_2022_P3(obskey="20231113UTa",imagetype='Map', Smoothing=False,
                    LatLims=[75,95],LonRng=20,CMpref=325,LonSys='1',
                    showbands=False,coef=[0.,0.],subproj='NEZ')

MJP.Map_Jup_Atm_2022_P3(obskey="20231113UTb",imagetype='Map', Smoothing=False,
                    LatLims=[75,95],LonRng=20,CMpref=325,LonSys='1',
                    showbands=False,coef=[0.,0.],subproj='NEZ')

MJP.Map_Jup_Atm_2022_P3(obskey="20231113UTa",imagetype='Map', Smoothing=False,
                    LatLims=[75,95],LonRng=20,CMpref=345,LonSys='1',
                    showbands=False,coef=[0.,0.],subproj='NEZ')

MJP.Map_Jup_Atm_2022_P3(obskey="20231113UTb",imagetype='Map', Smoothing=False,
                    LatLims=[75,95],LonRng=20,CMpref=345,LonSys='1',
                    showbands=False,coef=[0.,0.],subproj='NEZ')



MJP.Map_Jup_Atm_2022_P3(obskey="20231129UTa",imagetype='Map', Smoothing=False,
                    LatLims=[75,95],LonRng=20,CMpref=45,LonSys='1',
                    showbands=False,coef=[0.,0.],subproj='NEZ')
MJP.Map_Jup_Atm_2022_P3(obskey="20231129UTb",imagetype='Map', Smoothing=False,
                    LatLims=[75,95],LonRng=20,CMpref=45,LonSys='1',
                    showbands=False,coef=[0.,0.],subproj='NEZ')
MJP.Map_Jup_Atm_2022_P3(obskey="20231129UTc",imagetype='Map', Smoothing=False,
                    LatLims=[75,95],LonRng=20,CMpref=45,LonSys='1',
                    showbands=False,coef=[0.,0.],subproj='NEZ')
MJP.Map_Jup_Atm_2022_P3(obskey="20231206UTa",imagetype='Map', Smoothing=False,
                    LatLims=[75,95],LonRng=20,CMpref=45,LonSys='1',
                    showbands=False,coef=[0.,0.],subproj='NEZ')



MJP.Map_Jup_Atm_2022_P3(obskey="20231217UTa",imagetype='Map', Smoothing=False,
                    LatLims=[75,95],LonRng=20,CMpref=250,LonSys='1',
                    showbands=False,coef=[0.,0.],subproj='NEZ')
"""