# -*- coding: utf-8 -*-
"""
Created on Sun Jan 26 11:25:43 2025

@author: smhil
"""

def pipe_one_obs(obskeys,target="Jupiter",close=False):
    import sys
    import time
    ts=time.time()
    drive='c:'
    sys.path.append(drive+'/Astronomy/Python Play')
    sys.path.append(drive+'/Astronomy/Python Play/Util_P3')
    sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy')
    sys.path.append('./Services')

    import make_L2_abs_data as ML2
    import make_L3_env_data as ML3
    import image_array_new as IAN
    import L3_Jup_Map_Plot as L3P
    import pylab as pl

    
    for obskey in obskeys:
    
        ML2.make_l2_abs_data(obsdate=obskey,target=target,imagetype='Img',
                             mask=True,CH4shift=0.0,NH3shift=0.0)
        ML2.make_l2_abs_data(obsdate=obskey,target=target,imagetype='Map',
                             mask=True,CH4shift=0.0,NH3shift=0.0)
        
        ML3.make_L3_env_data(obsdate=obskey,target=target,imagetype='Img',
                             Smoothing=False,First=True)
        ML3.make_L3_env_data(obsdate=obskey,target=target,imagetype='Map',
                             Smoothing=False,First=True)
        
        IAN.image_array_new(obsdate=obskey,target=target,imagetype='Img',
                            contour=False)
        
        L3P.L3_Jup_Map_Plot(obskey=obskey,imagetype='Map',target=target,
                        Smoothing=False,LatLims=[45,135],LonRng=45,CMpref='subobs',
                        LonSys='2',showbands=False,coef=[0.,0.],subproj='',
                        figxy=[8.0,4.0],FiveMicron=False,ROI=False)
        if close:
            pl.close('all')
        
    elapsed=time.time()-ts
    print("Elapsed Time ",elapsed," seconds")