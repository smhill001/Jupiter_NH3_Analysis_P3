# -*- coding: utf-8 -*-
def extract_profile(pth,filename,LonCtr='CM',LonRng=45.):
    """
    Created on Mon Jul 31 19:54:53 2023
    
    @author: smhil
    
    PURPOSE: General utility to extract profiles from fits map files. 
    
    First version only does meridional plots of varying widths
    """
    import numpy as np
    from astropy.convolution import convolve, Box1DKernel
    from astropy.io import fits
    import RetrievalLibrary as RL
    import sys
    sys.path.append('./Services')
    import read_master_calibration
    import RetrievalLibrary as RL
                 
    hdulist=fits.open(pth+filename)
    hdulist.info()
    hdr=hdulist[0].header
    mapdata=hdulist[0].data
    hdulist.close()
    ###########################################################################
    # Get ephemeris data from file name date-time string and set lat and lon lims
    ###########################################################################             
    sec=str(int(str(filename[16:17]))*6)
    obsdate=(filename[0:10]+"_"+filename[11:13]+":"+filename[13:15]+":"+sec.zfill(2))
    eph=RL.get_WINJupos_ephem(obsdate)
    Real_CM2=float(eph[1].strip())
    CM2=Real_CM2#+delta_CM2
    LonLims=[360-int(CM2+LonRng),360-int(CM2-LonRng)]
    
    patch=RL.make_patch(mapdata,[0,180],LonLims,CM2,LonRng,pad=True)

    AvgMerid=np.flip(np.mean(patch[:,:],axis=1))
    StdMerid=np.flip(np.std(patch[:,:],axis=1))
    Lats=np.linspace(-89.5,89.5,180)
    
    return(Lats,AvgMerid,StdMerid)
