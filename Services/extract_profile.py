# -*- coding: utf-8 -*-
def extract_profile(pth,filename,LonCtr='CM',ProfileHalfWidth=45.,
                    profile="Meridional"):
    """
    Created on Mon Jul 31 19:54:53 2023
    
    @author: smhil
    
    PURPOSE: General utility to extract profiles from fits map files. 
    
    First version only extracts meridional profiles of varying widths from
        FITS data files
    Second version extracts meridional and zonal profiles of various widths
        centered either on the Sys. 2 CM or on the equator respectively
        11/18/2023
    Third version suggested improvements: 
        1) Flexible centering in latitude and longitude profiles
        2) Return FITS header information from files
    """
    import numpy as np
    from astropy.convolution import convolve, Box1DKernel
    from astropy.io import fits
    import RetrievalLibrary as RL
    import sys
    sys.path.append('./Services')
    import read_master_calibration
    import get_WINJupos_ephem as get_ephem
                 
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
    eph=get_ephem.get_WINJupos_ephem(obsdate)
    Real_CM2=float(eph[1].strip())
    CM2=Real_CM2+0.00001 #the addition of 0.00001 eliminates a literal edge case

    if profile=="Meridional":
        LonLims=[360-int(CM2+ProfileHalfWidth),360-int(CM2-ProfileHalfWidth)]
        patch=RL.make_patch(mapdata,[0,180],LonLims,CM2,ProfileHalfWidth,pad=True)
        AvgProf=np.flip(np.mean(patch[:,:],axis=1))
        StdProf=np.flip(np.std(patch[:,:],axis=1))
        Coords=np.linspace(-89.5,89.5,180)
    
    if profile=="Zonal":
        LonLims=[360-int(CM2+180),360-int(CM2-180)]
        LatLims=[90-ProfileHalfWidth,90+ProfileHalfWidth]
        patch=RL.make_patch(mapdata,LatLims,LonLims,CM2,180,pad=True)
        print("LonLims=",LonLims)
        print("patch.shape=",patch.shape)
        AvgProf=np.flip(np.mean(patch[:,:],axis=0))
        StdProf=np.flip(np.std(patch[:,:],axis=0))
        Coords=np.arange(-180,180,1)#+0.5

    print("AvgProfile.shape,Coords.shape=",AvgProf.shape,Coords.shape)

    return(Coords,AvgProf,StdProf,CM2)
