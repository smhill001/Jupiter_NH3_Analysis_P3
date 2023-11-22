# -*- coding: utf-8 -*-
def read_fits_map_L2_L3(obsdate="20220919UTa",imagetype="Map",Level="L3"):
    """
    Created on Mon Nov 20 08:42:28 2023
    
    @author: smhil
    """
    import sys
    drive='c:'
    sys.path.append(drive+'/Astronomy/Python Play')
    sys.path.append(drive+'/Astronomy/Python Play/Util_P3')
    sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy')
    sys.path.append('./Services')

    import os
    #from matplotlib.pyplot import imread
    import pylab as pl
    import numpy as np
    from imageio import imwrite
    from matplotlib.pyplot import imread
    #from numpy import inf
    #from re import search
    from astropy.io import fits
    #from astropy.convolution import Gaussian2DKernel
    #from astropy.convolution import convolve
    import RetrievalLibrary as RL
    sys.path.append('./Services')
    import get_L2_abs_data as GAOD
    import make_L3_env_data
    import get_WINJupos_ephem as WJ_ephem

    target="Jupiter"
    sourcedata=obsdate+"_"+imagetype
    sourcefiles=GAOD.get_L2_abs_data()
    pathRGB='c:/Astronomy/Projects/Planets/'+target+'/Imaging Data/'+obsdate[0:10]+'/'
    pathout='C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Analysis Data/map plots/'
    diagout='C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Analysis Data/Map Plots Diagnostic/'

    RGBfile=sourcefiles[sourcedata]['RGBfile']+"_CM2_L360_MAP-BARE.png"
    if RGBfile != 'NA':
        RGB=imread(pathRGB+RGBfile)
        RGBsec=str(int(str(RGBfile[16:17]))*6)
        RGBtime=(RGBfile[0:10]+"_"+RGBfile[11:13]+":"+RGBfile[13:15]+":"+RGBsec.zfill(2))
        eph=WJ_ephem.get_WINJupos_ephem(RGBtime)
        #time.sleep(5)
        RGB_CM1=float(eph[0].strip())
        RGB_CM2=float(eph[1].strip())
        RGB_CM3=float(eph[2].strip())

    if Level=="L2":
        CH4suffix="-Jupiter_620CH4AbsMap"
        NH3suffix="-Jupiter_647NH3AbsMap"
        pathFITS='C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Analysis Data/L2 FITS/'
    elif Level=="L3":
        CH4suffix="-Jupiter_PCloud_Sys2"
        NH3suffix="-Jupiter_fNH3_Sys2"
        pathFITS='C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Analysis Data/L3 FITS/'

        
    try:
        PCloudfile=sourcefiles[sourcedata]['CH4file']+CH4suffix+\
                sourcefiles[sourcedata]['Metadata']['Variation']+".fits"
        variation=sourcefiles[sourcedata]['Metadata']['Variation']
    except:
        PCloudfile=sourcefiles[sourcedata]['CH4file']+CH4suffix+".fits"
        variation=""
    try:
        fNH3file=sourcefiles[sourcedata]['NH3file']+NH3suffix+\
                sourcefiles[sourcedata]['Metadata']['Variation']+".fits"
        variation=sourcefiles[sourcedata]['Metadata']['Variation']
    except:
        fNH3file=sourcefiles[sourcedata]['NH3file']+NH3suffix+".fits"
        variation=""

    PCloudhdulist=fits.open(pathFITS+PCloudfile)
    PCloudhdulist.info()
    PCloudhdr=PCloudhdulist[0].header
    PClouddata=PCloudhdulist[0].data
    PCloudhdulist.close()
    
    fNH3hdulist=fits.open(pathFITS+fNH3file)
    fNH3hdulist.info()
    fNH3hdr=fNH3hdulist[0].header
    fNH3data=fNH3hdulist[0].data
    fNH3hdulist.close()
    #pl.imshow(fNH3data)
        
    return(PCloudhdr,PClouddata,fNH3hdr,fNH3data,RGB,RGB_CM2)

