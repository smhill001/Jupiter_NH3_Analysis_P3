# -*- coding: utf-8 -*-
def read_fits_map_L2_L3(obskey="20231026UTa",imagetype="Map",Level="L3",
                        LonSys='3',Smoothing=False):
    """
    Created on Mon Nov 20 08:42:28 2023
    Called by: Map_Jup_Atm_2022_P3, currently only for L3 data to plot
               maps of fNH3 and PCld
    
    @author: smhil
    """
    import sys
    drive='c:'
    sys.path.append(drive+'/Astronomy/Python Play')
    sys.path.append(drive+'/Astronomy/Python Play/Util_P3')
    sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy')
    sys.path.append('./Services')

    from matplotlib.pyplot import imread
    from astropy.io import fits
    sys.path.append('./Services')
    import get_obs_list as getlist

    import get_WINJupos_ephem as WJ_ephem
    import numpy as np

    target="Jupiter"
    sourcedata=obskey#+"_"+imagetype
    sourcefiles=getlist.get_obs_list()
    pathRGB='c:/Astronomy/Projects/Planets/'+target+'/Imaging Data/'+obskey[0:10]+'/'

    if imagetype=="Map":
        RGBfile=sourcefiles[sourcedata]['RGBfile']+"_CM2_L360_MAP-BARE.png"
    elif imagetype=="Img":
        RGBfile=sourcefiles[sourcedata]['RGBfile']+".png"

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
        CH4suffix="-Jupiter_"+imagetype+"_L2TCH4"
        NH3suffix="-Jupiter_"+imagetype+"_L2TNH3"
        pathFITS='C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Analysis Data/L2 FITS/'
    elif Level=="L3":
        CH4suffix="-Jupiter_"+imagetype+"_L3PCld_S0"
        NH3suffix="-Jupiter_"+imagetype+"_L3fNH3_S0"
        pathFITS='C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Analysis Data/L3 FITS/'
        
    try:
        PCloudfile=sourcefiles[sourcedata]['CH4file'][0:17]+CH4suffix+\
                sourcefiles[sourcedata]['Metadata']['Variation']+".fits"
        variation=sourcefiles[sourcedata]['Metadata']['Variation']
    except:
        PCloudfile=sourcefiles[sourcedata]['CH4file'][0:17]+CH4suffix+".fits"
        variation=""
    try:
        fNH3file=sourcefiles[sourcedata]['NH3file'][0:17]+NH3suffix+\
                sourcefiles[sourcedata]['Metadata']['Variation']+".fits"
        variation=sourcefiles[sourcedata]['Metadata']['Variation']
    except:
        fNH3file=sourcefiles[sourcedata]['NH3file'][0:17]+NH3suffix+".fits"
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
    fNH3sza=fNH3hdulist[1].data
    fNH3eza=fNH3hdulist[2].data
    fNH3hdulist.close()
    
    NH3time=fNH3hdr["DATE-OBS"]
    NH3_CM1=float(fNH3hdr["CM1"])
    NH3_CM2=float(fNH3hdr["CM2"])
    NH3_CM3=float(fNH3hdr["CM3"])
        
    CH4time=PCloudhdr["DATE-OBS"]
    CH4_CM1=float(PCloudhdr["CM1"])
    CH4_CM2=float(PCloudhdr["CM2"])
    CH4_CM3=float(PCloudhdr["CM3"])
    print("***********LonSys=",LonSys)
    if LonSys=='1':
        print("LonSys1=",LonSys)
        if imagetype=="Map":
            RGBroll=RGB_CM2-RGB_CM1
            RGB=np.roll(RGB,int(RGBroll),axis=1)

        NH3roll=NH3_CM3-NH3_CM1
        fNH3datar=np.roll(fNH3data,int(NH3roll),axis=1)
        
        CH4roll=CH4_CM3-CH4_CM1
        PClouddatar=np.roll(PClouddata,int(CH4roll),axis=1)
        
        fNH3szar=np.roll(fNH3sza,int(NH3roll),axis=1)
        fNH3ezar=np.roll(fNH3eza,int(NH3roll),axis=1)
        
        NH3CM=NH3_CM1
        CH4CM=CH4_CM1
        RGBCM=RGB_CM1       
        #CM=NH3_CM1
        #Real_CM=NH3_CM1

    if LonSys=='2':
        print("In LonSys==2")
        NH3roll=NH3_CM3-NH3_CM2
        fNH3datar=np.roll(fNH3data,int(NH3roll),axis=1)
        
        CH4roll=CH4_CM3-CH4_CM2
        print("&&&&&&&&&& CH4roll=",CH4roll)
        PClouddatar=np.roll(PClouddata,int(CH4roll),axis=1)
        
        fNH3szar=np.roll(fNH3sza,int(NH3roll),axis=1)
        fNH3ezar=np.roll(fNH3eza,int(NH3roll),axis=1)
        
        NH3CM=NH3_CM2
        CH4CM=CH4_CM2
        RGBCM=RGB_CM2
        
    if LonSys=='3':
        if imagetype=="Map":
            RGBroll=RGB_CM3-RGB_CM2
            RGB=np.roll(RGB,int(-RGBroll),axis=1)

        fNH3datar=fNH3data
        PClouddatar=PClouddata
        fNH3szar=fNH3sza
        fNH3ezar=fNH3eza
        
        NH3CM=NH3_CM3
        CH4CM=CH4_CM3
        RGBCM=RGB_CM3
        
    #pl.imshow(fNH3data)
    #figamf,axsamf=pl.subplots(3,1,figsize=(8.0,4.0), dpi=150, facecolor="white")
    #amfdata=(1.0/fNH3sza+1.0/fNH3eza)/2.0

    #axsamf[0].imshow(fNH3data)
    #axsamf[1].imshow(np.flipud(amfdata),vmin=-5.0,vmax=5.0)
    #axsamf[2].imshow(fNH3eza,vmin=-5.0,vmax=5.0)
    return(PCloudhdr,PClouddatar,
           fNH3hdr,fNH3datar,
           fNH3szar,fNH3ezar,
           RGB,RGBCM,RGBtime)

