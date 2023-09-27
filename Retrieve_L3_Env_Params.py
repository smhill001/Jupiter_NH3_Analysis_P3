# -*- coding: utf-8 -*-
def Retrieve_L3_Env_Params(obsdate="20220919UTa",target="Jupiter",
                             imagetype='Map',CalModel='VLT-Obs-Final',
                             Smoothing=True,LonSys='2'):
    """
    Created on Tue Aug 22 11:01:44 2023
    
    CALLED by AmmoniaMapsScript_P3.py->Retrieve_Jup_Atm_2022_P3
    
    @author: smhil
    """
    import sys
    drive='c:'
    sys.path.append(drive+'/Astronomy/Python Play')
    sys.path.append(drive+'/Astronomy/Python Play/Util_P3')
    sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy')
    sys.path.append('./Services')

    import os
    from matplotlib.pyplot import imread
    import numpy as np
    from astropy.io import fits
    from astropy.convolution import Gaussian2DKernel
    from astropy.convolution import convolve
    #import RetrievalLibrary as RL
    import read_master_calibration
    sys.path.append('./Services')
    import get_L2_abs_data as GAOD
    import time
    import get_WINJupos_ephem as WJ_ephem

    # Retrieve_Jup_Atm_2022_P3(obsdate="20221019UT",target="Jupiter")   
    
    amagat=2.69e24 #Lodschmits number. (cm-2)
    gravity=2228.0 #cm/s^2
    mean_mol_wt=3.85e-24 #cgs or SI !!!!WHY IS THIS 2.2 FOR MENDIKOA!!!!!
    fCH4=1.81e-3
    STP=1.01e6  #dyne/cm^2 [(g-cm/s^2)/cm^2]`
    
    ###########################################################################
    #  DATA FILES AND METADATA DICTIONARY
    ###########################################################################
    sourcedata=obsdate+"_"+imagetype
    #print("obsdate, sourcedata=",obsdate,sourcedata)
    sourcefiles=GAOD.get_L2_abs_data()
    calibration,K_eff=read_master_calibration.read_master_calibration()
    
    ###########################################################################
    # OPEN AND READ DATA FILES (FITS MAPS of NH3 and CH4 Absorption (Transmission?))
    ###########################################################################   
    print("sourcefiles[sourcedata]['Metadata']=",sourcefiles[sourcedata]['Metadata'])          
    print("sourcefiles[sourcedata]['Metadata']['Telescope']=",sourcefiles[sourcedata]['Metadata']['Telescope'])          
    K_eff_CH4620=K_eff['CH4_620'][sourcefiles[sourcedata]['Metadata']['Telescope']]#0.428
    K_eff_NH3647=K_eff['NH3_647'][sourcefiles[sourcedata]['Metadata']['Telescope']]#2.964

    pathRGB='c:/Astronomy/Projects/Planets/'+target+'/Imaging Data/'+obsdate[0:10]+'/'
    pathFITS='C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Analysis Data/L2 FITS/'
    pathout='C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Analysis Data/L3 FITS/'

    for a in['L2Cal']:#,'L3Cal']:
        #CH4file=sourcefiles[sourcedata]['CH4file']+"-Jupiter_620CH4AbsMap.fits"
        #NH3file=sourcefiles[sourcedata]['NH3file']+"-Jupiter_647NH3AbsMap.fits"
        RGBfile=sourcefiles[sourcedata]['RGBfile']+"_CM2_L360_MAP-BARE.png"

        try:
            CH4file=sourcefiles[sourcedata]['CH4file']+"-Jupiter_620CH4AbsMap"+\
                    sourcefiles[sourcedata]['Metadata']['Variation']+".fits"
            variation=sourcefiles[sourcedata]['Metadata']['Variation']
        except:
            CH4file=sourcefiles[sourcedata]['CH4file']+"-Jupiter_620CH4AbsMap.fits"
            variation=""
        try:
            NH3file=sourcefiles[sourcedata]['NH3file']+"-Jupiter_647NH3AbsMap"+\
                    sourcefiles[sourcedata]['Metadata']['Variation']+".fits"
            variation=sourcefiles[sourcedata]['Metadata']['Variation']
        except:
            NH3file=sourcefiles[sourcedata]['NH3file']+"-Jupiter_647NH3AbsMap.fits"
            variation=""
        
        print("########## pathFITS+CH4file=",pathFITS+CH4file)
        CH4hdulist=fits.open(pathFITS+CH4file)
        CH4hdulist.info()
        CH4hdr=CH4hdulist[0].header
        CH4data=CH4hdulist[0].data
        CH4hdulist.close()
        
        NH3hdulist=fits.open(pathFITS+NH3file)
        NH3hdulist.info()
        NH3hdr=NH3hdulist[0].header
        NH3data=NH3hdulist[0].data
        NH3hdulist.close()
    
        if RGBfile != 'NA':
            RGB=imread(pathRGB+RGBfile)
    
        ###########################################################################
        # Get ephemeris data from file name date-time string and set lat and lon lims
        ###########################################################################             
        RGBsec=str(int(str(RGBfile[16:17]))*6)
        RGBtime=(RGBfile[0:10]+"_"+RGBfile[11:13]+":"+RGBfile[13:15]+":"+RGBsec.zfill(2))
        eph=WJ_ephem.get_WINJupos_ephem(RGBtime)
        #time.sleep(5)
        RGB_CM1=float(eph[0].strip())
        RGB_CM2=float(eph[1].strip())
        RGB_CM3=float(eph[2].strip())
        
        NH3sec=str(int(str(NH3file[16:17]))*6)
        NH3time=(NH3file[0:10]+"_"+NH3file[11:13]+":"+NH3file[13:15]+":"+NH3sec.zfill(2))
        eph=WJ_ephem.get_WINJupos_ephem(NH3time)
        #time.sleep(5)
        NH3_CM1=float(eph[0].strip())
        NH3_CM2=float(eph[1].strip())
        NH3_CM3=float(eph[2].strip())
        
        CH4sec=str(int(str(CH4file[16:17]))*6)
        CH4time=(CH4file[0:10]+"_"+CH4file[11:13]+":"+CH4file[13:15]+":"+CH4sec.zfill(2))
        eph=WJ_ephem.get_WINJupos_ephem(CH4time)
        #time.sleep(5)
        CH4_CM1=float(eph[0].strip())
        CH4_CM2=float(eph[1].strip())
        CH4_CM3=float(eph[2].strip())
            
        if LonSys=='1':
            RGBroll=RGB_CM2-RGB_CM1
            RGB=np.roll(RGB,int(RGBroll),axis=1)
            NH3CM=NH3_CM2
            NH3roll=NH3_CM2-NH3_CM1
            NH3data=np.roll(NH3data,int(NH3roll),axis=1)
            CH4roll=CH4_CM2-CH4_CM1
            CH4data=np.roll(CH4data,int(CH4roll),axis=1)
            NH3CM=NH3_CM1
            CH4CM=CH4_CM1
            RGBCM=RGB_CM1
    
            #CM=NH3_CM1
            #Real_CM=NH3_CM1
        elif LonSys=='2':
            print("LonSys=",LonSys)
            NH3CM=NH3_CM2
            CH4CM=CH4_CM2
            RGBCM=RGB_CM2
            #CM=Real_CM2
            #Real_CM=Real_CM2
        elif LonSys=='3':
            RGBroll=RGB_CM2-RGB_CM3
            RGB=np.roll(RGB,int(RGBroll),axis=1)
            NH3roll=NH3_CM2-NH3_CM3
            NH3data=np.roll(NH3data,int(NH3roll),axis=1)
            CH4roll=CH4_CM2-CH4_CM3
            CH4data=np.roll(CH4data,int(CH4roll),axis=1)
            NH3CM=NH3_CM3
            CH4CM=CH4_CM3
            RGBCM=RGB_CM3
    
            #roll=Real_CM2-Real_CM3
            #CM=Real_CM3
            #Real_CM=Real_CM3
            #CH4data=np.roll(CH4data,int(roll),axis=1)
            #NH3data=np.roll(NH3data,int(roll),axis=1)
            #RGB=np.roll(RGB,int(roll),axis=1)
           
        ###########################################################################
        # Set calibration and compute physical parameters
        ###########################################################################             
        #CH4GlobalTrans=calibration[CalModel]['CH4GlobalTrans']
        #NH3GlobalTrans=calibration[CalModel]['NH3GlobalTrans']
        
        kernel = Gaussian2DKernel(1)
        if Smoothing:
            print('aaaaaaaaaaaaaaaaaaa=',a)
            """print('CH4GlobalTrans,NH3GlobalTrans=',CH4GlobalTrans,NH3GlobalTrans)
            if a=='L3Cal':
                CH4_tau=-np.log(convolve((CH4data*CH4GlobalTrans),kernel,boundary='extend'))
                hdu = fits.PrimaryHDU((CH4data*CH4GlobalTrans).astype(np.float32))
                hdul = fits.HDUList([hdu])
                try:
                    os.remove(pathout+'CH4L3Cal.fits')
                except: 
                    print("file doesn't exist")
                hdul.writeto(pathout+'CH4L3Cal.fits')
                hdul.close()
                
                NH3_tau=-np.log(convolve((NH3data*NH3GlobalTrans),kernel,boundary='extend'))
                hdu = fits.PrimaryHDU((NH3data*NH3GlobalTrans).astype(np.float32))
                hdul = fits.HDUList([hdu])
                try:
                    os.remove(pathout+'NH3L3Cal.fits')
                except: 
                    print("file doesn't exist")
                hdul.writeto(pathout+'NH3L3Cal.fits')
                hdul.close()
 
            if a=='L2Cal':
                CH4_tau=-np.log(convolve(CH4data,kernel,boundary='extend'))
                hdu = fits.PrimaryHDU((CH4data).astype(np.float32))
                hdul = fits.HDUList([hdu])
                try:
                    os.remove(pathout+'CH4L.fits')
                except: 
                    print("file doesn't exist")
                hdul.writeto(pathout+'CH4.fits')
                hdul.close()

                NH3_tau=-np.log(convolve(NH3data,kernel,boundary='extend'))
                hdu = fits.PrimaryHDU((NH3data).astype(np.float32))
                hdul = fits.HDUList([hdu])
                try:
                    os.remove(pathout+'NH3.fits')
                except: 
                    print("file doesn't exist")
                hdul.writeto(pathout+'NH3.fits')
                hdul.close()
            smthtitle="Smoothed"
        else: 
            """
        CH4_tau=-np.log(CH4data)
        NH3_tau=-np.log(NH3data)
        CH4_Ncol=1000*CH4_tau/K_eff_CH4620
        NH3_Ncol=1000*NH3_tau/K_eff_NH3647
    
        CH4_Cloud_Press=(CH4_Ncol)*amagat*gravity*mean_mol_wt/(fCH4*STP)
        NH3=(NH3_Ncol/1000.)*amagat*gravity*mean_mol_wt/(CH4_Cloud_Press*STP)
        ##!!!! WOW!!! I need to calculate fNH3 the EASY way also and compare!!!
        fNH3=(1.81e-3)*NH3_Ncol/CH4_Ncol #-> It works 7/20/2023
        
        for BUNIT in ['Mole Fraction','Cloud-top Press']:
            if BUNIT=='Mole Fraction':
                hdu = fits.PrimaryHDU(fNH3.astype(np.float32))
                hdul = fits.HDUList([hdu])
                comment="NH3 ppm"
                fn='fNH3'
                Real_CM1=NH3_CM1
                Real_CM2=NH3_CM2
                Real_CM3=NH3_CM3
                time=NH3time
                file=NH3file
            elif BUNIT=='Cloud-top Press':
                hdu = fits.PrimaryHDU((CH4_Cloud_Press/4.).astype(np.float32))
                hdul = fits.HDUList([hdu])
                comment="mb"
                fn='PCloud'
                Real_CM1=CH4_CM1
                Real_CM2=CH4_CM2
                Real_CM3=CH4_CM3
                time=CH4time
                file=CH4file
                
            hdul[0].header['BITPIX']=-32
            hdul[0].header['DATE-OBS']=time.replace('_','T')+'Z'
            #hdul[0].header['DATE']=NH3time
            hdul[0].header['AUTHOR']='Hill, S. M.'
        
            hdul[0].header['OBJECT']='Jupiter'
            hdul[0].header['TELESCOP']=sourcefiles[sourcedata]['Metadata']['Telescope']
            hdul[0].header['INSTRUME']=sourcefiles[sourcedata]['Metadata']['Camera']
            hdul[0].header['SEEING']=sourcefiles[sourcedata]['Metadata']['Seeing']
            hdul[0].header['TRANSPAR']=sourcefiles[sourcedata]['Metadata']['Transparency']
            hdul[0].header['BUNIT']=(BUNIT,comment)
            hdul[0].header['CALIBRA']=(CalModel,'Disk-Integrated Cal Ref')
            hdul[0].header['VERSION']=('TBD','TBD')
            hdul[0].header['CTYPE1']=('Sys. 2 Longitude','deg')
            hdul[0].header['CTYPE2']=('PG Latitude','deg')
            
            hdul[0].header['CM1']=(Real_CM1,'Sys. 1 Long. Central Meridian')
            hdul[0].header['CM2']=(Real_CM2,'Sys. 1 Long. Central Meridian')
            hdul[0].header['CM3']=(Real_CM3,'Sys. 1 Long. Central Meridian')
            hdul[0].header['SMOOTH']=(Smoothing,'')
            hdul[0].header['KERNEL']=(1,'Gaussian')
            hdul[0].header['CH4ABSFL']=(CH4file,'Source file for CH4 Absorption')
            hdul[0].header['NH3ABSFL']=(NH3file,'Source file for NH3 Absorption')
            hdul[0].header['CONTXTFL']=(RGBfile,'Source file for RGB Context')
            
            #print(hdul[0].header[:])
            #fnout=path+'/'+NH3file[0:26]+'fNH3Abs647.fits'
            fnout=pathout+file[0:26]+fn+'_Sys'+LonSys+'.fits'
            if a=='L3Cal':
                fnout=pathout+file[0:26]+fn+'_Sys'+LonSys+'.fits'
            if a=='L2Cal':
                fnout=pathout+file[0:26]+fn+'_Sys'+LonSys+'.fits'

            try:
                os.remove(fnout)
            except: 
                print("file doesn't exist")
            hdul.writeto(fnout)
            hdul.close()
            
            outdict={"CH4":{"file":CH4file,"CM":CH4CM,"time":CH4time,"tau":CH4_tau,"Ncol":CH4_Ncol,"PCloud":CH4_Cloud_Press},
                     "NH3":{"file":NH3file,"CM":NH3CM,"time":NH3time,"tau":NH3_tau,"Ncol":NH3_Ncol,"fNH3":fNH3},
                     "RGB":{"file":RGBfile,"CM":RGBCM,"time":RGBtime,"RGB":RGB},
                     "pathFITS":pathFITS,"pathout":pathout}#,"Variation":variation}
        
    return(outdict)
