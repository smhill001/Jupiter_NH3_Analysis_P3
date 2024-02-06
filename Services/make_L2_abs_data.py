def make_l2_abs_data(obsdate="20231207UTc",target="Jupiter",
                     imagetype='Map'):
    """
    Created on Fri Sep 15 07:48:05 2023
    
    PURPOSE: Converts L1 normalized radiance RGB PNG maps into L2 transmission
             data for CH4 and NH3. Options are planned for
             maps vs images.
             The code can be called in batch mode by:
                 AmmoniaMapsScript_P3.py
                 
    UPDATES:
        2023-12-25: Made output FITS in sys 3 long instead of sys 2
                    and removed all CH4 channels other than 619nm to 
                    simplify the code. When it's desireable to address other
                    bands, e.g., 727 nm, 889 nm, that should be done
                    in new tailored codes rather than here.
    
    @author: smhil
    """
    import sys
    drive='c:'
    sys.path.append(drive+'/Astronomy/Python Play')
    sys.path.append(drive+'/Astronomy/Python Play/Util_P3')
    sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy')
    sys.path.append('./Services')

    import os
    import numpy as np
    from imageio import imwrite
    from numpy import inf
    from astropy.io import fits
    import get_obs_list as getlist
    import load_png as LP
    import get_WINJupos_ephem as WJ_ephem
    import read_master_calibration as RMC
    import make_sza_eza_planes as za
    import pylab as pl

    ###########################################################################
    #  DATA FILES AND METADATA DICTIONARY
    ###########################################################################
    sourcedata=obsdate#+"_"+imagetype
    sourcefiles=getlist.get_obs_list()
    
    if sourcefiles[sourcedata]["Telescope"]=="C11":
        CalModel="SCT-Obs-Final"
    elif sourcefiles[sourcedata]["Telescope"]=="VLT":
        CalModel="VLT-Obs-Final"

    calibration,K_eff=RMC.read_master_calibration()
    CH4GlobalTrans=calibration[CalModel]['CH4GlobalTrans']
    NH3GlobalTrans=calibration[CalModel]['NH3GlobalTrans']
    ###########################################################################
    # OBTAIN IMAGES TO DISPLAY, READ DATA, AND DETERMINE IMAGE ARRAY SIZE
    ###########################################################################             
    path='c:/Astronomy/Projects/Planets/'+target+'/Imaging Data/'+obsdate[0:10]+'/'
    if imagetype=="Img":
        NH3file=sourcefiles[sourcedata]['NH3file']+".png"
        CH4file=sourcefiles[sourcedata]['CH4file']+".png"
    elif imagetype=="Map":
        NH3file=sourcefiles[sourcedata]['NH3file']+"_CM2_L360_MAP-BARE.png"
        CH4file=sourcefiles[sourcedata]['CH4file']+"_CM2_L360_MAP-BARE.png"
    NH3_RGB=LP.load_png(path+NH3file)
    CH4_RGB=LP.load_png(path+CH4file)

    ###########################################################################
    # EPHEMERIS, COLOR SLOPE AND NH3 TRANSMISSION CALCULATIONS
    ###########################################################################
    NH3sec=str(int(str(NH3file[16:17]))*6) #COMPUTE FROM FRACTIONAL WINJUPOS MINUTE
    NH3time=(NH3file[0:10]+"_"+NH3file[11:13]+":"+NH3file[13:15]+":"
             +NH3sec.zfill(2))
    NH3eph=WJ_ephem.get_WINJupos_ephem(NH3time)
    NH3_CM1=float(NH3eph[0].strip())
    NH3_CM2=float(NH3eph[1].strip())
    NH3_CM3=float(NH3eph[2].strip())
    clrslp=(np.array(NH3_RGB[:,:,0]).astype(float)
            -np.array(NH3_RGB[:,:,2]).astype(float))/24.0 
    CNT647=15.0*clrslp+np.array(NH3_RGB[:,:,2])
    nh3abs=NH3GlobalTrans*np.array(NH3_RGB[:,:,1])/CNT647

    ###########################################################################
    # EPHEMERIS AND CH4 TRANSMISSION CALCULATIONS
    ###########################################################################
    CH4sec=str(int(str(CH4file[16:17]))*6) #COMPUTE FROM FRACTIONAL WINJUPOS MINUTE
    CH4time=(CH4file[0:10]+"_"+CH4file[11:13]+":"+CH4file[13:15]+":"
             +CH4sec.zfill(2))
    CH4eph=WJ_ephem.get_WINJupos_ephem(CH4time)
    CH4_CM1=float(CH4eph[0].strip())
    CH4_CM2=float(CH4eph[1].strip())
    CH4_CM3=float(CH4eph[2].strip()) 
    CNTSynth=np.array(CH4_RGB[:,:,2]) #632nm continuum
    ch4abs=CH4GlobalTrans*np.array(CH4_RGB[:,:,1])/CNTSynth
    ###########################################################################
    # Assuming the input (L2) FITS map files are in system 2 longitude 
    # coordinates, 'roll' the data array in longitude to System 3 Longitude
    # as the defacto standard. Once the L2 data are converted to system 3,
    # this step will no longer be necessasry.
    ###########################################################################
    #RGBroll=RGB_CM2-RGB_CM3
    #RGB=np.roll(RGB,int(RGBroll),axis=1)
    
    if imagetype=="Map":
        NH3roll=NH3_CM2-NH3_CM3
        NH3data=np.roll(nh3abs,int(NH3roll),axis=1)
        
        CLSLdata=np.roll(clrslp,int(NH3roll),axis=1)
        
        CH4roll=CH4_CM2-CH4_CM3
        CH4data=np.roll(ch4abs,int(CH4roll),axis=1)
        
        NH3CM=NH3_CM3
        CH4CM=CH4_CM3
    else:
        NH3data=nh3abs
        CLSLdata=clrslp
        CH4data=ch4abs

    #RGBCM=RGB_CM3
    eza,sza=za.make_sza_eza_planes(dateobs=NH3time.replace('_','T')+'Z')
 
   ###########################################################################
    # SET UP FILE PATH AND NAMES
    ###########################################################################
    pathout='C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Analysis Data/L2 FITS/'
    #fn=NH3file[0:26]+'L2CLSL'+imagetype
    #print("fn=",fn)
    """try:
        fnout=fnout+sourcefiles[sourcedata]['Variation']
    except:
        print('No Variation')"""
    
    ###########################################################################
    # Create FITS files and headers, including backplanes
    ###########################################################################
    for BUNIT in ['NH3 Trans','CH4 Trans','Clr Slp']:
        if BUNIT=='NH3 Trans':
            dataarray=NH3data
            hdu = fits.PrimaryHDU(NH3data.astype(np.float32))
            comment="NH3 Transmission"
            Real_CM1=NH3_CM1
            Real_CM2=NH3_CM2
            Real_CM3=NH3_CM3
            datetime=NH3time
            file=NH3file
            fn=file[0:25]+"_"+imagetype+'_L2TNH3'
            Range=[0.90,1.0] #scaled range for PNG file
        elif BUNIT=='CH4 Trans':
            dataarray=CH4data
            hdu = fits.PrimaryHDU(CH4data.astype(np.float32))
            comment="CH4 Transmission"
            Real_CM1=CH4_CM1
            Real_CM2=CH4_CM2
            Real_CM3=CH4_CM3
            datetime=CH4time
            file=CH4file
            fn=file[0:25]+"_"+imagetype+'_L2TCH4'
            Range=[0.80,1.0] #scaled range for PNG file
        elif BUNIT=='Clr Slp':
            dataarray=CLSLdata
            hdu = fits.PrimaryHDU(CLSLdata.astype(np.float32))
            comment="Color Slope (632 to 656 nm)"
            Real_CM1=NH3_CM1
            Real_CM2=NH3_CM2
            Real_CM3=NH3_CM3
            datetime=NH3time
            file=NH3file
            fn=file[0:25]+"_"+imagetype+'_L2CLSL'
            Range=[0.0,100.0] #scaled range for PNG file

            
        szadata=fits.ImageHDU(sza)
        ezadata=fits.ImageHDU(eza)
        hdul = fits.HDUList([hdu,szadata,ezadata])

        ###########################################################################
        # BUILD FITS FILE AND HEADER AND WRITE IT FOR COLOR SLOPE
        ###########################################################################       
        hdul[0].header['BITPIX']=-32
        hdul[0].header['DATE-OBS']=NH3time.replace('_','T')+'Z'
        hdul[0].header['AUTHOR']='Hill, S. M.'
        hdul[0].header['FILENAME']=fn
        
        hdul[0].header['OBJECT']='Jupiter'
        hdul[0].header['TELESCOP']=sourcefiles[sourcedata]['Telescope']
        hdul[0].header['INSTRUME']=sourcefiles[sourcedata]['Camera']
        hdul[0].header['SEEING']=sourcefiles[sourcedata]['Seeing']
        hdul[0].header['TRANSPAR']=sourcefiles[sourcedata]['Transparency']
        hdul[0].header['BUNIT']=(BUNIT,comment)
        hdul[0].header['CALIBRA']=(CalModel,'Disk-Integrated Cal Ref')
        hdul[0].header['VERSION']=('TBD','TBD')
        hdul[0].header['CTYPE1']=('Sys. 3 Longitude','deg')
        hdul[0].header['CTYPE2']=('PG Latitude','deg')
        
        hdul[0].header['CM1']=(Real_CM1,'Sys. 1 Long. Central Meridian')
        hdul[0].header['CM2']=(Real_CM2,'Sys. 1 Long. Central Meridian')
        hdul[0].header['CM3']=(Real_CM3,'Sys. 1 Long. Central Meridian')
        hdul[0].header['SOURCEFL']=(file,'Source file')
        
        fnout=pathout+fn

        try:
            os.remove(fnout+'.fits')
        except: 
            print("file doesn't exist")
        hdul.writeto(fnout+'.fits')
        hdul.close()
        
        #######################################################################
        # WRITE Scaled PNG File
        ######################################################################
        scl_arr = ((dataarray - Range[0]) * (1/(Range[1] - Range[0]) * 65534.9999)).astype('uint16')#
        imwrite(fnout+'.png', scl_arr)#.astype(np.uint16))
