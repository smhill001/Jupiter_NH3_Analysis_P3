# -*- coding: utf-8 -*-

def make_l2_abs_data(obsdate="20220919UTa",target="Jupiter",
                     imagetype='Map'):
    """
    Created on Fri Sep 15 07:48:05 2023
    
    @author: smhil
    """
    import sys
    drive='c:'
    sys.path.append(drive+'/Astronomy/Python Play')
    sys.path.append(drive+'/Astronomy/Python Play/Util_P3')
    sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy')

    import os
    import numpy as np
    from imageio import imwrite
    from numpy import inf
    from astropy.io import fits
    sys.path.append('./Services')
    import get_L1_img_data as getL1
    import load_png as LP
    import get_WINJupos_ephem as WJ_ephem
    import read_master_calibration as RMC
    
    for a in ['L2Cal']:
        ###########################################################################
        #  DATA FILES AND METADATA DICTIONARY
        ###########################################################################
        sourcedata=obsdate+"_"+imagetype
        sourcefiles=getL1.get_L1_img_data()
        calibration,K_eff=RMC.read_master_calibration()
        
        if sourcefiles[sourcedata]["Metadata"]["Telescope"]=="C11":
            CalModel="SCT-Obs-Final"
        elif sourcefiles[sourcedata]["Metadata"]["Telescope"]=="VLT":
            CalModel="VLT-Obs-Final"

        CH4GlobalTrans=calibration[CalModel]['CH4GlobalTrans']
        NH3GlobalTrans=calibration[CalModel]['NH3GlobalTrans']
        ###########################################################################
        # OBTAIN IMAGES TO DISPLAY, READ DATA, AND DETERMINE IMAGE ARRAY SIZE
        ###########################################################################             
        path='c:/Astronomy/Projects/Planets/'+target+'/Imaging Data/'+obsdate[0:10]+'/'
        NH3file=sourcefiles[sourcedata]['NH3file']
        CH4file=sourcefiles[sourcedata]['CH4file']
        
        NH3_RGB=LP.load_png(path+NH3file)
        ###########################################################################
        #!!CREATE MASK. WHY?? - Removed Mask on 11/17/2023. May be desired
        #   for images (removing limb artifacts) as opposed to maps
        ###########################################################################  
        indices=(NH3_RGB>0.002*NH3_RGB.max())
        mask=np.ones(NH3_RGB.shape)
        #mask[indices]=1.
        ###########################################################################
        # PARSE WINJUPOS TIME AND GET EPHEMERIS
        ###########################################################################
        sec=str(int(str(NH3file[16:17]))*6) #COMPUTE FROM FRACTIONAL WINJUPOS MINUTE
        NH3time=(NH3file[0:10]+"_"+NH3file[11:13]+":"+NH3file[13:15]+":"
                 +sec.zfill(2))
        eph=WJ_ephem.get_WINJupos_ephem(NH3time)
        ###########################################################################
        # COMPUTE COLOR SLOPE BY PIXEL ASSUMING CHANNEL 0 IS R(656) AND CHANNEL 2
        #   IS B(632) AND THERE IS NO BULK COLOR SLOPE. THEN COMPUTE THE EFFECTIVE
        #   CONTINUUM AT 647NM (NH3) AND THE CONTINUUM DIVIDED NH3 IMAGE.
        ###########################################################################
        clrslp=(np.array(NH3_RGB[:,:,0]).astype(float)
                -np.array(NH3_RGB[:,:,2]).astype(float))/24.0 
        CNT647=15.0*clrslp+np.array(NH3_RGB[:,:,2])
        nh3abs=NH3GlobalTrans*np.array(NH3_RGB[:,:,1])/CNT647
        ###########################################################################
        # SET UP FILE PATH AND NAMES
        ###########################################################################
        pathout='C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Analysis Data/L2 FITS/'
        fn=NH3file[0:26]+'647ClrSlp'+imagetype
        print("fn=",fn)
        fnout=pathout+fn
        try:
            fnout=fnout+sourcefiles[sourcedata]['Metadata']['Variation']
        except:
            print('No Variation')

        ###########################################################################
        # GET METADATA
        ###########################################################################
        NH3sec=str(int(str(NH3file[16:17]))*6)
        NH3time=(NH3file[0:10]+"_"+NH3file[11:13]+":"+NH3file[13:15]+":"+NH3sec.zfill(2))
        eph=WJ_ephem.get_WINJupos_ephem(NH3time)
        NH3_CM1=float(eph[0].strip())
        NH3_CM2=float(eph[1].strip())
        NH3_CM3=float(eph[2].strip())
        
        ###########################################################################
        # BUILD FITS FILE AND HEADER AND WRITE IT FOR COLOR SLOPE
        ###########################################################################
        hdu = fits.PrimaryHDU(clrslp.astype(np.float32))
        hdul = fits.HDUList([hdu])
        hdul[0].header['BITPIX']=-32
        hdul[0].header['AUTHOR']='Hill, S. M.'
        hdul[0].header['FILENAME']=fn
        hdul[0].header['OBJECT']='Jupiter'
        hdul[0].header['TELESCOP']=sourcefiles[sourcedata]['Metadata']['Telescope']
        hdul[0].header['INSTRUME']=sourcefiles[sourcedata]['Metadata']['Camera']
        hdul[0].header['SEEING']=sourcefiles[sourcedata]['Metadata']['Seeing']
        hdul[0].header['TRANSPAR']=sourcefiles[sourcedata]['Metadata']['Transparency']
        hdul[0].header['CALIBRA']=(CalModel,'Disk-Integrated Cal Ref')
        hdul[0].header['VERSION']=('TBD','TBD')
        hdul[0].header['CTYPE1']=('Sys. 2 Longitude','deg')
        hdul[0].header['CTYPE2']=('PG Latitude','deg')
        
        hdul[0].header['DATE-OBS']=NH3time.replace('_','T')+'Z'
        hdul[0].header['CM1']=(NH3_CM1,'Sys. 1 Long. Central Meridian')
        hdul[0].header['CM2']=(NH3_CM2,'Sys. 1 Long. Central Meridian')
        hdul[0].header['CM3']=(NH3_CM3,'Sys. 1 Long. Central Meridian')
        hdul[0].header['BUNIT']=("ClrSlp","Non-dimensional normalized color slope")
        hdul[0].header['SOURCEFL']=(NH3file,'Source file for NH3 color slope')

        try:
            os.remove(fnout+'.fits')
        except: 
            print("file doesn't exist")
        hdul.writeto(fnout+'.fits')
        hdul.close()
    
        clrslp16bit = np.nan_to_num((32000+500*(clrslp)).astype(np.uint16))
        imwrite(fnout+'.png', clrslp16bit.astype(np.uint16))
        #  !!!!OR, FOR SEPARATION OF FUNCTION, CREATE A STANDALONG FILE-WRITING FUNCTION   
        fn=NH3file[0:25]+'_647NH3Abs'+imagetype
        fnout=pathout+fn
        try:
            fnout=fnout+sourcefiles[sourcedata]['Metadata']['Variation']
        except:
            print('No Variation')
        
        ###########################################################################
        # BUILD FITS FILE AND HEADER AND WRITE IT FOR NH3 ABSORPTION
        ###########################################################################
        maskednh3=nh3abs*mask[:,:,1]
        #normnh3=maskednh3/maskednh3[maskednh3>0].mean()
        normnh3=maskednh3#/maskednh3[maskednh3>0].mean()
        hdu = fits.PrimaryHDU(normnh3.astype(np.float32))
        hdul = fits.HDUList([hdu])
        hdul[0].header['BITPIX']=-32
        hdul[0].header['AUTHOR']='Hill, S. M.'
        hdul[0].header['FILENAME']=fn
        hdul[0].header['OBJECT']='Jupiter'
        hdul[0].header['TELESCOP']=sourcefiles[sourcedata]['Metadata']['Telescope']
        hdul[0].header['INSTRUME']=sourcefiles[sourcedata]['Metadata']['Camera']
        hdul[0].header['SEEING']=sourcefiles[sourcedata]['Metadata']['Seeing']
        hdul[0].header['TRANSPAR']=sourcefiles[sourcedata]['Metadata']['Transparency']
        hdul[0].header['CALIBRA']=(CalModel,'Disk-Integrated Cal Ref')
        hdul[0].header['VERSION']=('TBD','TBD')
        hdul[0].header['CTYPE1']=('Sys. 2 Longitude','deg')
        hdul[0].header['CTYPE2']=('PG Latitude','deg')
        
        hdul[0].header['DATE-OBS']=NH3time.replace('_','T')+'Z'
        hdul[0].header['CM1']=(NH3_CM1,'Sys. 1 Long. Central Meridian')
        hdul[0].header['CM2']=(NH3_CM2,'Sys. 1 Long. Central Meridian')
        hdul[0].header['CM3']=(NH3_CM3,'Sys. 1 Long. Central Meridian')
        hdul[0].header['BUNIT']=("Trans","NH3 dimensionless transmission")
        hdul[0].header['SOURCEFL']=(NH3file,'Source file for NH3 Absorption')

        try:
            os.remove(fnout+'.fits')
        except: 
            print("file doesn't exist")
        hdul.writeto(fnout+'.fits')
        hdul.close()
        #nh3abs16bit = np.nan_to_num(((5.*65535.*(normnh3 - 0.9))*mask[:,:,1]).astype(np.uint16))
        normnh3scaled=np.nan_to_num(((5.*65535.*(normnh3*mask[:,:,1] - 0.9))))
        normnh3scaled[normnh3scaled<=0.]=0.0
        nh3abs16bit = normnh3scaled.astype(np.uint16)
        imwrite(fnout+'.png', nh3abs16bit)#.astype(np.uint16))
        
    
        ###########################################################################
        # COMPUTE METHANE IMAGE PANELS ALONG WITH COLOR SLOPE (VALID AT NH3 EPOCH)
        #   UP TO THREE SETS OF METHANE PANELS ARE COMPUTED (620NM, 730NM, 889NM)
        #   COULD ADD >1000NM, BUT THAT WOULD BE FOR A LATER DATE
        ###########################################################################
        #
        if CH4file != ['NA']: #NEED LOOPED CALLS HERE BASED ON NUMBER OF FILES AND LOOK UP FOR LABELS
            for ifile in range(0,len(CH4file)):
                CH4fl=CH4file[ifile]
                CH4_RGB=LP.load_png(path+CH4fl)
                print(CH4file[ifile])
                if "G620" in CH4fl:
                    wvstr='620'
                    CH4labels=['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))']
                elif "G730" in CH4fl:
                    wvstr='730'
                    CH4labels=['Synth. Continuum @ 730nm','730nm (CH4)','730/Cont. (CH4))']
                elif "G889" in CH4fl:
                    wvstr='889'
                    CH4labels=['Continuum @ 940nm','889nm (CH4)','889/940 (CH4))']
        
        ###########################################################################
        # PARSE WINJUPOS TIME AND GET EPHEMERIS
        ###########################################################################
                sec=str(int(str(CH4fl[16:17]))*6)
                CH4time=(CH4fl[0:10]+"_"+CH4fl[11:13]+":"+CH4fl[13:15]+":"+sec.zfill(2))
                eph=WJ_ephem.get_WINJupos_ephem(CH4time)
                MethaneHeader=CH4time[11:19]+" UT; CM1"+eph[0]+"; CM2"+eph[1]+"; CM3"+eph[2]+"; Alt"+eph[3]
                ###########################################################################
                # COMPUTE COLOR SLOPE BY PIXEL 
                ###########################################################################
                meta={'620nm (CH4)':{'dwave':-12,'maxmin':[0.9,1.1]},
                      '730nm (CH4)':{'dwave':98,'maxmin':[0.8,1.2]},
                      '889nm (CH4)':{'dwave':0,'maxmin':[0.5,2.0]}}
    
                if CH4labels[1]=='889nm (CH4)':
                    CNTSynth=np.array(CH4_RGB[:,:,0]) #940nm continuum
                    ch4abs=(np.array(CH4_RGB[:,:,1])+0.0001)/(CNTSynth+0.0001)
                elif CH4labels[1]=='730nm (CH4)':# or CH4labels[1]=='620nm (CH4)':
                    clrslp=(np.array(CH4_RGB[:,:,0]).astype(float)-np.array(CH4_RGB[:,:,2]).astype(float))/24.0 
                    CNTSynth=meta[CH4labels[1]]['dwave']*clrslp+np.array(CH4_RGB[:,:,2])
                elif CH4labels[1]=='620nm (CH4)': ###ADDITION FOR PIECEWISE CONTINUOUS RATHER THAN EXTRAPOLATED
                    CNTSynth=np.array(CH4_RGB[:,:,2]) #632nm continuum
    
                maxmin=meta[CH4labels[1]]['maxmin']
    
                ch4abs=CH4GlobalTrans*(np.array(CH4_RGB[:,:,1])+0.0001)/(CNTSynth+0.0001)
        
                ch4abs[ch4abs == inf] = 0
                ch4abs[ch4abs == -inf] = 0
                ch4abs[ch4abs == np.nan] = 0
                
                #nh3abs16bit = np.nan_to_num(((0.6667*65535*(ch4abs-0.5))*mask[:,:,0]).astype(np.uint16))
                #ch4abs16bit = np.nan_to_num(((5.*65535.*(ch4abs)-0.1)*mask[:,:,1]).astype(np.uint16))
                #imwrite(path+'/'+CH4file[0][0:26]+'CH4'+wvstr+'AbsPython.png', ch4abs16bit.astype(np.uint16))
                ###!!!!FIX CH4 FILE TIME IN NAMING!!!!
    
                if "G620" in CH4fl:
                    fn=CH4fl[0:25]+'_620CH4Abs'+imagetype
                    fnout=pathout+fn
                    """
                    if a=='L2Cal':
                        fnout=pathout+NH3file[0:26]+'_620CH4Abs'+a+'_'+imagetype
                    elif a=='L3Cal':
                        fnout=pathout+NH3file[0:26]+'_620CH4Abs'+a+'_'+imagetype
                        """
                    try:
                        fnout=fnout+sourcefiles[sourcedata]['Metadata']['Variation']
                    except:
                        print('No Variation')
                    ###########################################################
                    maskedch4=ch4abs*mask[:,:,1]
                    #normch4=maskedch4/maskednh3[maskedch4>0].mean()
                    ###########################################################
                    normch4=maskedch4#/maskednh3[maskedch4>0].mean()
                    
                    #hdu = fits.PrimaryHDU(normch4.astype(np.float32))
                    hdu = fits.PrimaryHDU(normch4.astype(np.float32))
                    hdul = fits.HDUList([hdu])
                    hdul[0].header['BITPIX']=-32
                    hdul[0].header['AUTHOR']='Hill, S. M.'
                    hdul[0].header['FILENAME']=fn
                    hdul[0].header['OBJECT']='Jupiter'
                    hdul[0].header['TELESCOP']=sourcefiles[sourcedata]['Metadata']['Telescope']
                    hdul[0].header['INSTRUME']=sourcefiles[sourcedata]['Metadata']['Camera']
                    hdul[0].header['SEEING']=sourcefiles[sourcedata]['Metadata']['Seeing']
                    hdul[0].header['TRANSPAR']=sourcefiles[sourcedata]['Metadata']['Transparency']
                    hdul[0].header['CALIBRA']=("TEST",'Disk-Integrated Cal Ref')
                    hdul[0].header['VERSION']=('TBD','TBD')
                    hdul[0].header['CTYPE1']=('Sys. 2 Longitude','deg')
                    hdul[0].header['CTYPE2']=('PG Latitude','deg')
                    
                    hdul[0].header['DATE-OBS']=CH4time.replace('_','T')+'Z'
                    hdul[0].header['CM1']=(NH3_CM1,'Sys. 1 Long. Central Meridian')
                    hdul[0].header['CM2']=(NH3_CM2,'Sys. 1 Long. Central Meridian')
                    hdul[0].header['CM3']=(NH3_CM3,'Sys. 1 Long. Central Meridian')
                    hdul[0].header['BUNIT']=("Trans","CH4 dimensionless transmission")
                    hdul[0].header['SOURCEFL']=(CH4file[0],'Source file for NH3 Absorption')

                    try:
                        os.remove(fnout+'.fits')
                    except: 
                        print("file doesn't exist")
                        
                    print("############# fnout=",fnout+".fits")
    
                    hdul.writeto(fnout+'.fits')
                    hdul.close()
    
                    normch4scaled=np.nan_to_num(((5.*65535.*(normch4*mask[:,:,1] - 0.9))))
                    normch4scaled[normch4scaled<=0.]=0.0
                    ch4abs16bit = normch4scaled.astype(np.uint16)
                    imwrite(fnout+'.png', ch4abs16bit)#.astype(np.uint16))