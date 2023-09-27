# -*- coding: utf-8 -*-

def Images_Maps_Arrays_P3(obsdate="20220919UTa",target="Jupiter",
                                     imagetype='Map'):
    """
    Created on Wed Oct 06 14:36:01 2021

    @author: Steven Hill

    PURPOSE:    Read a white balanced RGB (656,647,632) image and compute the
                color slope and NH3 absorption images.

    EXAMPLE:    Jupiter_NH3_CH4_Images_Arrays_P3(obsdate="20210905UT")
                
    """
    import sys
    drive='c:'
    sys.path.append(drive+'/Astronomy/Python Play')
    sys.path.append(drive+'/Astronomy/Python Play/Util_P3')
    sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy')

    import os
    from matplotlib.pyplot import imread
    import pylab as pl
    import numpy as np
    from imageio import imwrite
    from numpy import inf
    from astropy.io import fits
    import RetrievalLibrary as RL
    sys.path.append('./Services')
    import get_L1_img_data as getL1
    import load_png as LP
    import get_WINJupos_ephem as WJ_ephem
 
    ###########################################################################
    #  DATA FILES AND METADATA DICTIONARY
    #    !!!!SHOULD MAKE THIS A DATA OBJECT!!!!
    ###########################################################################
    sourcedata=obsdate+"_"+imagetype
    sourcefiles=getL1.get_L1_img_data()
    ###########################################################################
    # OBTAIN IMAGES TO DISPLAY, READ DATA, AND DETERMINE IMAGE ARRAY SIZE
    ###########################################################################             
    path='c:/Astronomy/Projects/Planets/'+target+'/Imaging Data/'+obsdate[0:10]+'/'
    NH3file=sourcefiles[sourcedata]['NH3file']
    CH4file=sourcefiles[sourcedata]['CH4file']
    NUVfile=sourcefiles[sourcedata]['Context']['NUVfile']
    RGBfile=sourcefiles[sourcedata]['Context']['RGBfile']

    print(NUVfile,RGBfile)
    nny=1
    #print(nny)
    if NH3file != 'NA':
        NH3_RGB=LP.load_png(path+NH3file)
        nny=nny+1
        
    if CH4file != ['NA']:
        nny=nny+len(CH4file)
        print("in CH4")
        
    if NUVfile != 'NA':
        NUV_RGB=LP.load_png(path+NUVfile)
    else:
        NUV_RGB=[]
        
    if RGBfile != 'NA':
        RGB_RGB=imread(path+RGBfile)
    else:
        RGB_RGB=[]
    #print(nny)
    
    if NUVfile !='NA' or RGBfile != 'NA':
        nny=nny+1

    ###########################################################################
    #!!CREATE MASK. WHY??
    ###########################################################################  
    indices=(NH3_RGB>0.2*NH3_RGB.max())
    mask=np.zeros(NH3_RGB.shape)
    mask[indices]=1.

    ###########################################################################
    # PARSE WINJUPOS TIME AND GET EPHEMERIS
    ###########################################################################
    sec=str(int(str(NH3file[16:17]))*6) #COMPUTE FROM FRACTIONAL WINJUPOS MINUTE
    NH3time=(NH3file[0:10]+"_"+NH3file[11:13]+":"+NH3file[13:15]+":"
             +sec.zfill(2))
    eph=WJ_ephem.get_WINJupos_ephem(NH3time)
    AmmoniaHeader=NH3time[11:19]+" UT; CM1"+eph[0]+"; CM2"+eph[1]+"; CM3"\
                    +eph[2]+"; Alt"+eph[3]

    ###########################################################################
    # COMPUTE COLOR SLOPE BY PIXEL ASSUMING CHANNEL 0 IS R(656) AND CHANNEL 2
    #   IS B(632) AND THERE IS NO BULK COLOR SLOPE. THEN COMPUTE THE EFFECTIVE
    #   CONTINUUM AT 647NM (NH3) AND THE CONTINUUM DIVIDED NH3 IMAGE.
    ###########################################################################
    clrslp=(np.array(NH3_RGB[:,:,0]).astype(float)
            -np.array(NH3_RGB[:,:,2]).astype(float))/24.0 
    CNT647=15.0*clrslp+np.array(NH3_RGB[:,:,2])
    nh3abs=np.array(NH3_RGB[:,:,1])/CNT647

    ###########################################################################
    # SET UP FILE PATH AND NAMES
    ###########################################################################
    pathout='C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Analysis Data/L2 FITS/'
    fnout=pathout+NH3file[0:26]+'ClrSlp'+imagetype
        
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
    hdul[0].header['DATE-OBS']=NH3time.replace('_','T')+'Z'
    #hdul[0].header['DATE']=NH3time
    hdul[0].header['AUTHOR']='Hill, S. M.'

    hdul[0].header['OBJECT']='Jupiter'
    hdul[0].header['TELESCOP']=sourcefiles[sourcedata]['Metadata']['Telescope']
    hdul[0].header['INSTRUME']=sourcefiles[sourcedata]['Metadata']['Camera']
    hdul[0].header['SEEING']=sourcefiles[sourcedata]['Metadata']['Seeing']
    hdul[0].header['TRANSPAR']=sourcefiles[sourcedata]['Metadata']['Transparency']
    hdul[0].header['BUNIT']=("TEST","comment")
    hdul[0].header['CALIBRA']=("TEST",'Disk-Integrated Cal Ref')
    hdul[0].header['VERSION']=('TBD','TBD')
    hdul[0].header['CTYPE1']=('Sys. 2 Longitude','deg')
    hdul[0].header['CTYPE2']=('PG Latitude','deg')
    
    hdul[0].header['CM1']=(NH3_CM1,'Sys. 1 Long. Central Meridian')
    hdul[0].header['CM2']=(NH3_CM2,'Sys. 1 Long. Central Meridian')
    hdul[0].header['CM3']=(NH3_CM3,'Sys. 1 Long. Central Meridian')
    hdul[0].header['SMOOTH']=('Smoothing','')
    hdul[0].header['KERNEL']=(1,'Gaussian')
    #hdul[0].header['CH4ABSFL']=(CH4file,'Source file for CH4 Absorption')
    hdul[0].header['NH3ABSFL']=(NH3file,'Source file for NH3 Absorption')
    hdul[0].header['CONTXTFL']=(RGBfile,'Source file for RGB Context')
    
    try:
        os.remove(fnout)
    except: 
        print("file doesn't exist")
    hdul.writeto(fnout)
    hdul.close()

    print(hdul[0].header[:])

    clrslp16bit = np.nan_to_num((32000+500*(clrslp)).astype(np.uint16))
    imwrite(fnout+'.png', clrslp16bit.astype(np.uint16))
    #  !!!!OR, FOR SEPARATION OF FUNCTION, CREATE A STANDALONG FILE-WRITING FUNCTION   
    fnout=pathout+NH3file[0:25]+'_647NH3Abs'+imagetype
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
    hdul[0].header['DATE-OBS']=NH3time.replace('_','T')+'Z'
    #hdul[0].header['DATE']=NH3time
    hdul[0].header['AUTHOR']='Hill, S. M.'

    hdul[0].header['OBJECT']='Jupiter'
    hdul[0].header['TELESCOP']=sourcefiles[sourcedata]['Metadata']['Telescope']
    hdul[0].header['INSTRUME']=sourcefiles[sourcedata]['Metadata']['Camera']
    hdul[0].header['SEEING']=sourcefiles[sourcedata]['Metadata']['Seeing']
    hdul[0].header['TRANSPAR']=sourcefiles[sourcedata]['Metadata']['Transparency']
    hdul[0].header['BUNIT']=("TEST","comment")
    hdul[0].header['CALIBRA']=("TEST",'Disk-Integrated Cal Ref')
    hdul[0].header['VERSION']=('TBD','TBD')
    hdul[0].header['CTYPE1']=('Sys. 2 Longitude','deg')
    hdul[0].header['CTYPE2']=('PG Latitude','deg')
    
    hdul[0].header['CM1']=(NH3_CM1,'Sys. 1 Long. Central Meridian')
    hdul[0].header['CM2']=(NH3_CM2,'Sys. 1 Long. Central Meridian')
    hdul[0].header['CM3']=(NH3_CM3,'Sys. 1 Long. Central Meridian')
    hdul[0].header['SMOOTH']=('Smoothing','')
    hdul[0].header['KERNEL']=(1,'Gaussian')
    #hdul[0].header['CH4ABSFL']=(CH4file,'Source file for CH4 Absorption')
    hdul[0].header['NH3ABSFL']=(NH3file,'Source file for NH3 Absorption')
    hdul[0].header['CONTXTFL']=(RGBfile,'Source file for RGB Context')
    print(hdul[0].header[:])
    #fnout=path+'/'+NH3file[0:26]+'NH3Abs647.fits'
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
    # SET UP PLOT/CANVAS AND CREATE FIGURE HEADER
    ###########################################################################
        
    print(nny)
    if nny==3:
        ny,nx,dx,dy=3,3,6.0,6.0
    elif nny==4:
        ny,nx,dx,dy=4,3,5.0,10.0
    elif nny==5:
        ny,nx,dx,dy=5,3,5.0,13.0
    elif nny==6:
        ny,nx,dx,dy=6,3,5.0,10.0


    fig,ax=pl.subplots(ny,nx,figsize=(dx,dy), dpi=150, facecolor="black")#,
                          #sharey=True,sharex=True)
    fig.suptitle(NH3file[0:10],x=0.5,ha='center',color='w',fontsize=10)
    #fig.suptitle("Methane and Ammonia Experiment "+NH3file[0:10],
    #             x=0.5,ha='center',color='w',fontsize=12)
    ax[0,1].annotate("Steven Hill",[0.5,1.2],ha='center',xycoords='axes fraction',color='w',fontsize=8)
    ax[0,1].annotate("Denver, Colorado",[0.5,1.1],ha='center',xycoords='axes fraction',color='w',fontsize=8)
    
    ax[0,0].annotate(sourcefiles[sourcedata]['Metadata']['Telescope'],[0.07,1.3],
                 ha='left',xycoords='axes fraction',color='w',fontsize=8)
    ax[0,0].annotate("5.60m FL",[0.07,1.2],
                 ha='left',xycoords='axes fraction',color='w',fontsize=8)
    ax[0,0].annotate("ASI120MM",[0.07,1.1],
                 ha='left',xycoords='axes fraction',color='w',fontsize=8)
    ax[0,2].annotate("Seeing="+sourcefiles[sourcedata]['Metadata']['Seeing'],[0.95,1.3],
                 ha='right',xycoords='axes fraction',color='w',fontsize=8)
    ax[0,2].annotate("Transparency="+sourcefiles[sourcedata]['Metadata']['Transparency'],[0.95,1.2],
                 ha='right',xycoords='axes fraction',color='w',fontsize=8)
    
    ###########################################################################
    # COMPUTE AMMONIA IMAGE PANELS ALONG WITH COLOR SLOPE (VALID AT NH3 EPOCH)
    ###########################################################################
    temp=Ammonia_Panels(clrslp,nh3abs,CNT647,NH3_RGB,mask,AmmoniaHeader,ax) #INCLUDE "ROW" HERE BASED ON nx

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

            ch4abs=(np.array(CH4_RGB[:,:,1])+0.0001)/(CNTSynth+0.0001)
    
            ch4abs[ch4abs == inf] = 0
            ch4abs[ch4abs == -inf] = 0
            ch4abs[ch4abs == np.nan] = 0
            
            #nh3abs16bit = np.nan_to_num(((0.6667*65535*(ch4abs-0.5))*mask[:,:,0]).astype(np.uint16))
            #ch4abs16bit = np.nan_to_num(((5.*65535.*(ch4abs)-0.1)*mask[:,:,1]).astype(np.uint16))
            #imwrite(path+'/'+CH4file[0][0:26]+'CH4'+wvstr+'AbsPython.png', ch4abs16bit.astype(np.uint16))
            ###!!!!FIX CH4 FILE TIME IN NAMING!!!!

            if "G620" in CH4fl:
                fnout=pathout+CH4fl[0:25]+'_620CH4Abs'+imagetype
                try:
                    fnout=fnout+sourcefiles[sourcedata]['Metadata']['Variation']
                except:
                    print('No Variation')

                maskedch4=ch4abs*mask[:,:,1]
                #normch4=maskedch4/maskednh3[maskedch4>0].mean()
                normch4=maskedch4#/maskednh3[maskedch4>0].mean()
                
                hdu = fits.PrimaryHDU(normch4.astype(np.float32))
                hdul = fits.HDUList([hdu])
                hdul[0].header['BITPIX']=-32
                hdul[0].header['DATE-OBS']=NH3time.replace('_','T')+'Z'
                #hdul[0].header['DATE']=NH3time
                hdul[0].header['AUTHOR']='Hill, S. M.'
            
                hdul[0].header['OBJECT']='Jupiter'
                hdul[0].header['TELESCOP']=sourcefiles[sourcedata]['Metadata']['Telescope']
                hdul[0].header['INSTRUME']=sourcefiles[sourcedata]['Metadata']['Camera']
                hdul[0].header['SEEING']=sourcefiles[sourcedata]['Metadata']['Seeing']
                hdul[0].header['TRANSPAR']=sourcefiles[sourcedata]['Metadata']['Transparency']
                hdul[0].header['BUNIT']=("TEST","comment")
                hdul[0].header['CALIBRA']=("TEST",'Disk-Integrated Cal Ref')
                hdul[0].header['VERSION']=('TBD','TBD')
                hdul[0].header['CTYPE1']=('Sys. 2 Longitude','deg')
                hdul[0].header['CTYPE2']=('PG Latitude','deg')
                
                hdul[0].header['CM1']=(NH3_CM1,'Sys. 1 Long. Central Meridian')
                hdul[0].header['CM2']=(NH3_CM2,'Sys. 1 Long. Central Meridian')
                hdul[0].header['CM3']=(NH3_CM3,'Sys. 1 Long. Central Meridian')
                hdul[0].header['SMOOTH']=('Smoothing','')
                hdul[0].header['KERNEL']=(1,'Gaussian')
                #hdul[0].header['CH4ABSFL']=(CH4file,'Source file for CH4 Absorption')
                hdul[0].header['NH3ABSFL']=(NH3file,'Source file for NH3 Absorption')
                hdul[0].header['CONTXTFL']=(RGBfile,'Source file for RGB Context')
                print(hdul[0].header[:])                
                #print(CH4file[ifile][0:26])
                try:
                    os.remove(fnout+'.fits')
                except: 
                    print("file doesn't exist")
                hdul.writeto(fnout+'.fits')
                hdul.close()

                normch4scaled=np.nan_to_num(((5.*65535.*(normch4*mask[:,:,1] - 0.9))))
                normch4scaled[normch4scaled<=0.]=0.0
                ch4abs16bit = normch4scaled.astype(np.uint16)
                imwrite(fnout+'.png', ch4abs16bit)#.astype(np.uint16))


            temp=Methane_Panels(ch4abs,CNTSynth,CH4_RGB,mask,MethaneHeader,ax,
                                maxmin,ifile,CH4labels=CH4labels) #INCLUDE "ROW" HERE BASED ON nx
           
    # CONTEXT PANELS (3)
    if NUVfile != 'NA' or RGBfile != 'NA':
        try:
            Contextlabels = sourcefiles[sourcedata]['Contextlabels']
        except KeyError:
            Contextlabels=['380nm','RGB']

        tmp=Context_Panels(ax,ny,NUVFile=NUVfile,NUV_RGB=NUV_RGB,RGBFile=RGBfile,RGB_RGB=RGB_RGB,mask=[],
                           Contextlabels=Contextlabels) #INCLUDE "ROW" HERE BASED ON nx

    for i in range(0,3):
        for j in range(0,3):
            print(i,j)
            ax[i,j].axis('off')

    fig.subplots_adjust(left=0.00, bottom=0.05, right=1.0, top=0.90,
                wspace=0.001, hspace=0.10)
    pathmaparraybyfilter=\
        'C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Analysis Data/Map Array by Filter/'

    fig.savefig(pathmaparraybyfilter+obsdate+'-'+target+"_"+
                imagetype+"_SMHill.png",dpi=150,bbox_inches = 'tight')
    fig.savefig(pathmaparraybyfilter+obsdate+'-'+target+"_"+
                imagetype+"_SMHill.jpg",dpi=150,bbox_inches = 'tight')
    #imageio.imwrite(path+"PublicImageArray16bit"+obsdate+".png",fig.astype(np.uint16))
    
    #figx,axx=pl.subplots(1,1,figsize=(4,4), dpi=150, facecolor="black",
    #                      sharex=False)
    #ch4abs_scaled=(np.array(ch4abs)-1.0)/10.0+1.0
    #axx.imshow((np.array(nh3abs)/np.array(ch4abs_scaled))*np.array(mask[:,:,0]),'gist_gray')

###############################################################################
###############################################################################
    
def Ammonia_Panels(clrslp,nh3abs,CNT647,NH3_RGB,mask,AmmoniaHeader,ax,
                   NH3labels=['656nm (Red Cont.)','647nm (NH3)','632nm (Blue Cont.)']):
    ###########################################################################
    # PURPOSE: COMPUTE AND PLOT PANELS FOR NH3 ANALYSIS
    # INPUTS:  NH3 FILE NAME, NH3_RGB PNG FILE, mask, and axis system to plot 
    #          on
    # OUTPUTS: FLOATING POINT ARRAYS OF THE COLOR SLOPE AND THE CONTINUUM
    #          DIVIDED 647NM AMMONIA OBSERVATIONS
    # ASSUMES: INPUT IS A NORMALIZED (WHITE BALANCED) RGB ARRAY WITH
    #          R(656), G(647), AND B(632)
    # IMPROVEMENTS:
    #          !!!!SPLIT COMPUTATIONS FROM PLOTTING INTO SEPARATE FUNCTIONS
    #          !!!!COMBINE FILE WRITING WITH COMPUTATIONS
    #          !!!!BETTER PLOT SCALING
    import numpy as np
    #from astropy.io import fits
    #import pylab as pl
    #import sharpen as sharp
    ###########################################################################
    # PLOT 1ST ROW: CONTINUUM R(656), B(632), COLOR SLOPE AND ANNOTATE
    #   !!!!SCALING IS MANUAL HERE, PERHAPS IT SHOULD BE AUTOMATED?
    #   !!!!NOT SURE HOW SHARPENING IS FUNCTIONING HERE, IT DOESN'T SEEM TO 
    #       BE IMPORTED. MAYBE IT IS NATIVE TO imshow?
    ###########################################################################
    ax[0,0].imshow(sharpen(NH3_RGB[:,:,0]),'gist_gray',vmin=0.0,
                   vmax=np.nanmax(NH3_RGB[:,:,0])*1.4)
    ax[0,0].annotate(NH3labels[0],[0.5,0.03],color='white',ha='center',
                    xycoords='axes fraction',fontsize=8)
    ax[0,1].imshow(sharpen(NH3_RGB[:,:,2])*mask[:,:,1],'gist_gray',
                   vmin=0.0,vmax=np.nanmax(NH3_RGB[:,:,0])*1.4)
    ax[0,1].annotate(NH3labels[2],[0.5,0.03],color='white',ha='center',
                    xycoords='axes fraction',fontsize=8)
    # CUSTOM MASKING FOR COLOR SLOPE SINCE IT CAN BE POSITIVE AND NEGATIVE
    ax[0,2].imshow(clrslp*mask[:,:,1]+50*(mask[:,:,1]-1.0),'gist_gray',
                   vmin=-50.,vmax=50.)
    ax[0,2].annotate('656/632 (Color Slope)',[0.5,0.03],color='white',ha='center',
                    xycoords='axes fraction',fontsize=8)
    ###########################################################################
    # PLOT 2ND ROW: 647NM SYNTH CONTINUUM, 647NM OBSERVATION, CONTINUUM DIVIDED
    #   647NM NH3 OBSERVATION AND ANNOTATE
    #   !!!!SCALING IS MANUAL HERE, PERHAPS IT SHOULD BE AUTOMATED?
    #   !!!!NOT SURE HOW SHARPENING IS FUNCTIONING HERE, IT DOESN'T SEEM TO 
    #       BE IMPORTED. MAYBE IT IS NATIVE TO imshow?
    ###########################################################################       
    ax[1,0].imshow(sharpen(CNT647)*mask[:,:,1],'gist_gray',
                   vmin=0.0,vmax=np.nanmax(NH3_RGB[:,:,0])*1.4)
    ax[1,0].annotate('Synth. Continuum @ 647nm',[0.5,0.03],color='white',
                     ha='center',xycoords='axes fraction',fontsize=8)
    ax[1,1].imshow(sharpen(NH3_RGB[:,:,1]),'gist_gray',vmin=0.0,
                   vmax=np.nanmax(NH3_RGB[:,:,0])*1.4)
    ax[1,1].annotate(NH3labels[1],[0.5,0.03],color='white',ha='center',
                    xycoords='axes fraction',fontsize=8)
    ax[1,2].imshow(nh3abs*mask[:,:,1],'gist_gray',vmin=0.9,vmax=1.1)
    ax[1,2].annotate('647/Cont. (NH3)',[0.51,0.03],color='white',ha='center',
                    xycoords='axes fraction',fontsize=8)
    ###########################################################################
    # ANNOTATE CONTINUUM ROW AND AMMONIA ROW OF PANELS
    ###########################################################################
    ax[0,1].annotate("CONTINUUM",[-1.2,0.93],color='white',ha='left',
                    xycoords='axes fraction',fontsize=9)
    ax[0,1].annotate(AmmoniaHeader,[2.2,0.93],color='white',ha='right',
                    xycoords='axes fraction',fontsize=8)
    ax[1,1].annotate(r"AMMONIA (NH$_3$)",[-1.2,0.93],color='white',ha='left',
                    xycoords='axes fraction',fontsize=9)
    ax[1,1].annotate(AmmoniaHeader,[2.2,0.93],color='white',ha='right',
                    xycoords='axes fraction',fontsize=8)
    ax[0,1].set_zorder(1) #required so not blocked by the ax[0,2] image
    ax[1,1].set_zorder(1) #required so not blocked by the ax[0,2] image
    #see https://stackoverflow.com/questions/29735743/getting-text-to-display-in-front-of-subplot-images
    
    return clrslp,nh3abs
    
def Methane_Panels(ch4abs,CNTSynth,CH4_RGB,mask,MethaneHeader,ax,maxmin,ifile,
                   CH4labels=['656nm (Cont.)','889nm (CH4)','889/Cont. (CH4))']):
    import numpy as np
    #import sharpen as sharp
    
    ax[2+ifile,0].imshow(sharpen(CNTSynth)*mask[:,:,1],'gist_gray',vmin=0.0,vmax=np.nanmax(CNTSynth)*1.2)
    ax[2+ifile,1].imshow(sharpen(CH4_RGB[:,:,1])*mask[:,:,1],'gist_gray',vmin=0.0,vmax=np.nanmax(CH4_RGB[:,:,0])*1.2)

    ax[2+ifile,2].imshow(ch4abs*mask[:,:,1],'gist_gray',vmin=maxmin[0],vmax=maxmin[1])
    
    for icol in range(0,3):
        ax[2+ifile,icol].annotate(CH4labels[icol],[0.5,0.03],color='white',ha='center',
                        xycoords='axes fraction',fontsize=8)
    
    ax[2+ifile,1].annotate(r"METHANE (CH$_4$)",[-1.2,0.93],color='white',ha='left',
                    xycoords='axes fraction',fontsize=9)
    ax[2+ifile,1].annotate(MethaneHeader,[2.2,0.93],color='white',ha='right',
                    xycoords='axes fraction',fontsize=8)
    ax[2+ifile,1].set_zorder(1) #required so not blocked by the ax[2,2] image
    #see https://stackoverflow.com/questions/29735743/getting-text-to-display-in-front-of-subplot-images
    return ch4abs

def Context_Panels(ax,ny,NUVFile='NA',NUV_RGB=[],RGBFile='NA',RGB_RGB=[],mask=[],
                   Contextlabels=['380nm','RGB']):
    import numpy as np
    import get_WINJupos_ephem as WJ_ephem
    if NUVFile != 'NA': 
        sec=str(int(str(NUVFile[16:17]))*6)
        NUVtime=(NUVFile[0:10]+"_"+NUVFile[11:13]+":"+NUVFile[13:15]+":"+sec.zfill(2))
        eph=WJ_ephem.get_WINJupos_ephem(NUVtime)
        if np.array(mask).size != 0:
            ax[ny-1,0].imshow(np.array(NUV_RGB[:,:,0])*np.array(mask[:,:,0]),'gist_gray')
        else:
            ax[ny-1,0].imshow(NUV_RGB[:,:,0],'gist_gray')
        ax[ny-1,0].annotate(Contextlabels[0]+'; '+NUVtime[11:19]+" UT; Alt"+eph[3],[0.5,0.03],color='white',ha='center',
                        xycoords='axes fraction',fontsize=8)
        ax[ny-1,0].annotate("CM1"+eph[0]+"; CM2"+eph[1]+"; CM3"+eph[2],
                         [0.5,-0.05],color='white',ha='center',
                        xycoords='axes fraction',fontsize=7)
        ax[ny-1,1].imshow(np.zeros(NUV_RGB.shape))
        ax[ny-1,2].imshow(np.zeros(NUV_RGB.shape))
        
    print("RGBFile=",RGBFile)
    if RGBFile != 'NA': 
        sec=str(int(str(RGBFile[16:17]))*6)
        RGBtime=(RGBFile[0:10]+"_"+RGBFile[11:13]+":"+RGBFile[13:15]+":"+sec.zfill(2))
        eph=WJ_ephem.get_WINJupos_ephem(RGBtime)
        ax[ny-1,1].imshow(RGB_RGB)
        ax[ny-1,1].annotate(Contextlabels[1]+'; '+RGBtime[11:19]+" UT; Alt"+eph[3],[0.5,0.03],color='white',ha='center',
                        xycoords='axes fraction',fontsize=8)
        ax[ny-1,1].annotate("CM1"+eph[0]+"; CM2"+eph[1]+"; CM3"+eph[2],
                         [0.5,-0.05],color='white',ha='center',
                        xycoords='axes fraction',fontsize=7)
    
        ax[ny-1,0].imshow(np.zeros(RGB_RGB.shape))
        ax[ny-1,2].imshow(np.zeros(RGB_RGB.shape))
    
    ContextHeader="CONTEXT"
    ax[ny-1,1].annotate(ContextHeader,[-1.2,0.93],color='white',ha='left',
                    xycoords='axes fraction',fontsize=9)

    ax[ny-1,2].annotate("NOTES",[0.5,0.93],color='white',ha='center',
                     xycoords='axes fraction',fontsize=8)
    ax[ny-1,2].annotate("Ammonia and methane \n"+
                     "absorption are explored \n"+
                     "using in-band images \n"+
                     "divided by continuum \n"+
                     "images. Computed ratio \n"
                     "images are not sharpened \n"
                     "to minimize artifacts.",
                     [0.0,0.85],color='white',ha='left',va='top',
                     xycoords='axes fraction',fontsize=7)

    ax[ny-1,1].set_zorder(1) #required so not blocked by the ax[2,2] image
    #see https://stackoverflow.com/questions/29735743/getting-text-to-display-in-front-of-subplot-images 

def sharpen(image):
    from astropy.convolution import Gaussian2DKernel
    from astropy.convolution import convolve
    kernel = Gaussian2DKernel(x_stddev=5)
    blurred=convolve(image,kernel)
    tst=image+0.99*(image-blurred) #Need to revalidate that this is correct
    return tst