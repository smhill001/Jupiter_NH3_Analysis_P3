def image_array_new(obsdate="20240925UTa",target="Jupiter",
                     imagetype='Img',contour=False):
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
    sys.path.append('./Maps')
    import read_fits_map_L2_L3 as RFM

    ###########################################################################
    #  DATA FILES AND METADATA DICTIONARY
    ###########################################################################
    sourcedata=obsdate#+"_"+imagetype
    sourcefiles=getlist.get_obs_list(planet=target)
    
    if sourcefiles[sourcedata]["Telescope"]=="C11":
        CalModel="SCT-Obs-Final"
    elif sourcefiles[sourcedata]["Telescope"]=="VLT":
        CalModel="VLT-Obs-Final"

    calibration,K_eff=RMC.read_master_calibration()
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
    NH3eph=WJ_ephem.get_WINJupos_ephem(NH3time,planet=target)
    NH3_CM1=float(NH3eph[0].strip())
    NH3_CM2=float(NH3eph[1].strip())
    NH3_CM3=float(NH3eph[2].strip())

    ###########################################################################
    # EPHEMERIS AND CH4 TRANSMISSION CALCULATIONS
    ###########################################################################
    CH4sec=str(int(str(CH4file[16:17]))*6) #COMPUTE FROM FRACTIONAL WINJUPOS MINUTE
    CH4time=(CH4file[0:10]+"_"+CH4file[11:13]+":"+CH4file[13:15]+":"
             +CH4sec.zfill(2))
    CH4eph=WJ_ephem.get_WINJupos_ephem(CH4time,planet=target)
    CH4_CM1=float(CH4eph[0].strip())
    CH4_CM2=float(CH4eph[1].strip())
    CH4_CM3=float(CH4eph[2].strip()) 

    ###########################################################################
    # OPEN AND READ DATA FILES (FITS MAPS of NH3 and CH4 Transmission)
    ###########################################################################   
    pathRGB='c:/Astronomy/Projects/Planets/'+target+'/Imaging Data/'+obsdate[0:10]+'/'
    pathL2='C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Analysis Data/L2 FITS/'
    pathL3='C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Analysis Data/L3 FITS/'

    ###########################################################################
    # CH4 Transmission File name and read
    try:    #Set up to allow for parametric studies of different processing paths
        CH4file=sourcefiles[sourcedata]['CH4file'][0:18]+target+"_"+imagetype+"_L2TCH4"+\
                sourcefiles[sourcedata]['Variation']+".fits"
        variation=sourcefiles[sourcedata]['Variation']
    except:
        CH4file=sourcefiles[sourcedata]['CH4file'][0:18]+target+"_"+imagetype+"_L2TCH4.fits"
        variation=""

    CH4hdulist=fits.open(pathL2+CH4file)
    CH4hdulist.info()
    CH4hdr=CH4hdulist[0].header
    CH4data=CH4hdulist[0].data
    CH4hdulist.close()
    
    ###########################################################################
    # NH3 Transmission File name and read
    try:    #Set up to allow for parametric studies of different processing paths
        NH3file=sourcefiles[sourcedata]['NH3file'][0:18]+target+"_"+imagetype+"_L2TNH3"+\
                sourcefiles[sourcedata]['Variation']+".fits"
        variation=sourcefiles[sourcedata]['Variation']
    except:
        NH3file=sourcefiles[sourcedata]['NH3file'][0:18]+target+"_"+imagetype+"_L2TNH3.fits"
        variation=""
    
    NH3hdulist=fits.open(pathL2+NH3file)
    NH3hdulist.info()
    NH3hdr=NH3hdulist[0].header
    NH3data=NH3hdulist[0].data
    sza=NH3hdulist[1].data
    eza=NH3hdulist[2].data
    NH3hdulist.close()
    
    PCloudhdr,PClouddata,fNH3hdr,fNH3data,sza,eza,RGB,RGB_CM,RGBtime= \
                    RFM.read_fits_map_L2_L3(obskey=obsdate,target=target,
                                            imagetype="Img",Level="L3")
                    
    tx=np.linspace(0,65500.,5,endpoint=True)

    fig,axs=pl.subplots(3,3,figsize=(7.0,7.0), dpi=150, facecolor="black")
    axs[0,0].imshow(CH4_RGB[:,:,1],'gray')
    axs[0,0].set_title("a) CH4 Radiance",color='1.0',fontsize=9)
    axs[0,0].annotate(CH4time.replace("_","T")+"Z",[0.5,0.95],ha='center',
                  xycoords="axes fraction",fontsize=8,color='1.0')

    
    axs[1,0].imshow(NH3_RGB[:,:,1],'gray')
    axs[1,0].set_title("d) NH3 Radiance",color='1.0',fontsize=9)
    axs[1,0].annotate(NH3time.replace("_","T")+"Z",[0.5,0.95],ha='center',
                  xycoords="axes fraction",fontsize=8,color='1.0')

    print(CH4hdr["RANGE0"],CH4hdr["RANGE1"])
    temp = ((CH4data - CH4hdr["RANGE0"]) /(CH4hdr["RANGE1"] - CH4hdr["RANGE0"]) * 65534.9999)#
    temp[temp < 0] = 0
    CH4scl=temp.astype('uint16')
    show=axs[0,1].imshow(CH4scl,'gray')
    if contour:
        axs[0,1].contourf(CH4data,levels=[CH4hdr["RANGE1"],1.2],colors=['C0'],alpha=0.3)
    axs[0,1].set_title("b) CH4 Transmission",color='1.0',fontsize=9)
    axs[0,1].annotate(CH4hdr["DATE-OBS"],[0.5,0.95],ha='center',
                  xycoords="axes fraction",fontsize=8,color='1.0')
    cbar = pl.colorbar(show,ticks=tx,
               orientation='vertical',cmap='gray',
               ax=axs[0,1],fraction=0.046, pad=0.04)
    txlb=np.linspace(CH4hdr["RANGE0"],CH4hdr["RANGE1"],5,endpoint=True)
    cbar.ax.set_yticklabels(np.around(txlb,3))
    cbar.ax.tick_params(labelsize=6,colors='1.0')

    print(NH3hdr["RANGE0"],NH3hdr["RANGE1"])
    temp = ((NH3data - NH3hdr["RANGE0"]) * (1/(NH3hdr["RANGE1"] - NH3hdr["RANGE0"]) * 65534.9999))#
    temp[temp < 0] = 0
    NH3scl=temp.astype('uint16')
    show=axs[1,1].imshow(NH3scl,'gray')
    if contour:
        axs[1,1].contourf(NH3data,levels=[NH3hdr["RANGE1"],1.2],colors=['C0'],alpha=0.4)
    axs[1,1].set_title("e) NH3 Transmission",color='1.0',fontsize=9)
    axs[1,1].annotate(NH3hdr["DATE-OBS"],[0.5,0.95],ha='center',
                  xycoords="axes fraction",fontsize=8,color='1.0')
    cbar = pl.colorbar(show,ticks=tx,
               orientation='vertical',cmap='gray',
               ax=axs[1,1],fraction=0.046, pad=0.04)
    txlb=np.linspace(NH3hdr["RANGE0"],NH3hdr["RANGE1"],5,endpoint=True)
    cbar.ax.set_yticklabels(np.around(txlb,3))
    cbar.ax.tick_params(labelsize=6,colors='1.0')

    CH4tau=-np.log(CH4data)
    CH4tau=-np.log(CH4data)
    CH4tau[CH4tau>0.9]=0.0
    temp = ((CH4tau - (-np.log(CH4hdr["RANGE1"]))) * (1/(-np.log((CH4hdr["RANGE0"])) - -np.log(CH4hdr["RANGE1"])) * 65534.9999))#
    temp[temp < 0] = 0
    CH4scl=temp.astype('uint16')
    show=axs[0,2].imshow(CH4scl,'gray')
    if contour:
        axs[0,2].contourf(CH4data,levels=[CH4hdr["RANGE1"],1.2],colors=['C0'],alpha=0.4)
    axs[0,2].set_title("c) CH4 Opacity",color='1.0',fontsize=9)
    axs[0,2].annotate(CH4hdr["DATE-OBS"],[0.5,0.95],ha='center',
                  xycoords="axes fraction",fontsize=8,color='1.0')
    cbar = pl.colorbar(show,ticks=tx,
               orientation='vertical',cmap='gray',
               ax=axs[0,2],fraction=0.046, pad=0.04)
    txlb=np.linspace(-np.log(CH4hdr["RANGE1"]),-np.log(CH4hdr["RANGE0"]),5,endpoint=True)
    cbar.ax.set_yticklabels(np.around(txlb,3))
    cbar.ax.tick_params(labelsize=6,colors='1.0')

    NH3tau=-np.log(NH3data)
    NH3tau[NH3tau>0.9]=0.0
    temp = ((NH3tau - (-np.log(NH3hdr["RANGE1"]))) * (1/(-np.log((NH3hdr["RANGE0"])) - -np.log(NH3hdr["RANGE1"])) * 65534.9999))#
    temp[temp < 0] = 0
    NH3scl=temp.astype('uint16')
    show=axs[1,2].imshow(NH3scl,'gray')
    if contour:
        axs[1,2].contourf(NH3data,levels=[NH3hdr["RANGE1"],1.2],colors=['C0'],alpha=0.4)
    axs[1,2].set_title("f) NH3 Opacity",color='1.0',fontsize=9)
    axs[1,2].annotate(NH3hdr["DATE-OBS"],[0.5,0.95],ha='center',
                  xycoords="axes fraction",fontsize=8,color='1.0')
    cbar = pl.colorbar(show,ticks=tx,
               orientation='vertical',cmap='gray',
               ax=axs[1,2],fraction=0.046, pad=0.04)
    txlb=np.linspace(-np.log(NH3hdr["RANGE1"]),-np.log(NH3hdr["RANGE0"]),5,endpoint=True)
    cbar.ax.set_yticklabels(np.around(txlb,3))
    cbar.ax.tick_params(labelsize=6,colors='1.0')


    temp = ((PClouddata - PCloudhdr["RANGE0"]) * (1/(PCloudhdr["RANGE1"] - PCloudhdr["RANGE0"]) * 65534.9999))#
    temp[temp < 0] = 0
    Pclscl=temp.astype('uint16')
    Pclscl[Pclscl==0]=65535
    show=axs[2,0].imshow(Pclscl,'gray_r')
    axs[2,0].set_title("g) Effective Cloud Top (mb)",color='1.0',fontsize=9)
    axs[2,0].annotate(PCloudhdr["DATE-OBS"],[0.5,0.95],ha='center',
                  xycoords="axes fraction",fontsize=8,color='1.0')
    if contour:
        axs[2,0].contourf(CH4data,levels=[CH4hdr["RANGE1"],1.2],colors=['C0'],alpha=0.4)

    cbar = pl.colorbar(show,ticks=tx,
               orientation='vertical',cmap='gray',
               ax=axs[2,0],fraction=0.046, pad=0.04)
    txlb=np.linspace(PCloudhdr["RANGE0"],PCloudhdr["RANGE1"],5,endpoint=True)
    print(txlb)
    cbar.ax.set_yticklabels(np.around(txlb,3))
    cbar.ax.tick_params(labelsize=6,colors='1.0')
    
    axs[2,1].imshow(RGB)
    axs[2,1].set_title("h) RGB Context",color='1.0',fontsize=9)
    axs[2,1].annotate(RGBtime.replace("_","T")+"Z",[0.5,0.95],ha='center',
                  xycoords="axes fraction",fontsize=8,color='1.0')

    temp = ((fNH3data - fNH3hdr["RANGE0"]) * (1/(fNH3hdr["RANGE1"] - fNH3hdr["RANGE0"]) * 65534.9999))#
    temp[temp < 0] = 0
    NH3scl=temp.astype('uint16')
    show=axs[2,2].imshow(NH3scl,'gray')
    axs[2,2].set_title("i) NH3 Mole Fraction (ppm)",color='1.0',fontsize=9)
    axs[2,2].annotate(fNH3hdr["DATE-OBS"],[0.5,0.95],ha='center',
                  xycoords="axes fraction",fontsize=8,color='1.0')
    if contour:
        axs[2,2].contourf(CH4data,levels=[CH4hdr["RANGE1"],1.2],colors=['C0'],alpha=0.4)
        axs[2,2].contourf(NH3data,levels=[NH3hdr["RANGE1"],1.2],colors=['C0'],alpha=0.4)

    cbar = pl.colorbar(show,ticks=tx,
               orientation='vertical',cmap='gray',
               ax=axs[2,2],fraction=0.046, pad=0.04)
    txlb=np.linspace(fNH3hdr["RANGE0"],fNH3hdr["RANGE1"],5,endpoint=True)
    print(txlb)
    cbar.ax.set_yticklabels(np.around(txlb,3))
    cbar.ax.tick_params(labelsize=6,colors='1.0')
    
    pathout='C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Analysis Data/Img Plots Diagnostic/'

    fnout=obsdate
    #imwrite(pathout+fnout+'.png', scl_arr)#.astype(np.uint16))
    fig.suptitle(obsdate,color='1.0')
    
    fig.subplots_adjust(left=0.05, bottom=0.05, right=0.90, top=0.90,
                wspace=0.25, hspace=0.15)     

    fig.savefig(pathout+fnout+'.png',dpi=300,bbox_inches = 'tight')

    return(fig)
