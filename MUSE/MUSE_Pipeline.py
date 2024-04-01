def MUSE_Pipeline(date):
    """
    Created on Thu Jan 12 14:21:05 2023
    
    PROGRAM: ReadMUSE
    PURPOSE: Read MUSE fits file for Jupiter observation on 9/19/2022 and
                1) Create synthetic images for Ammonia abundance analysis
                    (620, 632, 647, and 656 nm)
                2) Create synthetic RGB images
                3) Create full disk radiance spectrum
                4) Create estimated albedo spectrum for comparison to Karkoschka
                5) Write filter files for synthetic filter bands
    @author: smhil
    """
    
    import sys
    drive='c:'
    sys.path.append(drive+'/Astronomy/Python Play')
    sys.path.append(drive+'/Astronomy/Python Play/Util_P3')
    sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy_P3')
    
    import os
    import pylab as pl
    import numpy as np
    from imageio import imwrite
    from astropy.io import fits
    import GeneralSpecUtils_P3 as GSU
    sys.path.append('./MUSE')
    import ReadMUSE as RM
    import MUSE_Spectrum as MS
    from copy import deepcopy
    import time
    start_time=time.time()

    ###########################################################################
    # Read MUSE data cube and filter metadata
    ###########################################################################
    MUSEhdr,MUSEhdr1,MUSEdata,MUSEzen,MUSEszen,wavelength,filterdata,path=RM.ReadMUSE(date)

    ###############################################################################
    # Compute spectrum, smoothed spectrum, albedo spectrum, and write filter files
    ###############################################################################
    temp=MS.MUSE_Spectrum(date,MUSEhdr,MUSEdata,wavelength,filterdata,path)

    ###############################################################################
    # Compute air mass factor as by pixel from emissiong and incidence angles
    ###############################################################################
    """
    fig,axs=pl.subplots(2,2,figsize=(7.0,4.5), dpi=150, facecolor="white",
                        sharey=True,sharex=True)      
    fig.suptitle("MUSE Ammonia Analysis Synthetic Band Images",ha='center')
    
    amf=1./np.cos(MUSEzen/180*np.pi)+1./np.cos(MUSEszen/180*np.pi)
    
    axs[0,0].imshow(amf,"gray",origin='lower')
    axs[0,0].set_title("AMF")
    
    hdu = fits.PrimaryHDU(amf.astype(np.float32))
    hdul = fits.HDUList([hdu])
    fnout=path+"Test0730amf.fits"
    try:
        os.remove(fnout)
    except: 
        print("file doesn't exist")
    hdul.writeto(fnout)
    hdul.close()
    amfscaled=np.nan_to_num(65535.*amf*0.1)
    amfscaled[amfscaled<=0.]=0.0
    amf16bit = np.flipud(amfscaled.astype(np.uint16))
    imwrite(path+'amf.png', amf16bit)
    """

    ###############################################################################
    # Load MUSE Image and Smooth if date=20220919UT
    ###############################################################################

    dateobs=MUSEhdr['DATE-OBS']
    fn=dateobs[0:10]+'-'+dateobs[11:13]+dateobs[14:16]+'_'+str(float(dateobs[17:19])/60.)[2]+'-Jupiter-'
    
    print(dateobs)
    print(fn)
    print(MUSEhdr['DATE-OBS'])
    print(MUSEhdr['OBJECT'])

    #filtpath='C:/Astronomy/Projects/Techniques/InstrumentPerformance-P3/Filters/'
    #filterdata['620']['FiltTrans']=np.loadtxt(filtpath+filterdata['620']['transfile'],usecols=range(2))


    ###########################################################################
    #!!!LOOP OVER MUSE_SPECTRUM COLUMN BY COLUMN TO BUILD SMOOTHED, FILTER-
    #!!! CONVOLVED IMAGE ARRAYS
    ###########################################################################
    if date=="20220919UT":
        kernel_size = 14
        kernel = np.ones(kernel_size) / kernel_size
        
        MUSE_convolved0=deepcopy(MUSEdata)
        MUSE_convolved=deepcopy(MUSEdata)
        for i in range(0,MUSEhdr1['NAXIS1']):
            for j in range(0,MUSEhdr1['NAXIS2']):
                MUSE_convolved0[:,i,j] = np.convolve(MUSEdata[:,i,j], kernel,
                                                     mode='same')
                MUSE_convolved[:,i,j] = np.convolve(MUSE_convolved0[:,i,j], 
                                                    kernel, mode='same')
                WaveGrid,SignalonGrid=GSU.uniform_wave_grid(wavelength,
                                                            MUSE_convolved[:,i,j],
                                                            Extend=False,Fine=False)

    elif date=="20220730UT":
        WaveGrid,SignalonGrid=GSU.uniform_wave_grid(wavelength,MUSESpec,Extend=False,Fine=False)        
    MuseSpecGrid=np.column_stack((WaveGrid,SignalonGrid))
    print("MuseSpecGrid.shape=",MuseSpecGrid.shape)


    F620_wvs=filterdata["620"]["wvs"]
    F620_idx=[np.argmin(abs(F620_wvs[0]-wavelength)),np.argmin(abs(F620_wvs[1]-wavelength))]
    MUSE620=np.mean(MUSEdata[F620_idx[0]:F620_idx[1],:,:],axis=0)
    for i in range(0,MUSEhdr1['NAXIS1']):
        for j in range(0,MUSEhdr1['NAXIS2']):
            MUSE620_ProductV=GSU.SpectrumMath(FilterTransmission620,MUSE620[:,i,j],"Multiply")
            MUSE620[i,j]=sum(MUSE620_ProductV[F620_idx[0]:F620_idx[1],1])/ \
                            sum(FilterTransmission620[F620_idx[0]:F620_idx[1],1])

            MUSE_convolved0[:,i,j] = np.convolve(MUSEdata[:,i,j], kernel, mode='same')
            MUSE_convolved[:,i,j] = np.convolve(MUSE_convolved0[:,i,j], kernel, mode='same')
    MUSE620scaled=np.nan_to_num(65535.*MUSE620*100.)
    MUSE620scaled[MUSE620scaled<=0.]=0.0
    MUSE620abs16bit = np.flipud(MUSE620scaled.astype(np.uint16))
    imwrite(path+fn+'MUSE620-10nm.png', MUSE620abs16bit)
    
    F632_wvs=filterdata["632"]["wvs"]
    F632_idx=[np.argmin(abs(F632_wvs[0]-wavelength)),np.argmin(abs(F632_wvs[1]-wavelength))]
    MUSE632=np.mean(MUSEdata[F632_idx[0]:F632_idx[1],:,:],axis=0)
    MUSE632scaled=np.nan_to_num(65535.*MUSE632*100.)
    MUSE632scaled[MUSE632scaled<=0.]=0.0
    MUSE632abs16bit = np.flipud(MUSE632scaled.astype(np.uint16))
    imwrite(path+fn+'MUSE632-10nm.png', MUSE632abs16bit)
    
    F647_wvs=filterdata["647"]["wvs"]
    F647_idx=[np.argmin(abs(F647_wvs[0]-wavelength)),np.argmin(abs(F647_wvs[1]-wavelength))]
    MUSE647=np.mean(MUSEdata[F647_idx[0]:F647_idx[1],:,:],axis=0)
    MUSE647scaled=np.nan_to_num(65535.*MUSE647*100.)
    MUSE647scaled[MUSE647scaled<=0.]=0.0
    MUSE647abs16bit = np.flipud(MUSE647scaled.astype(np.uint16))
    imwrite(path+fn+'MUSE647-10nm.png', MUSE647abs16bit)
    
    F656_wvsBLU=filterdata["656"]["wvsBLU"]
    F656_idxBLU=[np.argmin(abs(F656_wvsBLU[0]-wavelength)),np.argmin(abs(F656_wvsBLU[1]-wavelength))]
    MUSE656BLU=np.mean(MUSEdata[F656_idxBLU[0]:F656_idxBLU[1],:,:],axis=0)
    F656_wvsRED=filterdata["656"]["wvsRED"]
    F656_idxRED=[np.argmin(abs(F656_wvsRED[0]-wavelength)),np.argmin(abs(F656_wvsRED[1]-wavelength))]
    MUSE656RED=np.mean(MUSEdata[F656_idxRED[0]:F656_idxRED[1],:,:],axis=0)
    MUSE656=(MUSE656BLU+MUSE656RED)/2.
    
    MUSE656scaled=np.nan_to_num(65535.*MUSE656*100.)
    MUSE656scaled[MUSE656scaled<=0.]=0.0
    MUSE656abs16bit = np.flipud(MUSE656scaled.astype(np.uint16))
    imwrite(path+fn+'MUSE656-6nmSplit.png', MUSE656abs16bit)


    
    fig,axs=pl.subplots(2,2,figsize=(7.0,4.5), dpi=150, facecolor="white",
                        sharey=True,sharex=True)      
    fig.suptitle("MUSE Ammonia Analysis Synthetic Band Images",ha='center')
    
    axs[0,0].imshow(MUSE620,"gray",origin='lower')
    axs[0,0].set_title("615-625 nm")
    axs[0,1].imshow(MUSE632,"gray",origin='lower')
    axs[0,1].set_title("627-637 nm")
    axs[1,0].imshow(MUSE647,"gray",origin='lower')
    axs[1,0].set_title("642-652 nm")
    axs[1,1].imshow(MUSE656,"gray",origin='lower')
    axs[1,1].set_title("651-654, 658-661 nm")
    
    F450_wvs=[475.,480.]
    #F450_wvs=[471.,480.]
    F450_idx=[np.argmin(abs(F450_wvs[0]-wavelength)),np.argmin(abs(F450_wvs[1]-wavelength))]
    MUSE450=np.mean(MUSEdata[F450_idx[0]:F450_idx[1],:,:],axis=0)
    MUSE450scaled=np.nan_to_num(65535.*MUSE450*100.)
    MUSE450scaled[MUSE450scaled<=0.]=0.0
    MUSE450abs16bit = np.flipud(MUSE450scaled.astype(np.uint16))
    imwrite(path+fn+'MUSE450.png', MUSE450abs16bit)
    
    F550_wvs=[530.,580.]
    F550_idx=[np.argmin(abs(F550_wvs[0]-wavelength)),np.argmin(abs(F550_wvs[1]-wavelength))]
    MUSE550=np.mean(MUSEdata[F550_idx[0]:F550_idx[1],:,:],axis=0)
    MUSE550scaled=np.nan_to_num(65535.*MUSE550*100.)
    MUSE550scaled[MUSE550scaled<=0.]=0.0
    MUSE550abs16bit = np.flipud(MUSE550scaled.astype(np.uint16))
    imwrite(path+fn+'MUSE550.png', MUSE550abs16bit)
    
    F650_wvs=[630.,680.]
    F650_idx=[np.argmin(abs(F650_wvs[0]-wavelength)),np.argmin(abs(F650_wvs[1]-wavelength))]
    MUSE650=np.mean(MUSEdata[F650_idx[0]:F650_idx[1],:,:],axis=0)
    MUSE650scaled=np.nan_to_num(65535.*MUSE650*100.)
    MUSE650scaled[MUSE450scaled<=0.]=0.0
    MUSE650abs16bit = np.flipud(MUSE650scaled.astype(np.uint16))
    imwrite(path+fn+'MUSE650.png', MUSE650abs16bit)
    
    fig1,axs1=pl.subplots(2,2,figsize=(7.0,4.5), dpi=150, facecolor="white",
                        sharey=True,sharex=True)      
    fig1.suptitle("MUSE RGB Synthetic Band Images",ha='center')
    axs1[0,0].imshow(MUSE450,"gray",origin='lower')
    axs1[0,0].set_title("471-480 nm")
    axs1[0,1].imshow(MUSE550,"gray",origin='lower')
    axs1[0,1].set_title("530-580 nm")
    axs1[1,0].imshow(MUSE650,"gray",origin='lower')
    axs1[1,0].set_title("630-680 nm")
    
    print("--- %s seconds ---" % (time.time() - start_time))