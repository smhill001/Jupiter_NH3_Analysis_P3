def MUSE_Spectrum(date):
    """
    Created on Thu Jan 12 14:21:05 2023
    
    PROGRAM: MUSE_Spectrum
    PURPOSE: Read MUSE fits file for Jupiter observations on 9/19/2022 and
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
    import ReadMUSE as RM
    
    MUSEhdr,MUSEdata,MUSEzen,MUSEszen,wavelength,filtdata,path=RM.ReadMUSE(date)

    ###########################################################################
    # Compute spectrum and smooth
    ###########################################################################
    dt=MUSEhdr["DATE-OBS"][0:10]+"-"+MUSEhdr["DATE-OBS"][11:13]+\
       MUSEhdr["DATE-OBS"][14:16]+"_"+str(int(10*(int(MUSEhdr["DATE-OBS"][17:19])/60.)))
    print(dt)
    MUSESpec=np.nanmean(MUSEdata[:,:,:],axis=(1,2))
    print("MUSESpec.shape=",MUSESpec.shape)#,MUSESpec)
    
    kernel_size = 4
    kernel = np.ones(kernel_size) / kernel_size
    data_convolved = np.convolve(MUSESpec, kernel, mode='same')
    print("*******wavelength,data_convolved=",wavelength.shape,data_convolved.shape)
    
    ###########################################################################
    # Regrid spectra to 0.5 nm; use smoothed spectrum for 2022-09-19
    ###########################################################################
    if date=="20220919UT":
        WaveGrid,SignalonGrid=GSU.uniform_wave_grid(wavelength,data_convolved,Extend=False,Fine=False)
    elif date=="20220730UT":
        WaveGrid,SignalonGrid=GSU.uniform_wave_grid(wavelength,MUSESpec,Extend=False,Fine=False)        
    MuseSpecGrid=np.column_stack((WaveGrid,SignalonGrid))
    print("MuseSpecGrid.shape=",MuseSpecGrid.shape)
    
    ###########################################################################
    # Load Solar Spectrum Reference and compute Albedo Spectrum
    ###########################################################################
    RefPath="c:/Astronomy/Python Play/SPLibraries/SpectralReferenceFiles/ReferenceLibrary/"
    G2V=np.loadtxt(RefPath+"g2v.dat", dtype=float, usecols=(0,1))
    print("G2V.shape=",G2V.shape)
    AlbedoSpec=GSU.SpectrumMath(MuseSpecGrid,G2V,"Divide")
    
    ###########################################################################
    # Write output spectral data
    ###########################################################################
    pathout="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/MUSE/"
    np.savetxt(pathout+dt+'_MUSE_Albedo.txt',AlbedoSpec,delimiter=" ",
               fmt="%10.3F %10.7F")
    
    SmoothedSpec=np.column_stack((WaveGrid,SignalonGrid))
    np.savetxt(pathout+dt+'_MUSE_Spec.txt',SmoothedSpec,delimiter=" ",
               fmt="%10.3F %10.7F") 
    
    ###########################################################################
    # Compute mean signal 
    ###########################################################################    
    F620_idx=[np.argmin(abs(filtdata['620']['wvs'][0]-WaveGrid)),np.argmin(abs(filtdata['620']['wvs'][1]-WaveGrid))]
    MUSE620Albedo=np.mean(AlbedoSpec[F620_idx[0]:F620_idx[1],1],axis=0)
    print("MUSE620Albedo=",MUSE620Albedo)
    
    F632_idx=[np.argmin(abs(filtdata['632']['wvs'][0]-WaveGrid)),np.argmin(abs(filtdata['632']['wvs'][1]-WaveGrid))]
    MUSE632Albedo=np.mean(AlbedoSpec[F632_idx[0]:F632_idx[1],1],axis=0)
    print("MUSE632Albedo=",MUSE632Albedo)
    
    F647_idx=[np.argmin(abs(filtdata['647']['wvs'][0]-WaveGrid)),np.argmin(abs(filtdata['647']['wvs'][1]-WaveGrid))]
    MUSE647Albedo=np.mean(AlbedoSpec[F647_idx[0]:F647_idx[1],1],axis=0)
    print("MUSE647Albedo=",MUSE647Albedo)
    
    F656_idxRED=[np.argmin(abs(filtdata['656']['wvsRED'][0]-WaveGrid)),np.argmin(abs(filtdata['656']['wvsRED'][1]-WaveGrid))]
    MUSE656AlbedoRED=np.mean(AlbedoSpec[F656_idxRED[0]:F656_idxRED[1],1],axis=0)
    F656_idxBLU=[np.argmin(abs(filtdata['656']['wvsBLU'][0]-WaveGrid)),np.argmin(abs(filtdata['656']['wvsBLU'][1]-WaveGrid))]
    MUSE656AlbedoBLU=np.mean(AlbedoSpec[F656_idxBLU[0]:F656_idxBLU[1],1],axis=0)
    MUSE656Albedo=(MUSE656AlbedoRED+MUSE656AlbedoBLU)/2.
    print("MUSE656Albedo=",MUSE656Albedo)
    
    MUSERelSlope=(MUSE656Albedo-MUSE632Albedo)/(656.-632.)
    MUSE647Slope=MUSE632Albedo+MUSERelSlope*15.
    MUSE620Slope=MUSE632Albedo+MUSERelSlope*(-12.)
    MUSENH3Abs=MUSE647Albedo/MUSE647Slope
    MUSECH4Abs=MUSE620Albedo/MUSE620Slope
    
    print("")
    print("NH3 Abs=",MUSENH3Abs)
    print("CH4 Abs=",MUSECH4Abs)
    
    fig2,axs2=pl.subplots(1,1,figsize=(6.0,4.0), dpi=150, facecolor="white",
                        sharey=True,sharex=True)      
    axs2.plot(wavelength,MUSESpec,label='Raw MUSE Spectrum')
    axs2.plot(WaveGrid,SignalonGrid,label='Smoothed, Regridded to 0.5 nm)')
    axs2.plot(AlbedoSpec[:,0],AlbedoSpec[:,1],label="Albedo (G2V Reference: Pickles (1998))")
    axs2.set_xlim(450.,950.)
    axs2.set_ylim(0,0.001)
    axs2.set_xlabel("Wavelength (nm)")
    axs2.set_title("VLT-MUSE Spectrum "+dt)
    axs2.legend(fontsize=8)
    axs2.annotate("NH3 Abs="+str(MUSENH3Abs)[0:6],[0.01,0.95],
                  xycoords="axes fraction",fontsize=8)
    axs2.annotate("CH4 Abs="+str(MUSECH4Abs)[0:6],[0.01,0.90],
                  xycoords="axes fraction",fontsize=8)
    
    fig2.savefig(pathout+dt+'_MUSE_Spec.png')
