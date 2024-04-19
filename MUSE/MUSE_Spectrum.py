def MUSE_Spectrum(date):#,MUSEhdr,MUSEdata,wavelength,filterdata,path):
    """
    Created on Thu Jan 12 14:21:05 2023
    
    PROGRAM: MUSE_Spectrum
    PURPOSE: Read MUSE fits file for Jupiter observations on 9/19/2022 and
                1) tbd
    @author: smhil
    """
    
    import sys
    drive='c:'
    sys.path.append(drive+'/Astronomy/Python Play')
    sys.path.append(drive+'/Astronomy/Python Play/Util_P3')
    sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy_P3')
    sys.path.append(drive+'/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Molecular Absorption/code')
    
    import os
    import pylab as pl
    import numpy as np
    from imageio import imwrite
    from astropy.io import fits
    import GeneralSpecUtils_P3 as GSU
    import ReadMUSE as RM
    import copy
    import sys
    sys.path.append('./Molecular Absorption/code')
    import get_karkoschka_data as GKD
    import get_albedo_continua_crossections as gACC
    import NH3_Filter_Library_P3 as NFL
    import get_filter_base_dict as GFBD
    import planetmapper as pm
    
    ###########################################################################
    # Read MUSE data cube and filter metadata
    ###########################################################################
    MUSEhdr,MUSEhdr1,MUSEdata,MUSEzen,MUSEszen,wavelength,filterdata,path=RM.ReadMUSE(date)

    ###########################################################################
    # Compute spectrum and smooth
    ###########################################################################
    datetime=MUSEhdr["DATE-OBS"][0:10]+"-"+MUSEhdr["DATE-OBS"][11:13]+\
       MUSEhdr["DATE-OBS"][14:16]+"_"+str(int(10*(int(MUSEhdr["DATE-OBS"][17:19])/60.)))
    threshold=np.nanmax(MUSEdata[:,:,:])*0.001
    print("threshold=",threshold)
    MUSEdata[MUSEdata<threshold]=np.nan
    MUSESpec=np.nanmean(MUSEdata[:,:,:],axis=(1,2))


    mask=copy.deepcopy(MUSEdata)
    mask[mask>=threshold]=1.0
    count=np.nansum(mask[:,:,:],axis=(1,2)) 
    MUSESum=np.nansum(MUSEdata[:,:,:],axis=(1,2))    
    
    if date=="20220919UT":
        kernel_size = 14
        Gkernel_size=1
    elif date=="20220730UT":
        kernel_size = 1
        Gkernel_size=5
    kernel = np.ones(kernel_size) / kernel_size
    print(kernel.shape)
    data_convolved0 = np.convolve(MUSESpec, kernel, mode='same')
    data_convolved = np.convolve(data_convolved0, kernel, mode='same')
    
    ###########################################################################
    # Regrid smoothed spectra to 0.5 nm
    ###########################################################################
    WaveGrid,SignalonGrid=GSU.uniform_wave_grid(wavelength,data_convolved,Extend=False,Fine=False)
    MuseSpecGrid=np.column_stack((WaveGrid,SignalonGrid))
    print("MuseSpecGrid.shape=",MuseSpecGrid.shape)
    
    ###########################################################################
    # Load Solar Spectrum Reference and compute Albedo Spectrum
    ###########################################################################
    RefPath="c:/Astronomy/Python Play/SPLibraries_P3/SpectralReferenceFiles/ReferenceLibrary/"
    #ss="Pickles"
    ss="Kurucz"
    if ss=="Pickles":
        G2V=np.loadtxt(RefPath+"g2v.dat", dtype=float, usecols=(0,1))
        G2Vgrid=G2V
    elif ss=="Kurucz":
        G2V=np.loadtxt(RefPath+"combined_chance_kurucz_no_HDR.txt", dtype=float, usecols=(0,1))
    
        print("G2V.shape=",G2V.shape)
        print("G2V[0,0:5]=",G2V[0,0:5])
        print("G2V[0:5,0]=",G2V[0:5,0])
        
        print()
        print("G2V[1,0:5]=",G2V[1,1:5])
        
        G2V[:,0]=G2V[:,0]*1000.
        G2V[:,1]=G2V[:,1]/(10000*(4.*np.pi*(5.2*1.496e11)**2)) #convert from M**-2 to cm**-2
        G2V[:,1]=G2V[:,1]/(2.*np.pi)
        WaveGrid2,SignalonGrid2=GSU.uniform_wave_grid(G2V[:,0],G2V[:,1],Extend=False,Fine=False)

        Gkernel = np.ones(Gkernel_size) / Gkernel_size
        print(Gkernel.shape)
        temp0 = np.convolve(SignalonGrid2, Gkernel, mode='same')
        temp1 = np.convolve(temp0, Gkernel, mode='same')

        G2Vgrid=np.column_stack((WaveGrid2,temp1))
    #fig,axs=pl.subplots(1,1,figsize=(6.0,4.0), dpi=150, facecolor="white",
    #                    sharey=True,sharex=True)      
    #axs.plot(G2V[:,0],G2V[:,1])

    AlbedoSpec=GSU.SpectrumMath(MuseSpecGrid,G2Vgrid,"Divide")
    
    ###########################################################################
    # Write output spectral data
    ###########################################################################
    pathout="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/MUSE/data output/"
    np.savetxt(pathout+datetime+'_MUSE_Albedo.txt',AlbedoSpec,delimiter=" ",
               fmt="%10.3F %10.7F")
    
    SmoothedSpec=np.column_stack((WaveGrid,SignalonGrid))
    np.savetxt(pathout+datetime+'_MUSE_Spec.txt',SmoothedSpec,delimiter=" ",
               fmt="%10.3F %10.7F") 
    
    ###########################################################################
    # Compute mean signal
    ###########################################################################    
    F620_idx=[np.argmin(abs(filterdata['620']['wvs'][0]-WaveGrid)),np.argmin(abs(filterdata['620']['wvs'][1]-WaveGrid))]
    MUSE620Albedo=np.mean(AlbedoSpec[F620_idx[0]:F620_idx[1],1],axis=0)
    MUSE620Radiance=np.mean(MuseSpecGrid[F620_idx[0]:F620_idx[1],1],axis=0)
    print("MUSE620Albedo=",MUSE620Albedo)
    print("MUSE620Radiance=",MUSE620Radiance)
    
    F632_idx=[np.argmin(abs(filterdata['632']['wvs'][0]-WaveGrid)),np.argmin(abs(filterdata['632']['wvs'][1]-WaveGrid))]
    MUSE632Albedo=np.mean(AlbedoSpec[F632_idx[0]:F632_idx[1],1],axis=0)
    MUSE632Radiance=np.mean(MuseSpecGrid[F632_idx[0]:F632_idx[1],1],axis=0)
    print("MUSE632Albedo=",MUSE632Albedo)
    print("MUSE632Radiance=",MUSE632Radiance)
    
    F647_idx=[np.argmin(abs(filterdata['647']['wvs'][0]-WaveGrid)),np.argmin(abs(filterdata['647']['wvs'][1]-WaveGrid))]
    MUSE647Albedo=np.mean(AlbedoSpec[F647_idx[0]:F647_idx[1],1],axis=0)
    MUSE647Radiance=np.mean(MuseSpecGrid[F647_idx[0]:F647_idx[1],1],axis=0)
    print("MUSE647Albedo=",MUSE647Albedo)
    print("MUSE647Radiance=",MUSE647Radiance)
    
    F656_idxRED=[np.argmin(abs(filterdata['656']['wvsRED'][0]-WaveGrid)),np.argmin(abs(filterdata['656']['wvsRED'][1]-WaveGrid))]
    MUSE656AlbedoRED=np.mean(AlbedoSpec[F656_idxRED[0]:F656_idxRED[1],1],axis=0)

    F656_idxBLU=[np.argmin(abs(filterdata['656']['wvsBLU'][0]-WaveGrid)),np.argmin(abs(filterdata['656']['wvsBLU'][1]-WaveGrid))]
    MUSE656AlbedoBLU=np.mean(AlbedoSpec[F656_idxBLU[0]:F656_idxBLU[1],1],axis=0)
    MUSE656Albedo=(MUSE656AlbedoRED+MUSE656AlbedoBLU)/2.
    #MUSE656Radiance=(MUSE656AlbedoRED+MUSE656AlbedoBLU)/2.
    #print("MUSE656Albedo=",MUSE656Albedo)
    
    MUSERelSlope=(MUSE656Albedo-MUSE632Albedo)/(656.-632.)
    MUSE647cont=MUSE632Albedo+MUSERelSlope*15.
    #MUSE620cont=MUSE632Albedo+MUSERelSlope*(-12.)
    MUSE620cont=MUSE632Albedo
    MUSENH3Abs=MUSE647Albedo/MUSE647cont
    MUSECH4Abs=MUSE620Albedo/MUSE620cont
    
    print("")
    print("NH3 Abs=",MUSENH3Abs)
    print("CH4 Abs=",MUSECH4Abs)

    ###########################################################################
    # Begin Filter Convolution caluculations of transmission
    ###########################################################################
    path,filterdict=GFBD.get_filter_base_dict("SCT")
    kark=GKD.get_karkoschka_data()

    FilterTransmission620=filterdict['620']['FiltTrans']
    MUSE620_ProductV=GSU.SpectrumMath(FilterTransmission620,AlbedoSpec,"Multiply")
    MUSE620_ProductK=GSU.SpectrumMath(FilterTransmission620,kark,"Multiply")

    AvgTrans620V=sum(MUSE620_ProductV[F620_idx[0]:F620_idx[1],1])/ \
                    sum(FilterTransmission620[F620_idx[0]:F620_idx[1],1])
    AvgTrans620K=sum(MUSE620_ProductK[F620_idx[0]:F620_idx[1],1])/ \
                    sum(FilterTransmission620[F620_idx[0]:F620_idx[1],1])
    print("AvgTrans620=",AvgTrans620V,AvgTrans620K)                  
    
    FilterTransmission632=filterdict['632']['FiltTrans']
    MUSE632_ProductV=GSU.SpectrumMath(FilterTransmission632,AlbedoSpec,"Multiply")
    MUSE632_ProductK=GSU.SpectrumMath(FilterTransmission632,kark,"Multiply")
    AvgTrans632V=sum(MUSE632_ProductV[F632_idx[0]:F632_idx[1],1])/ \
                    sum(FilterTransmission632[F632_idx[0]:F632_idx[1],1])
    AvgTrans632K=sum(MUSE632_ProductK[F632_idx[0]:F632_idx[1],1])/ \
                    sum(FilterTransmission632[F632_idx[0]:F632_idx[1],1])
    print("AvgTrans632=",AvgTrans632V,AvgTrans632K)  
    print()

    CH4TransV=AvgTrans620V/AvgTrans632V
    CH4TransK=AvgTrans620K/AvgTrans632K
    print("CH4 Trans=",CH4TransV,CH4TransK)           
    print()
    FilterTransmission647=filterdict['647']['FiltTrans']
    MUSE647_ProductV=GSU.SpectrumMath(FilterTransmission647,AlbedoSpec,"Multiply")
    MUSE647_ProductK=GSU.SpectrumMath(FilterTransmission647,kark,"Multiply")
    AvgTrans647V=sum(MUSE647_ProductV[F647_idx[0]:F647_idx[1],1])/ \
                    sum(FilterTransmission647[F647_idx[0]:F647_idx[1],1])
    AvgTrans647K=sum(MUSE647_ProductK[F647_idx[0]:F647_idx[1],1])/ \
                    sum(FilterTransmission647[F647_idx[0]:F647_idx[1],1])
    print("AvgTrans647=",AvgTrans647V,AvgTrans647K)  
                
    FilterTransmission656=filterdict['656']['FiltTrans']
    MUSE656_ProductV=GSU.SpectrumMath(FilterTransmission656,AlbedoSpec,"Multiply")
    MUSE656_ProductK=GSU.SpectrumMath(FilterTransmission656,kark,"Multiply")
    AvgTrans656V=sum(MUSE656_ProductV[F656_idxBLU[0]:F656_idxRED[1],1])/ \
                    sum(FilterTransmission656[F656_idxBLU[0]:F656_idxRED[1],1])
    AvgTrans656K=sum(MUSE656_ProductK[F656_idxBLU[0]:F656_idxRED[1],1])/ \
                    sum(FilterTransmission656[F656_idxBLU[0]:F656_idxRED[1],1])
    print("AvgTrans656=",AvgTrans656V,AvgTrans656K) 
    print()

    NH3RelSlopeV=(AvgTrans656V-AvgTrans632V)/(656.-632.)
    NH3contV=AvgTrans632V+NH3RelSlopeV*15.
    NH3RelSlopeK=(AvgTrans656K-AvgTrans632K)/(656.-632.)
    NH3contK=AvgTrans632K+NH3RelSlopeK*15.

    NH3TransV=AvgTrans647V/NH3contV
    NH3TransK=AvgTrans647K/NH3contK
    print("NH3 Trans=",NH3TransV,NH3TransK)    
    print()

    ###########################################################################
    # Again, but with Karkoschka instead of VLT
    ###########################################################################     
    ###########################################################################
    # Do Plots
    ###########################################################################     
    scale=0.63
    fig2,axs2=pl.subplots(1,1,figsize=(6.0,6.0), dpi=150, facecolor="white",
                        sharey=True,sharex=True)      
    #axs2.plot(G2V[:,0],G2V[:,1],label=ss+" Solar Flux (W/cm**2/um")
    #axs2.plot(wavelength,MUSESpec,label='Raw MUSE Spectrum')
    #axs2.plot(WaveGrid,SignalonGrid,label='Smoothed, Regridded to 0.5 nm)')
    axs2.set_xlim(600.,680.)
    axs2.set_ylim(0,.0015)
    axs2.set_xlabel("Wavelength (nm)")
    axs2.set_ylabel("Radiance (W/m^2/sr/um)")
    axs2.set_title("VLT-MUSE Spectrum "+datetime)
    axs2.legend(fontsize=8,loc="lower right")
    axs2.annotate("NH3 Abs="+str(NH3TransV)[0:6],[0.01,0.95],
                  xycoords="axes fraction",fontsize=8)
    axs2.annotate("CH4 Abs="+str(CH4TransV)[0:6],[0.01,0.90],
                  xycoords="axes fraction",fontsize=8)
    
    axs3=axs2.twinx()
    axs3.plot(AlbedoSpec[:,0],AlbedoSpec[:,1]*scale,color='C3',
              label=ss+" Albedo")
    axs3.plot(kark[:,0],kark[:,1],color='C4',
              label="Karkoschka")
    axs3.plot(FilterTransmission620[:,0],FilterTransmission620[:,1],label='620')
    axs3.plot(FilterTransmission632[:,0],FilterTransmission632[:,1],label='632')
    axs3.plot(FilterTransmission647[:,0],FilterTransmission647[:,1],label='647')
    axs3.plot(FilterTransmission656[:,0],FilterTransmission656[:,1],label='656')
    axs3.set_ylim(0,1.)
    axs3.set_ylabel("Albedo")
    axs3.legend(fontsize=8,loc="upper right")

    fig2.savefig(pathout+datetime+'_MUSE_Spec.png')
    """
    MUSEzen[MUSEzen==-999]=np.nan
    lambert=np.cos(MUSEzen*np.pi/180)
    mean=np.nanmean(lambert)
    print(mean)
    MUSEszen[MUSEszen==-999]=np.nan
    lamberts=np.cos(MUSEzen*np.pi/180)
    means=np.nanmean(lamberts)
    print(means)
    """
    Jup=pm.Body("Jupiter",MUSEhdr["DATE-OBS"])
    print("Jupiter Dist=",Jup.target_distance)
    JupSun=pm.Body("Jupiter",MUSEhdr["DATE-OBS"],"Sun")
    print("Jupiter Dist=",JupSun.target_distance)

    fig4,axs4=pl.subplots(1,1,figsize=(6.0,3.5), dpi=150, facecolor="white",
                        sharey=True,sharex=True)
    axs4.set_xlim(600.,680.)
    axs4.set_ylim(0,1)
    axs4.set_xlabel("Wavelength (nm)")
    axs4.plot(FilterTransmission620[:,0],FilterTransmission620[:,1]/np.nanmax(FilterTransmission620[:,1]),label='620')
    axs4.plot(FilterTransmission632[:,0],FilterTransmission632[:,1]/np.nanmax(FilterTransmission632[:,1]),label='632')
    axs4.plot(FilterTransmission647[:,0],FilterTransmission647[:,1]/np.nanmax(FilterTransmission647[:,1]),label='647')
    axs4.plot(FilterTransmission656[:,0],FilterTransmission656[:,1]/np.nanmax(FilterTransmission656[:,1]),label='656')
    axs4.set_ylim(0,1.)
    axs4.set_ylabel("Normalized Transmission")
    axs4.legend(fontsize=8,loc="upper right")
