# -*- coding: utf-8 -*-
"""
Created on Sun Nov  6 16:47:21 2022

@author: smhil
"""

def Retrieve_Jup_Atm_2022_P3(obsdate="20221009UTa",target="Jupiter",imagetype='Map'):
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
    from re import search
    from astropy.io import fits
    from astropy.convolution import Gaussian2DKernel
    from astropy.convolution import convolve
    from astropy.io import fits

    # Retrieve_Jup_Atm_2022_P3(obsdate="20221019UT",target="Jupiter")   
    
    amagat=2.69e24 #Lodschmits number. (cm-2)
    gravity=2228.0 #cm/s^2
    mean_mol_wt=3.85e-24 #cgs or SI !!!!WHY IS THIS 2.2 FOR MENDIKOA!!!!!
    fCH4=1.81e-3
    STP=1.01e6  #dyne/cm^2 [(g-cm/s^2)/cm^2]
    
    K_eff_CH4620=0.428
    K_eff_NH3647=2.964

    ###########################################################################
    #  DATA FILES AND METADATA DICTIONARY
    #    !!!!SHOULD MAKE THIS A DATA OBJECT!!!!
    #    !!!!DOUBLE CHECK THAT FITS FILE TIME TAGS ARE ACCURATE BETWEEN CH4 AND NH3
    ###########################################################################
    sourcedata=obsdate+"_"+imagetype
    sourcefiles={'20220810UTa_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                                          'Seeing':'6/10','Transparency':'7/10'}, 
                               'CH4file':'2022-08-10-1013_0-Jupiter_620CH4AbsMap.fits',
                               'NH3file':'2022-08-10-1013_0-Jupiter_647NH3AbsMap.fits',
                               'RGBfile':'2022-08-10-1030_0-Jupiter_WV-R(AllRED)GB-RGB-WhtBal-Wavelets-Str_CM2_L360_MAP-BARE.png'},
                 '20220818UT':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'6/10','Transparency':'7/10'}, 
                               'CH4file':'2022-08-18-0733_4-Jupiter_CH4Abs620.fits',
                               'NH3file':'2022-08-18-0733_4-Jupiter_NH3Abs647.fits',
                               'RGBfile':'2022-08-18-0745_4-Jupiter_AllRED-WV-RGB-WhtBal-Wavelets_CM2_L360_MAP-BARE.png'},
                 '20220818UT_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                                         'Seeing':'6/10','Transparency':'7/10'}, 
                                            'CH4file':'2022-08-18-0733_4-Jupiter_620CH4AbsMap.fits',
                                            'NH3file':'2022-08-18-0733_4-Jupiter_647NH3AbsMap.fits',
                                            'RGBfile':'2022-08-18-0745_4-Jupiter_AllRED-WV-RGB-WhtBal-Wavelets_CM2_L360_MAP-BARE.png'},
                 '20220828UTa':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                                         'Seeing':'6/10','Transparency':'7/10'}, 
                                            'CH4file':'2022-08-28-0608_2-Jupiter_CH4Abs620.fits',
                                            'NH3file':'2022-08-28-0608_2-Jupiter_NH3Abs647.fits',
                                            'RGBfile':'NA'},
                 '20220828UTa_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                                         'Seeing':'6/10','Transparency':'7/10'}, 
                                            'CH4file':'2022-08-28-0608_2-Jupiter_620CH4AbsMap.fits',
                                            'NH3file':'2022-08-28-0608_2-Jupiter_647NH3AbsMap.fits',
                                            'RGBfile':'2022-08-28-0601_3-Jupiter-RGB-JamesWillinghan-j220828a1_CM2_L360_MAP-BARE.png'},
                 '20220830UTa':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                                         'Seeing':'6/10','Transparency':'7/10'}, 
                                            'CH4file':'2022-08-30-0559_1-Jupiter_CH4Abs620.fits',
                                            'NH3file':'2022-08-30-0559_1-Jupiter_NH3Abs647.fits'},
                 '20220901UTa':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                                         'Seeing':'6/10','Transparency':'7/10'}, 
                                            'CH4file':'2022-09-01-0604_9-Jupiter_CH4Abs620.fits',
                                            'NH3file':'2022-09-01-0604_9-Jupiter_NH3Abs647.fits',
                                            'RGBfile':'NA'},
                 '20220904UTa_Img':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                                         'Seeing':'6/10','Transparency':'7/10'}, 
                                            'CH4file':'2022-09-04-0638_9-Jupiter_620CH4AbsImg.fits',
                                            'NH3file':'2022-09-04-0638_2-Jupiter_647NH3AbsImg.fits',
                                            'RGBfile':'2022-09-04-1644_7-Jupiter-Yamane-j220904j4.jpg'},
                 '20220904UTa_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                                         'Seeing':'6/10','Transparency':'7/10'}, 
                                            'CH4file':'2022-09-04-0638_9-Jupiter_620CH4AbsMap.fits',
                                            'NH3file':'2022-09-04-0638_2-Jupiter_647NH3AbsMap.fits',
                                            'RGBfile':'2022-09-04-1644_7-Jupiter-Yamane-j220904j4_CM2_L360_MAP-BARE.png'},
                 '20220905UTa_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                                         'Seeing':'6/10','Transparency':'7/10'}, 
                                            'CH4file':'2022-09-05-0559_1-Jupiter_620CH4AbsMap.fits',
                                            'NH3file':'2022-09-05-0559_1-Jupiter_647NH3AbsMap.fits',
                                            'RGBfile':'2022-09-05-1559_0-Jupiter_j220905l1_Michael_Wong_CM2_L360_MAP-BARE.png'},
                 '20220912UTa':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                                         'Seeing':'6/10','Transparency':'7/10'}, 
                                            'CH4file':'2022-09-12-0533_4-Jupiter_CH4Abs620.fits',
                                            'NH3file':'2022-09-12-0533_4-Jupiter_NH3Abs647.fits',
                                            'RGBfile':'2022-09-12-0532_3-Jupiter_WV-R685G550B450-RGB-WhtBal-ClrSmth-Wavelets_CM2_L360_MAP-BARE.png'},
                 '20220912UTa_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                                         'Seeing':'6/10','Transparency':'7/10'}, 
                                            'CH4file':'2022-09-12-0533_4-Jupiter_620CH4AbsMap.fits',
                                            'NH3file':'2022-09-12-0533_4-Jupiter_647NH3AbsMap.fits',
                                            'RGBfile':'2022-09-12-0532_3-Jupiter_WV-R685G550B450-RGB-WhtBal-ClrSmth-Wavelets_CM2_L360_MAP-BARE.png'},
                 '20220919UTa':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                                         'Seeing':'6/10','Transparency':'7/10'}, 
                                            'CH4file':'2022-09-19-0453_4-Jupiter_CH4Abs620.fits',
                                            'NH3file':'2022-09-19-0453_4-Jupiter_NH3Abs647.fits',
                                            'RGBfile':'2022-09-19-0518_7-Jupiter_WV-R685G550B450-RGB-WhtBal-Wavelets_CM2_L360_MAP-BARE.png'},
                 '20220919UTa_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                                         'Seeing':'6/10','Transparency':'7/10'}, 
                                            'CH4file':'2022-09-19-0453_4-Jupiter_620CH4AbsMap.fits',
                                            'NH3file':'2022-09-19-0453_4-Jupiter_647NH3AbsMap.fits',
                                            'RGBfile':'2022-09-19-0518_7-Jupiter_WV-R685G550B450-RGB-WhtBal-Wavelets_CM2_L360_MAP-BARE.png'},
                 '20221009UTa_Img':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                                         'Seeing':'6/10','Transparency':'7/10'}, 
                                            'CH4file':'2022-10-09-0401_5-Jupiter_620CH4AbsImg.fits',
                                            'NH3file':'2022-10-09-0401_5-Jupiter_647NH3AbsImg.fits',
                                            'RGBfile':'2022-10-09-0339_0-Jupiter_WV-R685G550B450-RGB-WhtBal-Wavelets_CM2_L360_MAP-BARE.png'},
                 '20221009UTa_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                                         'Seeing':'6/10','Transparency':'7/10'}, 
                                            'CH4file':'2022-10-09-0401_5-Jupiter_620CH4AbsMap.fits',
                                            'NH3file':'2022-10-09-0401_5-Jupiter_647NH3AbsMap.fits',
                                            'RGBfile':'2022-10-09-0339_0-Jupiter_WV-R685G550B450-RGB-WhtBal-Wavelets_CM2_L360_MAP-BARE.png'},
                 '20221009UTb':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                                         'Seeing':'6/10','Transparency':'7/10'}, 
                                            'CH4file':'2022-10-09-0524_5-Jupiter_CH4Abs620.fits',
                                            'NH3file':'2022-10-09-0524_5-Jupiter_NH3Abs647.fits',
                                            'RGBfile':'2022-10-09-0542_8-Jupiter_NoWV-R685G550B450-RGB-WhtBal-Str0to160-Wavelets_CM2_L360_MAP-BARE.png'},
                 '20221009UTb_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                                         'Seeing':'6/10','Transparency':'7/10'}, 
                                            'CH4file':'2022-10-09-0524_5-Jupiter_620CH4AbsMap.fits',
                                            'NH3file':'2022-10-09-0524_5-Jupiter_647NH3AbsMap.fits',
                                            'RGBfile':'2022-10-09-0542_8-Jupiter_NoWV-R685G550B450-RGB-WhtBal-Str0to160-Wavelets_CM2_L360_MAP-BARE.png'},
                 '20221013UTa':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                                         'Seeing':'6/10','Transparency':'7/10'}, 
                                            'CH4file':'2022-10-13-0345_5-Jupiter-CH4Abs620.fits',
                                            'NH3file':'2022-10-13-0345_5-Jupiter-NH3Abs647.fits',
                                            'RGBfile':'2022-10-13-0402_0-Jupiter-WV-R685G550B450-RGB-WhtBal-Wavelets_CM2_L360_MAP-BARE.png'},
                 '20221013UT_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                                         'Seeing':'6/10','Transparency':'7/10'}, 
                                            'CH4file':'2022-10-13-0345_5-Jupiter-620CH4AbsMap.fits',
                                            'NH3file':'2022-10-13-0345_5-Jupiter-647NH3AbsMap.fits',
                                            'RGBfile':'2022-10-13-0402_0-Jupiter-WV-R685G550B450-RGB-WhtBal-Wavelets_CM2_L360_MAP-BARE.png'},
                 '20221019UT_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                                         'Seeing':'6/10','Transparency':'7/10'}, 
                                            'CH4file':'2022-10-19-0342_4-Jupiter-620CH4AbsMap.fits',
                                            'NH3file':'2022-10-19-0342_4-Jupiter-647NH3AbsMap.fits',
                                            'RGBfile':'2022-10-19-0358_2-Jupiter-WV-R685G550B450-RGB-WhtBal-Wavelets_CM2_L360_MAP-BARE.png'},
                 '20221020UT_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                                         'Seeing':'6/10','Transparency':'7/10'}, 
                                            'CH4file':'2022-10-20-0440_4-Jupiter-620CH4AbsMap.fits',
                                            'NH3file':'2022-10-20-0440_4-Jupiter-647NH3AbsMap.fits',
                                            'RGBfile':'2022-10-20-0422_6-Jupiter-WV-R685G550B450-RGB-WhtBal-Wavelets_CM2_L360_MAP-BARE.png'},
                 '20221021UT_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                                         'Seeing':'6/10','Transparency':'7/10'}, 
                                            'CH4file':'2022-10-21-0358_6-Jupiter-620CH4AbsMap.fits',
                                            'NH3file':'2022-10-21-0358_6-Jupiter-647NH3AbsMap.fits',
                                            'RGBfile':'2022-10-21-0342_1-Jupiter-WV-R685G550B450-RGB-WhtBal-Wavelets_CM2_L360_MAP-BARE.png'}}
    ###########################################################################
    # OBTAIN IMAGES TO DISPLAY AND DETERMINE IMAGE ARRAY SIZE
    ###########################################################################             
    path='c:/Astronomy/Projects/Planets/'+target+'/Imaging Data/'+obsdate[0:10]+'/'
    CH4file=sourcefiles[sourcedata]['CH4file']
    NH3file=sourcefiles[sourcedata]['NH3file']
    RGBfile=sourcefiles[sourcedata]['RGBfile']

    sec=str(int(str(NH3file[16:17]))*6)
    NH3time=(NH3file[0:10]+"_"+NH3file[11:13]+":"+NH3file[13:15]+":"+sec.zfill(2))
    eph=get_WINJupos_ephem(NH3time)
    CM2=float(eph[1].strip())
    print("********** CM2=",CM2)
    #LatLims=[45,135]
    #LonLims=[360-int(CM2+45),360-int(CM2-45)]
    LatLims=[30,150]
    LonLims=[360-int(CM2+60),360-int(CM2-60)]



    print("CH4file=",CH4file)
    CH4hdulist=fits.open(path+CH4file)
    CH4hdulist.info()
    CH4hdr=CH4hdulist[0].header
    CH4data=CH4hdulist[0].data
    CH4hdulist.close()
    
    NH3hdulist=fits.open(path+NH3file)
    NH3hdulist.info()
    NH3hdr=NH3hdulist[0].header
    NH3data=NH3hdulist[0].data
    NH3hdulist.close()

    if RGBfile != 'NA':
        RGB=imread(path+RGBfile)

    CH4GlobalTrans=0.858 #just a test number for now
    NH3GlobalTrans=0.928
    #CH4GlobalTrans=0.858 #just a test number for now
    #NH3GlobalTrans=0.961

    fig,axs=pl.subplots(2,3,figsize=(7.0,4.5), dpi=150, facecolor="white",
                        sharey=True,sharex=True)
    #fig.suptitle(obsdate+", CM2="+str(int(CM2)),x=0.5,ha='center',color='k')
    fig.suptitle(NH3time.replace("_"," ")+", CM2="+str(int(CM2)),x=0.5,ha='center',color='k')

    if imagetype=='Map':
        for iPlot in range(0,6):
                # Set up
            i=int(iPlot/3)                           #Plot row
            j=np.mod(iPlot,3)                   #Plot column
            axs[i,j].grid(linewidth=0.2)
            axs[i,j].ylim=[-45.,45.]
            axs[i,j].xlim=[360-LonLims[0],360-LonLims[1]]
            axs[i,j].set_xticks(np.linspace(450,0,31), minor=False)
            xticklabels=np.array(np.mod(np.linspace(450,0,31),360))
            axs[i,j].set_xticklabels(xticklabels.astype(int))
            axs[i,j].set_yticks(np.linspace(-45,45,7), minor=False)
            axs[i,j].tick_params(axis='both', which='major', labelsize=7)

            axs[i,j].set_adjustable('box') 

    kernel = Gaussian2DKernel(1)
    CH4_tau=-np.log(CH4data*CH4GlobalTrans)
    if imagetype=='Map':
        CH4_tau_patch=make_patch(CH4_tau,LatLims,LonLims,CM2)
        show=axs[0,0].imshow(CH4_tau_patch, "gist_heat", origin='upper',vmin=0.10,vmax=0.20,  
                   extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],
                           90-LatLims[0]],#vmin=0,vmax=1.2,
                           aspect="equal")
        temp=make_contours_CH4_patch(axs[0,0],convolve(CH4_tau_patch,kernel),LatLims,LonLims,
                               lvls=np.linspace(0.15,0.2,num=5,endpoint=True),frmt='%3.2f')

        cbar = pl.colorbar(show, ticks=[0.10,0.12,0.14,0.16,0.18,0.20], 
                   orientation='vertical',cmap='gist_heat',
                   ax=axs[0,0],fraction=0.046, pad=0.04)
        cbar.ax.tick_params(labelsize=6,color="k")#if iSession >1:

        
    else:
        show=axs[0,0].imshow(CH4_tau, "gist_heat", origin='upper',vmin=0.10,vmax=0.20,
                           aspect="equal")
        temp=make_contours_CH4(axs[0,0],np.flip(convolve(CH4_tau,kernel),axis=0),
                               lvls=np.linspace(0.15,0.20,num=5,endpoint=True))
    
    axs[0,0].set_title("CH4 620nm Optical Depth",fontsize=8)

    CH4_Ncol=CH4_tau/K_eff_CH4620
    if imagetype=='Map':
        CH4_Ncol_patch=make_patch(CH4_Ncol,LatLims,LonLims,CM2)
        show=axs[0,1].imshow(CH4_Ncol_patch*1000., "gist_heat", origin='upper',vmin=200,vmax=500,  
                   extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],
                           90-LatLims[0]],#vmin=0,vmax=1.2,
                           aspect="equal")
        temp=make_contours_CH4_patch(axs[0,1],convolve(CH4_Ncol_patch,kernel)*1000.,LatLims,LonLims,
                               lvls=np.linspace(200.,500.,num=5,endpoint=True),frmt='%3.0f')
        
        cbar = pl.colorbar(show, ticks=[200,300,400,500], 
                   orientation='vertical',cmap='gist_heat',
                   ax=axs[0,1],fraction=0.046, pad=0.04)
        cbar.ax.tick_params(labelsize=6,color="k")#if iSession >1:

    else:
        show=axs[0,1].imshow(CH4_Ncol*1000., "gist_heat", origin='upper',vmin=200,vmax=500,
                           aspect="equal")
        temp=make_contours_CH4(axs[0,1],np.flip(convolve(CH4_Ncol,kernel)*1000.,axis=0),
                               lvls=np.linspace(200.,500.,num=5,endpoint=True),frmt='%3.0f')
        
    axs[0,1].set_title("N[CH4] (m-atm)",fontsize=8)

    CH4_Cloud_Press=CH4_Ncol*amagat*gravity*mean_mol_wt/(fCH4*STP)
    
    if RGBfile == 'NA':
        show=axs[0,2].imshow(CH4_Cloud_Press, "gist_heat", origin='upper',vmin=2.0,vmax=6.0,
                           aspect="equal")
        temp=make_contours_CH4(axs[0,2],CH4_Cloud_Press,
                               lvls=np.linspace(2.0,6.0,num=5,endpoint=True))
        axs[0,2].set_title("Scattering Pressure (Pa)",fontsize=8)
    else:
        if imagetype=='Map':
            RGB_patch=make_patch(RGB,LatLims,LonLims,CM2)
            print("@@@@@@@ RGB_patch.shape()=",RGB_patch.shape)
            show=axs[0,2].imshow(RGB_patch, vmin=0,vmax=2^16,  
                       extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],
                               90-LatLims[0]],#vmin=0,vmax=1.2,
                               aspect="equal")
        else:    
            show=axs[0,2].imshow(RGB, aspect="equal",vmin=0,vmax=2^16)
            axs[0,2].set_title("Scattering Pressure (Pa)",fontsize=8)



    NH3_tau=-np.log(NH3data*NH3GlobalTrans)
    if imagetype=='Map':
        NH3_tau_patch=make_patch(NH3_tau,LatLims,LonLims,CM2)
        show=axs[1,0].imshow(NH3_tau_patch, "gist_heat", origin='upper',vmin=0.02,vmax=0.12,  
                   extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],
                           90-LatLims[0]],#vmin=0,vmax=1.2,
                           aspect="equal")
        temp=make_contours_CH4_patch(axs[1,0],convolve(NH3_tau_patch,kernel),LatLims,LonLims,
                               lvls=np.linspace(0.02,0.12,num=5,endpoint=True),frmt='%3.0f')
        
        cbar = pl.colorbar(show, ticks=[0.02,0.04,0.06,0.08,0.10,0.12], 
                   orientation='vertical',cmap='gist_heat',
                   ax=axs[1,0],fraction=0.046, pad=0.04)
        cbar.ax.tick_params(labelsize=6,color="k")#if iSession >1:

    else:
        show=axs[1,0].imshow(NH3_tau, "gist_heat", origin='upper',vmin=0.02,vmax=0.12,  
                           aspect="equal")
        temp=make_contours_CH4(axs[1,0],np.flip(convolve(NH3_tau,kernel),axis=0),
                               lvls=np.linspace(0.02,0.12,num=5,endpoint=True))
        
    axs[1,0].set_title("NH3 647nm Optical Depth",fontsize=8)

    NH3_Ncol=NH3_tau/K_eff_NH3647
    if imagetype=='Map':
        NH3_Ncol_patch=make_patch(NH3_Ncol,LatLims,LonLims,CM2)
        show=axs[1,1].imshow(NH3_Ncol_patch*1000., "gist_heat", origin='upper',vmin=10,vmax=40,  
                   extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],
                           90-LatLims[0]],#vmin=0,vmax=1.2,
                           aspect="equal")
        temp=make_contours_CH4_patch(axs[1,1],convolve(NH3_Ncol_patch,kernel)*1000.,LatLims,LonLims,
                              lvls=np.linspace(20.,40.,num=5,endpoint=True),frmt='%3.0f')
        
        cbar = pl.colorbar(show, ticks=[10,20,30,40], 
                           orientation='vertical',cmap='gist_heat',
                           ax=axs[1,1],fraction=0.046, pad=0.04)
        cbar.ax.tick_params(labelsize=6,color="k")#if iSession >1:
        
    else:    
        show=axs[1,1].imshow(NH3_Ncol*1000., "gist_heat", origin='upper',vmin=10,vmax=40,
                           aspect="equal")
        temp=make_contours_CH4(axs[1,1],np.flip(convolve(NH3_Ncol,kernel)*1000.,axis=0),
                               lvls=np.linspace(20.,40.,num=5,endpoint=True),frmt='%3.0f')
    
    axs[1,1].set_title("N[NH3] (m-atm)",fontsize=8)
    

    fNH3=NH3_Ncol*amagat*gravity*mean_mol_wt/(CH4_Cloud_Press*STP)
    if imagetype=='Map':
        fNH3_patch=make_patch(fNH3,LatLims,LonLims,CM2)
        show=axs[1,2].imshow(fNH3_patch*1e6, "jet", origin='upper',vmin=100,vmax=160,  
                   extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],
                           90-LatLims[0]],#vmin=0,vmax=1.
                           aspect="equal")
        temp=make_contours_CH4_patch(axs[1,2],convolve(fNH3_patch,kernel)*1.0e6,LatLims,LonLims,
                               lvls=np.linspace(0.00010*1.0e6,0.00015*1.0e6,num=5,endpoint=True),frmt='%3.0f')
        
        cbar = pl.colorbar(show, ticks=[100,120,140,160], 
                           orientation='vertical',cmap='gist_heat',
                           ax=axs[1,2],fraction=0.046, pad=0.04)
        cbar.ax.tick_params(labelsize=6,color='C7')#if iSession >1:


    else:    
        show=axs[1,2].imshow(fNH3*1e6, "jet", origin='upper',vmin=100,vmax=160,
                           aspect="equal")
        temp=make_contours_CH4(axs[1,2],np.flip(convolve(fNH3*1.0e6,kernel),axis=0),
                               lvls=np.linspace(0.00010*1e6,0.00015*1e6,num=5,endpoint=True),frmt='%3.0f')
    
    axs[1,2].set_title("Ammonia Abundance Index (ppm)",fontsize=8)


    if RGBfile != 'NA':
        if imagetype=='Map':
            
            temp=make_contours_CH4_patch(axs[0,2],convolve(fNH3_patch,kernel)*1.0e6,LatLims,LonLims,
                                   lvls=np.linspace(0.00010*1.0e6,0.00015*1.0e6,num=5,endpoint=True),frmt='%3.0f')
        else:
            temp=make_contours_CH4(axs[0,2],np.flip(convolve(fNH3,kernel),axis=0),
                                   lvls=np.linspace(0.00010,0.00015,num=5,endpoint=True))

    axs[0,0].set_ylabel("Planetographic Latitude (deg)",fontsize=8)
    axs[1,0].set_ylabel("Planetographic Latitude (deg)",fontsize=8)
    axs[1,0].set_xlabel("Sys. 2 Longitude (deg)",fontsize=8)
    axs[1,1].set_xlabel("Sys. 2 Longitude (deg)",fontsize=8)
    axs[1,2].set_xlabel("Sys. 2 Longitude (deg)",fontsize=8)

    fig.subplots_adjust(left=0.10, bottom=0.08, right=0.94, top=0.90,
                wspace=0.25, hspace=0.05)     

  
    fig.savefig(path+"/"+obsdate+"-Jupiter-Retrieval"+"_CMII_"+
               str(CM2)+"-Map.png",dpi=300)

    hdu = fits.PrimaryHDU(fNH3.astype(np.float32))
    hdul = fits.HDUList([hdu])
    hdul[0].header['BITPIX']=-32
    print(hdul[0].header[:])
    fnout=path+'/'+NH3file[0:26]+'fNH3Abs647.fits'
    try:
        os.remove(fnout)
    except: 
        print("file doesn't exist")
    hdul.writeto(fnout)
    hdul.close()
    #nh3abs16bit = np.nan_to_num(((5.*65535.*(normnh3 - 0.9))*mask[:,:,1]).astype(np.uint16))
    fNH3scaled=np.nan_to_num(((5000.*65535.*fNH3 - 0.9)))
    fNH3scaled[fNH3scaled<=0.]=0.0
    fNH3abs16bit = fNH3scaled.astype(np.uint16)
    imwrite(path+'/'+NH3file[0:26]+'fNH3Python.png', fNH3abs16bit)#.astype(np.uint16))

    ###########################################################################
    ## Just RGB and Abundance
    ###########################################################################
    
    fig1,axs1=pl.subplots(1,2,figsize=(8.0,4.0), dpi=150, facecolor="white",
                        sharey=True,sharex=True)
    fig1.suptitle(NH3time.replace("_"," ")+", CM2="+str(int(CM2)),x=0.5,ha='center',color='k')

    for ix in range(0,1):
        axs1[ix].grid(linewidth=0.2)
        axs1[ix].ylim=[-45.,45.]
        axs1[ix].xlim=[360-LonLims[0],360-LonLims[1]]
        axs1[ix].set_xticks(np.linspace(450,0,31), minor=False)
        xticklabels=np.array(np.mod(np.linspace(450,0,31),360))
        axs1[ix].set_xticklabels(xticklabels.astype(int))
        axs1[ix].set_yticks(np.linspace(-45,45,7), minor=False)
        axs1[ix].tick_params(axis='both', which='major', labelsize=9)

        axs1[ix].set_adjustable('box') 


    show=axs1[0].imshow(fNH3_patch*1e6, "jet", origin='upper',vmin=100,vmax=160,  
               extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],
                       90-LatLims[0]],#vmin=0,vmax=1.
                       aspect="equal")
    temp=make_contours_CH4_patch(axs1[0],convolve(fNH3_patch,kernel)*1.0e6,LatLims,LonLims,
                           lvls=np.linspace(0.00010*1.0e6,0.00015*1.0e6,num=5,endpoint=True),frmt='%3.0f')
    
    cbar = pl.colorbar(show, ticks=[100,120,140,160], 
                       orientation='vertical',cmap='gist_heat',
                       ax=axs1[0],fraction=0.046, pad=0.04)
    cbar.ax.tick_params(labelsize=8,color='k')#if iSession >1:
    axs1[0].set_title("Ammonia Abundance Index (ppm)",fontsize=10)

    gamma=1.3
    show=axs1[1].imshow(np.power(np.array(RGB_patch).astype(float),1.3),   
               extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],
                       90-LatLims[0]],#vmin=0,vmax=1.2,
                       aspect="equal")
    temp=make_contours_CH4_patch(axs1[1],convolve(fNH3_patch,kernel)*1.0e6,LatLims,LonLims,
                           lvls=np.linspace(0.00010*1.0e6,0.00015*1.0e6,num=5,endpoint=True),frmt='%3.0f')
    box = axs1[1].get_position()
    
    
    axs1[1].tick_params(axis='both', which='major', labelsize=9)
    axs1[1].set_title("RGB Context Image",fontsize=10)

    axs1[0].set_ylabel("Planetographic Latitude (deg)",fontsize=10)
    axs1[0].set_xlabel("Sys. 2 Longitude (deg)",fontsize=10)
    axs1[1].set_xlabel("Sys. 2 Longitude (deg)",fontsize=10)

    fig1.subplots_adjust(left=0.10, bottom=0.03, right=0.98, top=0.95,
                wspace=0.25, hspace=0.05)     
    axs1[1].set_position([box.x0+0.03, box.y0-0.01, box.width * 1.015, box.height * 1.015])

    fig.savefig(path+"/"+obsdate+"-Jupiter-Retrieval-NH3-RGB-only"+"_CMII_"+
              str(CM2)+"-Map.png",dpi=300)

    return()

def make_contours_CH4(ax,CH4Abs_conv,lvls=[0.71,0.73,0.75,0.77,0.79],frmt='%3.1e'):
    """
    PURPOSE: Overlay countours of NH3 absorption data on Jovian maps.
             Specifically designed for equivalent widths with mean values of
             ~0.55nm
    """
    cs=ax.contour(CH4Abs_conv,origin='upper', 
                  colors=['w','w','w','w','w'], alpha=0.5,levels=lvls,
                  linewidths=[0.5,0.5,1.0,0.5,0.5],
                  linestyles=['dashed','dashed','solid','dashed','dashed'])
    #ax.clabel(cs,[19.0,19.5,20.0,20.5,21.0],inline=True,fmt='%2.1f',fontsize=8)
    ax.clabel(cs,lvls,inline=True,fmt=frmt,fontsize=8)
    
def make_contours_CH4_patch(ax,CH4Abs_conv,LatLims,LonLims,lvls=[0.71,0.73,0.75,0.77,0.79],frmt='%3.1e'):
    """
    PURPOSE: Overlay countours of NH3 absorption data on Jovian maps.
             Specifically designed for equivalent widths with mean values of
             ~0.55nm
    """
    cs=ax.contour(CH4Abs_conv,origin='upper', 
                  extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],90-LatLims[0]],
                  colors=['w','w','w','w','w'], alpha=0.5,levels=lvls,
                  linewidths=[0.5,0.5,1.0,0.5,0.5],
                  linestyles=['dashed','dashed','solid','dashed','dashed'])
    #ax.clabel(cs,[19.0,19.5,20.0,20.5,21.0],inline=True,fmt='%2.1f',fontsize=8)
    ax.clabel(cs,lvls,inline=True,fmt=frmt,fontsize=8)

def load_png(file_path):
    """
    Purpose: Properly load a 48-bit PNG file
    Read from KITTI .png file
    Args:
        file_path string: file path(absolute)
    Returns:
        data (numpy.array): data of image in (Height, Width, 3) layout
    
    FROM: https://www.programcreek.com/python/example/98900/png.Reader
    """
    import png
    import numpy as np

    flow_object = png.Reader(filename=file_path)
    flow_direct = flow_object.asDirect()
    flow_data = list(flow_direct[2])
    (w, h) = flow_direct[3]['size']

    flow = np.zeros((h, w, 3), dtype=np.float64)
    for i in range(len(flow_data)):
        flow[i, :, 0] = flow_data[i][0::3]
        flow[i, :, 1] = flow_data[i][1::3]
        flow[i, :, 2] = flow_data[i][2::3]

    return flow.astype(np.uint16) 
    
def get_WINJupos_ephem(dateobs):
    # Example call: get_WINJupos_ephem('2021-09-05_04:09:00')
    import win32com.shell.shell as shell
    import time
    #shell.ShellExecuteEx(lpVerb='runas', lpFile='cmd.exe', lpParameters='/c '+commands) #run as admin
    ###########################################################################
    # WRITE *.BAT FILE WITH COMMAND SCRIPT FOR WINJUPOS
    ###########################################################################   
    batfile=open("WINJupos_CM.bat",'w')
    Line1='cd "\Program Files\WinJUPOS 12.1.1"\r\n'
    #Line2='WinJUPOS.x64.exe Jupiter /GetCM:2021-09-05_04:09:00 /GeoLong:-104.9 /GeoLat:39.7 /GetAlt >"c:\Astronomy\Projects\SAS 2021 Ammonia\Jupiter_NH3_Analysis_P3\cm.txt"\r\n'
    Line2='WinJUPOS.x64.exe Jupiter /GetCM:'+dateobs+' /GeoLong:-104.9 /GeoLat:39.7 /GetAlt >"c:\Astronomy\Projects\SAS 2021 Ammonia\Jupiter_NH3_Analysis_P3\cm.txt"\r\n'
    #batfile.writelines(["Line1\r\n","Line2\r\n"])
    batfile.writelines([Line1,Line2])
    batfile.close()
    ###########################################################################
    # EXECUTE *.BAT COMMAND FILE FOR WINJUPOS AND WAIT TO READ RESULT FILE
    ###########################################################################
    commands = "WINJupos_CM.bat"  
    shell.ShellExecuteEx(lpFile='cmd.exe', lpParameters='/c '+commands)
    time.sleep(1)
    ephemfile=open("C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/cm.txt",'r')
    Lines=ephemfile.readlines()
    ephemfile.close()
    LineString=str(Lines[0])
    print(LineString)
    ###########################################################################
    # PARSE OUTPUT FILE AND CREATE STRING ARRAY EPH FOR CM1, CM2, CM3, AND ALT
    ###########################################################################
    start=[i for i, letter in enumerate(LineString) if letter == "="]
    end=[i for i, letter in enumerate(LineString) if letter == "°"]
    eph=[]
    for i in range(0,3):
        temp=LineString[int(start[i])+1:int(end[i])]
        #print("CM"+str(i+1)+" = "+LineString[int(start[i])+1:int(end[i])])#,Linestring[start[1]:end[1]],Linestring[start[2]:end[2]])
        print("CM"+str(i+1)+" = "+temp)#,Linestring[start[1]:end[1]],Linestring[start[2]:end[2]])
        eph.extend([temp])        
    temp=LineString[int(start[3])+1:int(end[3])]
    print("Alt =  "+temp)
    eph.extend([temp])        
    return eph

def make_patch(Map,LatLims,LonLims,CM2deg):
    """
    Purpose: Make a map patch and handle the case where the data overlap
             the map edges. This is designed for a map with Jovian longitude
             conventions that with the left boundary at 360 ascending from
             the right boundary at 0. In WinJUPOS, the actual map setting
             shows the left boundary at zero, which is of course, also 360.
    """
    import numpy as np
    patch=np.copy(Map[LatLims[0]:LatLims[1],LonLims[0]:LonLims[1]])
    if CM2deg<45:
        patch=np.concatenate((np.copy(Map[LatLims[0]:LatLims[1],LonLims[0]-1:360]),
                              np.copy(Map[LatLims[0]:LatLims[1],0:LonLims[1]-360])),axis=1)
    if CM2deg>315:
        patch=np.concatenate((np.copy(Map[LatLims[0]:LatLims[1],360+LonLims[0]:360]),
                              np.copy(Map[LatLims[0]:LatLims[1],0:LonLims[1]])),axis=1)
    return patch
