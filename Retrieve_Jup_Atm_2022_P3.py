# -*- coding: utf-8 -*-
"""
Created on Sun Nov  6 16:47:21 2022

@author: smhil
"""

def Retrieve_Jup_Atm_2022_P3(obsdate="20221009UTa",target="Jupiter",
                             imagetype='Map',CalModel='SCT-Obs-Final',Smoothing=True,
                             LatLims=[45,135],LonRng=45,delta_CM2=0,showbands=False):
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
    import RetrievalLibrary as RL

    # Retrieve_Jup_Atm_2022_P3(obsdate="20221019UT",target="Jupiter")   
    
    amagat=2.69e24 #Lodschmits number. (cm-2)
    gravity=2228.0 #cm/s^2
    mean_mol_wt=3.85e-24 #cgs or SI !!!!WHY IS THIS 2.2 FOR MENDIKOA!!!!!
    fCH4=1.81e-3
    STP=1.01e6  #dyne/cm^2 [(g-cm/s^2)/cm^2]`
    
    ###########################################################################
    #  DATA FILES AND METADATA DICTIONARY
    #    !!!!SHOULD MAKE THIS A DATA OBJECT!!!!
    #    !!!!DOUBLE CHECK THAT FITS FILE TIME TAGS ARE ACCURATE BETWEEN CH4 AND NH3
    #    !!!!Make sure that variable map boundaries can be input!!!!
    ###########################################################################
    sourcedata=obsdate+"_"+imagetype
    sourcefiles={'20220810UTa_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                                          'Seeing':'6/10','Transparency':'7/10'}, 
                               'CH4file':'2022-08-10-1013_0-Jupiter_620CH4AbsMap.fits',
                               'NH3file':'2022-08-10-1013_0-Jupiter_647NH3AbsMap.fits',
                               'RGBfile':'2022-08-10-1030_0-Jupiter_WV-R(AllRED)GB-RGB-WhtBal-Wavelets-Str_CM2_L360_MAP-BARE.png'},
                 '20220812UT_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                                         'Seeing':'6/10','Transparency':'7/10'}, 
                                            'CH4file':'2022-08-12-1025_5-Jupiter_620CH4AbsMap.fits',
                                            'NH3file':'2022-08-12-1025_5-Jupiter_647NH3AbsMap.fits',
                                            'RGBfile':'2022-08-12-1033_8-Jupiter_WV-R(AllRed)GB-RGB-WhtBal-Wavelets_CM2_L360_MAP-BARE.png'},
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
                 '20220830UT_Img':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                                         'Seeing':'6/10','Transparency':'7/10'}, 
                                            'CH4file':'2022-08-30-0559_1-Jupiter_CH4Abs620.fits',
                                            'NH3file':'2022-08-30-0559_1-Jupiter_NH3Abs647.fits'},
                 '20220830UT_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                                         'Seeing':'6/10','Transparency':'7/10'}, 
                                            'CH4file':'2022-08-30-0559_2-Jupiter_620CH4AbsMap.fits',
                                            'NH3file':'2022-08-30-0559_1-Jupiter_647NH3AbsMap.fits',
                                            'RGBfile':'2022-08-28-1429_7-Jupiter-RGB-Arakawa-j220828e1_CM2_L360_MAP-BARE.png'},
                 #!!! 7/16/2023: SOME SORT OF PROBLEM WITH CONTOUR LEVELS PREVENTS THIS FROM RUNNING; NEED TO RERUN
                 '20220901UT_Img':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                                         'Seeing':'6/10','Transparency':'7/10'}, 
                                            'CH4file':'2022-09-01-0604_9-Jupiter_CH4Abs620.fits',
                                            'NH3file':'2022-09-01-0604_9-Jupiter_NH3Abs647.fits',
                                            'RGBfile':'NA'},
                 '20220901UT_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                                         'Seeing':'6/10','Transparency':'7/10'}, 
                                            'CH4file':'2022-09-01-0604_9-Jupiter_620CH4AbsMap.fits',
                                            'NH3file':'2022-09-01-0604_9-Jupiter_647NH3AbsMap.fits',
                                            'RGBfile':'2022-09-01-0618_0-Jupiter-Rivera-j220901l1_CM2_L360_MAP-BARE.png'},
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
                 '20220905UTaSmth_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                                         'Seeing':'6/10','Transparency':'7/10'}, 
                                            'CH4file':'2022-09-05-0559_2-Jupiter_620CH4AbsMap.fits',
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
                 '20220913UT_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                                         'Seeing':'6/10','Transparency':'7/10'}, 
                                            'CH4file':'2022-09-13-0457_4-Jupiter_620CH4AbsMap.fits',
                                            'NH3file':'2022-09-13-0457_4-Jupiter_647NH3AbsMap.fits',
                                            'RGBfile':'2022-09-13-0455_6-Jupiter_WV-R685G550B450-RGB-WhtBal-ClrSmth-Wavelets_CM2_L360_MAP-BARE.png'},
                 '20220919UTa':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                                         'Seeing':'6/10','Transparency':'7/10'}, 
                                            'CH4file':'2022-09-19-0453_4-Jupiter_CH4Abs620.fits',
                                            'NH3file':'2022-09-19-0453_4-Jupiter_NH3Abs647.fits',
                                            'RGBfile':'2022-09-19-0518_7-Jupiter_WV-R685G550B450-RGB-WhtBal-Wavelets_CM2_L360_MAP-BARE.png'},
                 '20220919UTa_Map':{'Metadata':{'Telescope':'VLT','FL':'5600mm','Camera':'ASI120MM',
                                                         'Seeing':'6/10','Transparency':'7/10'}, 
                                            'CH4file':'2022-09-19-0352_3-Jupiter-620CH4AbsMap.fits',
                                            'NH3file':'2022-09-19-0352_3-Jupiter-647NH3AbsMap.fits',
                                            'RGBfile':'2022-09-19-0352_3-Jupiter-MUSE-R650G550B450-RGB-WhtBal-Filtered-Wavelets_CM2_L360_MAP-BARE.png'},
                 '20220919UTb_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                                         'Seeing':'6/10','Transparency':'7/10'}, 
                                            'CH4file':'2022-09-19-0453_4-Jupiter_620CH4AbsMap.fits',
                                            'NH3file':'2022-09-19-0453_4-Jupiter_647NH3AbsMap.fits',
                                            'RGBfile':'2022-09-19-0518_7-Jupiter_WV-R685G550B450-RGB-WhtBal-Wavelets_CM2_L360_MAP-BARE.png'},
                 '20220925UTa_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                                         'Seeing':'6/10','Transparency':'7/10'}, 
                                            'CH4file':'2022-09-25-0615_4-Jupiter_620CH4AbsMap.fits',
                                            'NH3file':'2022-09-25-0615_6-Jupiter_647NH3AbsMap.fits',
                                            'RGBfile':'2022-09-25-0546_6-Jupiter_WV-R685G550B450-RGB-WhtBal-Wavelets_CM2_L360_MAP-BARE.png'},
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
#                                            'CH4file':'2022-10-21-0358_6-Jupiter-620CH4AbsMap-Sharp.fits',
#                                            'NH3file':'2022-10-21-0358_6-Jupiter-647NH3AbsMap-Gauss2.fits',
                                            'CH4file':'2022-10-21-0358_6-Jupiter-620CH4AbsMap.fits',
                                            'NH3file':'2022-10-21-0358_6-Jupiter-647NH3AbsMap.fits',
                                            'RGBfile':'2022-10-21-0342_1-Jupiter-WV-R685G550B450-RGB-WhtBal-Wavelets_CM2_L360_MAP-BARE.png'},
                 '20230113UT_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                                         'Seeing':'8/10','Transparency':'7/10'}, 
                                            'CH4file':'2023-01-13-0046_2-Jupiter-620CH4AbsMap.fits',
                                            'NH3file':'2023-01-13-0046_2-Jupiter-647NH3AbsMap.fits',
                                            'RGBfile':'2023-01-13-0046_0-Jupiter-WV-R685G550B450-RGB-WhtBal-Wavelets_CM2_L360_MAP-BARE.png'}}

    calibration={'AGU 2022':{'CH4GlobalTrans':0.858,'NH3GlobalTrans':0.928},
                 'Model 1':{'CH4GlobalTrans':0.880,'NH3GlobalTrans':0.960},
                 'Model 2':{'CH4GlobalTrans':0.878,'NH3GlobalTrans':0.940},
                 'Observed':{'CH4GlobalTrans':0.910,'NH3GlobalTrans':0.972},
                 'VLT-MUSE':{'CH4GlobalTrans':0.861,'NH3GlobalTrans':0.961},
                 'SCT-Obs-Final':{'CH4GlobalTrans':0.920,'NH3GlobalTrans':0.972},
                 'VLT-Obs-Final':{'CH4GlobalTrans':0.893,'NH3GlobalTrans':0.962}}

    K_eff={'CH4_620':{'C11':0.427,'VLT':0.454},
           'NH3_647':{'C11':2.955,'VLT':3.129}}
    ###########################################################################
    # OPEN AND READ DATA FILES (FITS MAPS of NH3 and CH4 Absorption (Transmission?))
    ###########################################################################             
    K_eff_CH4620=K_eff['CH4_620'][sourcefiles[sourcedata]['Metadata']['Telescope']]#0.428
    K_eff_NH3647=K_eff['NH3_647'][sourcefiles[sourcedata]['Metadata']['Telescope']]#2.964

    path='c:/Astronomy/Projects/Planets/'+target+'/Imaging Data/'+obsdate[0:10]+'/'
    CH4file=sourcefiles[sourcedata]['CH4file']
    NH3file=sourcefiles[sourcedata]['NH3file']
    RGBfile=sourcefiles[sourcedata]['RGBfile']

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

    ###########################################################################
    # Get ephemeris data from file name date-time string and set lat and lon lims
    ###########################################################################             
    sec=str(int(str(NH3file[16:17]))*6)
    NH3time=(NH3file[0:10]+"_"+NH3file[11:13]+":"+NH3file[13:15]+":"+sec.zfill(2))
    eph=RL.get_WINJupos_ephem(NH3time)
    Real_CM2=float(eph[1].strip())
    CM2=Real_CM2+delta_CM2
    LonLims=[360-int(CM2+LonRng),360-int(CM2-LonRng)]

    ###LatLims=[45,135]
    ###LonLims=[360-int(CM2+45),360-int(CM2-45)]
    #LatLims=[30,150]
    #LonLims=[360-int(CM2+60),360-int(CM2-60)]
    #LatLims=[70,130]
    #LonLims=[360-int(25+30),360-int(25-30)]

    ###########################################################################
    # Set calibration and compute physical parameters
    ###########################################################################             
    CH4GlobalTrans=calibration[CalModel]['CH4GlobalTrans']
    NH3GlobalTrans=calibration[CalModel]['NH3GlobalTrans']
    
    kernel = Gaussian2DKernel(1)
    if Smoothing:
        CH4_tau=-np.log(convolve((CH4data*CH4GlobalTrans),kernel,boundary='extend'))
        NH3_tau=-np.log(convolve((NH3data*NH3GlobalTrans),kernel,boundary='extend'))
        smthtitle="Smoothed"
    else: 
        CH4_tau=-np.log(CH4data*CH4GlobalTrans)
        NH3_tau=-np.log(NH3data*NH3GlobalTrans)
        smthtitle="Unsmoothed"

    CH4_Ncol=1000*CH4_tau/K_eff_CH4620
    NH3_Ncol=1000*NH3_tau/K_eff_NH3647

    CH4_Cloud_Press=(CH4_Ncol/1000.)*amagat*gravity*mean_mol_wt/(fCH4*STP)
    fNH3=(NH3_Ncol/1000.)*amagat*gravity*mean_mol_wt/(CH4_Cloud_Press*STP)
    ##!!!! WOW!!! I need to calculate fNH3 the EASY way also and compare!!!
    #fNH3=NH3_Ncol/CH4_Ncol
    ###########################################################################
    # Set up figure and axes for plots
    ###########################################################################             

    fig,axs=pl.subplots(2,3,figsize=(7.0,4.5), dpi=150, facecolor="white",
                        sharey=True,sharex=True)      
    fig.suptitle(NH3time.replace("_"," ")+", CM2="+str(int(Real_CM2))+", Calibration = "+CalModel+", Data "+smthtitle,x=0.5,ha='center',color='k')

    if imagetype=='Map':
        for iPlot in range(0,6):
            i=int(iPlot/3)                      #Plot row
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

    ###########################################################################
    # Plot CH4 Optical Depth
    ###########################################################################
    if imagetype=='Map':
        plot_patch(CH4_tau,LatLims,LonLims,CM2,LonRng,"gist_heat",axs[0,0],'%3.3f')
    else:
        show=axs[0,0].imshow(CH4_tau, "gist_heat", origin='upper',vmin=0.10,vmax=0.20,
                           aspect="equal")
        temp=RL.make_contours_CH4(axs[0,0],np.flip(convolve(CH4_tau,kernel),axis=0),
                               lvls=np.linspace(0.15,0.20,num=5,endpoint=True))
    
    axs[0,0].set_title("CH4 620nm Optical Depth",fontsize=8)

    ###########################################################################
    # Plot CH4 Column Density (m-atm)
    ###########################################################################
    if imagetype=='Map':
        plot_patch(CH4_Ncol,LatLims,LonLims,CM2,LonRng,"gist_heat",axs[0,1],'%3.0f')
    else:
        show=axs[0,1].imshow(CH4_Ncol*1000., "gist_heat", origin='upper',vmin=200,vmax=500,
                           aspect="equal")
        temp=RL.make_contours_CH4(axs[0,1],np.flip(convolve(CH4_Ncol,kernel)*1000.,axis=0),
                               lvls=np.linspace(200.,500.,num=5,endpoint=True),frmt='%3.0f')
        
    axs[0,1].set_title("N[CH4] (m-atm)",fontsize=8)

    ###########################################################################
    # Plot CH4 Cloud Top Pressure (mb?) OR RGB Image
    ###########################################################################    
    if RGBfile == 'NA':
        show=axs[0,2].imshow(CH4_Cloud_Press, "gist_heat", origin='upper',vmin=2.0,vmax=6.0,
                           aspect="equal")
        temp=RL.make_contours_CH4(axs[0,2],CH4_Cloud_Press,
                               lvls=np.linspace(2.0,6.0,num=5,endpoint=True))
        axs[0,2].set_title("Scattering Pressure (Pa)",fontsize=8)
    else:
        if imagetype=='Map':
            RGB_patch=make_patch_RGB(RGB,LatLims,LonLims,CM2,LonRng)
            print("@@@@@@@ RGB_patch.shape()=",RGB_patch.shape)
            show=axs[0,2].imshow(RGB_patch, vmin=0,vmax=2^16,  
                       extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],
                               90-LatLims[0]],#vmin=0,vmax=1.2,
                               aspect="equal")
        else:    
            show=axs[0,2].imshow(RGB, aspect="equal",vmin=0,vmax=2^16)
            axs[0,2].set_title("Scattering Pressure (Pa)",fontsize=8)

    ###########################################################################
    # Plot NH3 Optical Depth
    ###########################################################################
    if imagetype=='Map':
        plot_patch(NH3_tau,LatLims,LonLims,CM2,LonRng,"gist_heat",axs[1,0],'%3.3f')
    else:
        show=axs[1,0].imshow(NH3_tau, "gist_heat", origin='upper',vmin=0.02,vmax=0.12,  
                           aspect="equal")
        temp=RL.make_contours_CH4(axs[1,0],np.flip(convolve(NH3_tau,kernel),axis=0),
                               lvls=np.linspace(0.02,0.12,num=5,endpoint=True))
        
    axs[1,0].set_title("NH3 647nm Optical Depth",fontsize=8)

    ###########################################################################
    # Plot NH3 Column Density (km-atm)
    ###########################################################################
    if imagetype=='Map':
        plot_patch(NH3_Ncol,LatLims,LonLims,CM2,LonRng,"gist_heat",axs[1,1],'%3.1f')
    else:    
        show=axs[1,1].imshow(NH3_Ncol*1000., "gist_heat", origin='upper',vmin=10,vmax=40,
                           aspect="equal")
        temp=RL.make_contours_CH4(axs[1,1],np.flip(convolve(NH3_Ncol,kernel)*1000.,axis=0),
                               lvls=np.linspace(20.,40.,num=5,endpoint=True),frmt='%3.0f')
    
    axs[1,1].set_title("N[NH3] (m-atm)",fontsize=8)
    
    ###########################################################################
    # Plot NH3 Mole Fraction - Ammonia Absorption Index
    ###########################################################################
    if imagetype=='Map':
        fNH3_patch_mb,vn_fNH3,vx_fNH3,tx_fNH3=plot_patch(fNH3*1e6,LatLims,LonLims,CM2,LonRng,"jet",axs[1,2],'%3.0f')

    else:    
        show=axs[1,2].imshow(fNH3*1e6, "jet", origin='upper',vmin=100,vmax=160,
                           aspect="equal")
        temp=RL.make_contours_CH4(axs[1,2],np.flip(convolve(fNH3*1.0e6,kernel),axis=0),
                               lvls=np.linspace(0.00010*1e6,0.00015*1e6,num=5,endpoint=True),frmt='%3.0f')
    
    axs[1,2].set_title("Ammonia Abundance Index (ppm)",fontsize=8)

    if RGBfile != 'NA':
        if imagetype=='Map':
            temp=RL.make_contours_CH4_patch(axs[0,2],fNH3_patch_mb,LatLims,LonLims,
                                   lvls=tx_fNH3,frmt='%3.0f',clr='k')
        else:
            temp=RL.make_contours_CH4(axs[0,2],np.flip(convolve(fNH3,kernel),axis=0),
                                   lvls=np.linspace(0.00010,0.00015,num=5,endpoint=True))

    axs[0,0].set_ylabel("Planetographic Latitude (deg)",fontsize=8)
    axs[1,0].set_ylabel("Planetographic Latitude (deg)",fontsize=8)
    axs[1,0].set_xlabel("Sys. 2 Longitude (deg)",fontsize=8)
    axs[1,1].set_xlabel("Sys. 2 Longitude (deg)",fontsize=8)
    axs[1,2].set_xlabel("Sys. 2 Longitude (deg)",fontsize=8)

    fig.subplots_adjust(left=0.10, bottom=0.08, right=0.94, top=0.90,
                wspace=0.25, hspace=0.05)     

    fig.savefig(path+"/"+obsdate+"-Jupiter-Retrieval"+"-CMII_"+
              str(Real_CM2)+"-"+CalModel+"-"+smthtitle+"-Map.png",dpi=300)

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
    #fig1.suptitle(NH3time.replace("_"," ")+", CM2="+str(int(CM2)),x=0.5,ha='center',color='k')
    fig1.suptitle(NH3time.replace("_"," ")+", CM2="+str(int(Real_CM2))+", Calibration = "+CalModel+", Data "+smthtitle,x=0.5,ha='center',color='k')

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


    show=axs1[0].imshow(fNH3_patch_mb, "jet", origin='upper',vmin=vn_fNH3,vmax=vx_fNH3,  
               extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],
                       90-LatLims[0]],#vmin=0,vmax=1.
                       aspect="equal")
    temp=RL.make_contours_CH4_patch(axs1[0],fNH3_patch_mb,LatLims,LonLims,
                           lvls=tx_fNH3,frmt='%3.0f',clr='k')
    
    cbar = pl.colorbar(show, ticks=tx_fNH3, 
                       orientation='vertical',cmap='gist_heat',
                       ax=axs1[0],fraction=0.046, pad=0.04)
    cbar.ax.set_yticklabels(np.around(tx_fNH3,1))
    cbar.ax.tick_params(labelsize=8,color='k')#if iSession >1:
        
    axs1[0].set_title(r'$\bar{f_c}(NH3) (ppm)$',fontsize=10)

    gamma=1.3
    RGB4Display=np.power(np.array(RGB_patch).astype(float),gamma)
    RGB4Display=RGB4Display/RGB4Display.max()
    show=axs1[1].imshow(RGB4Display,
               extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],
                       90-LatLims[0]],
                       aspect="equal")
    temp=RL.make_contours_CH4_patch(axs1[1],fNH3_patch_mb,LatLims,LonLims,
                           tx_fNH3,frmt='%3.0f',clr='k')
    box = axs1[1].get_position()
    
    belt={"SSTB":[-39.6,-36.2],
          "STB":[-32.4,-27.1],
          "SEB":[-19.7,-7.2],
          "NEB":[6.9,17.4],
          "NTB":[24.2,31.4],
          "NNTB":[35.4,39.6]}
    
    zone={"STZ":[-36.2,-32.4],
          "STrZ":[-27.1,-19.7],
          "EZ":[-7.2,6.9],
          "NTrZ":[17.4,24.2],
          "NTZ":[31.4,35.4]}

    if showbands:
        for zb in belt:
            print(zb,belt[zb])
            axs1[0].fill_between([360-LonLims[0],360-LonLims[1]],[belt[zb][0],belt[zb][0]],[belt[zb][1],belt[zb][1]],
                                    color="0.5",alpha=0.4)
            axs1[1].fill_between([360-LonLims[0],360-LonLims[1]],[belt[zb][0],belt[zb][0]],[belt[zb][1],belt[zb][1]],
                                    color="0.5",alpha=0.4)
        #axs1[1].annotate(zb,xy=[np.mean(belt[zb]),51],ha="center")
    #for zb in zone:
        #axs1[1].annotate(zb,xy=[np.mean(zone[zb]),51],ha="center")

    axs1[1].tick_params(axis='both', which='major', labelsize=9)
    axs1[1].set_title("RGB Context Image",fontsize=10)

    axs1[0].set_ylabel("Planetographic Latitude (deg)",fontsize=10)
    axs1[0].set_xlabel("Sys. 2 Longitude (deg)",fontsize=10)
    axs1[1].set_xlabel("Sys. 2 Longitude (deg)",fontsize=10)

    fig1.subplots_adjust(left=0.10, bottom=0.03, right=0.98, top=0.95,
                wspace=0.25, hspace=0.05)     
    axs1[1].set_position([box.x0+0.03, box.y0-0.01, box.width * 1.015, box.height * 1.015])

    fig1.savefig(path+"/"+obsdate+"-Jupiter-Retrieval-NH3_RGB_only"+"-CMII_"+
              str(Real_CM2)+"-"+CalModel+"-"+smthtitle+"-Map.png",dpi=300)
    
    ###########################################################################
    ## Just RGB and Cloud Pressure
    ###########################################################################
    #!!!!!!!!!! FIX THIS !!!!!!!!!!!!!
    fig2,axs2=pl.subplots(1,2,figsize=(8.0,4.0), dpi=150, facecolor="white",
                        sharey=True,sharex=True)
    #fig1.suptitle(NH3time.replace("_"," ")+", CM2="+str(int(CM2)),x=0.5,ha='center',color='k')
    fig2.suptitle(NH3time.replace("_"," ")+", CM2="+str(int(Real_CM2))+", Calibration = "+CalModel+", Data "+smthtitle,x=0.5,ha='center',color='k')

    for ix in range(0,1):
        axs2[ix].grid(linewidth=0.2)
        axs2[ix].ylim=[-45.,45.]
        axs2[ix].xlim=[360-LonLims[0],360-LonLims[1]]
        axs2[ix].set_xticks(np.linspace(450,0,31), minor=False)
        xticklabels=np.array(np.mod(np.linspace(450,0,31),360))
        axs2[ix].set_xticklabels(xticklabels.astype(int))
        axs2[ix].set_yticks(np.linspace(-45,45,7), minor=False)
        axs2[ix].tick_params(axis='both', which='major', labelsize=9)

        axs2[ix].set_adjustable('box') 

    Pcloud_patch,vn,vx,tx=plot_patch(CH4_Cloud_Press/4.0,LatLims,LonLims,CM2,LonRng,"jet",
                                     axs2[0],'%3.2f',cont=False,cbar_reverse=True)

    temp=RL.make_contours_CH4_patch(axs2[0],Pcloud_patch,LatLims,LonLims,
                           lvls=tx,frmt='%3.2f',clr='k')
    #temp=RL.make_contours_CH4_patch(axs2[0],fNH3_patch_mb,LatLims,LonLims,
    #                       lvls=tx_fNH3,frmt='%3.0f',clr='r')
    
    axs2[0].set_title("Cloud Top Pressure (bars)",fontsize=10)

    gamma=1.3
    RGB4Display=np.power(np.array(RGB_patch).astype(float),gamma)
    RGB4Display=RGB4Display/RGB4Display.max()
    show=axs2[1].imshow(RGB4Display,
               extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],
                       90-LatLims[0]],
                       aspect="equal")
    temp=RL.make_contours_CH4_patch(axs2[1],Pcloud_patch,LatLims,LonLims,
                           tx,frmt='%3.2f',clr='k')
    #temp=RL.make_contours_CH4_patch(axs2[1],fNH3_patch_mb,LatLims,LonLims,
    #                       lvls=tx_fNH3,frmt='%3.0f',clr='r')
    box = axs2[1].get_position()
    
    belt={"SSTB":[-39.6,-36.2],
          "STB":[-32.4,-27.1],
          "SEB":[-19.7,-7.2],
          "NEB":[6.9,17.4],
          "NTB":[24.2,31.4],
          "NNTB":[35.4,39.6]}
    
    zone={"STZ":[-36.2,-32.4],
          "STrZ":[-27.1,-19.7],
          "EZ":[-7.2,6.9],
          "NTrZ":[17.4,24.2],
          "NTZ":[31.4,35.4]}

    if showbands:
        for zb in belt:
            print(zb,belt[zb])
            axs2[0].fill_between([360-LonLims[0],360-LonLims[1]],[belt[zb][0],belt[zb][0]],[belt[zb][1],belt[zb][1]],
                                    color="0.5",alpha=0.4)
            axs2[1].fill_between([360-LonLims[0],360-LonLims[1]],[belt[zb][0],belt[zb][0]],[belt[zb][1],belt[zb][1]],
                                    color="0.5",alpha=0.4)
        #axs1[1].annotate(zb,xy=[np.mean(belt[zb]),51],ha="center")
    #for zb in zone:
        #axs1[1].annotate(zb,xy=[np.mean(zone[zb]),51],ha="center")

    axs2[1].tick_params(axis='both', which='major', labelsize=9)
    axs2[1].set_title("RGB Context Image",fontsize=10)

    axs2[0].set_ylabel("Planetographic Latitude (deg)",fontsize=10)
    axs2[0].set_xlabel("Sys. 2 Longitude (deg)",fontsize=10)
    axs2[1].set_xlabel("Sys. 2 Longitude (deg)",fontsize=10)

    fig2.subplots_adjust(left=0.10, bottom=0.03, right=0.98, top=0.95,
                wspace=0.25, hspace=0.05)     
    axs2[1].set_position([box.x0+0.03, box.y0-0.01, box.width * 1.015, box.height * 1.015])

    fig2.savefig(path+"/"+obsdate+"-Jupiter-Retrieval-Pcloud_only"+"-CMII_"+
              str(Real_CM2)+"-"+CalModel+"-"+smthtitle+"-Map.png",dpi=300)



    plot_scatter(Pcloud_patch,fNH3_patch_mb,path,obsdate,Real_CM2,LatLims)


    return()


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
    

def make_patch_RGB(Map,LatLims,LonLims,CM2deg,LonRng,pad=True):
    """
    Purpose: Make a map patch and handle the case where the data overlap
             the map edges. This is designed for a map with Jovian longitude
             conventions that with the left boundary at 360 ascending from
             the right boundary at 0. In WinJUPOS, the actual map setting
             shows the left boundary at zero, which is of course, also 360.
    """
    import numpy as np
    patch=np.copy(Map[LatLims[0]:LatLims[1],LonLims[0]:LonLims[1],:])
    print("####################### Patch RGB shape",patch.shape)
    if CM2deg<LonRng:
        patch=np.concatenate((np.copy(Map[LatLims[0]:LatLims[1],LonLims[0]-1:360,:]),
                              np.copy(Map[LatLims[0]:LatLims[1],0:LonLims[1]-360,:])),axis=1)
    if CM2deg>360.-LonRng:
        patch=np.concatenate((np.copy(Map[LatLims[0]:LatLims[1],360+LonLims[0]:360,:]),
                              np.copy(Map[LatLims[0]:LatLims[1],0:LonLims[1],:])),axis=1)
    print("####################### Patch RGB shape",patch.shape)
    #if pad:
    #    patch_pad=np.pad(patch,5,mode='reflect')

    return patch

def plot_patch(fullmap,LatLims,LonLims,CM2,LonRng,colorscale,axis,frmt,cont=True,cbar_reverse=False):
    import numpy as np
    import pylab as pl
    import RetrievalLibrary as RL

    patch=RL.make_patch(fullmap,LatLims,LonLims,CM2,LonRng)
    vn,vx,n=0.10,0.20,6
    #vn=np.min(patch)
    #vx=np.max(patch)
    
    vn=np.mean(patch)-3.0*np.std(patch)
    vx=np.mean(patch)+3.0*np.std(patch)
    
    print(np.mean(patch),vn,vx)
    tx=np.linspace(vn,vx,n,endpoint=True)
    show=axis.imshow(patch, colorscale, origin='upper',vmin=vn,vmax=vx,  
               extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],
                       90-LatLims[0]],
                       aspect="equal")
    if cont:
        temp=RL.make_contours_CH4_patch(axis,patch,LatLims,LonLims,
                           lvls=tx,frmt=frmt,clr='k')

    cbar = pl.colorbar(show, ticks=tx, 
               orientation='vertical',cmap='gist_heat',
               ax=axis,fraction=0.046, pad=0.04)
    cbar.ax.set_yticklabels(np.around(tx,3))
    cbar.ax.tick_params(labelsize=6,color="k")#if iSession >1:
    if cbar_reverse:
        cbar.ax.invert_yaxis()

    return patch,vn,vx,tx

def plot_scatter(patch1,patch2,filepath,obsdate,Real_CM2,LatLims):
    import pylab as pl
    import numpy as np
    import copy
    
    BZ={"SSTB":[-39.6,-36.2],
          "STZ":[-36.2,-32.4],
          "STB":[-32.4,-27.1],
          "STrZ":[-27.1,-19.7],
          "SEB":[-19.7,-7.2],
          "SEZ":[-7.2,0.0],
          "NEZ":[0.0,6.9],
          "NEB":[6.9,17.4],
          "NTrZ":[17.4,24.2],
          "NTB":[24.2,31.4],
          "NTZ":[31.4,35.4],
          "NNTB":[35.4,39.6]}

    BZind=copy.deepcopy(BZ)   
    BZkeys=BZ.keys()
    patch1=patch1*1000.

    figcor,axscor=pl.subplots(1,1,figsize=(6.0,4.5), dpi=150, facecolor="white",
                        sharey=True,sharex=True)          

    for key in BZ.keys():
        print(key,BZ[key],[90,90]-np.array(BZ[key]),LatLims)
        BZind[key][0]=int(90-BZ[key][0])-LatLims[0]
        BZind[key][1]=int(90-BZ[key][1])-LatLims[0]
        
        if BZind[key][0]<0:
            BZind[key][0]=0
        if BZind[key][1]<0:
            BZind[key][1]=0
        if BZind[key][0]>(LatLims[1]-LatLims[0]):
            BZind[key][0]=LatLims[1]-LatLims[0]
        if BZind[key][1]>(LatLims[1]-LatLims[0]):
            BZind[key][1]=LatLims[1]-LatLims[0]
        
        if BZind[key][0]==BZind[key][1]:
            print("do nothing")
        else:
            axscor.scatter(patch1[BZind[key][1]:BZind[key][0],:],
                           patch2[BZind[key][1]:BZind[key][0],:],
                           marker="o",s=5.0,
                           alpha=0.8,label=key)

    #print(patch1.shape,patch2.shape)
    #axscor.scatter(patch1,patch2,marker="o",s=0.5,edgecolor='C0',alpha=0.7,label="Other")
    """axscor.scatter(patch1[0:6,:],patch2[0:6,:],marker="o",s=5.0,color='C1',alpha=0.8,label="NTB")
    axscor.scatter(patch1[7:13,:],patch2[7:13,:],marker="o",s=5.0,color='C2',alpha=0.8,label="NTrZ")
    axscor.scatter(patch1[14:22,:],patch2[14:22,:],marker="o",s=5.0,color='C3',alpha=0.8,label="NEB")
    axscor.scatter(patch1[23:29,:],patch2[23:29,:],marker="o",s=5.0,color='C4',alpha=0.8,label="NEZ")
    axscor.scatter(patch1[30:36,:],patch2[30:36,:],marker="o",s=5.0,color='C5',alpha=0.8,label="SEZ")
    axscor.scatter(patch1[37:45,:],patch2[37:45,:],marker="o",s=5.0,color='C6',alpha=0.8,label="SEB")"""
    #axscor.set_title("Scatter Plot")
     
    axscor.grid(linewidth=0.2)
    axscor.set_xlim(500.,900.)
    axscor.set_ylim(50.,140)
    #axscor.set_xticks(np.linspace(-45.,45.,7), minor=False)
    #axscor.set_yticks(np.linspace(0.0,1.0,6), minor=False)
    #axscor.tick_params(axis='both', which='major', labelsize=8)
    axscor.set_xlabel("Effective Cloud-top Pressure (mb)",fontsize=10)
    axscor.set_ylabel("Ammonia Mole Fraction (ppm)",fontsize=10)
    axscor.legend(fontsize=9,loc=1)
    
    figcor.subplots_adjust(left=0.12, right=0.97, top=0.93, bottom=0.11)
    figcor.savefig(filepath+"/"+obsdate+"-Jupiter-Retrieval-Scatter"+"-CMII_"+str(Real_CM2)+".png",dpi=300)
  

