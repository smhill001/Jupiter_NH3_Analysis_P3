# -*- coding: utf-8 -*-

def GRS_NH3_Comparison_2022(LatLims=[70,130],CM2=20,LonRng=30,target='GRSSCT'):
    """
    Created on Wed Nov 23 09:07:27 2022
    Example for paper profile plot with +/-20 deg long:
        GRS_NH3_Comparison_2022(LatLims=[45,135],CM2=25,LonRng=20,target='GRS')
    @author: smhil
    
    Call used for 2023 Paper:
        GRS_NH3_Comparison_2022(LatLims=[90,130],CM2=25,LonRng=20,target='GRSCM')
        
    """
    import sys
    drive='c:'
    sys.path.append(drive+'/Astronomy/Python Play')
    sys.path.append(drive+'/Astronomy/Python Play/Util_P3')
    sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy')
    sys.path.append('./Services')
    sys.path.append('./Profiles')

    import os
    from matplotlib.pyplot import imread
    import pylab as pl
    import numpy as np
    from numpy import inf
    from re import search
    from astropy.io import fits
    from astropy.convolution import Gaussian2DKernel
    from astropy.convolution import convolve
    import plot_TEXES_Groups_P3 as PTG
    import RetrievalLibrary as RL
    import get_obs_list as getlist

    import get_WINJupos_ephem

    ###########################################################################
    #  DATA FILES AND METADATA DICTIONARY
    #    !!!!SHOULD MAKE THIS A DATA OBJECT!!!!
    #    !!!!DOUBLE CHECK THAT FITS FILE TIME TAGS ARE ACCURATE BETWEEN CH4 AND NH3
    ###########################################################################
    #    sourcedata=obsdate+"_"+imagetype
    #!!!!Need to clarify GRSDates set vs GRS SCT+VLT or SCT only or VLT only
    sourcefiles=getlist.get_obs_list()

    GRSVLTDates=['20220730UT','20220919UTa']
    GRSSCTDates=['20220810UTa','20220818UTa','20220828UTa','20220904UTa',
                 '20220919UTb','20221013UTa','20221020UTa','20230113UTa']
    GRSCMDates=['20220730UT','20220818UT','20220828UT','20220904UT',
                '20220919UTa','20221013UT','20221020UT']
    
    GRSFilesSys2={'20220730UT':{'fNH3file':'2022-07-30-0729_8-Jupiter_fNH3_Sys2.fits',
                               'RGBfile':'2022-07-30-0729_8-Jupiter_R650G550B480-MUSE-RGB-Str_CM2_L360_MAP-BARE.png'},
                  '20220810UT':{'fNH3file':'2022-08-10-1013_0-Jupiter_fNH3Abs647.fits',
                                'RGBfile':'2022-08-10-1030_0-Jupiter_WV-R(AllRED)GB-RGB-WhtBal-Wavelets-Str_CM2_L360_MAP-BARE.png'},
                 '20220818UT':{'fNH3file':'2022-08-18-0733_4-Jupiter_fNH3Abs647.fits',
                               'RGBfile':'2022-08-18-0745_4-Jupiter_AllRED-WV-RGB-WhtBal-Wavelets_CM2_L360_MAP-BARE.png'},
                 '20220828UT':{'fNH3file':'2022-08-28-0608_2-Jupiter_fNH3Abs647.fits',
                               'RGBfile':'2022-08-28-0601_3-Jupiter-RGB-JamesWillinghan-j220828a1_CM2_L360_MAP-BARE.png'},
                 '20220904UT':{'fNH3file':'2022-09-04-0638_2-Jupiter_fNH3Abs647.fits',
                               'RGBfile':'2022-09-04-1644_7-Jupiter-Yamane-j220904j4_CM2_L360_MAP-BARE.png'},
                 '20220919UTa':{'fNH3file':'2022-09-19-0352_3-Jupiter-fNH3Abs647.fits',
                               'RGBfile':'2022-09-19-0518_7-Jupiter_WV-R685G550B450-RGB-WhtBal-Wavelets_CM2_L360_MAP-BARE.png'},
                 '20220919UTb':{'fNH3file':'2022-09-19-0453_4-Jupiter_fNH3Abs647.fits',
                               'RGBfile':'2022-09-19-0518_7-Jupiter_WV-R685G550B450-RGB-WhtBal-Wavelets_CM2_L360_MAP-BARE.png'},
                 '20221013UT':{'fNH3file':'2022-10-13-0345_5-Jupiter-fNH3Abs647.fits',
                               'RGBfile':'2022-10-13-0402_0-Jupiter-WV-R685G550B450-RGB-WhtBal-Wavelets_CM2_L360_MAP-BARE.png'},
                 '20221020UT':{'fNH3file':'2022-10-20-0440_4-Jupiter-fNH3Abs647.fits',
                               'RGBfile':'2022-10-20-0422_6-Jupiter-WV-R685G550B450-RGB-WhtBal-Wavelets_CM2_L360_MAP-BARE.png'},
                 '20230113UT':{'fNH3file':'2023-01-13-0046_2-Jupiter-fNH3Abs647.fits',
                               'RGBfile':'2023-01-13-0046_0-Jupiter-WV-R685G550B450-RGB-WhtBal-Wavelets_CM2_L360_MAP-BARE.png'}}

    OvalBADates=['20220812UT','20220901UT','20220925UT']
    OvalBAFilesSys2={'20220812UT':{'fNH3file':'2022-08-12-1025_5-Jupiter_fNH3Abs647.fits',
                               'RGBfile':'2022-08-12-1033_8-Jupiter_WV-R(AllRed)GB-RGB-WhtBal-Wavelets_CM2_L360_MAP-BARE.png'},
                 '20220901UT':{'fNH3file':'2022-09-01-0604_9-Jupiter_fNH3Abs647.fits',
                               'RGBfile':'2022-09-01-0618_0-Jupiter-Rivera-j220901l1_CM2_L360_MAP-BARE.png'},
                 '20220925UT':{'fNH3file':'2022-09-25-0615_6-Jupiter_fNH3Abs647.fits',
                               'RGBfile':'2022-09-25-0546_6-Jupiter_WV-R685G550B450-RGB-WhtBal-Wavelets_CM2_L360_MAP-BARE.png'}}

    WS6Dates=['20220901UT','20220913UT','20221009UTb','20221019UT']
    WS6FilesSys2={'20220901UT':{'fNH3file':'2022-09-01-0604_9-Jupiter_fNH3Abs647.fits',
                               'RGBfile':'2022-09-01-0618_0-Jupiter-Rivera-j220901l1_CM2_L360_MAP-BARE.png'},
                 '20220913UT':{'fNH3file':'2022-09-13-0457_4-Jupiter_fNH3Abs647.fits',
                               'RGBfile':'2022-09-13-0455_6-Jupiter_WV-R685G550B450-RGB-WhtBal-ClrSmth-Wavelets_CM2_L360_MAP-BARE.png'},
                 '20221009UTb':{'fNH3file':'2022-10-09-0524_5-Jupiter_fNH3Abs647.fits',
                               'RGBfile':'2022-10-09-0542_8-Jupiter_NoWV-R685G550B450-RGB-WhtBal-Str0to160-Wavelets_CM2_L360_MAP-BARE.png'},
                 '20221019UT':{'fNH3file':'2022-10-19-0342_4-Jupiter-fNH3Abs647.fits',
                               'RGBfile':'2022-10-19-0358_2-Jupiter-WV-R685G550B450-RGB-WhtBal-Wavelets_CM2_L360_MAP-BARE.png'}}
 
    SysII_147Dates=['20220905UT','20221009UTa','20221021UT']
    SysII_147FilesSys2={'20220905UT':{'fNH3file':'2022-09-05-0559_1-Jupiter_fNH3Abs647.fits',
                               'RGBfile':'2022-09-05-1559_0-Jupiter_j220905l1_Michael_Wong_CM2_L360_MAP-BARE.png'},
                 '20221009UTa':{'fNH3file':'2022-10-09-0401_5-Jupiter_fNH3Abs647.fits',
                               'RGBfile':'2022-10-09-0339_0-Jupiter_WV-R685G550B450-RGB-WhtBal-Wavelets_CM2_L360_MAP-BARE.png'},
                 '20221021UT':{'fNH3file':'2022-10-21-0358_6-Jupiter-fNH3Abs647.fits',
                               'RGBfile':'2022-10-21-0342_1-Jupiter-WV-R685G550B450-RGB-WhtBal-Wavelets_CM2_L360_MAP-BARE.png'}}

    DarkFeatDates=['20220810UT','20220905UT','20220912UT','20220919UT','20221019UT','20221021UT']
    
    if target=='GRSSCT':  #CM2~22
        Dates=GRSSCTDates
        sourcefiles=GRSFilesSys2
    if target=='GRSVLT':  #CM2~22
        Dates=GRSVLTDates
        sourcefiles=GRSFilesSys2
    if target=='GRSCM':  #CM2~22
        Dates=GRSCMDates
        sourcefiles=GRSFilesSys2
    elif target=="OvalBA":  #CM2~270
        Dates=OvalBADates
        sourcefiles=OvalBAFilesSys2
    elif target=="WS6":
        Dates=WS6Dates
        sourcefiles=WS6FilesSys2
    elif target=="SysII_147":  #CM2=147
        Dates=SysII_147Dates
        sourcefiles=SysII_147FilesSys2
    
    
    NDates=len(Dates)

    fig,axs=pl.subplots(2,NDates,figsize=(12.0,4.0), dpi=150, facecolor="white",
                        sharey=True,sharex=True)
    fig.suptitle('2022 GRS NH3',x=0.5,ha='center',color='k')
    LonLims=[360-int(CM2+LonRng),360-int(CM2-LonRng)]
    target='Jupiter'
    kernel = Gaussian2DKernel(1)
    #clrmap='gist_heat'
    clrmap='jet'
    
    MeridEWArray=np.zeros((LatLims[1]-LatLims[0],len(Dates)))
    MeridEWArrayStd=np.zeros((LatLims[1]-LatLims[0],len(Dates)))
    statsarray=np.zeros((NDates,4))

    first=True
    firstRGB=True
    for iPlot in range(0,2*NDates):
            # Set up
        i=int(iPlot/NDates)                           #Plot row
        j=np.mod(iPlot,NDates)                   #Plot column
        axs[i,j].grid(linewidth=0.2)
        axs[i,j].ylim=[-45.,45.]
        axs[i,j].xlim=[360-LonLims[0],360-LonLims[1]]
        print(360-LonLims[0],360-LonLims[1])
        print(LatLims)
        axs[i,j].set_xticks(np.linspace(450,0,31), minor=False)
        xticklabels=np.array(np.mod(np.linspace(450,0,31),360))
        axs[i,j].set_xticklabels(xticklabels.astype(int))
        axs[i,j].set_yticks(np.linspace(-45,45,7), minor=False)
        axs[i,j].tick_params(axis='both', which='major', labelsize=7)

        axs[i,j].set_adjustable('box') 
        #path='c:/Astronomy/Projects/Planets/'+target+'/Imaging Data/'+Dates[j][0:10]+'/'
        path='c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Analysis Data/L3 FITS/'
        print(Dates[j])

        try:    #Set up to allow for parametric studies of different processing paths
            NH3file=sourcefiles[target]['NH3file'][0:17]+"-Jupiter_Map_L2TNH3"+\
                    sourcefiles[target]['Variation']+".fits"
            variation=sourcefiles[target]['Variation']
        except:
            NH3file=sourcefiles[target]['NH3file'][0:17]+"-Jupiter_Map_L2TNH3.fits"
            variation=""
        
        NH3hdulist=fits.open(path+NH3file)
        NH3hdulist.info()
        NH3hdr=NH3hdulist[0].header
        NH3data=NH3hdulist[0].data
        sza=NH3hdulist[1].data
        eza=NH3hdulist[2].data
        NH3hdulist.close()
        #!!!WORK IN PROGRESS - MODERNIZING TO USE get_obs_list etc.
        """
        fNH3hdulist=fits.open(path+sourcefiles[Dates[j]]['fNH3file'])
        fNH3hdulist.info()
        fNH3hdr=fNH3hdulist[0].header
        fNH3data=fNH3hdulist[0].data
        fNH3hdulist.close()
        """
        fNH3_patch=RL.make_patch(fNH3data,LatLims,LonLims,CM2,LonRng)*1.0e6
        vn=np.mean(fNH3_patch)-3.0*np.std(fNH3_patch)
        vx=np.mean(fNH3_patch)+3.0*np.std(fNH3_patch)
        n=6
        tx=np.linspace(vn,vx,n,endpoint=True)
        print('$$$$$$$$$$ tx=',tx,vn,vx,n)
        
        stats=RL.patchstats(fNH3_patch)
        statsarray[j,:]=stats
        
        fNH3_patch_smooth=RL.make_patch(convolve(fNH3data,kernel,boundary='extend'),LatLims,LonLims,CM2,LonRng)*1.0e6
        if first:
            fNH3_patch_avg=np.zeros(fNH3_patch_smooth.shape)
        if i>>0:
            print()
            print(vn,vx)

            show=axs[1,j].imshow(fNH3_patch_smooth, clrmap, origin='upper',vmin=vn,vmax=vx,  
                       extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],
                               90-LatLims[0]],#vmin=0,vmax=1.2,
                               aspect="equal")
            temp=RL.make_contours_CH4_patch(axs[1,j],fNH3_patch_smooth,LatLims,LonLims,
                                   lvls=tx,frmt='%3.0f',clr='k')
            axs[1,j].set_xlabel("Sys. 2 Longitude (deg)",fontsize=9)
            
            sfile=sourcefiles[Dates[j]]['fNH3file']
            sec=str(int(str(sfile[16:17]))*6)
            sfiletime=(sfile[0:10]+"_"+sfile[11:13]+":"+sfile[13:15]+":"+sec.zfill(2))
            eph=get_WINJupos_ephem.get_WINJupos_ephem(sfiletime)
            ObsCM2=float(eph[1].strip())


            #MeridEW=np.flip(np.mean(fNH3_patch_smooth[:,:],axis=1),axis=0)
            #MeridEWerror=np.flip(np.std(fNH3_patch_smooth[:,:],axis=1),axis=0)
            print("############ fNH3_patch.shape",fNH3_patch.shape)
            MeridEW=np.flip(np.mean(fNH3_patch[:,:],axis=1),axis=0)
            MeridEWerror=np.flip(np.std(fNH3_patch[:,:],axis=1),axis=0)
            Lats=np.linspace(-44.5,44.5,90)
            #print DateCounter; Date
            MeridEWArray[:,j]=MeridEW[:]

            MeridEWArrayStd[:,j]=MeridEWerror[:]

            axs[1,j].set_title("CM2 = "+str(ObsCM2),fontsize=9)
            fNH3_patch_avg=fNH3_patch_avg+fNH3_patch_smooth/NDates
            first=False
            
        RGB=imread(path+sourcefiles[Dates[j]]['RGBfile'])
        RGB_patch=RL.make_patch(RGB,LatLims,LonLims,CM2,LonRng)
        if firstRGB:
            RGB_patch_avg=np.zeros(RGB_patch.shape)
        if i==0:
            show=axs[0,j].imshow(RGB_patch, vmin=0,vmax=2^16,  
                       extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],
                               90-LatLims[0]],
                               aspect="equal")
            temp=RL.make_contours_CH4_patch(axs[0,j],fNH3_patch_smooth,LatLims,LonLims,
                                   lvls=tx,frmt='%3.0f',clr='k')
            RGB_patch_avg=RGB_patch_avg+RGB_patch/NDates
            firstRGB=False
            axs[0,j].set_title(Dates[j],fontsize=9)

    axs[0,0].set_ylabel("Planetographic Latitude (deg)",fontsize=8)
    axs[1,0].set_ylabel("Planetographic Latitude (deg)",fontsize=8)


    fig.subplots_adjust(left=0.05, bottom=0.09, right=0.98, top=0.90,
                wspace=0.10, hspace=0.10)  
    
    pathout='c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/'

    fig.savefig(pathout+"GRS_NH3_Comparison_2022.png",dpi=300)
    for j in range(0,NDates):
        print(statsarray[j,:],(statsarray[j,3]-statsarray[j,2])/(2.*statsarray[j,1]))

    ###########################################################################
    # Make GRS 2022 average plot 
    ###########################################################################


    figavg,axsavg=pl.subplots(1,2,figsize=(8.0,4.0), dpi=150, facecolor="white",
                        sharey=True,sharex=True)
    #fig1.suptitle(NH3time.replace("_"," ")+", CM2="+str(int(CM2)),x=0.5,ha='center',color='k')
    figavg.suptitle('Average 2022 SCT GRS Observations',x=0.5,ha='center',color='k')

    for ix in range(0,1):
        axsavg[ix].grid(linewidth=0.2)
        axsavg[ix].ylim=[-45.,45.]
        axsavg[ix].xlim=[360-LonLims[0],360-LonLims[1]]
        axsavg[ix].set_xticks(np.linspace(450,0,31), minor=False)
        xticklabels=np.array(np.mod(np.linspace(450,0,31),360))
        axsavg[ix].set_xticklabels(xticklabels.astype(int))
        axsavg[ix].set_yticks(np.linspace(-45,45,7), minor=False)
        axsavg[ix].tick_params(axis='both', which='major', labelsize=10)

        axsavg[ix].set_adjustable('box') 

    vn=np.mean(fNH3_patch_avg)-3.0*np.std(fNH3_patch_avg)
    vx=np.mean(fNH3_patch_avg)+3.0*np.std(fNH3_patch_avg)
    n=6
    tx=np.linspace(vn,vx,n,endpoint=True)


    show=axsavg[0].imshow(fNH3_patch_avg, "jet", origin='upper',vmin=vn,vmax=vx,  
               extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],
                       90-LatLims[0]],#vmin=0,vmax=1.
                       aspect="equal")
    temp=RL.make_contours_CH4_patch(axsavg[0],fNH3_patch_avg,LatLims,LonLims,
                           lvls=tx,frmt='%3.0f',clr='k')
    
    cbar = pl.colorbar(show, ticks=tx, 
                       orientation='vertical',cmap='gist_heat',
                       ax=axsavg[0],fraction=0.046, pad=0.04)
    cbar.ax.set_yticklabels(np.around(tx,1))

    cbar.ax.tick_params(labelsize=8,color='k')#if iSession >1:
    axsavg[0].set_title(r'$\bar{f_c} (NH_3)$ (ppm)',fontsize=10)

    gamma=1.3
    RGB4Display_avg=np.power(np.array(RGB_patch_avg).astype(float),gamma)
    RGB4Display_avg=RGB4Display_avg/RGB4Display_avg.max()

    show=axsavg[1].imshow(RGB4Display_avg,   
               extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],
                       90-LatLims[0]],#vmin=0,vmax=1.2,
                       aspect="equal")
    temp=RL.make_contours_CH4_patch(axsavg[1],fNH3_patch_avg,LatLims,LonLims,
                           tx,frmt='%3.0f',clr='k')
    box = axsavg[1].get_position()
    
    
    #axsavg[1].tick_params(axis='both', which='major', labelsize=9)
    axsavg[1].set_title("RGB Context Image",fontsize=10)

    axsavg[0].set_ylabel("Planetographic Latitude (deg)",fontsize=10)
    axsavg[0].set_xlabel("Sys. 2 Longitude (deg)",fontsize=10)
    axsavg[1].set_xlabel("Sys. 2 Longitude (deg)",fontsize=10)

    figavg.subplots_adjust(left=0.10, bottom=0.03, right=0.98, top=0.95,
                wspace=0.25, hspace=0.05)     
    axsavg[1].set_position([box.x0+0.03, box.y0-0.01, box.width * 1.015, box.height * 1.015])

    #figavg.savefig(path+"/"+obsdate+"-Jupiter-Retrieval-NH3_RGB_only"+"-CMII_"+
    #          str(CM2)+"-"+CalModel+"-"+smthtitle+"-Map.png",dpi=300)
    figavg.savefig(pathout+"GRS_NH3_Comparison_AVG_2022.png",dpi=300)
    
    ###########################################################################
    # Make Fletcher Plot
    ###########################################################################
    
    figfletch,axsfletch=pl.subplots(1,2,figsize=(8.0,4.0), dpi=150, facecolor="white",
                        sharey=True,sharex=True)
    figfletch.suptitle('Gemini/TEXES observations in 2014',x=0.5,ha='center',color='k')

    print(fNH3_patch_avg.shape,RGB_patch_avg.shape)
    print(RGB_patch_avg.min(),RGB_patch_avg.max())
    for ix in range(0,1):
        axsfletch[ix].grid(linewidth=0.2)
        #axsavg[i,j].ylim=[-45.,45.]
        #axsavg[i,j].xlim=[360-LonLims[0],360-LonLims[1]]
        print(360-LonLims[0],360-LonLims[1])
        print(LatLims)
        axsfletch[ix].set_xticks(np.linspace(450,0,31), minor=False)
        xticklabels=np.array(np.mod(np.linspace(450,0,31),360))
        axsfletch[ix].set_xticklabels(xticklabels.astype(int))

        axsfletch[ix].set_yticks(np.linspace(-45,45,7), minor=False)
        axsfletch[ix].tick_params(axis='both', which='major', labelsize=10)

    FletcherRGBfile='C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Fletcher-2016-RGB-Screenshot 2022-11-26 232747.jpeg'
    FletcherRGB=imread(FletcherRGBfile)
    FletcherNH3file='C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Fletcher-2016-AmmoniaVMR-Screenshot 2022-11-26 232643.jpeg'
    FletcherNH3=imread(FletcherNH3file)
    show=axsfletch[1].imshow(FletcherRGB,  
               extent=[148,88,90-LatLims[1],
                       90-LatLims[0]],#vmin=0,vmax=1.2,
                       aspect="equal")
    axsfletch[1].set_title("RGB Context  \n[Fletcher et al., 2016] (adapted)",fontsize=9,wrap=True)

    show=axsfletch[0].imshow(fNH3_patch_avg, clrmap, origin='upper',vmin=120,vmax=150,  
               extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],
                       90-LatLims[0]],#vmin=0,vmax=1.2,
                       aspect="equal")
    cbar = pl.colorbar(show, ticks=[122,130.5,139.5,148], 
                       orientation='vertical',cmap='jet',
                       ax=axsfletch[0],fraction=0.046, pad=0.04)
    cbar.ax.tick_params(labelsize=8,color='k')#if iSession >1:
    cbar.ax.set_yticklabels(['5', '10', '15','20'],color="k")
    box = axsfletch[1].get_position()

    statsavg=RL.patchstats(fNH3_patch_avg)
    print()
    print(statsavg,(statsavg[3]-statsavg[2])/(2.*statsavg[1]))
    print()
    print(np.nanmean(statsarray[:,0]),np.nanstd(statsarray[:,0]))
    
    #####Fletcher########
    #!! 8/4/2023 - there's some problem with the Fletcher data scaling properly
    #              in the map extent for this figure.
    #show=axsfletch[0].imshow(FletcherNH3,  
    #           extent=[148,88,90-LatLims[1],
    #                   90-LatLims[0]],#vmin=0,vmax=1.2,
    #                   aspect="equal")
    show=axsfletch[0].imshow(FletcherNH3,  
               extent=[135,85,90-LatLims[1],
                       90-LatLims[0]],#vmin=0,vmax=1.2,
                       aspect="equal")

    axsfletch[0].set_title("Ammonia Abundance at 500mb (ppm) \n[Fletcher et al., 2016] (adapted)",fontsize=9,wrap=True)

    axsfletch[0].set_ylabel("Planetographic Latitude (deg)",fontsize=10)
    axsfletch[0].set_xlabel("Sys. 3 Longitude (deg)",fontsize=10)
    axsfletch[1].set_xlabel("Sys. 3 Longitude (deg)",fontsize=10)
    

    figfletch.subplots_adjust(left=0.10, bottom=0.03, right=0.98, top=0.95,
                wspace=0.25, hspace=0.05)     
    #axsfletch[1].set_position([box.x0+0.0, box.y0-0., box.width * 1., box.height * 1.])
    axsfletch[1].set_position([box.x0+0.03, box.y0-0.01, box.width * 1.015, box.height * 1.015])

    figfletch.savefig(pathout+"GRS_NH3_Comparison_Fletcher.png",dpi=300)
    
    ###########################################################################
    # MAKE PROFILE AVERAGE PLOT
    ###########################################################################

    #figprof,axsprof=pl.subplots(1,1,figsize=(6.0,4.0), dpi=150, facecolor="white")
    figprof,axsprof=pl.subplots(1,1,figsize=(6.0,6.0), dpi=150, facecolor="white")
    figprof.suptitle('Ammonia Abundance Meridional Profile',x=0.5,ha='center',color='k')

    AvgMeridEW=np.mean(MeridEWArray[:,:],axis=1)
    StdMeridEW=np.std(MeridEWArray[:,:],axis=1)

    MeridEW=np.flip(np.mean(fNH3_patch_avg[:,:],axis=1),axis=0)
    MeridEWerror=np.flip(np.std(fNH3_patch_avg[:,:],axis=1),axis=0)
    #Lats=np.linspace(-44.5,44.5,90)
    Lats=np.linspace(90.0-LatLims[1],90.0-LatLims[0],LatLims[1]-LatLims[0]) #!!!! Problem! doesn't work for 
    #!!!! anything other than 45-> -45 (45,135)

    #print('^^^^^^^^^',Lats)
    #print('^^^^^^^^^',MeridEW)
    axsprof.plot(Lats,AvgMeridEW,color='k',
                 label=r'This Work, Celestron 11 Seasonal Avg., $\bar{f_c} (NH_3)$') 
    axsprof.fill_between(Lats, AvgMeridEW-StdMeridEW, AvgMeridEW+StdMeridEW,color='k',alpha=.1)

    plevel=0.752910
    PTG.plot_TEXES_Groups(axsprof,clr='C2',prs=plevel,mult=1000000.)
    plevel=0.657540
    PTG.plot_TEXES_Groups(axsprof,clr='r',prs=plevel,mult=1000000.)
    plevel=0.574240
    PTG.plot_TEXES_Groups(axsprof,clr='C4',prs=plevel,mult=1000000.)
    
    path20220919='C:/Astronomy/Projects/Planets/Jupiter/Imaging Data/20220919UT/'
    VLTMUSEProfile=np.loadtxt(path20220919+
                              #'Profile of 2022-09-19-0352_3-Jupiter-fNH3Abs647-VLT-final.csv',
                              'Profile of 2022-09-19-0352_3-Jupiter-fNH3Abs647-VLT-final-2deg-wide.csv',
                              usecols=range(2),delimiter=",")
    axsprof.plot(VLTMUSEProfile[:,0]-45.,VLTMUSEProfile[:,1]*1e6,color='k',
                 linestyle="dashed",label=r'This Work, VLT-MUSE, $\bar{f_c} (NH_3)$') 
      
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

    for zb in belt:
        print(zb,belt[zb])
        axsprof.fill_between([belt[zb][0],belt[zb][1]],np.array([0.,0.]),np.array([1000.,1000.]),
                                color="0.5",alpha=0.2)
        axsprof.annotate(zb,xy=[np.mean(belt[zb]),1],ha="center")
    for zb in zone:
        axsprof.annotate(zb,xy=[np.mean(zone[zb]),1],ha="center")

    #for j in range(0,7):
    #    axsprof.plot(Lats,MeridEWArray[:,j])
    #    axsprof.fill_between(Lats, MeridEWArray[:,j]-MeridEWArrayStd[:,j],MeridEWArray[:,j]+MeridEWArrayStd[:,j],alpha=.2)
    axsprof.set_xlim(-30.,30.)
    axsprof.set_xlabel("Planetographic Latitude (deg)",fontsize=10)
    axsprof.set_ylim(0.,300.)
    axsprof.set_ylabel("Ammonia Abundance (ppm)",fontsize=10)
    axsprof.legend(loc=1,ncol=2, borderaxespad=0.,prop={'size':8})
    axsprof.grid(linewidth=0.2)
    axsprof.tick_params(axis='both', which='major', labelsize=8)

    figprof.subplots_adjust(left=0.10, bottom=0.12, right=0.98, top=0.92)  
    figprof.savefig(pathout+"GRS_NH3_AbundanceProfile_AVG_2022.png",dpi=300)

    #axsprof[0].fill_between(Lats, MeridEW-MeridEWerror, MeridEW+MeridEWerror,alpha=.2)
    #axsprof[0].set_title(Date+", CM2="+str(int(CM2deg))+", ASI120MM",fontsize=12)
    #axsprof[0].set_xlabel("Planetographic Latitude (deg)",fontsize=8)


