# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 12:49:14 2018

@author: Steven Hill

PURPOSE:    This code creates arrays of map patches to evaluate NH3 absorption
            compared to visible RGB imagery, visible panchromatic reflectivity,
            889CH4 band reflectivity, 380NUV band reflectivity, and local 
            color slope.

HISTORY:    This code evolved from PlanetMapTest.py committed on 2/21/2019 at
            https://github.com/smhill001/PlanetMaps/blob/master/PlanetMapTest.py.
            That code provided simple equatorial and polar map views for any planet
            as well as a bit of limited patch mapping.
            
VERSION:    2.0 (12/31/2021) - Produces CMOS and CCD maps and given date as long
            as the file naming and data format conventions are followed.
            - Produces zonal and meridional 645nm EW profiles for each single 
              NH3 observation session.
            - Produces (landscape only) map patch arrays with countour overlays of 
              NH3 absorption EW. Maps patche extent is hardcoded to +/-45 deg
              of the equator and +/-45 deg planetographic latitude from the 
              central meridion of the NH3 absorption map. Maps are computed in
              System II longitude.
            - Produces cumulative (average) profiles with standard deviations
              for a campaign of sessions with overlays of data from Teifel and Fletcher.
            - Deprecates need for an external configuration file (mapconfig.txt).
            
            1.0 - Produces all CMOS maps for 2020 and 2021 and all CCD maps
            for 2020. 
            - Has limited or unverified capability to produce cumulative 
              patches for meridional analysis.
            - Does plots Cartopy longitudes, not correct Jovian longitudes
              from WinJUPOS
            - Poorly commented and probably has extraneous code
            - Dependency on F:\Astronomy\Python Play\PlanetMapsMapConfig.txt
            
EXAMPLE:    AmmoniaMaps(DateSelection=["2021-09-23"])

"""

def CH4_EW_Map_Test_P3(coords='map',cont=True,zonecorr=[0,0],DateSelection='All',orientation='Landscape'):
    ###########################################################################
    # Set up environment and paths
    #
    import sys
    drive='C:'
    sys.path.append(drive+'/Astronomy/Python Play')
    sys.path.append(drive+'/Astronomy/Python Play/Util_P3')
    sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy_P3')
    #sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy')
    sys.path.append(drive+'/Astronomy/Python Play/SPLibraries_P3')
    #import scipy.ndimage as nd
    from matplotlib.pyplot import imread
    import pylab as pl
    import numpy as np
    from astropy.convolution import Gaussian2DKernel
    from astropy.convolution import convolve
    from datetime import datetime
    import ephem
    import EWLibV006_P3 as EWL
    import plot_TEXES_Groups_P3 as PTG

    ###########################################################################
    #Initialization section
    drive="C:"
    path="/Astronomy/Projects/Planets/Jupiter/Imaging Data/Mapping/"
    
    amagat=2.69e24 #Lodschmits number?
    gravity=2228.0
    mean_mol_wt=3.85e-24
    fCH4=1.81e-3
    STP=1.01e6

    observer = ephem.Observer()
    #Location from Google Maps at 483 S Oneida Way, Denver 80224
    observer.lon = ephem.degrees('-104.907985')
    observer.lat = ephem.degrees('39.708200')
    
    CM2=NewGetMaps()     
    print("############################")
    print("len(CM2_L360)=",len(CM2))

    #Set up labels for the six kinds of maps
    PlotTypes=["NH3Abs","CH4Abs","ClrSlp"]
    
    #Create empty arrays for the aggregation of NH3 absorption patch data, to
    #  be used as input to ProfileComparison.py
    MeridEWArray=np.zeros((90,len(DateSelection)))
    ZoneEWArray=np.zeros((90,len(DateSelection)))
    ###########################################################################
    #Begin looping over observing sessions
    DateCounter=0
    for Date in DateSelection:
        #Compute Ephermeris for NH3 observation CM2
        MapsforDate=[k for k in CM2 if Date in k]              
        NH3Abs_fn=[k for k in MapsforDate if "NH3Abs" in k]        
        strdate=NH3Abs_fn[0][0:15]
        print(strdate)
        dates=[datetime.strptime(strdate,"%Y-%m-%d-%H%M")]
        date_list_datetime,elev_list,airmass_list,CMI_list,CMII_list,\
            Io_vis_list,Europa_vis_list, \
            Ganymede_vis_list,Callisto_vis_list= \
            EWL.JupiterEphemLists(dates,observer)
        print(CMII_list[0]*180./np.pi+45,CMII_list[0]-45*180./np.pi)
        LatLims=[45,135]
        CM2deg=CMII_list[0]*180./np.pi
        print("CM2deg=",CM2deg)
        LonLims=[360-int(CM2deg+45),360-int(CM2deg-45)]
        print("LonLims=",LonLims)
        
        #Left pad CM2 strings 
        if int(CM2deg)<10:
            CM2str="00"+str(int(CM2deg))
        elif int(CM2deg)<100:
            CM2str="0"+str(int(CM2deg))
        else:
            CM2str=str(int(CM2deg))
            
        print("Date=",Date)     
        ##### Compute Jupiter ephemeris
        #MapsforDate=[k for k in CM2 if Date in k]              
        #######################################################################
        # Create CH4Abs patch and scale to EW(nm) for contour computations
        CH4Abs_fn=[k for k in MapsforDate if "CH4Abs" in k]        
        tmp=load_png(drive+path+CH4Abs_fn[0])
        CH4Abs_map=tmp[:,:,0]
        CH4Abs_patch=make_patch(CH4Abs_map,LatLims,LonLims,CM2deg)
        CH4Abs_patch_EW=(1.0-0.104*(CH4Abs_patch*2.0/((2.0**16.)-1.0)+0.0))*20.2/0.896
        
        # Make smoothed patch for NH3Abs contours
        kernel = Gaussian2DKernel(1)
        print(kernel)
        #Empirical longitudinal limb darkening correction
        if zonecorr[:]==[0,0]:
            ZoneLimbDarkCorr=np.ones(90)
        else:
            ZoneLimbDarkCorr=np.mean(CH4Abs_patch_EW[45-zonecorr[0]:45-zonecorr[1],:],axis=0) ###To flip or not to flip in long?
        ZoneLimbDarkCorr=ZoneLimbDarkCorr/np.max(ZoneLimbDarkCorr)
        print(ZoneLimbDarkCorr.shape,CH4Abs_patch_EW.shape)
        CH4Abs_patch_EW_corr=CH4Abs_patch_EW/ZoneLimbDarkCorr[None,:]
        CH4Abs_patch_EW_conv = convolve(CH4Abs_patch_EW_corr, kernel)
        
        CH4Abs_patch_EW_corr_WaveNum=(CH4Abs_patch_EW_corr*((1./(889.0*1e-8))/889.0))/10.
        NCH4_lin=3.0*(CH4Abs_patch_EW_corr_WaveNum/2.2)/5753. #band strength from Hill,2015 work
        #Random factor of 3 above to 'compensate' for band saturation.
        NCH4_lin_conv = convolve(NCH4_lin, kernel)
        Pcloud=NCH4_lin*amagat*gravity*mean_mol_wt/(fCH4*STP)
        Pcloud_conv = convolve(Pcloud, kernel)

        #######################################################################
        # Create NH3Abs patch and scale to EW(nm) for contour computations
        tmp=load_png(drive+path+NH3Abs_fn[0])
        NH3Abs_map=tmp[:,:,0]
        NH3Abs_patch=make_patch(NH3Abs_map,LatLims,LonLims,CM2deg)
        NH3Abs_patch_EW=(1.0-0.961*(NH3Abs_patch*0.2/((2.0**16.)-1.0)+0.9))*0.55/0.039
        
        # Make smoothed patch for NH3Abs contours
        kernel = Gaussian2DKernel(1)
        print(kernel)
        #Empirical longitudinal limb darkening correction
        if zonecorr[:]==[0,0]:
            ZoneLimbDarkCorr=np.ones(90)
        else:
            ZoneLimbDarkCorr=np.mean(NH3Abs_patch_EW[45-zonecorr[0]:45-zonecorr[1],:],axis=0) ###To flip or not to flip in long?
        ZoneLimbDarkCorr=ZoneLimbDarkCorr/np.max(ZoneLimbDarkCorr)
        print(ZoneLimbDarkCorr.shape,NH3Abs_patch_EW.shape)
        NH3Abs_patch_EW_corr=NH3Abs_patch_EW/ZoneLimbDarkCorr[None,:]
        NH3Abs_conv = convolve(NH3Abs_patch_EW_corr, kernel)

        NH3Abs_patch_EW_corr_WaveNum=(NH3Abs_patch_EW_corr*((1./(646.0*1e-8))/646.0))/10.
        NNH3_lin=(NH3Abs_patch_EW_corr_WaveNum/2.2)/630. #band strength from Lutz&Owen 1980
        NNH3_lin_conv = convolve(NNH3_lin, kernel)


        fNH3=NNH3_lin*amagat*gravity*mean_mol_wt/(Pcloud*STP)

        #######################################################################
        # Compute and plot single-session meridional and zonal profiles from 
        #   unadjusted patch. Also store the profiles for creation of average
        #   profile plots across all sessions in a give run.
        MeridEW=np.flip(np.mean(CH4Abs_patch_EW[:,:],axis=1),axis=0)
        MeridEWerror=np.flip(np.std(CH4Abs_patch_EW[:,:],axis=1),axis=0)
        Lats=np.linspace(-44.5,44.5,90)
        #print DateCounter; Date
        MeridEWArray[:,DateCounter]=MeridEW[:]
        
        ZoneEW=np.mean(CH4Abs_patch_EW[:,:],axis=0) 
        ZoneEWerror=np.std(CH4Abs_patch_EW[:,:],axis=0)
        Lons=np.linspace(-44.5,44.5,90)
        ZoneEWArray[:,DateCounter]=ZoneEW[:]
        
        figprof,axsprof=pl.subplots(2,1,figsize=(6.0,6.0), dpi=150, facecolor="white")
        figprof.suptitle("Ammonia Absorption Profile",x=0.5,ha='center',color='k')

        axsprof[0].plot(Lats,MeridEW)
        axsprof[0].fill_between(Lats, MeridEW-MeridEWerror, MeridEW+MeridEWerror,alpha=.2)
        axsprof[0].set_title(Date+", CM2="+str(int(CM2deg))+", ASI120MM",fontsize=12)
        axsprof[0].set_xlabel("Planetographic Latitude (deg)",fontsize=8)

        axsprof[1].plot(Lons,ZoneEW)
        axsprof[1].fill_between(Lons, ZoneEW-ZoneEWerror, ZoneEW+ZoneEWerror,alpha=.2)
        axsprof[1].set_xlabel("Longitude from Meridian (deg)",fontsize=8)

        for iprof in range(0,2):
            axsprof[iprof].grid(linewidth=0.2)
            axsprof[iprof].xlim=[-45.,45.]
            axsprof[iprof].ylim=[18.,22.]
            axsprof[iprof].set_xticks(np.linspace(-45.,45.,7), minor=False)
            axsprof[iprof].set_yticks(np.linspace(18.,22.,9), minor=False)
            axsprof[iprof].tick_params(axis='both', which='major', labelsize=7)
            axsprof[iprof].set_ylabel("Equivalent Width (nm)",fontsize=8)

        figprof.savefig(drive+path+"NH3 Map Plots/"+Date+"-Jupiter-NH3"+"_CMII_"+
                   CM2str+"-Profile.png",dpi=300)
        
        #######################################################################
        #Begin Looping over individual patches or bands
        fig,axs=pl.subplots(2,3,figsize=(6.0,4.5), dpi=150, facecolor="white",
                            sharey=True,sharex=True)
        fig.suptitle(Date+", CM2="+str(int(CM2deg)),x=0.5,ha='center',color='k')

        figcor,axscor=pl.subplots(2,2,figsize=(6.0,4.5), dpi=150, facecolor="white",
                            sharey=True,sharex=True)
        figcor.suptitle(Date+", CM2="+str(int(CM2deg)),x=0.5,ha='center',color='k')
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
            
        ###############################################################
        # NH3 Plots
        #
        mapshow=axs[0,0].imshow(NH3Abs_patch_EW_corr, "gist_heat", origin='upper',  
                   extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],
                           90-LatLims[0]],vmin=0,vmax=1.2,
                           aspect="equal")
        axs[0,0].set_title("NH3Abs EW(nm)",fontsize=8)

        mapshow=axs[0,1].imshow(NNH3_lin, "gist_heat", origin='upper',  
                   extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],
                           90-LatLims[0]],vmin=0.005,vmax=0.015,
                           aspect="equal")
        axs[0,1].set_title("NH3 (km-amagat)",fontsize=8)

        mapshow=axs[0,2].imshow(fNH3, "gist_heat", origin='upper',  
                   extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],
                           90-LatLims[0]],vmin=0.0002,vmax=0.0005,
                           aspect="equal")
        axs[0,2].set_title("fNH3 (vmf)",fontsize=8)
        
        if cont:
            temp=make_contours_CH4(axs[0,0],NH3Abs_conv,LatLims,LonLims,
                                   lvls=[0.4,0.5,0.6,0.7,0.8])
            temp=make_contours_CH4(axs[0,1],NNH3_lin,LatLims,LonLims,
                                   lvls=[0.0050,0.0075,0.0100,0.0125,0.0150])
            temp=make_contours_CH4(axs[0,2],fNH3,LatLims,LonLims,
                                   lvls=[0.00025,0.00030,0.00035,0.0004])
        
        ###############################################################
        # CH4 Plots
        #
        mapshow=axs[1,0].imshow(CH4Abs_patch_EW, "gist_heat", origin='upper',  
                   extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],
                           90-LatLims[0]],vmin=19.0,vmax=21.0,
                           aspect="equal")
        axs[1,0].set_title("CH4Abs EW(nm)",fontsize=8)

        mapshow=axs[1,1].imshow(NCH4_lin, "gist_heat", origin='upper',  
                   extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],
                           90-LatLims[0]],vmin=0.059,vmax=0.063,
                           aspect="equal")
        axs[1,1].set_title("NCH4(km-amagat)",fontsize=8)

        mapshow=axs[1,2].imshow(Pcloud, "gist_heat", origin='upper',  
                   extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],
                           90-LatLims[0]],vmin=0.72,vmax=0.82,
                           aspect="equal")
        axs[1,2].set_title("Cloud Top Pressure (Pa)",fontsize=8)
        if cont:
            temp=make_contours_CH4(axs[1,0],CH4Abs_patch_EW_conv,LatLims,LonLims,
                                   lvls=[19.0,19.5,20.0,20.5,21.0])
            temp=make_contours_CH4(axs[1,1],NCH4_lin,LatLims,LonLims,
                                   lvls=[0.059,0.060,0.061,0.062,0.063])
            temp=make_contours_CH4(axs[1,2],Pcloud_conv,LatLims,LonLims,
                                   lvls=[0.72,0.74,0.76,0.78,0.80])
            
        DateCounter=DateCounter+1
        fig.subplots_adjust(left=0.10, bottom=0.08, right=0.98, top=0.90,
                    wspace=0.10, hspace=0.10)            
        fig.savefig(drive+path+"NH3 Map Plots/"+Date+"-Jupiter-NH3"+"_CMII_"+
                   CM2str+"-Map.png",dpi=300)
        
        figcor.subplots_adjust(left=0.10, bottom=0.10, right=0.98, top=0.90,
                    wspace=0.10, hspace=0.20)            
        figcor.savefig(drive+path+"NH3 Map Plots/"+Date+"-Jupiter-NH3"+"_CMII_"+
                   CM2str+"-Correlation.png",dpi=300)

    ###########################################################################
    # Compute and plot cumulative profiles from campaign sessions
    #
    AvgMeridEW=np.mean(MeridEWArray[:,:],axis=1)
    StdMeridEW=np.std(MeridEWArray[:,:],axis=1)
    AvgZoneEW=np.mean(ZoneEWArray[:,:],axis=1)
    StdZoneEW=np.std(ZoneEWArray[:,:],axis=1)
    figavgprof,axsavgprof=pl.subplots(2,1,figsize=(6.0,6.0), dpi=150, facecolor="white")
    figavgprof.suptitle("Average Ammonia Absorption Profile",x=0.5,ha='center',color='k')

    # Meridional profile & reference data
    axsavgprof[0].plot(Lats,AvgMeridEW,label='This Work',color='C0')
    axsavgprof[0].fill_between(Lats, AvgMeridEW-StdMeridEW, AvgMeridEW+StdMeridEW,color='C0',alpha=.2)
    axsavgprof[0].set_xlabel("Planetographic Latitude (deg)",fontsize=8)

    #PTG.plot_Teifel(axsavgprof[0],clr='C0')
    axsavgprof[0].legend(fontsize=7,loc=2)

    ax2 = axsavgprof[0].twinx()  # instantiate a second axes that shares the same x-axis
    ax2.ticklabel_format(axis='y',style='sci',scilimits=(0,1))
    ax2.tick_params(axis='y', which='major', labelsize=7)
    PTG.plot_TEXES_Groups(ax2,clr='C2')

    ax2.set_ylabel("NH3 Mole Fraction at 440mb",fontsize=8)
    ax2.yaxis.label.set_color('C2')
    ax2.legend(fontsize=7,loc=1)

    # Zonal Profile
    axsavgprof[1].plot(Lons,AvgZoneEW)
    axsavgprof[1].fill_between(Lons, AvgZoneEW-StdZoneEW, AvgZoneEW+StdZoneEW,alpha=.2)
    axsavgprof[1].set_xlabel("Longitude (deg)",fontsize=8)
    
    # Plot layout details and labeling
    for iavgprof in range(0,2):
        axsavgprof[iavgprof].grid(linewidth=0.2)
        axsavgprof[iavgprof].xlim=[-45.,45.]
        axsavgprof[iavgprof].ylim=[18.,22.]
        axsavgprof[iavgprof].set_xticks(np.linspace(-45.,45.,7), minor=False)
        axsavgprof[iavgprof].set_yticks(np.linspace(18.,22.,9), minor=False)
        axsavgprof[iavgprof].tick_params(axis='both', which='major', labelsize=7)
        axsavgprof[iavgprof].set_ylabel("Equivalent Width (nm)",fontsize=8)
        
    figavgprof.savefig(drive+path+"NH3 Map Plots/Jupiter-NH3_CMII_AvgProfile.png",dpi=300)
    return 0


def NewGetMaps():
    import sys
    drive='c:'
    sys.path.append(drive+'/Astronomy/Python Play')
    sys.path.append(drive+'/Astronomy/Python Play/Util_P3')
    sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy_P3')

    import os

    path='c:/Astronomy/Projects/Planets/Jupiter/Imaging Data/Mapping/'
    fnlist = os.listdir(path)
    #print fnlist
    pnglist=[k for k in fnlist if '.png' in k]
    #print pnglist
    CM2_L360=[k for k in pnglist if 'CM2_L360_MAP-BARE' in k]
    #print wavelets
    CM2_L360.sort()
    
    filelistCM2=CM2_L360
    
    return filelistCM2

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

def make_contours_CH4(ax,CH4Abs_conv,LatLims,LonLims,lvls=[0.71,0.73,0.75,0.77,0.79]):
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
    ax.clabel(cs,lvls,inline=True,fmt='%3.1e',fontsize=8)

