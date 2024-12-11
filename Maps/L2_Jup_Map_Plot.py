def L2_Jup_Map_Plot(obskey="20240925UTa",imagetype='Map',target="Jupiter",
                        Smoothing=False,LatLims=[45,135],LonRng=45,
                        CMpref='subobs',LonSys='2',showbands=False,
                        coef=[-0.05,-0.05],subproj='',figxy=[8.0,4.0],FiveMicron=False):
    """
    Created on Sun Nov  6 16:47:21 2022
    
    PURPOSE: Create maps of environmental parameters paired with RGB context
             maps. Based on Retrieve_Jup_Atm_2022_P3, which ALSO performed
             the calibration phase. So now I've separated that module into 
             a calibration module, make_L3_env_data.py and this plotting
             module.
             
    EXAMPLES:
        Map_Jup_Atm_P3(obskey="20240925UTa",imagetype='Map',target="Jupiter",
                                Smoothing=False,LatLims=[45,135],LonRng=45,
                                CMpref='subobs',LonSys='2',showbands=False,
                                coef=[0.,0.],subproj='',figxy=[8.0,4.0],
                                FiveMicron=False)
        
        Map_Jup_Atm_P3(obskey="20240730UTa",imagetype='Map',target="Jupiter",
                                Smoothing=False,LatLims=[45,135],LonRng=45,
                                CMpref='subobs',LonSys='2',showbands=False,
                                coef=[0.,0.],subproj='',figxy=[8.0,4.0],
                                FiveMicron=True)
    
    @author: smhil
    """    
    import sys
    drive='c:'
    sys.path.append(drive+'/Astronomy/Python Play')
    sys.path.append(drive+'/Astronomy/Python Play/Util_P3')
    sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy')
    sys.path.append('./Services')

    import os
    import pylab as pl
    import numpy as np
    from imageio import imwrite
    from astropy.io import fits
    import RetrievalLibrary as RL
    sys.path.append('./Maps')
    import read_fits_map_L2_L3 as RFM
    import plot_patch as PP
    import make_patch_RGB as MPRGB
    import copy
    import map_and_context as mac
    import map_and_scatter as mas

    fNH3low=0.0
    fNH3high=0.07
    PCldlow=0.1
    PCldhigh=0.15
    micronlow=0.5
    micronhigh=3.5
    
    if (not FiveMicron) or FiveMicron=="png":
        PCldhdr,PClddata,fNH3hdr,fNH3data,sza,eza,RGB,RGB_CM,RGBtime= \
                        RFM.read_fits_map_L2_L3(obskey=obskey,LonSys=LonSys,
                                                imagetype="Map",Level="L2",
                                                target=target,FiveMicron=FiveMicron)
    elif FiveMicron=="fits":
        PCldhdr,PClddata,fNH3hdr,fNH3data,sza,eza,RGB,RGB_CM,RGBtime,micronhdr,microndatar= \
                        RFM.read_fits_map_L2_L3(obskey=obskey,LonSys=LonSys,
                                                imagetype="Map",Level="L2",
                                                target=target,FiveMicron=FiveMicron)
                    
    pathmapplots='C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Analysis Data/L3 Plots/'+subproj+'/'
    if not os.path.exists(pathmapplots):
        os.makedirs(pathmapplots)
    ###########################################################################
    # Special for limb correction
    ###########################################################################             
    amfdata=(1.0/sza+1.0/eza)/2.0
    #figamf,axsamf=pl.subplots(figsize=(8.0,4.0), dpi=150, facecolor="white")
    #axsamf.imshow(amfdata,vmin=-5.,vmax=5.)
    if obskey=="20220730UTa":
        pathFITS='C:/Astronomy/Projects/Planets/Jupiter/Imaging Data/Mapping/'
        amf=fits.open(pathFITS+"2022-07-30-amf_CM2_L360_MAP-BARE.fit")
        amf.info()
        amfhdr=amf[0].header
        amfdata=5.*amf[0].data/65535.
        amf.close()
    elif obskey=="20220919UTa":                
        pathFITS='C:/Astronomy/Projects/Planets/Jupiter/Imaging Data/Mapping/'
        amf=fits.open(pathFITS+"2022-09-19-amf_CM2_L360_MAP-BARE.fit")
        amf.info()
        amfhdr=amf[0].header
        amfdata=5.*amf[0].data/65535.
        amf.close()

    #pl.imshow(amfdata)
    ###########################################################################
    # Set up figure and axes for plots
    ###########################################################################             
    if CMpref=='subobs':
        fNH3PlotCM=fNH3hdr["CM"+LonSys]
        PCldPlotCM=PCldhdr["CM"+LonSys]
    else:
        fNH3PlotCM=CMpref
        PCldPlotCM=CMpref
    NH3LonLims=[360-int(fNH3PlotCM+LonRng),360-int(fNH3PlotCM-LonRng)]
    print("#######fNH3PlotCM=",fNH3PlotCM)
    print("fNH3PlotCM+LonRng,fNH3PlotCM-LonRng=",fNH3PlotCM+LonRng,fNH3PlotCM-LonRng)
    print("#######NH3LonLims=",NH3LonLims)
    print("#######360-NH3LonLims=",360-np.array(NH3LonLims))
    if Smoothing:
        smthtitle="Smoothed"
    else: 
        smthtitle="Unsmoothed"
    CalModel=fNH3hdr['CALIBRA']

    ###########################################################################
    ## Just RGB and Abundance
    ###########################################################################
    fNH3_patch_mb,TestfNH3,tx_fNH3,fnNH3=mac.map_and_context(-np.log(fNH3data),fNH3hdr,
                                                       RGB,RGBtime,
                                                       LonSys,LatLims,NH3LonLims,
                                                       LonRng,fNH3PlotCM,
                                                       amfdata,coef[0],fNH3low,fNH3high,
                                                       showbands,FiveMicron,figxy,
                                                       "jet",pathmapplots,Level='L2',
                                                       suptitle="Ammonia Opacity",
                                                       cbar_title="Ammonia Opacity")

    ###########################################################################
    ## Just RGB and Cloud Pressure
    ###########################################################################
    PCld_patch,TestPCld,tx_PCld,fnPCld=mac.map_and_context(-np.log(PClddata),PCldhdr,
                                                        RGB,RGBtime,
                                                        LonSys,LatLims,NH3LonLims,
                                                        LonRng,fNH3PlotCM,
                                                        amfdata,coef[1],PCldlow,PCldhigh,
                                                        showbands,FiveMicron,figxy,
                                                        "jet",pathmapplots,Level='L2',
                                                        suptitle="Methane Opacity",
                                                        cbar_title="Methane Opacity")

    ###########################################################################
    ## Just RGB and 5 micron
    ###########################################################################
    micron_patch,Testmicron,tx_micron,fn5um=mac.map_and_context(np.log10(microndatar),micronhdr,
                                                          RGB,RGBtime,
                                                          LonSys,LatLims,NH3LonLims,
                                                          LonRng,fNH3PlotCM,
                                                          amfdata,0.0,micronlow,micronhigh,
                                                          showbands,FiveMicron,figxy,
                                                          "gist_heat",pathmapplots,Level='L2',
                                                          suptitle="5 micron Radiance (log10)",
                                                          cbar_title="Log10(5um radiance)")

       
    ###########################################################################
    ## Compute Scatter Plot (PCloud vs fNH3)
    ###########################################################################
    mas.map_and_scatter(fNH3_patch_mb,PCld_patch,-np.log(PClddata),fNH3hdr,LonSys,
                        LatLims,NH3LonLims,LonRng,PCldPlotCM,fnNH3,
                        coef[0],tx_fNH3,fNH3low,fNH3high,PCldlow,PCldhigh,
                        figxy,"gray",pathmapplots,"Methane Opacity & Ammonia Opacity (contours)",
                        "Methane vs Ammonia Opacity",Level='L2',cbar_rev=False,
                        axis_inv=False,cbar_title="Methane Opacity")
    
    ###########################################################################
    ## Compute Scatter Plot (PCloud vs 5um radiance)
    ###########################################################################
    mas.map_and_scatter(PCld_patch,micron_patch,np.log10(microndatar),fNH3hdr,LonSys,
                        LatLims,NH3LonLims,LonRng,PCldPlotCM,fnNH3,
                        coef[1],tx_PCld,PCldlow,PCldhigh,micronlow,micronhigh,
                        figxy,"gist_heat",pathmapplots,"5um Radiance & Methane Opacity (contours)",
                        "Methane Opacity vs 5um Radiance",FiveMicron=True,Level='L2',cbar_rev=False,swap_xy=True,
                        axis_inv=False,cbar_title="Log10(5um Radiance)")

    ###########################################################################
    ## Compute Scatter Plot (fNH3 vs 5um radiance)
    ###########################################################################
    mas.map_and_scatter(fNH3_patch_mb,micron_patch,np.log10(microndatar),fNH3hdr,LonSys,
                        LatLims,NH3LonLims,LonRng,PCldPlotCM,fnNH3,
                        0.0,tx_fNH3,fNH3low,fNH3high,micronlow,micronhigh,
                        figxy,"gist_heat",pathmapplots,"5um Radiance & Ammonia Opacity (contours)",
                        "Ammonia Opacity vs 5um Radiance",FiveMicron=True,Level='L2',cbar_rev=False,swap_xy=True,
                        axis_inv=False,cbar_title="Log10(5um Radiance)")
    
    
    
    #return(fig1,axs1,fig2,axs2,fig3,axs3)


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

