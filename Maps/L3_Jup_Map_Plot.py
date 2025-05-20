def L3_Jup_Map_Plot(obskey="20250106UTf",imagetype='Map',target="Jupiter",
                        Smoothing=False,LatLims=[45,135],LonRng=45,
                        CMpref='subobs',LonSys='2',showbands=False,
                        coef=[0.,0.],subproj='',figxy=[8.0,4.0],FiveMicron=False,
                        ROI=False,ctbls=["terrain_r","Blues"]):
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
    import numpy as np
    from astropy.io import fits
    sys.path.append('./Maps')
    import read_fits_map_L2_L3 as RFM
    import make_patch as MP
    import map_and_context as mac
    import map_and_scatter as mas

    if ctbls[0]=="jet":
        fNH3low=60
        fNH3high=160
        PCldlow=1200
        PCldhigh=2000
    if ctbls[0]=="terrain_r":
        fNH3low=60
        fNH3high=160
        PCldlow=1600
        PCldhigh=2200
        
    micronlow=0.5
    micronhigh=3.5

    if (not FiveMicron) or FiveMicron=="png":
        PCldhdr,PClddata,fNH3hdr,fNH3data,sza,eza,RGB,RGB_CM,RGBtime= \
                        RFM.read_fits_map_L2_L3(obskey=obskey,LonSys=LonSys,
                                                imagetype="Map",Level="L3",
                                                target=target,FiveMicron=FiveMicron)
    elif FiveMicron=="fits":
        PCldhdr,PClddata,fNH3hdr,fNH3data,sza,eza,RGB,RGB_CM,RGBtime,micronhdr,microndatar= \
                        RFM.read_fits_map_L2_L3(obskey=obskey,LonSys=LonSys,
                                                imagetype="Map",Level="L3",
                                                target=target,FiveMicron=FiveMicron)
                    
    pathmapplots='C:/Astronomy/Projects/SAS 2021 Ammonia/Data/L3 Plots/'+subproj+'/'
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
        if FiveMicron:
            micronPlotCM=micronhdr["CM"+LonSys]
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
    amfpatch=MP.make_patch(amfdata,LatLims,NH3LonLims,fNH3PlotCM,LonRng,pad=True)

    ###########################################################################
    ## Just RGB and Abundance
    ###########################################################################
    fNH3_patch_mb,TestfNH3,tx_fNH3,fnNH3=mac.map_and_context(fNH3data,fNH3hdr,
                                                       RGB,RGBtime,
                                                       LonSys,LatLims,NH3LonLims,
                                                       LonRng,fNH3PlotCM,
                                                       amfdata,coef[0],fNH3low,fNH3high,
                                                       showbands,FiveMicron,figxy,
                                                       ctbls[0],pathmapplots,Level='L3',
                                                       suptitle="Ammonia Mole Fraction",
                                                       cbar_title="Ammonia Mole Fraction (ppm)",
                                                       ROI=ROI)
    

    ###########################################################################
    ## Just RGB and Cloud Pressure
    ###########################################################################
    PCld_patch,TestPCld,tx_PCld,fnPCld=mac.map_and_context(PClddata,PCldhdr,
                                                        RGB,RGBtime,
                                                        LonSys,LatLims,NH3LonLims,
                                                        LonRng,PCldPlotCM,
                                                        amfdata,coef[1],PCldlow,PCldhigh,
                                                        showbands,FiveMicron,figxy,
                                                        ctbls[1],pathmapplots,Level='L3',
                                                        suptitle="Cloud Top Pressure",
                                                        cbar_rev=True,
                                                        cbar_title="Cloud Top Pressure (mb)",
                                                        ROI=ROI)

    ###########################################################################
    ## Just RGB and 5 micron
    ###########################################################################
    if FiveMicron:
        micron_patch,Testmicron,tx_micron,fn5um=mac.map_and_context(np.log10(microndatar),micronhdr,
                                                          RGB,RGBtime,
                                                          LonSys,LatLims,NH3LonLims,
                                                          LonRng,PCldPlotCM,
                                                          amfdata,0.0,micronlow,micronhigh,
                                                          showbands,FiveMicron,figxy,
                                                          "gist_heat",pathmapplots,Level='L3',
                                                          suptitle="5 micron Radiance (log10)",
                                                          cbar_title="Log10(5um radiance)")

   
    ###########################################################################
    ## Compute Band or ROI Scatter Plot (PCloud vs fNH3)
    ###########################################################################
    dateobs,roilabel,mean1,stdv1,mean2,stdv2,meanamf=\
        mas.map_and_scatter(fNH3_patch_mb,PCld_patch,PClddata,fNH3hdr,LonSys,
        LatLims,NH3LonLims,LonRng,PCldPlotCM,fnNH3,
        coef[0],tx_fNH3,fNH3low,fNH3high,PCldlow,PCldhigh,
        figxy,"Blues",pathmapplots,"PCloud & fNH3 (contours)",
        "PCloud vs fNH3",Level='L3',cbar_rev=True,cbar_title="Cloud-top Pressure (mb)",
        axis_inv=True,ROI=ROI,amfpatch=amfpatch)
    
    ###########################################################################
    ## Compute Scatter Plot (PCloud vs 5um radiance)
    ###########################################################################
    if FiveMicron:
        mas.map_and_scatter(PCld_patch,micron_patch,np.log10(microndatar),fNH3hdr,LonSys,
                        LatLims,NH3LonLims,LonRng,PCldPlotCM,fnNH3,
                        coef[1],tx_PCld,PCldlow,PCldhigh,micronlow,micronhigh,
                        figxy,"gist_heat",pathmapplots,"5um Radiance & PCloud (contours)",
                        "PCloud vs 5um Radiance",Level='L3',FiveMicron=True,cbar_rev=False,swap_xy=True,
                        axis_inv=True,cbar_title="Log10(5um radiance)")

    ###########################################################################
    ## Compute Scatter Plot (fNH3 vs 5um radiance)
    ###########################################################################
    if FiveMicron:
        mas.map_and_scatter(fNH3_patch_mb,micron_patch,np.log10(microndatar),fNH3hdr,LonSys,
                        LatLims,NH3LonLims,LonRng,PCldPlotCM,fnNH3,
                        0.0,tx_fNH3,fNH3low,fNH3high,micronlow,micronhigh,
                        figxy,"gist_heat",pathmapplots,"fNH3 vs 5um Radiance",
                        "5um Radiance & fNH3 (contours)",Level='L3',FiveMicron=True,
                        cbar_rev=False,swap_xy=False,
                        axis_inv=True,cbar_title="Log10(5um radiance)")
   
    #return(fig1,axs1,fig2,axs2,fig3,axs3)
    return(dateobs,roilabel,mean1,stdv1,mean2,stdv2,meanamf)


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

