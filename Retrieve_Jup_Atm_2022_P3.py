# -*- coding: utf-8 -*-


def Retrieve_Jup_Atm_2022_P3(obsdate="20220919UTa",target="Jupiter",
                             imagetype='Map',CalModel='VLT-Obs-Final',
                             Smoothing=True,LatLims=[45,135],LonRng=45,
                             delta_CM2=0,LonSys='2',showbands=False):
    """
    Created on Sun Nov  6 16:47:21 2022
    
    #!!! Need to figure out how best to manage longitude limits on plots,
    #!!! e.g., where to use NH3LonLims vs something else
    
    @author: smhil
    """    
    import sys
    drive='c:'
    sys.path.append(drive+'/Astronomy/Python Play')
    sys.path.append(drive+'/Astronomy/Python Play/Util_P3')
    sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy')
    sys.path.append('./Services')

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
    import RetrievalLibrary as RL
    import read_master_calibration
    sys.path.append('./Services')
    import get_abs_obs_data as GAOD
    import Retrieve_L3_Env_Params as RL3
    
    ###########################################################################
    # Retrieve_Jup_Atm_2022_P3(obsdate="20221019UT",target="Jupiter")   
    #CH4_tau,CH4_Ncol,CH4_Cloud_Press,NH3_tau,NH3_Ncol,fNH3,RGBfile,RGB,CM,NH3time,NH3file,path,pathfNH3=\
    #    RL3.Retrieve_L3_Env_Params(obsdate=obsdate,target="Jupiter",
    #                           imagetype='Map',CalModel='SCT-Obs-Final',
    #                           Smoothing=True,LonSys=LonSys)
    R=RL3.Retrieve_L3_Env_Params(obsdate=obsdate,target="Jupiter",
                               imagetype='Map',CalModel='SCT-Obs-Final',
                               Smoothing=True,LonSys=LonSys)
    ###########################################################################
    # Set up figure and axes for plots
    ###########################################################################             
    #CM2=Real_CM2+delta_CM2
    NH3LonLims=[360-int(R["NH3"]["CM"]+LonRng),360-int(R["NH3"]["CM"]-LonRng)]
    if Smoothing:
        smthtitle="Smoothed"
    else: 
        smthtitle="Unsmoothed"
        
    fig,axs=pl.subplots(2,3,figsize=(7.0,4.5), dpi=150, facecolor="white",
                        sharey=True,sharex=True)      
    fig.suptitle(R["NH3"]["time"].replace("_"," ")+", CM"+LonSys+"="
                 +str(int(R["NH3"]["CM"]))+", Calibration = "+CalModel+", Data "
                 +smthtitle,x=0.5,ha='center',color='k')

    if imagetype=='Map':
        for iPlot in range(0,6):
            i=int(iPlot/3)                      #Plot row
            j=np.mod(iPlot,3)                   #Plot column
            axs[i,j].grid(linewidth=0.2)
            axs[i,j].ylim=[-45.,45.]
            axs[i,j].xlim=[360-NH3LonLims[0],360-NH3LonLims[1]]
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
        plot_patch(R["CH4"]["tau"],LatLims,NH3LonLims,R["NH3"]["CM"],LonRng,"gist_heat",axs[0,0],
                   '%3.3f',vn=0.05,vx=0.15,n=6)
    else:
        show=axs[0,0].imshow(R["CH4"]["tau"], "gist_heat", origin='upper',vmin=0.10,vmax=0.20,
                           aspect="equal")
        #temp=RL.make_contours_CH4(axs[0,0],np.flip(convolve(CH4_tau,kernel),axis=0),
        #                       lvls=np.linspace(0.15,0.20,num=5,endpoint=True))
    
    axs[0,0].set_title("CH4 620nm Optical Depth",fontsize=8)

    ###########################################################################
    # Plot CH4 Column Density (m-atm)
    ###########################################################################
    if imagetype=='Map':
        plot_patch(R["CH4"]["Ncol"],LatLims,NH3LonLims,R["CH4"]["CM"],LonRng,"gist_heat",axs[0,1],
                   '%3.0f',vn=150,vx=350,n=6)
    else:
        show=axs[0,1].imshow(R["CH4"]["Ncol"]*1000., "gist_heat", 
                             origin='upper',vmin=200,vmax=500,aspect="equal")
        #temp=RL.make_contours_CH4(axs[0,1],np.flip(convolve(CH4_Ncol,kernel)*1000.,axis=0),
        #                       lvls=np.linspace(200.,500.,num=5,endpoint=True),frmt='%3.0f')
        
    axs[0,1].set_title("N[CH4] (m-atm)",fontsize=8)

    ###########################################################################
    # Plot CH4 Cloud Top Pressure (mb?) OR RGB Image
    ###########################################################################    
    if R["RGB"]["file"] == 'NA':
        show=axs[0,2].imshow(R["CH4"]["PCloud"], "gist_heat", origin='upper',vmin=2.0,vmax=6.0,
                           aspect="equal")
        temp=RL.make_contours_CH4(axs[0,2],R["CH4"]["PCloud"],
                               lvls=np.linspace(2.0,6.0,num=5,endpoint=True))
        axs[0,2].set_title("Scattering Pressure (Pa)",fontsize=8)
    else:
        if imagetype=='Map':
            RGB_patch=make_patch_RGB(R["RGB"]["RGB"],LatLims,NH3LonLims,R["RGB"]["CM"],LonRng)
            print("@@@@@@@ RGB_patch.shape()=",RGB_patch.shape)
            show=axs[0,2].imshow(RGB_patch, vmin=0,vmax=2^16,  
                       extent=[360-NH3LonLims[0],360-NH3LonLims[1],90-LatLims[1],
                               90-LatLims[0]],#vmin=0,vmax=1.2,
                               aspect="equal")
        else:    
            show=axs[0,2].imshow(R["RGB"]["RGB"], aspect="equal",vmin=0,vmax=2^16)
            axs[0,2].set_title("Scattering Pressure (Pa)",fontsize=8)

    ###########################################################################
    # Plot NH3 Optical Depth
    ###########################################################################
    if imagetype=='Map':
        plot_patch(R["NH3"]["tau"],LatLims,NH3LonLims,R["NH3"]["CM"],LonRng,
                   "gist_heat",axs[1,0],'%3.3f',vn=0.02,vx=0.05,n=6)
    else:
        show=axs[1,0].imshow(R["NH3"]["tau"], "gist_heat", origin='upper',
                             vmin=0.02,vmax=0.12,aspect="equal")
        #temp=RL.make_contours_CH4(axs[1,0],np.flip(convolve(NH3_tau,kernel),axis=0),
        #                       lvls=np.linspace(0.02,0.12,num=5,endpoint=True))
        
    axs[1,0].set_title("NH3 647nm Optical Depth",fontsize=8)

    ###########################################################################
    # Plot NH3 Column Density (km-atm)
    ###########################################################################
    if imagetype=='Map':
        plot_patch(R["NH3"]["Ncol"],LatLims,NH3LonLims,R["NH3"]["CM"],LonRng,"gist_heat",axs[1,1],
                   '%3.1f',vn=5,vx=20,n=6)
    else:    
        show=axs[1,1].imshow(R["NH3"]["Ncol"]*1000., "gist_heat", origin='upper',
                             vmin=10,vmax=40,aspect="equal")
        #temp=RL.make_contours_CH4(axs[1,1],np.flip(convolve(NH3_Ncol,kernel)*1000.,axis=0),
        #                       lvls=np.linspace(20.,40.,num=5,endpoint=True),frmt='%3.0f')
    
    axs[1,1].set_title("N[NH3] (m-atm)",fontsize=8)
    
    ###########################################################################
    # Plot NH3 Mole Fraction - Ammonia Absorption Index
    ###########################################################################
    if imagetype=='Map':
        fNH3_patch_mb,vn_fNH3,vx_fNH3,tx_fNH3=\
            plot_patch(R["NH3"]["fNH3"]*1e6,LatLims,NH3LonLims,R["NH3"]["CM"],
                       LonRng,"jet",axs[1,2],'%3.0f',vn=50,vx=150,n=6)

    else:    
        show=axs[1,2].imshow(R["NH3"]["fNH3"], "jet", origin='upper',
                             vmin=100,vmax=160,aspect="equal")
        #temp=RL.make_contours_CH4(axs[1,2],np.flip(convolve(fNH3*1.0e6,kernel),axis=0),
        #                       lvls=np.linspace(0.00010*1e6,0.00015*1e6,num=5,endpoint=True),frmt='%3.0f')
    
    axs[1,2].set_title("Ammonia Abundance Index (ppm)",fontsize=8)

    if R["RGB"]["file"] != 'NA':
        if imagetype=='Map':
            temp=RL.make_contours_CH4_patch(axs[0,2],fNH3_patch_mb,LatLims,NH3LonLims,
                                   lvls=tx_fNH3,frmt='%3.0f',clr='k')
        #else:
            #temp=RL.make_contours_CH4(axs[0,2],np.flip(convolve(fNH3,kernel),axis=0),
            #                       lvls=np.linspace(0.00010,0.00015,num=5,endpoint=True))

    axs[0,0].set_ylabel("Planetographic Latitude (deg)",fontsize=8)
    axs[1,0].set_ylabel("Planetographic Latitude (deg)",fontsize=8)
    axs[1,0].set_xlabel("Sys. "+LonSys+" Longitude (deg)",fontsize=8)
    axs[1,1].set_xlabel("Sys. "+LonSys+" Longitude (deg)",fontsize=8)
    axs[1,2].set_xlabel("Sys. "+LonSys+" Longitude (deg)",fontsize=8)

    fig.subplots_adjust(left=0.10, bottom=0.08, right=0.94, top=0.90,
                wspace=0.25, hspace=0.05)     
    pathout="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Analysis Data/Map Plots Diagnostic/"
    fig.savefig(pathout+obsdate+"-Jupiter-Retrieval"+"-CMII_"+
              str(R["NH3"]["CM"])+"-"+CalModel+"-"+smthtitle+"-Map.png",dpi=300)

    #nh3abs16bit = np.nan_to_num(((5.*65535.*(normnh3 - 0.9))*mask[:,:,1]).astype(np.uint16))
    fNH3scaled=np.nan_to_num(((5000.*65535.*R["NH3"]["fNH3"] - 0.9)))
    fNH3scaled[fNH3scaled<=0.]=0.0
    fNH3abs16bit = fNH3scaled.astype(np.uint16)
    fnout=R["pathfNH3"]+R["NH3"]["file"][0:26]+'fNH3_Sys'+LonSys+'.png'
    imwrite(fnout, fNH3abs16bit)#.astype(np.uint16))

    ###########################################################################
    ## Just RGB and Abundance
    ###########################################################################
    
    fig1,axs1=pl.subplots(1,2,figsize=(8.0,4.0), dpi=150, facecolor="white",
                        sharey=True,sharex=True)
    #fig1.suptitle(NH3time.replace("_"," ")+", CM2="+str(int(CM2)),x=0.5,ha='center',color='k')
    fig1.suptitle(R["NH3"]["time"].replace("_"," ")+", CM"+LonSys+"="
                  +str(int(R["NH3"]["CM"]))+", Calibration = "+CalModel+", Data "+smthtitle,x=0.5,ha='center',color='k')

    for ix in range(0,1):
        axs1[ix].grid(linewidth=0.2)
        axs1[ix].ylim=[-45.,45.]
        axs1[ix].xlim=[360-NH3LonLims[0],360-NH3LonLims[1]]
        axs1[ix].set_xticks(np.linspace(450,0,31), minor=False)
        xticklabels=np.array(np.mod(np.linspace(450,0,31),360))
        axs1[ix].set_xticklabels(xticklabels.astype(int))
        axs1[ix].set_yticks(np.linspace(-45,45,7), minor=False)
        axs1[ix].tick_params(axis='both', which='major', labelsize=9)

        axs1[ix].set_adjustable('box') 


    show=axs1[0].imshow(fNH3_patch_mb, "jet", origin='upper',vmin=vn_fNH3,vmax=vx_fNH3,  
               extent=[360-NH3LonLims[0],360-NH3LonLims[1],90-LatLims[1],
                       90-LatLims[0]],#vmin=0,vmax=1.
                       aspect="equal")
    temp=RL.make_contours_CH4_patch(axs1[0],fNH3_patch_mb,LatLims,NH3LonLims,
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
               extent=[360-NH3LonLims[0],360-NH3LonLims[1],90-LatLims[1],
                       90-LatLims[0]],
                       aspect="equal")
    temp=RL.make_contours_CH4_patch(axs1[1],fNH3_patch_mb,LatLims,NH3LonLims,
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
            axs1[0].fill_between([360-NH3LonLims[0],360-NH3LonLims[1]],[belt[zb][0],belt[zb][0]],[belt[zb][1],belt[zb][1]],
                                    color="0.5",alpha=0.4)
            axs1[1].fill_between([360-NH3LonLims[0],360-NH3LonLims[1]],[belt[zb][0],belt[zb][0]],[belt[zb][1],belt[zb][1]],
                                    color="0.5",alpha=0.4)
        #axs1[1].annotate(zb,xy=[np.mean(belt[zb]),51],ha="center")
    #for zb in zone:
        #axs1[1].annotate(zb,xy=[np.mean(zone[zb]),51],ha="center")

    axs1[1].tick_params(axis='both', which='major', labelsize=9)
    axs1[1].set_title("RGB Context Image",fontsize=10)

    axs1[0].set_ylabel("Planetographic Latitude (deg)",fontsize=10)
    axs1[0].set_xlabel("Sys. "+LonSys+" Longitude (deg)",fontsize=10)
    axs1[1].set_xlabel("Sys. "+LonSys+" Longitude (deg)",fontsize=10)

    fig1.subplots_adjust(left=0.10, bottom=0.03, right=0.98, top=0.95,
                wspace=0.25, hspace=0.05)     
    axs1[1].set_position([box.x0+0.03, box.y0-0.01, box.width * 1.015, box.height * 1.015])
    pathmapplots='C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Analysis Data/map plots/'

    fig1.savefig(pathmapplots+obsdate+"-Jupiter-Retrieval-NH3_RGB_only"+"-CMII_"+
              str(R["NH3"]["CM"])+"-"+CalModel+"-"+smthtitle+"-Map.png",dpi=300)
    
    ###########################################################################
    ## Just RGB and Cloud Pressure
    ###########################################################################
    #!!!!!!!!!! FIX THIS !!!!!!!!!!!!!
    fig2,axs2=pl.subplots(1,2,figsize=(8.0,4.0), dpi=150, facecolor="white",
                        sharey=True,sharex=True)
    fig2.suptitle(R["NH3"]["time"].replace("_"," ")+", CM"+LonSys+"="
                  +str(int(R["NH3"]["CM"]))+", Calibration = "+CalModel+", Data "+smthtitle,x=0.5,ha='center',color='k')

    for ix in range(0,1):
        axs2[ix].grid(linewidth=0.2)
        axs2[ix].ylim=[-45.,45.]
        axs2[ix].xlim=[360-NH3LonLims[0],360-NH3LonLims[1]]
        axs2[ix].set_xticks(np.linspace(450,0,31), minor=False)
        xticklabels=np.array(np.mod(np.linspace(450,0,31),360))
        axs2[ix].set_xticklabels(xticklabels.astype(int))
        axs2[ix].set_yticks(np.linspace(-45,45,7), minor=False)
        axs2[ix].tick_params(axis='both', which='major', labelsize=9)

        axs2[ix].set_adjustable('box') 

    Pcloud_patch,vn,vx,tx=plot_patch(R["CH4"]["PCloud"]/4.0,LatLims,NH3LonLims,
                                     R["CH4"]["CM"],LonRng,"jet",
                                     axs2[0],'%3.2f',cont=False,
                                     cbar_reverse=True,vn=400,vx=900,n=6)

    temp=RL.make_contours_CH4_patch(axs2[0],Pcloud_patch,LatLims,NH3LonLims,
                           lvls=tx,frmt='%3.2f',clr='k')
    #temp=RL.make_contours_CH4_patch(axs2[0],fNH3_patch_mb,LatLims,NH3LonLims,
    #                       lvls=tx_fNH3,frmt='%3.0f',clr='r')
    
    axs2[0].set_title("Cloud Top Pressure (bars)",fontsize=10)

    gamma=1.3
    RGB4Display=np.power(np.array(RGB_patch).astype(float),gamma)
    RGB4Display=RGB4Display/RGB4Display.max()
    show=axs2[1].imshow(RGB4Display,
               extent=[360-NH3LonLims[0],360-NH3LonLims[1],90-LatLims[1],
                       90-LatLims[0]],
                       aspect="equal")
    temp=RL.make_contours_CH4_patch(axs2[1],Pcloud_patch,LatLims,NH3LonLims,
                           tx,frmt='%3.2f',clr='k')
    #temp=RL.make_contours_CH4_patch(axs2[1],fNH3_patch_mb,LatLims,NH3LonLims,
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
            axs2[0].fill_between([360-NH3LonLims[0],360-NH3LonLims[1]],[belt[zb][0],belt[zb][0]],[belt[zb][1],belt[zb][1]],
                                    color="0.5",alpha=0.4)
            axs2[1].fill_between([360-NH3LonLims[0],360-NH3LonLims[1]],[belt[zb][0],belt[zb][0]],[belt[zb][1],belt[zb][1]],
                                    color="0.5",alpha=0.4)
        #axs1[1].annotate(zb,xy=[np.mean(belt[zb]),51],ha="center")
    #for zb in zone:
        #axs1[1].annotate(zb,xy=[np.mean(zone[zb]),51],ha="center")

    axs2[1].tick_params(axis='both', which='major', labelsize=9)
    axs2[1].set_title("RGB Context Image",fontsize=10)

    axs2[0].set_ylabel("Planetographic Latitude (deg)",fontsize=10)
    axs2[0].set_xlabel("Sys. "+LonSys+" Longitude (deg)",fontsize=10)
    axs2[1].set_xlabel("Sys. "+LonSys+" Longitude (deg)",fontsize=10)

    fig2.subplots_adjust(left=0.10, bottom=0.03, right=0.98, top=0.95,
                wspace=0.25, hspace=0.05)     
    axs2[1].set_position([box.x0+0.03, box.y0-0.01, box.width * 1.015, box.height * 1.015])
    pathmapplots='C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Analysis Data/map plots/'

    fig2.savefig(pathmapplots+obsdate+"-Jupiter-Retrieval-Pcloud_only"+"-CMII_"+
              str(R["CH4"]["CM"])+"-"+CalModel+"-"+smthtitle+"-Map.png",dpi=300)


    pathScatter='C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Analysis Data/Scatter Plots/'

    plot_scatter(Pcloud_patch,fNH3_patch_mb,pathScatter,obsdate,R["NH3"]["CM"],
                 LatLims)

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
    #print("####################### Patch RGB shape",patch.shape)
    if CM2deg<LonRng:
        patch=np.concatenate((np.copy(Map[LatLims[0]:LatLims[1],LonLims[0]-1:360,:]),
                              np.copy(Map[LatLims[0]:LatLims[1],0:LonLims[1]-360,:])),axis=1)
    if CM2deg>360.-LonRng:
        patch=np.concatenate((np.copy(Map[LatLims[0]:LatLims[1],360+LonLims[0]:360,:]),
                              np.copy(Map[LatLims[0]:LatLims[1],0:LonLims[1],:])),axis=1)
    #print("####################### Patch RGB shape",patch.shape)
    #if pad:
    #    patch_pad=np.pad(patch,5,mode='reflect')

    return patch

def plot_patch(fullmap,LatLims,LonLims,CM2,LonRng,colorscale,axis,frmt,
               cont=True,cbar_reverse=False,vn=0.10,vx=0.20,n=6):
    import numpy as np
    import pylab as pl
    import RetrievalLibrary as RL

    patch=RL.make_patch(fullmap,LatLims,LonLims,CM2,LonRng)
    np.nan_to_num(patch, copy=False, nan=0.0, posinf=0.0, neginf=0.0)
    #vn=np.mean(patch)-3.0*np.std(patch)
    #vx=np.mean(patch)+3.0*np.std(patch)
    tx=np.linspace(vn,vx,n,endpoint=True)
    
    print(np.mean(patch),vn,vx)

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
    
    print("000000000: ",patch1.shape,patch2.shape)
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
    #patch1=patch1*1000.

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
  

