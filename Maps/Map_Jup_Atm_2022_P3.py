def Map_Jup_Atm_2022_P3(obskey="20231103UTa",imagetype='Map',
                        Smoothing=False,LatLims=[45,135],LonRng=45,
                        CMpref='subobs',LonSys='2',showbands=False,
                        coef=[0.,0.],subproj=''):
    """
    Created on Sun Nov  6 16:47:21 2022
    
    PURPOSE: Create maps of environmental parameters paired with RGB context
             maps. Based on Retrieve_Jup_Atm_2022_P3, which ALSO performed
             the calibration phase. So now I've separated that module into 
             a calibration module, make_L3_env_data.py and this plotting
             module.
    
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
    
    target="Jupiter"
    PCloudhdr,PClouddata,fNH3hdr,fNH3data,sza,eza,RGB,RGB_CM= \
                    RFM.read_fits_map_L2_L3(obskey=obskey,LonSys=LonSys,
                                            imagetype="Map",Level="L3")
                    
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
        PCldPlotCM=PCloudhdr["CM"+LonSys]
    else:
        fNH3PlotCM=CMpref
        PCldPlotCM=CMpref
    NH3LonLims=[360-int(fNH3PlotCM+LonRng),360-int(fNH3PlotCM-LonRng)]
    if Smoothing:
        smthtitle="Smoothed"
    else: 
        smthtitle="Unsmoothed"
    CalModel=fNH3hdr['CALIBRA']

    ###########################################################################
    ## Just RGB and Abundance
    ###########################################################################
    fig1,axs1=pl.subplots(1,2,figsize=(8.0,4.0), dpi=150, facecolor="white",
                        sharey=True,sharex=True)
    fig1.suptitle(fNH3hdr["DATE-OBS"].replace("_"," ")+", CM"+LonSys+"="
                  +str(int(fNH3PlotCM))+", Calibration = "+CalModel+", Data "+smthtitle,x=0.5,ha='center',color='k')

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

    TestfNH3=fNH3data*amfdata**coef[0]
    
    fNH3_patch_mb,vn,vx,tx_fNH3=plot_patch(TestfNH3*1e6,LatLims,NH3LonLims,
                                     fNH3PlotCM,LonRng,"jet",
                                     axs1[0],'%3.2f',cont=False,n=6,vn=50,vx=150)

    temp=RL.make_contours_CH4_patch(axs1[0],fNH3_patch_mb,LatLims,NH3LonLims,
                           lvls=tx_fNH3,frmt='%3.0f',clr='k')
    
    axs1[0].set_title(r'$\bar{f_c}(NH3) (ppm)$',fontsize=10)

    gamma=1.3

    #Logic in RGB_patch depends on LonLims and CM being consistent
    RGB_patch=make_patch_RGB(RGB,LatLims,NH3LonLims,fNH3PlotCM,LonRng)
    
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
            #print(zb,belt[zb])
            axs1[0].fill_between([360-NH3LonLims[0],360-NH3LonLims[1]],[belt[zb][0],belt[zb][0]],[belt[zb][1],belt[zb][1]],
                                    color="0.5",alpha=0.2)
            axs1[1].fill_between([360-NH3LonLims[0],360-NH3LonLims[1]],[belt[zb][0],belt[zb][0]],[belt[zb][1],belt[zb][1]],
                                    color="0.8",alpha=0.1)
        #axs1[1].annotate(zb,xy=[np.mean(belt[zb]),51],ha="center")
    #for zb in zone:
        #axs1[1].annotate(zb,xy=[np.mean(zone[zb]),51],ha="center")

    axs1[1].tick_params(axis='both', which='major', labelsize=9)
    axs1[1].set_title("RGB Context Image",fontsize=10)

    axs1[0].set_ylabel("Planetographic Latitude (deg)",fontsize=10)
    axs1[0].set_xlabel("Sys. "+LonSys+" Longitude (deg)",fontsize=10)
    axs1[1].set_xlabel("Sys. "+LonSys+" Longitude (deg)",fontsize=10)
    axs1[1].grid(linewidth=0.2)

    fig1.subplots_adjust(left=0.10, bottom=0.03, right=0.98, top=0.95,
                wspace=0.25, hspace=0.05)     
    axs1[1].set_position([box.x0+0.03, box.y0-0.01, box.width * 1.015, box.height * 1.015])
    """
    try:
        fout=pathmapplots+obskey+"-Jupiter-Retrieval-NH3_RGB_only"+"-CMII_"+\
                  str(fNH3PlotCM)+"-"+CalModel+"-"+smthtitle+"-Map-"#+\
                  R["Variation"]+".png"
    except:
        fout=pathmapplots+obskey+"-Jupiter-Retrieval-NH3_RGB_only"+"-CMII_"+\
                  str(fNH3PlotCM)+"-"+CalModel+"-"+smthtitle+"-Map.png"
    """
    if coef[0]==0.0:
        correction='_C0'
    else:
        correction='_C1'
    fn=fNH3hdr["FILENAME"][:-5]+correction+'_Sys'+LonSys+'_N'+\
                str(90-LatLims[0])+'-S'+str(LatLims[1]-90)+\
                '_Lon'+str(NH3LonLims[0]).zfill(3)+'-'+str(NH3LonLims[1]).zfill(3)+'.png'
    fig1.savefig(pathmapplots+fn,dpi=300)
    
    ###########################################################################
    ## Just RGB and Cloud Pressure
    ###########################################################################
    #!!!!!!!!!! FIX THIS !!!!!!!!!!!!!
    fig2,axs2=pl.subplots(1,2,figsize=(8.0,4.0), dpi=150, facecolor="white",
                        sharey=True,sharex=True)
    fig2.suptitle(fNH3hdr["DATE-OBS"].replace("_"," ")+", CM"+LonSys+"="
                  +str(int(PCldPlotCM))+", Calibration = "+CalModel+", Data "+smthtitle,x=0.5,ha='center',color='k')

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

    #Pcloud_patch,vn,vx,tx=plot_patch(PClouddata,LatLims,NH3LonLims,
    #                                 PCldPlotCM,LonRng,"jet",
    #                                 axs2[0],'%3.2f',cont=False,
    #                                 cbar_reverse=True,vn=400,vx=900,n=6)
    TestPCloud=PClouddata*amfdata**coef[1]
    Pcloud_patch,vn,vx,tx=plot_patch(TestPCloud,LatLims,NH3LonLims,
                                     PCldPlotCM,LonRng,"jet",
                                     axs2[0],'%3.2f',cont=False,
                                     cbar_reverse=True,vn=400,vx=900,n=6)

    """##########TEST CODE
    hdu = fits.PrimaryHDU((R["CH4"]["PCloud"]).astype(np.float32))
    hdul = fits.HDUList([hdu])
    try:
        os.remove(R["pathFITS"]+'PCloud.fits')
    except: 
        print("file doesn't exist")
    hdul.writeto(R["pathFITS"]+'PCloud.fits')
    hdul.close()
    ##########TEST CODE"""

    temp=RL.make_contours_CH4_patch(axs2[0],Pcloud_patch,LatLims,NH3LonLims,
                           lvls=tx,frmt='%3.2f',clr='k')
    #temp=RL.make_contours_CH4_patch(axs2[0],fNH3_patch_mb,LatLims,NH3LonLims,
    #                       lvls=tx_fNH3,frmt='%3.0f',clr='r')
    
    axs2[0].set_title("Cloud Top Pressure (mb)",fontsize=10)

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
            #print(zb,belt[zb])
            axs2[0].fill_between([360-NH3LonLims[0],360-NH3LonLims[1]],[belt[zb][0],belt[zb][0]],[belt[zb][1],belt[zb][1]],
                                    color="0.5",alpha=0.2)
            axs2[1].fill_between([360-NH3LonLims[0],360-NH3LonLims[1]],[belt[zb][0],belt[zb][0]],[belt[zb][1],belt[zb][1]],
                                    color="0.8",alpha=0.1)
        #axs1[1].annotate(zb,xy=[np.mean(belt[zb]),51],ha="center")
    #for zb in zone:
        #axs1[1].annotate(zb,xy=[np.mean(zone[zb]),51],ha="center")

    axs2[1].tick_params(axis='both', which='major', labelsize=9)
    axs2[1].set_title("RGB Context Image",fontsize=10)

    axs2[0].set_ylabel("Planetographic Latitude (deg)",fontsize=10)
    axs2[0].set_xlabel("Sys. "+LonSys+" Longitude (deg)",fontsize=10)
    axs2[1].set_xlabel("Sys. "+LonSys+" Longitude (deg)",fontsize=10)
    axs2[1].grid(linewidth=0.2)

    fig2.subplots_adjust(left=0.10, bottom=0.03, right=0.98, top=0.95,
                wspace=0.25, hspace=0.05)     
    axs2[1].set_position([box.x0+0.03, box.y0-0.01, box.width * 1.015, box.height * 1.015])
    
    if coef[0]==0.0:
        correction='_C0'
    else:
        correction='_C1'
    fn=PCloudhdr["FILENAME"][:-5]+correction+'_Sys'+LonSys+'_N'+\
                str(90-LatLims[0])+'-S'+str(LatLims[1]-90)+\
                '_Lon'+str(NH3LonLims[0]).zfill(3)+'-'+str(NH3LonLims[1]).zfill(3)+'.png'

    fig2.savefig(pathmapplots+fn,dpi=300)

    fn=fNH3hdr["FILENAME"][:-12]+'Scat_'+fNH3hdr["FILENAME"][-7:-5]+\
                correction+'_Sys'+LonSys+'_N'+\
                str(90-LatLims[0])+'-S'+str(LatLims[1]-90)+\
                '_Lon'+str(NH3LonLims[0]).zfill(3)+'-'+str(NH3LonLims[1]).zfill(3)+'.png'
    plot_scatter(Pcloud_patch,fNH3_patch_mb,pathmapplots,fn,obskey,fNH3PlotCM,
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
    if CM2deg<LonRng:
        patch=np.concatenate((np.copy(Map[LatLims[0]:LatLims[1],LonLims[0]-1:360,:]),
                              np.copy(Map[LatLims[0]:LatLims[1],0:LonLims[1]-360,:])),axis=1)
    if CM2deg>360-LonRng:
        patch=np.concatenate((np.copy(Map[LatLims[0]:LatLims[1],360+LonLims[0]:360,:]),
                              np.copy(Map[LatLims[0]:LatLims[1],0:LonLims[1],:])),axis=1)
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
    
    #print(np.mean(patch),vn,vx)

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

def plot_scatter(patch1,patch2,filepath,fn,obskey,Real_CM2,LatLims):
    import pylab as pl
    import numpy as np
    import copy
    
    #print("000000000: ",patch1.shape,patch2.shape)
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

    figcor,axscor=pl.subplots(1,1,figsize=(6.0,4.), dpi=150, facecolor="white",
                        sharey=True,sharex=True)          

    for key in BZ.keys():
        #print(key,BZ[key],[90,90]-np.array(BZ[key]),LatLims)
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
                           marker="o",s=2.0,
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
    axscor.set_xlim(700.,1100.)
    axscor.set_ylim(70.,160)
    #axscor.set_xticks(np.linspace(-45.,45.,7), minor=False)
    #axscor.set_yticks(np.linspace(0.0,1.0,6), minor=False)
    #axscor.tick_params(axis='both', which='major', labelsize=8)
    axscor.set_xlabel("Effective Cloud-top Pressure (mb)",fontsize=10)
    axscor.set_ylabel("Ammonia Mole Fraction (ppm)",fontsize=10)
    axscor.set_title("Distribution of Ammonia vs Cloud-Top Pressure")
    axscor.legend(fontsize=9)
    
    figcor.subplots_adjust(left=0.12, right=0.97, top=0.93, bottom=0.11)
    figcor.savefig(filepath+fn,dpi=300)
  

