def map_and_context(mapdata,maphdr,RGB,RGBtime,LonSys,LatLims,LonLims,LonRng,PlotCM,
                    amfdata,coef,low,high,showbands,FiveMicron,figxy,ct,pathout,
                    Level='L3',suptitle="Test",cbar_rev=False,cbar_title="Test",
                    ROI=False):
    """
    PURPOSE:    To create a pair of plots, the right one representing a mapped
                data set, e.g., ammonia abundance or cloud pressure, and the 
                right one showing an RGB or IGB context image.

    Parameters
    ----------
    mapdata : TYPE
        DESCRIPTION.
    maphdr : TYPE
        DESCRIPTION.
    RGB : TYPE
        DESCRIPTION.
    RGBtime : TYPE
        DESCRIPTION.
    LonSys : TYPE
        DESCRIPTION.
    LatLims : TYPE
        DESCRIPTION.
    LonLims : TYPE
        DESCRIPTION.
    LonRng : TYPE
        DESCRIPTION.
    PlotCM : TYPE
        DESCRIPTION.
    amfdata : TYPE
        DESCRIPTION.
    coef : TYPE
        DESCRIPTION.
    low : TYPE
        DESCRIPTION.
    high : TYPE
        DESCRIPTION.
    showbands : TYPE
        DESCRIPTION.
    FiveMicron : TYPE
        DESCRIPTION.
    figxy : TYPE
        DESCRIPTION.
    ct : TYPE
        DESCRIPTION.
    pathout : TYPE
        DESCRIPTION.
    Level : TYPE, optional
        DESCRIPTION. The default is 'L3'.
    suptitle : TYPE, optional
        DESCRIPTION. The default is "Test".
    cbar_rev : TYPE, optional
        DESCRIPTION. The default is False.
    cbar_title : TYPE, optional
        DESCRIPTION. The default is "Test".
    ROI : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    None.

    """
    import sys
    drive='c:'
    sys.path.append(drive+'/Astronomy/Python Play')
    sys.path.append(drive+'/Astronomy/Python Play/Util_P3')
    sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy')
    sys.path.append('./Services')

    import pylab as pl
    import numpy as np
    sys.path.append('./Maps')
    import plot_patch as PP
    import make_patch_RGB as MPRGB
    import make_patch as MP
    import plot_contours_on_patch as PC

    ###########################################################################
    ## Just RGB and Abundance
    ###########################################################################
    #ct="jet"
    context_title=RGBtime
    
    fig1,axs1=pl.subplots(1,2,figsize=(figxy[0],figxy[1]), dpi=150, facecolor="white",
                        sharey=True,sharex=True)
    fig1.suptitle(suptitle,x=0.5,ha='center',color='k')

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

    data_patch=MP.make_patch(mapdata,LatLims,LonLims,PlotCM,LonRng)
    data_patch,vn,vx,tx=PP.plot_patch(data_patch,LatLims,LonLims,
                                     PlotCM,LonRng,ct,axs1[0],'%3.2f',
                                     n=6,vn=low,vx=high,
                                     cbar_title=cbar_title,cbar_reverse=cbar_rev)
    
    temp=PC.plot_contours_on_patch(axs1[0],data_patch,LatLims,LonLims,
                           lvls=tx,frmt='%3.0f',clr='k')
        
    if maphdr["TELESCOP"]=="NASA IRTF":
        axs1[0].set_title("IRTF 5um Radiance (Log10(arb. units))",fontsize=10)
    else:
        if maphdr["BUNIT"]=="Cloud-top Press":
            axs1[0].set_title("Cloud Top Pressure (mb)",fontsize=10)        
        if maphdr["BUNIT"]=="Mole Fraction":
            axs1[0].set_title(r'$\bar{f_c}(NH3) (ppm)$',fontsize=10)

    axs1[0].set_title(maphdr["DATE-OBS"],fontsize=10)

    gamma=1.3

    #Logic in RGB_patch depends on LonLims and CM being consistent
    #RGB_patch=MPRGB.make_patch_RGB(RGB,LatLims,LonLims,PlotCM,LonRng)
    RGB_patch=MP.make_patch(RGB,LatLims,LonLims,PlotCM,LonRng)
    
    RGB4Display=np.power(np.array(RGB_patch).astype(float),gamma)
    RGB4Display=RGB4Display/RGB4Display.max()
    show=axs1[1].imshow(RGB4Display,
               extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],
                       90-LatLims[0]],
                       aspect="equal")
    temp=PC.plot_contours_on_patch(axs1[1],data_patch,LatLims,LonLims,
                           tx,frmt='%3.0f',clr='k')

    if ROI:
        for R in ROI:
            for iax in range(0,2):
                axs1[iax].plot(np.array([ROI[R][2]+ROI[R][3],ROI[R][2]-ROI[R][3],
                              ROI[R][2]-ROI[R][3],ROI[R][2]+ROI[R][3],
                              ROI[R][2]+ROI[R][3]]),
                              90.-np.array([ROI[R][0],ROI[R][0],ROI[R][1],
                              ROI[R][1],ROI[R][0]]))

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
            axs1[0].fill_between([360-LonLims[0],360-LonLims[1]],[belt[zb][0],belt[zb][0]],[belt[zb][1],belt[zb][1]],
                                    color="0.5",alpha=0.2)
            axs1[1].fill_between([360-LonLims[0],360-LonLims[1]],[belt[zb][0],belt[zb][0]],[belt[zb][1],belt[zb][1]],
                                    color="0.8",alpha=0.1)
        #axs1[1].annotate(zb,xy=[np.mean(belt[zb]),51],ha="center")
    #for zb in zone:
        #axs1[1].annotate(zb,xy=[np.mean(zone[zb]),51],ha="center")
    
    axs1[1].tick_params(axis='both', which='major', labelsize=9)
    axs1[1].set_title(context_title,fontsize=10)

    axs1[0].set_ylabel("Planetographic Latitude (deg)",fontsize=10)
    axs1[0].set_xlabel("Sys. "+LonSys+" Longitude (deg)",fontsize=10)
    axs1[1].set_xlabel("Sys. "+LonSys+" Longitude (deg)",fontsize=10)
    axs1[1].grid(linewidth=0.2)

    fig1.subplots_adjust(left=0.10, bottom=0.03, right=0.98, top=0.95,
                wspace=0.25, hspace=0.05)     
    axs1[1].set_position([box.x0+0.03, box.y0-0.01, box.width * 1.015, box.height * 1.015])

    if coef==0.0:
        correction='_C0'
    else:
        correction='_C1'
    #print("360-LonLims[0],360-LonLims[1]=",360-LonLims[0],360-LonLims[1])
    if FiveMicron:
        fnskeleton=correction+'_Sys'+LonSys+'_N'+\
                    str(90-LatLims[0])+'-S'+str(LatLims[1]-90)+\
                    '_Lon'+str(np.mod(360-LonLims[1],360)).zfill(3)+'-'+\
                        str(np.mod(360-LonLims[0],360)).zfill(3)+'_5micron.png'
    else:
        fnskeleton=correction+'_Sys'+LonSys+'_N'+\
                    str(90-LatLims[0])+'-S'+str(LatLims[1]-90)+\
                    '_Lon'+str(np.mod(360-LonLims[1],360)).zfill(3)+'-'+\
                        str(np.mod(360-LonLims[0],360)).zfill(3)+'.png'

    ###!!!! Temporary fix for the fact that L2 header FILENAME doesn't contain
    ###!!!! the file extension and L3 header FILENAME does!
    if Level=='L2':
        fnout=maphdr["FILENAME"]+fnskeleton
    elif Level=='L3':
        fnout=maphdr["FILENAME"][:-5]+fnskeleton
        
    fig1.savefig(pathout+fnout,dpi=300)
    
    return(data_patch,mapdata,tx,fnout,RGB4Display)