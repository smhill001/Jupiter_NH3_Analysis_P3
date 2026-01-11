def map_and_scatter(patchx,patchy,mapydata,mapyhdr,LonSys,
                    LatLims,LonLims,LonRng,PlotCM,fnout,
                    coef,txin,xlow,xhigh,ylow,yhigh,figxy,
                    ct,pathout,Ltitle,Rtitle,Level='L3',FiveMicron=False,
                    cbar_rev=False,swap_xy=False,axis_inv=False,cbar_title="Test",
                    suptitle="Test",ROI=False):
    """
    PURPOSE:    Makes a pair of plots, the left one is a patch map of one data
                set overlayed by another patch map data set. The right plot is
                a scatter plot of the two patch maps.
                
    CALLS:      plot_map_scatter
                plot_roi_scatter

    Parameters
    ----------
    patchx : TYPE
        DESCRIPTION.
    patchy : TYPE
        DESCRIPTION.
    mapydata : TYPE
        DESCRIPTION.
    mapyhdr : TYPE
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
    fnout : TYPE
        DESCRIPTION.
    coef : TYPE
        DESCRIPTION.
    txin : TYPE
        DESCRIPTION.
    xlow : TYPE
        DESCRIPTION.
    xhigh : TYPE
        DESCRIPTION.
    ylow : TYPE
        DESCRIPTION.
    yhigh : TYPE
        DESCRIPTION.
    figxy : TYPE
        DESCRIPTION.
    ct : TYPE
        DESCRIPTION.
    pathout : TYPE
        DESCRIPTION.
    Ltitle : TYPE
        DESCRIPTION.
    Rtitle : TYPE
        DESCRIPTION.
    Level : TYPE, optional
        DESCRIPTION. The default is 'L3'.
    FiveMicron : TYPE, optional
        DESCRIPTION. The default is False.
    cbar_rev : TYPE, optional
        DESCRIPTION. The default is False.
    swap_xy : TYPE, optional
        DESCRIPTION. The default is False.
    axis_inv : TYPE, optional
        DESCRIPTION. The default is False.
    cbar_title : TYPE, optional
        DESCRIPTION. The default is "Test".
    suptitle : TYPE, optional
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
    import make_patch as MP
    import plot_contours_on_patch as PC
    import copy
    import plot_map_scatter as pms
    import plot_roi_scatter as prs

    ###########################################################################
    ## Compute Scatter Plot (PCloud vs fNH3)
    ###########################################################################
    fig3,axs3=pl.subplots(1,2,figsize=(figxy[0],figxy[1]), dpi=150, facecolor="white")
    #fig3.suptitle(suptitle,x=0.5,ha='center',color='k')
    fig3.suptitle(mapyhdr["DATE-OBS"].replace("_"," ")+", CM"+LonSys+"="
                  +str(int(PlotCM)),x=0.5,ha='center',color='k')

    axs3[0].grid(linewidth=0.2)
    axs3[0].ylim=[-45.,45.]
    axs3[0].xlim=[360-LonLims[0],360-LonLims[1]]
    axs3[0].set_xticks(np.linspace(450,0,31), minor=False)
    xticklabels=np.array(np.mod(np.linspace(450,0,31),360))
    axs3[0].set_xticklabels(xticklabels.astype(int))
    axs3[0].set_yticks(np.linspace(-45,45,7), minor=False)
    axs3[0].tick_params(axis='both', which='major', labelsize=9)
    axs3[0].set_ylabel("Planetographic Latitude (deg)",fontsize=10)
    axs3[0].set_xlabel("Sys. "+LonSys+" Longitude (deg)",fontsize=10)
    axs3[0].set_title(Ltitle,fontsize=10)

    axs3[0].set_adjustable('box') 
    axs3[1].set_adjustable('box') 

    #Pcloud_patch,vn,vx,tx=PP.plot_patch(PClouddata,LatLims,NH3LonLims,
    #                                 PCldPlotCM,LonRng,"jet",
    #                                 axs2[0],'%3.2f',cont=False,
    #                                 cbar_reverse=True,vn=400,vx=900,n=6)
    Testy_patch=MP.make_patch(mapydata,LatLims,LonLims,PlotCM,LonRng)
    Testy_patch,vn,vx,tx=PP.plot_patch(Testy_patch,LatLims,LonLims,
                                     PlotCM,LonRng,ct,
                                     axs3[0],'%3.2f',
                                     cbar_reverse=cbar_rev,vn=ylow,vx=yhigh,n=6,
                                     cbar_title=cbar_title)
    temp=PC.plot_contours_on_patch(axs3[0],patchx,LatLims,LonLims,
                           txin,frmt='%3.0f',clr='k')

    if coef==0.0:
        correction='_C0'
    else:
        correction='_C1'
    
    if Level=='L2':
        fnout=fnout.replace('TNH3',Rtitle)
    elif Level=='L3':
        fnout=fnout.replace('fNH3',Rtitle)
    print(patchx.shape,patchy.shape)
    
    if swap_xy and not ROI:
        roilabel,mean1,stdv1,mean2,stdv2,BZ=pms.plot_map_scatter(patchx,patchy,PlotCM,
                 LatLims,axs3[1],xlow,xhigh,ylow,yhigh,FiveMicron,axis_inv=axis_inv)
    if not swap_xy and not ROI:     
        roilabel,mean1,stdv1,mean2,stdv2,BZ=pms.plot_map_scatter(patchy,patchx,PlotCM,
                 LatLims,axs3[1],ylow,yhigh,xlow,xhigh,FiveMicron,axis_inv=axis_inv)
        
    if swap_xy and ROI:
        print("Calling ROI")
        roilabel,mean1,stdv1,mean2,stdv2=prs.plot_roi_scatter(patchx,patchy,PlotCM,
                 LatLims,LonLims,axs3[1],xlow,xhigh,ylow,yhigh,FiveMicron,
                 axis_inv=axis_inv,ROI=ROI,amfpatch=amfpatch)
    if not swap_xy and ROI:    
        print("Calling ROI")
        roilabel,mean1,stdv1,mean2,stdv2,meanamf=prs.plot_roi_scatter(patchy,patchx,PlotCM,
                 LatLims,LonLims,axs3[1],ylow,yhigh,xlow,xhigh,FiveMicron,
                 axis_inv=axis_inv,ROI=ROI,amfpatch=amfpatch)
        
        
    axs3[1].tick_params(axis='both', which='major', labelsize=9)
    axs3[1].set_title(Rtitle,fontsize=10)
    
    if "PCloud" in Rtitle and "5um" in Rtitle:
        axs3[1].set_ylabel("Cloud Top Pressure (mb)",fontsize=10)        
    if "fNH3" in Rtitle and "5um" in Rtitle:
        axs3[1].set_ylabel("Ammonia Mole Fraction (ppm)",fontsize=10)
    if "Methane" and "Ammonia" in Rtitle:
        axs3[1].set_xlabel("Ammonia Opacity")
        axs3[1].set_ylabel("Methane Opacity")
    if "Methane" in Rtitle and "5um" in Rtitle:
        axs3[1].set_ylabel("Methane Opacity",fontsize=10)        
    if "Ammonia" in Rtitle and "5um" in Rtitle:
        axs3[1].set_ylabel("Ammonia Opacity",fontsize=10)        
        axs3[1].set_xlabel("5um Radiance (Log10(arb. units)",fontsize=10)        
        
    if ROI:
        for R in ROI:
            axs3[0].plot(np.array([ROI[R][2]+ROI[R][3],ROI[R][2]-ROI[R][3],
                          ROI[R][2]-ROI[R][3],ROI[R][2]+ROI[R][3],
                          ROI[R][2]+ROI[R][3]]),
                          90.-np.array([ROI[R][0],ROI[R][0],ROI[R][1],
                          ROI[R][1],ROI[R][0]]))
    else:
        BZind=copy.deepcopy(BZ)   
        BZkeys=BZ.keys()
        BZind=copy.deepcopy(BZ)   
        BZkeys=BZ.keys()
        #patch1=patch1*1000.
    
        #figcor,axscor=pl.subplots(1,1,figsize=(6.0,4.), dpi=150, facecolor="white",
        #                    sharey=True,sharex=True)          
    
        clrind=0
        for key in BZ.keys():
            #print(key,BZ[key],[90,90]-np.array(BZ[key]),LatLims)
            #######################################################################
            # Compute the axis fraction and plot vertical bars using axvspan
            BZind[key][0]=1.-((90-BZ[key][0])-LatLims[0])/(LatLims[1]-LatLims[0])
            BZind[key][1]=1.-((90-BZ[key][1])-LatLims[0])/(LatLims[1]-LatLims[0])
            #print(key,BZind[key])
                
            #print(key,BZ[key],[90,90]-np.array(BZ[key]),LatLims)
            clr='C'+str(clrind)
            #print(clr)
            if BZind[key][0]>1.0 or BZind[key][1]<0.0:
                print("do nothing")
            else:
                axs3[0].axvspan(360-LonLims[0],360-LonLims[0]-1,
                                ymin=BZind[key][1],ymax=BZind[key][0],alpha=1.0,
                                color=clr)
                axs3[0].axvspan(360-LonLims[1],360-LonLims[1]+1,
                                ymin=BZind[key][1],ymax=BZind[key][0],alpha=1.0,
                                color=clr)
                
                clrind=clrind+1

    box = axs3[1].get_position()
    axs3[1].set_position([box.x0+0.03, box.y0-0.01, box.width * 0.5, box.height * 1.015])    
    fig3.subplots_adjust(left=0.12, right=0.97, top=0.83, bottom=0.18, 
                         wspace=0.4)
    fig3.savefig(pathout+fnout,dpi=300)
    
    dateobs=mapyhdr["DATE-OBS"]
    if not ROI:
        meanamf=False

    return(dateobs,roilabel,mean1,stdv1,mean2,stdv2,meanamf)