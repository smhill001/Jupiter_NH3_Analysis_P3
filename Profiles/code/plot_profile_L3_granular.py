# -*- coding: utf-8 -*-
"""
Created on Tue Apr  1 07:48:52 2025

@author: smhil
"""

def plot_profile_L3_granular(ax_spag,ax_eta,reference,
                             ProfileHalfWidth=45.,LatPlotLims=[45,135],
                             ZonePlotHalfWidth=45.,param="fNH3",
                             profile="Meridional",clr='C0',width=1.0,
                             style='solid',smooth=False,colat=90):

    """
    PURPOSE:
    Plot swaths of L3 environmental parameters as either meridional or zonal
        profiles.
        Plots include spaghetti plots of reference batch atasets of observing 
        sessions, e.g., CMOS 2022, along with a summary plot consisting of 
        averages and standard deviations of each batch data set.
    
    Parameters
    ----------
    ax : Axis object
        Plot axes on which to plot the profile.
    reference : String
        Index key for Batch List of FITS files
    ProfileHalfWidth : Integer, degrees, optional
        Latitude halfwidth to be averaged for a Zonal profile 
        Longitude halfwidth to be averagedfor a meridional profile. 
        The default is 45.
    LatPlotLims : List, integer, optional
        Plot limits in CO-Latitude for the X-axis if meridional plot is made.
        The default is [45,135].
    ZonePlotHalfWidth : Integer, optional
        Halfwidth of the X-axis plot limits centered on the system II 
        central meridion. The default is 45.
    band : String, optional
        Molecular absorption band to the plotted.
        The default is "CH4".
    profile : String, optional
        Zonal or Meridional Profile.
        The default is "Meridional".
    clr : String, optional
        Line color. The default is 'C0'.
    width : Float, optional
        Line width. The default is 1.0.
    style : TYPE, optional
        Line style. The default is 'solid'.
    smooth : Boolean, optional
        Smooths only the SUMMARY profile with a 3 degree boxcar kernal. 
        The default is False.

    Returns
    -------
    None.

    """
    import numpy as np
    from astropy.convolution import convolve, Box1DKernel
    import sys
    sys.path.append('./Services')
    import extract_profile as EP
    import get_obs_list as GOL
    import get_batch_lists as GBL

    ###########################################################################
    # Get sourcefile list, batch list, and set input path
    ###########################################################################
    sourcefiles=GOL.get_obs_list()
    DataSets=GBL.get_batch_lists()
    pth="C:/Astronomy/Projects/SAS 2021 Ammonia/Data/L3 FITS/"

    if param=="PCld":
        paramkey='CH4file'
        factor=1.0
    elif param=="fNH3":
        paramkey='NH3file'
        factor=1.0
    
    if profile=="Meridional":
        xlabel="Planetographic Latitude (deg)"
    elif profile=="Zonal":
        xlabel="Longitude from CM (deg)"

    suffix={"PCld":"-Jupiter_Map_L3PCld_S0.fits",
            "fNH3":"-Jupiter_Map_L3fNH3_S0.fits"}

    ###########################################################################
    # LOOP OVER DATA (OBSERVATIONS) IN DATASETS
    ###########################################################################
    First=True
    Num=0

    for obskey in DataSets[reference]:
        print("*******obskey=",obskey)
        dataset=obskey#[0:10].replace('-','')+'UT'+obskey[10]#+"_Map"
        print("******dataset=",dataset)
        print(paramkey)
        print(sourcefiles[dataset][paramkey][0:17])
        print(suffix[param])
        try:
            file=sourcefiles[dataset][paramkey][0:17]+suffix[param]+\
                    sourcefiles[dataset]['Metadata']['Variation']
            variation=sourcefiles[dataset]['Metadata']['Variation']
        except:
            file=sourcefiles[dataset][paramkey][0:17]+suffix[param]
            variation=""
        print("file=",file)
        
        #######################################################################
        # Extract a profile out of a Level 2 or 3 FITS map file
        #######################################################################
        Lats,AvgProf,StdProf,CM3,amfAvgProf,patch,amfpatch=EP.extract_profile(pth,file,
                                                ProfileHalfWidth=ProfileHalfWidth,
                                                profile=profile,colat=colat)
        if First:
            AvgArr=AvgProf
            AvgamfArr=amfAvgProf
            if profile=="Meridional":
                AvgArr=np.reshape(AvgArr,(1,180))
                AvgamfArr=np.reshape(AvgamfArr,(1,180))
            elif profile=="Zonal":
                AvgArr=np.reshape(AvgArr,(1,360))
                AvgamfArr=np.reshape(AvgamfArr,(1,360))
            First=False
        else:
            AvgArr=np.vstack((AvgArr,AvgProf))
            AvgamfArr=np.vstack((AvgamfArr,amfAvgProf))

        #######################################################################
        # Increment counter and plot a line of spaghetti
        #######################################################################
        Num=Num+1                
        ax_spag.plot(Lats,AvgProf*factor,linewidth=0.5,label=file[5:17])
        ax_eta.scatter(amfpatch,patch,s=0.3,label=file[5:17],alpha=0.5)

    AvgPro=np.mean(AvgArr[:,:],axis=0)*factor
    AvgStd=np.std(AvgArr[:,:],axis=0)*factor
    Avgamf=np.mean(AvgamfArr[:,:],axis=0)*factor
    
    if profile=="Meridional":
        ax_spag.set_xlim(90-LatPlotLims[1],90-LatPlotLims[0])
        ax_spag.set_xlabel("Planetographic Latitude (deg)",fontsize=10)
        ax_eta.set_xlim(1,3)
        ax_eta.set_xlabel("Planetographic Latitude (deg)",fontsize=10)
    elif profile=="Zonal":
        ax_spag.set_xlim(-ZonePlotHalfWidth,ZonePlotHalfWidth)
        ax_spag.set_xlabel("Longitude from CM (deg)",fontsize=10)
        ax_eta.set_xlim(1,3)
        ax_eta.set_xlabel("Longitude from CM (deg)",fontsize=10)

    if param=="PCld":
        ax_spag.set_ylim(1200,2200)
        ax_spag.set_ylabel("Cloud Pressure (mb)",fontsize=10)
        ax_spag.invert_yaxis()
        ax_eta.set_ylim(1200,2200)
        ax_eta.set_ylabel("Cloud Pressure (mb)",fontsize=10)
        ax_eta.invert_yaxis()
    elif param=="fNH3":
        ax_spag.set_ylim(0,200)
        ax_spag.set_ylabel("Ammonia Abundance (ppm)",fontsize=10)
        ax_eta.set_ylim(0,200)
        ax_eta.set_ylabel("Ammonia Abundance (ppm)",fontsize=10)
    #axsspaghetti.legend(fontsize=6,ncol=5)
    ax_spag.set_title(reference)
    ax_spag.grid(linewidth=0.2)
    ax_eta.set_title(reference)
    ax_eta.grid(linewidth=0.2)
    
    ax_spag.annotate("ProfileHalfWidth="+str(ProfileHalfWidth),(0.01,0.01),
                          xycoords='subfigure fraction',size=8)
    ax_spag.annotate("Smoothing="+str(smooth),(0.01,0.03),
                          xycoords='subfigure fraction',size=8)

    #path="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/"
    #figspaghetti.savefig(path+"Profiles/output/Profile_"+reference+"_"+param+"_"+profile+".png",dpi=300)

    ###########################################################################
    # Compute average and standard deviation (unweighted) from all the 
    #   spaghetti profiles.
    ###########################################################################
    AvgPro=np.mean(AvgArr[:,:],axis=0)
    AvgStd=np.std(AvgArr[:,:],axis=0)
    Avgamf=np.mean(AvgamfArr[:,:],axis=0)

    if smooth:
        AvgProSmth=convolve(AvgPro, Box1DKernel(3))
        AvgStdSmth=convolve(AvgStd, Box1DKernel(3))
        OutPro,OutStd=AvgProSmth,AvgStdSmth
    else:
        OutPro,OutStd,Outamf=AvgPro,AvgStd,Avgamf

    return(Lats,OutPro,OutStd,Outamf,Num)
        