# -*- coding: utf-8 -*-
"""
Created on Sun Dec 05 08:48:16 2021

@author: Steven Hill
"""

def plot_TEXES_Groups(ax,clr="C2",prs=0.43798,mult=1.0):
    """
    PURPOSE:    This code reads and plots the NH3 mole fraction data from
                Fletcher et al., 2016 at 440mb pressure. It loops over the
                seven groups of observations and computes an average and
                standard deviation for plotting.
    """
    import sys
    import matplotlib.pyplot as pl
    import scipy
    import numpy as np
    import matplotlib.pyplot as pl
    sys.path.append('./Photometry/code')
    sys.path.append('c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Photometry/code/')

    import ComputeNetRateJupiter_P3 as CNRJ

    #print('***************** mult=',mult)
    pth="c:/Astronomy/Projects/SAS 2021 Ammonia/GeminiTEXES2017/ZonalResults/"
    pressure = scipy.fromfile(file=pth+"zmean_g1_retnh3_pressure.txt", dtype=float, count=-1, sep=" ")
    data=np.zeros((7,181))
    #print("data=",data.shape)

    ind=np.where(np.abs(pressure-prs)/pressure<0.01)    #Pressure index
    PL=np.ndarray.flatten(np.array(ind))[0]             #Pressure
    #print('&&&&&&&&&&&&&&&&&&&&&&&&',ind)
    #print('&&&',np.ndarray.flatten(np.array(ind))[0],'***PL=',PL)
    
    PL=np.ndarray.flatten(np.array(ind))[0]
    for i in range(1,8):
        latc = scipy.fromfile(file=pth+"zmean_g"+str(i)+"_retnh3_lat.txt", dtype=float, count=-1, sep=" ")
        latg=Centric_to_Graphic(latc)
        latsize=len(latg)
        #print("latsize=========",latsize)
        Start=latsize*PL
        End=latsize*(PL+1)
        """if i==4:
            Start=59*PL
            End=59*(PL+1)
            latc = scipy.fromfile(file=pth+"zmean_g4_retnh3_lat.txt", dtype=float, count=-1, sep=" ")
            latg=Centric_to_Graphic(latc)
        else:
            Start=61*PL
            End=61*(PL+1)
            latc = scipy.fromfile(file=pth+"zmean_g1_retnh3_lat.txt", dtype=float, count=-1, sep=" ")
            latg=Centric_to_Graphic(latc)"""
        #print("Start,End=",Start,End)
        #print(i)
        tmp = scipy.fromfile(file=pth+"zmean_g"+str(i)+"_retnh3_data.txt", dtype=float, count=-1, sep=" ")
        dat=tmp[Start:End]
        latgrid,tmpsig=CNRJ.uniform_lat_grid(latg,dat,Fine=True)

        #print("dat=",tmpsig.shape, tmpsig)
        data[i-1,:]=tmpsig
    #print(Start,End,pressure[PL])
    scaled_data_mean=np.mean(data,axis=0)*mult#*8.0e4
    scaled_data_std=np.std(data,axis=0)*mult#*8.0e4
    if ax!='None':
        ax.plot(latgrid,scaled_data_mean,linewidth=1.0,
                label='Fletcher etal, 2020 ('+str(int(prs*1000.))+'mb)',color=clr)
        ax.fill_between(latgrid, scaled_data_mean-scaled_data_std, scaled_data_mean+scaled_data_std,
                        color=clr,alpha=0.08)
    #ax.plot(latgrid,data[3,:]*mult)
    #for i in range(1,8):
    #    ax.plot(latgrid,data[i-1,:]*mult)
    return(scaled_data_mean,scaled_data_std)
    
def plot_Teifel(ax,clr='C0',width=1.5):
    """
    PURPOSE:    This code reads and plots the zonally averaged NH3 absorption
                EW at 645nm from data scanned from Teifel et al., 2018 
                figure 7.
    """
    import numpy as np
    from numpy import genfromtxt
    Lats=[-24.0,-12.5,0.0,12.5,24.0]
    EWs=np.array([5.92,6.78,6.75,6.35,5.38])*0.1
    #ax.scatter(Lats,EWs)
    pth="c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/"
    Teifel = np.array(genfromtxt(pth+"Teifel2018-Fig7.txt", delimiter=','))
    ax.plot(Teifel[4:24,0],Teifel[4:24,1]*0.1,label='Teifel etal, 2018',
            linewidth=width,color=clr)

def plot_Historical(ax,reference,clr='C0'):
    """
    PURPOSE:    This code reads and plots the zonally averaged NH3 absorption
                EW at 645nm from data scanned from Teifel et al., 2018 
                figure 7.
    """
    import matplotlib.pyplot as pl
    import numpy as np
    data={"Vdov2021":{"Region":["SPR","STB","STrZ","SEB","EZ","NEB","NTrZ","NTB","NPR","GRS"],
                      "Center_pgLat":[-45.00,-29.75,-23.40,-13.45,-0.15,15.60,20.80,27.8,45.00,-23.00],
                      "645EW":[4.2,5.1,5.7,6.0,6.2,5.5,4.4,4.7,4.6,4.7]},
          "Teif2018":{"Region":["STrZ","SEB","EZ","NEB","NTrZ"],
                            "Center_pgLat":[-23.40,-13.45,-0.15,15.60,20.8],
                            "645EW":[5.92,6.78,6.75,6.35,5.38]},
          "More1991":{"Region":["SPR","STB","STrZ","SEB","EZ","NEB","NTB","NPR"],
                            "Center_pgLat":[-45.00,-29.75,-23.40,-13.45,-0.15,20.80,27.8,45],
                            "645EW":[5.6,5.7,7.8,9.8,7.7,4.9,5.7,7.2]},
          "L&O1980":{"Region":["STrZ","SEB","EZ","NEB","GRS"],
                            "Center_pgLat":[-23.40,-13.45,-0.15,15.60,-23.0],
                            "645EW":[11.00,12.90,11.7,8.3,9.3]}}
    ax.scatter(data[reference]["Center_pgLat"],np.array(data[reference]["645EW"])*0.1,label=reference,color=clr)


def plot_profile_L2(ax1,ax2,ax3,ax4,reference,ProfileHalfWidth=45.,LatPlotLims=[45,135],
                    ZonePlotHalfWidth=45.,band="CH4",profile="Meridional",
                    clr='C0',width=1.0,style='solid',smooth=False):
    """
    PURPOSE:
    Plot L2 molecular absorption profiles as equivalent widths
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
    import read_master_calibration
    import extract_profile as EP
    import get_obs_list as GOL
    import pylab as pl
    import get_batch_lists as GBL

    ###########################################################################
    # Get sourcefile list, batch list, and set input path
    ###########################################################################
    sourcefiles=GOL.get_obs_list()
    DataSets=GBL.get_batch_lists()
    pth="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Analysis Data/L2 FITS/"

    ###########################################################################
    # Set linear fit for conversion of transmission to equivalent width
    #!!! This should be simplified to use only the single master calibration
    #!!! value derived from VLT. It could be incorporated into an external
    #!!! module
    ###########################################################################
    Trans2EW={"VLTXUSE":{"EW_slope":{"NH3":-12.15343491,"CH4":-12.74115968},
                     "EW_const":{"NH3":12.15195684,"CH4":12.73323535}},
              #"CMOS":{"EW_slope":{"NH3":-12.87087479,"CH4":-13.52007859},
              #       "EW_const":{"NH3":12.86940675,"CH4":13.51315382}},
              #"VLTMUSE":{"EW_slope":{"NH3":-12.87087479,"CH4":-13.52007859},
              #       "EW_const":{"NH3":12.86940675,"CH4":13.51315382}}}
              "CMOS":{"EW_slope":{"NH3":-12.04629662,"CH4":-12.86742542},
                     "EW_const":{"NH3":12.04540329,"CH4":12.86268728}},
              "VLTMUSE":{"EW_slope":{"NH3":-12.04629662,"CH4":-12.86742542},
                     "EW_const":{"NH3":12.04540329,"CH4":12.86268728}}}
   
    suffix={"CH4":"-Jupiter_Map_L2TCH4",
            "NH3":"-Jupiter_Map_L2TNH3"}
    
    First=True
    Num=0
    for obskey in DataSets[reference]:
        print("*******obskey=",obskey)
        dataset=obskey[0:10].replace('-','')+'UT'+obskey[10]#+"_Map"
        print("******dataset=",dataset)

        try:
            file=sourcefiles[dataset][band+'file'][0:17]+suffix[band]+\
                    sourcefiles[dataset]['Metadata']['Variation']
            variation=sourcefiles[dataset]['Metadata']['Variation']
        except:
            file=sourcefiles[dataset][band+'file'][0:17]+suffix[band]
            variation=""
        
        print("file=",file)
        #######################################################################
        # Extract a profile out of a Level 2 or 3 FITS map file
        #######################################################################
        Lats,AvgProf,StdProf,CM2,amfAvgProf=EP.extract_profile(pth,file+".fits",
                                                ProfileHalfWidth=ProfileHalfWidth,
                                                profile=profile)

        print("@@@@@@@@@@@@reference[5:]]=",reference[5:],Trans2EW[reference[5:]])
        
        ###!!!!!!!!!!NEED DIAGNOSTICS HERE!!!
        print("Trans2EW[reference[5:]]['EW_slope'][band]=",Trans2EW[reference[5:]]["EW_slope"][band])
        print("Trans2EW[reference[5:]]['EW_const'][band]=",Trans2EW[reference[5:]]["EW_const"][band])
        print("Avg(AvgProf)=",np.mean(AvgProf))
        Avg_EW=Trans2EW[reference[5:]]["EW_slope"][band]*AvgProf+\
            Trans2EW[reference[5:]]["EW_const"][band]
        Std_EW=Trans2EW[reference[5:]]["EW_slope"][band]*AvgProf+\
            Trans2EW[reference[5:]]["EW_const"][band]

        if First:
            AvgArr=Avg_EW
            AvgamfArr=amfAvgProf
            if profile=="Meridional":
                AvgArr=np.reshape(AvgArr,(1,180))
                AvgamfArr=np.reshape(AvgamfArr,(1,180))
            elif profile=="Zonal":
                AvgArr=np.reshape(AvgArr,(1,360))
                AvgamfArr=np.reshape(AvgamfArr,(1,360))
            First=False
        else:
            AvgArr=np.vstack((AvgArr,Avg_EW))
            AvgamfArr=np.vstack((AvgamfArr,amfAvgProf))

        #######################################################################
        # Increment counter and plot a line of spaghetti
        #######################################################################
        Num=Num+1                
        ax3.plot(Lats,Avg_EW,linewidth=0.5,label=file[5:15])
        ax4.scatter(amfAvgProf,Avg_EW,s=2,label=file[5:15])
        
    ###########################################################################
    # Configure plot layout parameters, labels, and titles
    ###########################################################################
    if profile=="Meridional":
        xlabel="Planetographic Latitude (deg)"
    elif profile=="Zonal":
        xlabel="Longitude from CM (deg)"
        
    ax3.tick_params(axis='both', which='major', labelsize=7)
    ax3.grid(linewidth=0.2)
    ax3.set_title(reference,fontsize=9)
    if reference=="2023 CMOS":
        ax3.legend(fontsize=5,ncol=4,loc="upper right",bbox_to_anchor=(2.0, 1.0))
    else:
        ax3.legend(fontsize=5,ncol=4)
    
    ax4.tick_params(axis='both', which='major', labelsize=7)
    ax4.grid(linewidth=0.2)
    ax4.set_title(reference,fontsize=9)
    if reference=="2023 CMOS":
        ax4.legend(fontsize=5,ncol=4,loc="upper right",bbox_to_anchor=(2.0, 1.0))
    else:
        ax4.legend(fontsize=5,ncol=4)

    
    path="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/"
    #figspaghetti.savefig(path+"Profiles/output/Profile_"+reference+
    #                     "_"+band+"_"+profile+"_Absorption.png",dpi=300)  

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

    ax1.plot(Lats,OutPro,color=clr,linewidth=width,linestyle=style,
            label=reference+' (Avg. '+str(Num)+')')  
    ax1.fill_between(Lats, OutPro-OutStd, OutPro+OutStd,
                    color=clr,alpha=.1)

    print(AvgamfArr.shape,OutPro.shape)
    ax2.scatter(Avgamf,OutPro,s=2,label=reference+' (Avg. '+str(Num)+')')


    return(Lats,OutPro,OutStd,Outamf)

def plot_profile_L3(ax1,ax2,ax3,ax4,reference,ProfileHalfWidth=45.,LatPlotLims=[45,135],
                    ZonePlotHalfWidth=45.,param="fNH3",profile="Meridional",
                    clr='C0',width=1.0,style='solid',smooth=True):

    #(ax,reference,LatLims=[45,135],LonRng=45.,profile="Meridional",
    #                 clr='C0',width=1.0,
    #                 style='solid',smooth=False,param='PCloud'):
    """
    PURPOSE:
    Plot L2 molecular absorption profiles as equivalent widths
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
    from astropy.io import fits
    import RetrievalLibrary as RL
    import pylab as pl
    import sys
    sys.path.append('./Services')
    #import read_master_calibration
    import RetrievalLibrary as RL
    import extract_profile as EP
    import get_obs_list as GOL
    import get_batch_lists as GBL

    ###########################################################################
    # Get sourcefile list, batch list, and set input path
    ###########################################################################
    sourcefiles=GOL.get_obs_list()
    DataSets=GBL.get_batch_lists()
    pth="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Analysis Data/L3 FITS/"

    if param=="PCld":
        paramkey='CH4file'
        factor=1.0
    elif param=="fNH3":
        paramkey='NH3file'
        factor=1.0
    
    #pth=data[reference]["path"]
    if profile=="Meridional":
        xlabel="Planetographic Latitude (deg)"
    elif profile=="Zonal":
        xlabel="Longitude from CM (deg)"

    suffix={"PCld":"-Jupiter_Map_L3PCld_S0.fits",
            "fNH3":"-Jupiter_Map_L3fNH3_S0.fits"}

    #figspaghetti,axsspaghetti=pl.subplots(1,1,figsize=(6.0,4.0), dpi=150,
    #                                      facecolor="white")
    #figamfscat,axsamfscat=pl.subplots(1,1,figsize=(6.0,4.0), dpi=150, 
    #                                      facecolor="white")

    First=True
    Num=0

    for obskey in DataSets[reference]:
        print("*******obskey=",obskey)
        dataset=obskey[0:10].replace('-','')+'UT'+obskey[10]#+"_Map"
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
        Lats,AvgProf,StdProf,CM3,amfAvgProf=EP.extract_profile(pth,file,
                                                ProfileHalfWidth=ProfileHalfWidth,
                                                profile=profile)
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
        ax3.plot(Lats,AvgProf*factor,linewidth=0.5,label=file[5:17])
        ax4.scatter(amfAvgProf,AvgProf,s=5,label=file[5:17])


    AvgPro=np.mean(AvgArr[:,:],axis=0)*factor
    AvgStd=np.std(AvgArr[:,:],axis=0)*factor
    Avgamf=np.mean(AvgamfArr[:,:],axis=0)*factor

    if profile=="Meridional":
        ax3.set_xlim(90-LatPlotLims[1],90-LatPlotLims[0])
        ax3.set_xlabel("Planetographic Latitude (deg)",fontsize=10)
        ax4.set_xlim(1,3)
        ax4.set_xlabel("Planetographic Latitude (deg)",fontsize=10)
    elif profile=="Zonal":
        ax3.set_xlim(-ZonePlotHalfWidth,ZonePlotHalfWidth)
        ax3.set_xlabel("Longitude from CM (deg)",fontsize=10)
        ax4.set_xlim(1,3)
        ax4.set_xlabel("Longitude from CM (deg)",fontsize=10)

    if param=="PCld":
        ax3.set_ylim(0,1100)
        ax3.set_ylabel("Effective Cloud-Top Pressure (mb)",fontsize=10)
        ax3.invert_yaxis()
        ax4.set_ylim(0,1100)
        ax4.set_ylabel("Effective Cloud-Top Pressure (mb)",fontsize=10)
        ax4.invert_yaxis()
    elif param=="fNH3":
        ax3.set_ylim(0,200)
        ax3.set_ylabel("Column-Average Ammonia Abundance (ppm)",fontsize=10)
        ax4.set_ylim(0,200)
        ax4.set_ylabel("Column-Average Ammonia Abundance (ppm)",fontsize=10)
    #axsspaghetti.legend(fontsize=6,ncol=5)
    ax3.set_title(reference)
    ax3.grid(linewidth=0.2)

    ax3.annotate("ProfileHalfWidth="+str(ProfileHalfWidth),(0.01,0.01),
                          xycoords='subfigure fraction',size=8)
    ax3.annotate("Smoothing="+str(smooth),(0.01,0.03),
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
        
    ax1.plot(Lats,OutPro,color=clr,linewidth=width,linestyle=style,
            label=reference+' (Avg. '+str(Num)+')')  
    ax1.fill_between(Lats, OutPro-OutStd, OutPro+OutStd,
                    color=clr,alpha=.1)
    ax2.scatter(Outamf,OutPro,s=5,label=reference+' (Avg. '+str(Num)+')')


    return(Lats,OutPro,OutStd,Outamf)
        
def Centric_to_Graphic(Latc):
    #Formula used is from Simon and Beebe, 1996
    import numpy as np
    Req=71492.0
    Rp=66845.0
    Latg=Latc
    for i in range(len(Latc)):
        Latg[i]=np.arctan(((Req/Rp)**2)*np.tan(Latc[i]*np.pi/180.))*180.0/np.pi
    #print(Latg)
    return Latg