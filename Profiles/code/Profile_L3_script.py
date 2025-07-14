def Profile_L3_script(figxy=[6.0,6.0],LatPlotLims=[45,135],ZonePlotHalfWidth=60,
                      bands=["STrZ","SEB","SEZ","NEZ","NEB","NTrZ"],
                      colors=["C0","C1","C2","C3","C4","C5"]):
    """
    PURPOSE:    Compute changes over time in meridional profiles, plotted on
                a scatter plot.
                
    CALLS: Profile_L3 and plot_profile_scatter

    Parameters
    ----------
    figxy : TYPE, optional
        DESCRIPTION. The default is [6.0,6.0].
    bands : TYPE, optional
        DESCRIPTION. The default is ["STrZ","SEB","SEZ","NEZ","NEB","NTrZ"].
    colors : TYPE, optional
        DESCRIPTION. The default is ["C0","C1","C2","C3","C4","C5"].

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
    
    import os
    import pylab as pl
    import numpy as np
    import matplotlib.dates as mdates
    from datetime import datetime

    import plot_profile_scatter as pps
    import Profile_L3_multi as PL3M
    
    ###########################################################################
    # COMPUTE CLOUD PRESSURE AND AMMONIA PROFILES FOR FOUR OBSERVATION SERIES
    ###########################################################################
    CH4SCT22,CH4VLT22,CH4SCT23,CH4SCT24=PL3M.Profile_L3_multi(param="PCld",
                                                       profile="Meridional",
                                                       ProfileHalfWidth=2,
                                                       LatPlotLims=LatPlotLims,
                                                       ZonePlotHalfWidth=ZonePlotHalfWidth,
                                                       smooth=False,inset=False)
    
    NH3SCT22,NH3VLT22,NH3SCT23,NH3SCT24=PL3M.Profile_L3_multi(param="fNH3",
                                                       profile="Meridional",
                                                       ProfileHalfWidth=2,
                                                       LatPlotLims=LatPlotLims,
                                                       ZonePlotHalfWidth=ZonePlotHalfWidth,
                                                       smooth=False,inset=False)
    
    fig3,axs3=pl.subplots(1,1,figsize=(figxy[0],figxy[1]), dpi=150, facecolor="white")
    fig3.suptitle("Cloud-Top Pressure vs Ammonia Abundance: 2022-2024",
                  fontsize=14,x=0.5,ha='center',color='k')

    ###########################################################################
    # SET UP AND MAKE SCATTER PLOT OF SCT PROFILE DATA BY MERIDIONAL BAND
    ###########################################################################
    # Set plot parameters
    fNH3low=60
    fNH3high=140
    PCldlow=1400
    PCldhigh=2000
    #micronlow=0.5
    #micronhigh=3.5

    axs3.tick_params(axis='both', which='major', labelsize=10)
    axs3.set_title("Cloud-Top Pressure vs Ammonia Abundance: 2022-2024",fontsize=12)
    
    pps.plot_profile_scatter(CH4SCT22["Pro"],NH3SCT22["Pro"],CH4SCT22["Lats"],axs3,PCldlow,PCldhigh,
                     fNH3low,fNH3high,False,bands,colors,marker='o',Level="L3",
                     leg=True,axis_inv=True)
    pps.plot_profile_scatter(CH4SCT23["Pro"],NH3SCT23["Pro"],CH4SCT23["Lats"],axs3,PCldlow,PCldhigh,
                     fNH3low,fNH3high,False,bands,colors,marker='^',Level="L3",
                     leg=False,axis_inv=True)
    pps.plot_profile_scatter(CH4SCT24["Pro"],NH3SCT24["Pro"],CH4SCT24["Lats"],axs3,PCldlow,PCldhigh,
                     fNH3low,fNH3high,False,bands,colors,marker='s',Level="L3",
                     leg=False,axis_inv=True)
    
    axs3.scatter(62,1410,marker='o',c='k',s=50)
    axs3.scatter(62,1440,marker='^',c='k',s=50)
    axs3.scatter(62,1470,marker='s',c='k',s=50)
    
    axs3.annotate('2022',xy=(65,1410), xycoords='data',xytext=(64,1410),fontsize=9,verticalalignment='center_baseline')
    axs3.annotate('2023',xy=(65,1440), xycoords='data',xytext=(64,1440),fontsize=9,verticalalignment='center_baseline')
    axs3.annotate('2024',xy=(65,1470), xycoords='data',xytext=(64,1470),fontsize=9,verticalalignment='center_baseline')
    path="/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Profiles/output/"

    fig3.subplots_adjust(left=0.14, bottom=0.10, right=0.94, top=0.92)  

    fig3.savefig(path+"Clouds+Ammonia 2022-2024.png",dpi=300)
