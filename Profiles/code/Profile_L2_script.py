def Profile_L2_script(figxy=[6.0,6.0],bands=["STrZ","SEB","SEZ","NEZ","NEB","NTrZ"],
                      colors=["C0","C1","C2","C3","C4","C5"]):

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
    import Profile_L2 as PL2

    fNH3low=0.2
    fNH3high=0.8
    PCldlow=1.2
    PCldhigh=1.7
    micronlow=0.5
    micronhigh=3.5

    CH4SCT22,CH4VLT22,CH4SCT23,CH4VLT23,CH4SCT24=PL2.Profile_L2(band="CH4",
                                                       profile="Meridional",
                                                       ProfileHalfWidth=1,
                                                       LatPlotLims=[45,135],
                                                       ZonePlotHalfWidth=60,
                                                       smooth=False)
    
    NH3SCT21,NH3VLT21,NH3SCT22,NH3VLT22,NH3SCT23,NH3VLT23,NH3SCT24=PL2.Profile_L2(band="NH3",
                                                       profile="Meridional",
                                                       ProfileHalfWidth=1,
                                                       LatPlotLims=[45,135],
                                                       ZonePlotHalfWidth=60,
                                                       smooth=False)
    
    fig3,axs3=pl.subplots(1,1,figsize=(figxy[0],figxy[1]), dpi=150, facecolor="white")
    #fig3.suptitle(suptitle,x=0.5,ha='center',color='k')
    fig3.suptitle("Molecular Absorption: 2022-2024",
                  fontsize=14,x=0.5,ha='center',color='k')

    axs3.tick_params(axis='both', which='major', labelsize=10)
    
    pps.plot_profile_scatter(CH4SCT22["Pro"],NH3SCT22["Pro"],CH4SCT22["Lats"],axs3,PCldlow,PCldhigh,
                     fNH3low,fNH3high,False,bands,colors,marker='o',Level="L2",leg=True,axis_inv=False)
    pps.plot_profile_scatter(CH4SCT23["Pro"],NH3SCT23["Pro"],CH4SCT23["Lats"],axs3,PCldlow,PCldhigh,
                     fNH3low,fNH3high,False,bands,colors,marker='^',Level="L2",leg=False,axis_inv=False)
    pps.plot_profile_scatter(CH4SCT24["Pro"],NH3SCT24["Pro"],CH4SCT24["Lats"],axs3,PCldlow,PCldhigh,
                     fNH3low,fNH3high,False,bands,colors,marker='s',Level="L2",leg=False,axis_inv=False)

    axs3.scatter(0.22,1.68,marker='o',c='k',s=50)
    axs3.scatter(0.22,1.66,marker='^',c='k',s=50)
    axs3.scatter(0.22,1.64,marker='s',c='k',s=50)
    
    axs3.annotate('2022',xy=(0.22,1.68), xycoords='data',xytext=(0.23,1.68),fontsize=9,verticalalignment='center_baseline')
    axs3.annotate('2023',xy=(0.22,1.66), xycoords='data',xytext=(0.23,1.66),fontsize=9,verticalalignment='center_baseline')
    axs3.annotate('2024',xy=(0.22,1.64), xycoords='data',xytext=(0.23,1.64),fontsize=9,verticalalignment='center_baseline')
    path="/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Profiles/output/"


    fig3.subplots_adjust(left=0.14, bottom=0.10, right=0.94, top=0.92)  

    fig3.savefig(path+"Absorption 2022-2024.png",dpi=300)
    