def SEB_Outbreak_Script(collection="20241129-20241129 SEB OB",
                        ctbls=["terrain_r","Blues"],localmax=False,close=False):

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


    import L3_Jup_Map_Plot as L3JMP
    import get_map_collection as gmc
    import MakeContiguousMap as MCM
    
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

    
    obslist,dummy=gmc.get_map_collection(collection)
    
    if collection=="20241115-20241115":
        ROI={"CORE":[103,106,324.0,2.0],
             "NE Plume":[99,102,320.0,2.0],
             "SW Plume":[105,108,329.0,2.0],
             "SEB Ref1":[102,104,334,2.0],
             "SEB Ref2":[100,108,311,2.0],
             "STrZ Ref":[111,115,323.0,3.0]}     
        obslist=["20241115UTa","20241115UTb","20241115UTc","20241115UTd","20241115UTe"]
        CM=[322]

    if collection=="20241129-20241129 SEB OB":
        ROI={"CORE":[102,105,323.0,2.0],
             "NE Plume":[98,101,320.0,2.0],
             "SW Plume":[104,107,329.0,2.0],
             "SEB Ref1":[102,104,334,2.0],
             "SEB Ref2":[100,108,311,2.0],
             "STrZ Ref":[111,115,323.0,3.0]}
        obslist=["20241129UTg","20241129UTh","20241129UTi",
                 "20241129UTj","20241129UTk","20241129UTl","20241129UTm"]
        CM=[322]
    
    if collection=="20241202-20241202 SEB OB":
        ROI={"CORE":[102,105,326.0,2.0],
             "NE Plume":[101,104,318.0,2.0],
             "SW Plume":[106,109,333.0,2.0],
             "SEB Ref1":[102,104,336,2.0],
             "SEB Ref2":[100,108,312,2.0],
             "STrZ Ref":[111,115,326.0,3.0]}
        obslist=["20241202UTc","20241202UTd","20241202UTe","20241202UTf",
                 "20241202UTg","20241202UTh"]
        CM=[325]
    
    tempfNH3=np.zeros([180,360])
    tempPCloud=np.zeros([180,360])
    
    First=True
    coefs=[0.75,0.45]
    #coefs=[0.0,0.0]
    for o in obslist:
        dateobs,roilabel,mean1,stdv1,mean2,stdv2,dummy=\
            L3JMP.L3_Jup_Map_Plot(obskey=o,imagetype='Map',target="Jupiter",
                        Smoothing=False,LatLims=[90,120],LonRng=15,CMpref=323,
                        LonSys='2',showbands=False,coef=coefs,
                        subproj='convection',figxy=[8.0,4.0],FiveMicron=False,ROI=ROI,
                        ctbls=["terrain_r","Blues"])
            
        print("*********",dateobs,roilabel,mean1,stdv1,mean2,stdv2)
        if First:
            datearr=[dateobs]
            roilabelarr=[roilabel]
            mean1arr=[mean1]
            stdv1arr=[stdv1]
            mean2arr=[mean2]
            stdv2arr=[stdv2]
        else:
            datearr.append(dateobs)
            roilabelarr.append(roilabel)
            mean1arr.append(mean1)
            stdv1arr.append(stdv1)
            mean2arr.append(mean2)
            stdv2arr.append(stdv2)
        First=False
            
    if close:
        pl.close('all')

    ###############################################################################
    #
    ###############################################################################
    fig,axs=pl.subplots(1,2,figsize=(8,4), dpi=150, facecolor="white")
    #fig3.suptitle(suptitle,x=0.5,ha='center',color='k')
    fig.suptitle(dateobs[0:10]+" Regions of Interest vs Time",x=0.5,ha='center',color='k')
    axs[0].set_title("Effective Cloud Top Pressure (mb)",fontsize=10)
    axs[1].set_title("Ammonia Mole Fraction (ppm)",fontsize=10)
    fig1,axs1=pl.subplots(1,2,figsize=(8,4), dpi=150, facecolor="white",sharey=True)
    #fig3.suptitle(suptitle,x=0.5,ha='center',color='k')
    fig1.suptitle(dateobs[0:10]+" ROI Cloud Top Pressure vs Ammonia Abundance",x=0.5,ha='center',color='k')
    axs1[0].set_title("Individual ROI Sample Averages",fontsize=10)
    axs1[1].set_title("ROIs Averaged over Time",fontsize=10)
    
    
    print()
    for k in range(0,len(datearr)):    
        datearr[k]=datetime.strptime(datearr[k], '%Y-%m-%dT%H:%M:%SZ')
    datearray=np.array(datearr)
    
    mean1array=np.array(mean1arr)
    mean2array=np.array(mean2arr)
    stdv1array=np.array(stdv1arr)
    stdv2array=np.array(stdv2arr)
    roilabelarray=np.array(roilabelarr)
    print(roilabelarray,mean1array,stdv1array,mean2array,stdv2array)

    meanmean1=np.mean(mean1array,axis=0)
    stdvstdv1=np.std(mean1array,axis=0)
    meanmean2=np.mean(mean2array,axis=0)
    stdvstdv2=np.std(mean2array,axis=0)
    
    print("!!!!!!!!!!!!!!!!!!!!!!!")
    print(obslist)
    print(mean1array.shape)
    
    for i in range(0,mean1array.shape[1]):
        axs[0].errorbar(datearray,mean1array[:,i],yerr=stdv1array[:,i],label=roilabelarray[0,i],
                        capsize=3.0,elinewidth=0.5,capthick=0.5)
        axs[1].errorbar(datearray,mean2array[:,i],yerr=stdv2array[:,i],label=roilabelarray[0,i],
                        capsize=3.0,elinewidth=0.5,capthick=0.5)

        axs1[0].errorbar(mean2array[:,i],mean1array[:,i],xerr=stdv2array[:,i],yerr=stdv1array[:,i],
                         label=roilabelarray[0,i],linewidth=0.0,
                         elinewidth=0.5,marker="o",markersize=3.0)
        axs1[1].errorbar(meanmean2[i],meanmean1[i],xerr=stdvstdv2[i],yerr=stdvstdv1[i],
                         label=roilabelarray[0,i],
                         elinewidth=0.5,marker="o",markersize=3.0)
        
    axs[0].legend(ncols=2,fontsize=8)
    axs1[0].legend(ncols=2,fontsize=8)
    
    axs[0].tick_params(axis='x', which='major', labelsize=6)
    axs[0].tick_params(axis='y', which='major', labelsize=8)
    
    axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[0].xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
    axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[1].xaxis.set_minor_formatter(mdates.DateFormatter("%H:%M"))
    
    axs[1].tick_params(axis='x', which='major', labelsize=6)
    axs[1].tick_params(axis='y', which='major', labelsize=8)
    axs1[0].tick_params(axis='both', which='major', labelsize=8)
    axs1[1].tick_params(axis='both', which='major', labelsize=8)
    
    axs[0].set_ylabel("Cloud Top Pressure (mb)",fontsize=10)        
    axs[1].set_ylabel("Ammonia Mole Fraction (ppm)",fontsize=10) 
    
    axs[0].set_xlabel("Date",fontsize=10)        
    axs[1].set_xlabel("Date",fontsize=10) 
           
    axs1[0].set_ylabel("Cloud Top Pressure (mb)",fontsize=10)        
    axs1[0].set_xlabel("Ammonia Mole Fraction (ppm)",fontsize=10) 
    axs1[1].set_xlabel("Ammonia Mole Fraction (ppm)",fontsize=10) 
    """
    if "PCloud" in Rtitle and "5um" in Rtitle:
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
    """
    
    
    axs[0].grid(linewidth=0.2)
    axs[0].set_ylim(PCldlow,PCldhigh)
    #axs[0].set_xlim(0,5)
    axs[0].invert_yaxis()
    
    axs[1].grid(linewidth=0.2)
    axs[1].set_ylim(fNH3low,fNH3high)
    #axs[1].set_xlim(0,5)
    
    axs1[0].grid(linewidth=0.2)
    axs1[0].set_ylim(PCldlow,PCldhigh)
    axs1[0].set_xlim(fNH3low,fNH3high)
    axs1[0].invert_yaxis()
    
    axs1[1].grid(linewidth=0.2)
    axs1[1].set_ylim(PCldlow,PCldhigh)
    axs1[1].set_xlim(fNH3low,fNH3high)
    axs1[1].invert_yaxis()
    
    #box = axs3[1].get_position()
    #axs3[1].set_position([box.x0+0.03, box.y0-0.01, box.width * 0.5, box.height * 1.015])    
    fig.subplots_adjust(left=0.08, right=0.97, top=0.85, bottom=0.12)
    fig1.subplots_adjust(left=0.08, right=0.97, top=0.85, bottom=0.12)
    #fig3.savefig(pathout+fnout,dpi=300)
    pathout='C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Studies/convection/'
    fig.savefig(pathout+collection+" ROI Time Series",dpi=300)
    fig1.savefig(pathout+collection+" ROI Scatter",dpi=300)

    
    ###############################################################################
    #!!DO BLENDED MAP SEPARATELY USING MAKECONTIGUOUSMAP.PY
    ###############################################################################
    
    MCM.MakeContiguousMap(False,False,False,collection=collection,LonSys='2',
                          FiveMicron=False,lats=[90,120],LonLims=[CM[0]-15,CM[0]+15],
                          figsz=[3.0,6.0],ROI=ROI,variance=True,proj="convection",
                          ctbls=["terrain_r","Blues"],localmax=localmax)
