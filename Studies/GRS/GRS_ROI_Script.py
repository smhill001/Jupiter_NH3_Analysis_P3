def GRS_ROI_Script(collection="2022 GRS",dateformat="%H:%M",ctbls=["terrain_r","Blues"]):

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
    import copy
    
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
    
    if collection=="2022 GRS VLT":
        ROI={"CORE":[110,114,22.0,3.0],
             "Chimney":[102,106,22,1.0],
             "Collar":[106,108,22.0,3.0],
             "Hook":[109,113,31.0,1.0],
             "SEBn":[108,110,7.0,3.0],
             "SEBs":[113,115,7.0,3.0],
             "SEBz":[110,113,7.0,3.0],
             "STrZ":[116,118,7.0,3.0],
             "Wake":[115,117,37.0,2.0],             
             "SEZ":[92,94,22,5.0]}
            
        obslist=["20220730UTa","20220919UTa"]
        CM=[22,25]

    if collection=="2022 GRS SCT":
        ROI={"CORE":[110,114,22.0,3.0],
             "Chimney":[102,106,22,1.0],
             "Collar":[106,108,22.0,3.0],
             "Hook":[109,113,31.0,1.0],
             "SEBn":[108,110,7.0,3.0],
             "SEBs":[113,115,7.0,3.0],
             "SEBz":[110,113,7.0,3.0],
             "STrZ":[116,118,7.0,3.0],
             "Wake":[115,117,37.0,2.0],             
             "SEZ":[92,94,22,5.0]}

        obslist=["20220810UTa","20220818UTa","20220828UTa",
                 "20220904UTa","20221013UTa","20221020UTa"]
        CM=[22.5,23,23,24,26,27]

    if collection=="2022 GRS":
        ROI={"CORE":[110,114,22.0,3.0],
             "Chimney":[102,106,22,1.0],
             "Collar":[106,108,22.0,3.0],
             "Hook":[109,113,31.0,1.0],
             "Wake":[115,117,37.0,2.0],             
             "SEBn":[108,110,7.0,3.0],
             "SEBs":[113,115,7.0,3.0],
             "SEBz":[110,113,7.0,3.0],
             "STrZ":[116,118,7.0,3.0],
             "SEZ":[92,94,22,5.0]}

        obslist=["20220730UTa","20220810UTa","20220818UTa","20220828UTa",
                 "20220904UTa","20220919UTa","20221013UTa","20221020UTa"]
        CM=[22,22.5,23,23,24,25,26,27]

    if collection=="2023 GRS VLT":
        ROI={"CORE":[110,114,42.0,3.0],
             "Chimney":[102,106,42,1.0],
             "Collar":[106,108,42.0,3.0],
             "Hook":[109,113,52.0,1.0],
             "Wake":[115,117,57.0,2.0],             
             "SEBn":[108,110,27.0,3.0],
             "SEBs":[113,115,27.0,3.0],
             "SEBz":[110,113,27.0,3.0],
             "STrZ":[116,118,27.0,3.0],
             "SEZ":[92,94,42,5.0]}

        obslist=["20230812UTa","20230815UTa",
                 "20230923UTa"]
        CM=[42,42,44]
        
    if collection=="2023 GRS SCT":
        ROI={"CORE":[110,114,42.0,3.0],
             "Chimney":[102,106,42,1.0],
             "Collar":[106,108,42.0,3.0],
             "Hook":[109,113,52.0,1.0],
             "Wake":[115,117,57.0,2.0],             
             "SEBn":[108,110,27.0,3.0],
             "SEBs":[113,115,27.0,3.0],
             "SEBz":[110,113,27.0,3.0],
             "STrZ":[116,118,27.0,3.0],
             "SEZ":[92,94,42,5.0]}

        obslist=["20230831UTa","20230924UTa","20231005UTa",
                 "20231017UTa","20231022UTa","20231103UTa",
                 "20231113UTa","20231207UTa"]
        CM=[42,44,44,44,45,45,45,46]
        
    if collection=="20241027 GRS":
        ROI={"CORE":[111,115,62.0,3.0],
             "Chimney":[102,106,62,1.0],
             "Collar":[106,108,62.0,3.0],
             "Hook":[109,113,72.0,1.0],
             "Wake":[115,117,77.0,2.0],             
             "SEBn":[108,110,47.0,3.0],
             "SEBs":[113,115,47.0,3.0],
             "SEBz":[110,113,47.0,3.0],
             "STrZ":[116,118,47.0,3.0],
             "SEZ":[92,94,62,5.0]}

        obslist=["20241027UTd","20241027UTe","20241027UTf","20241027UTg","20241027UTh"]
        CM=[62,62,62,62,62]
        
    if collection=="2024 GRS":
        ROI={"CORE":[111,115,62.0,3.0],
             "Chimney":[102,106,62,1.0],
             "Collar":[106,108,62.0,3.0],
             "Hook":[109,113,72.0,1.0],
             "SEBn":[108,110,47.0,3.0],
             "SEBs":[113,115,47.0,3.0],
             "SEBz":[110,113,47.0,3.0],
             "STrZ":[116,118,47.0,3.0],
             "Wake":[115,117,77.0,2.0],             
             "SEZ":[92,94,62,5.0]}

        obslist=["20241009UTc","20241027UTd","20241027UTf","20241027UTg"]
        CM=[62,62,62,62]
        
    if collection=="20241118 GRS":
        ROI={"CORE":[110,114,62.0,3.0],
             "Chimney":[100,104,62,1.0],
             "Collar":[105,107,62.0,3.0],
             "Hook":[109,113,72.0,1.0],
             "Wake":[115,117,77.0,2.0],             
             "SEBn":[108,110,47.0,3.0],
             "SEBs":[113,115,47.0,3.0],
             "SEBz":[110,113,47.0,3.0],
             "STrZ":[116,118,47.0,3.0],
             "SEZ":[92,94,62,5.0]}

        obslist=["20241118UTe","20241118UTf","20241118UTg",
                 "20241118UTh","20241118UTi"]
        CM=[62,62,62,62,62,62]

    if collection=="20241128 GRS":
        ROI={"CORE":[110,114,62.0,3.0],
             "Chimney":[100,104,62,1.0],
             "Collar":[105,107,62.0,3.0],
             "Hook":[109,113,72.0,1.0],
             "Wake":[115,117,77.0,2.0],             
             "SEBn":[108,110,47.0,3.0],
             "SEBs":[113,115,47.0,3.0],
             "SEBz":[110,113,47.0,3.0],
             "STrZ":[116,118,47.0,3.0],
             "SEZ":[92,94,62,5.0]}

        obslist=["20241128UTa","20241128UTb","20241128UTc",
                 #"20241128UTd",
                 "20241128UTe"]
        CM=[62,62,62,62,62,62]

    if collection=="20241203 GRS":
        ROI={"CORE":[110,114,62.0,3.0],
             "Chimney":[100,104,62,1.0],
             "Collar":[105,107,62.0,3.0],
             "Hook":[109,113,72.0,1.0],
             "Wake":[115,117,77.0,2.0],             
             "SEBn":[108,110,47.0,3.0],
             "SEBs":[113,115,47.0,3.0],
             "SEBz":[110,113,47.0,3.0],
             "STrZ":[116,118,47.0,3.0],
             "SEZ":[92,94,62,5.0]}

        obslist=["20241203UTa","20241203UTb","20241203UTc",
                 "20241203UTd","20241203UTe"]
        CM=[62,62,62,62,62,62]


    tempfNH3=np.zeros([180,360])
    tempPCloud=np.zeros([180,360])
    
    ###########################################################################
    # Loop over observations in obslist
    #  For each observation, return the mean and standard deviation for each
    #  ROI. Load these into arrays containing for each ROI.
    ###########################################################################
    First=True
    i=0
    for o in obslist:
        ROImod=copy.deepcopy(ROI)
        for key in ROI:
            ROImod[key][2]=ROI[key][2]+(CM[i]-CM[0])
        dateobs,roilabel,mean1,stdv1,mean2,stdv2=\
            L3JMP.L3_Jup_Map_Plot(obskey=o,imagetype='Map',target="Jupiter",
                        Smoothing=False,LatLims=[90,130],LonRng=20,CMpref=CM[i],
                        LonSys='2',showbands=False,coef=[0.75,0.45],
                        subproj='GRS',figxy=[8.0,4.0],FiveMicron=False,ROI=ROImod,
                        ctbls=ctbls)
            
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
        i=i+1
            
    ###############################################################################
    # Create two figures, one for time series of ROIs and the other PCloud vs fNH3
    ###############################################################################
    fig,axs=pl.subplots(1,2,figsize=(8,4), dpi=150, facecolor="white")
    #fig3.suptitle(suptitle,x=0.5,ha='center',color='k')
    fig.suptitle(collection+" Regions of Interest vs Time",x=0.5,ha='center',color='k')
    axs[0].set_title("Effective Cloud Top Pressure (mb)",fontsize=10)
    axs[1].set_title("Ammonia Mole Fraction (ppm)",fontsize=10)
    fig1,axs1=pl.subplots(1,2,figsize=(8,4), dpi=150, facecolor="white",sharey=True)
    #fig3.suptitle(suptitle,x=0.5,ha='center',color='k')
    fig1.suptitle(collection+" ROI Cloud Top Pressure vs Ammonia Abundance",x=0.5,ha='center',color='k')
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
    
    
    
    First=True
    for i in range(0,mean1array.shape[1]):
        if First:
            fNH3meanarr=[meanmean2[i]]
            PCldmeanarr=[meanmean1[i]]
        else:
            fNH3meanarr.append(meanmean2[i])
            PCldmeanarr.append(meanmean1[i])
        First=False
    fNH3meanarr=np.array(fNH3meanarr)
    PCldmeanarr=np.array(PCldmeanarr)
        
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
        
        
    axs[0].legend(ncols=4,fontsize=7,loc="upper center")
    axs1[0].legend(ncols=4,fontsize=7,loc="upper center")
    
    axs[0].tick_params(axis='x', which='major', labelsize=6)
    axs[0].tick_params(axis='y', which='major', labelsize=8)
    
    axs[0].xaxis.set_major_formatter(mdates.DateFormatter(dateformat))
    axs[0].xaxis.set_minor_formatter(mdates.DateFormatter(dateformat))
    axs[1].xaxis.set_major_formatter(mdates.DateFormatter(dateformat))
    axs[1].xaxis.set_minor_formatter(mdates.DateFormatter(dateformat))
    
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
    
    p=np.polyfit(fNH3meanarr[0:7], PCldmeanarr[0:7], 1)
    fit=p[0]*fNH3meanarr+p[1]
    axs1[1].plot(fNH3meanarr,fit)  
    
    #box = axs3[1].get_position()
    #axs3[1].set_position([box.x0+0.03, box.y0-0.01, box.width * 0.5, box.height * 1.015])    
    fig.subplots_adjust(left=0.08, right=0.97, top=0.85, bottom=0.12)
    fig1.subplots_adjust(left=0.08, right=0.97, top=0.85, bottom=0.12)
    #fig3.savefig(pathout+fnout,dpi=300)
    pathout='C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Studies/GRS/'
    fig.savefig(pathout+collection+" ROI Time Series",dpi=300)
    fig1.savefig(pathout+collection+" ROI Scatter",dpi=300)

    
    ###############################################################################
    #!!DO BLENDED MAP SEPARATELY USING MAKECONTIGUOUSMAP.PY
    ###############################################################################
    
    for key in ROI:
        ROImod[key][2]=ROI[key][2]+(np.mean(CM[:])-CM[0])

    MCM.MakeContiguousMap(False,False,False,collection=collection,LonSys='2',
                          FiveMicron=False,lats=[90,130],LonLims=[int(CM[0]-20),int(CM[0]+20)],
                          figsz=[3.0,6.0],ROI=ROImod,variance=True,proj="GRS",ctbls=['terrain_r','Blues'])

    print()
    print(fNH3meanarr,PCldmeanarr)
    print()
    print(p[0],p[1])
    
