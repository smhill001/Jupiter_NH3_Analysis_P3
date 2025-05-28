def ROI_Script(collection="20241202-20241202 NEDF 340",dateformat="%b:%d",
                     ctbls=["terrain_r","Blues"],
                     LatLims=[75,105],LonRng=20,LonSys='1',close=False,coefs=[0.75,0.45],
                     localmax=False,writecatalog=False):

    import sys
    drive='c:'
    sys.path.append(drive+'/Astronomy/Python Play')
    sys.path.append(drive+'/Astronomy/Python Play/Util_P3')
    sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy')
    sys.path.append('./Services')
    sys.path.append('../Studies')
    sys.path.append('../Studies/NEZ')
    
    import os
    import pylab as pl
    import numpy as np
    import matplotlib.dates as mdates
    from datetime import datetime

    import L3_Jup_Map_Plot as L3JMP
    import MakeContiguousMap as MCM
    import NEDF_ROI_collections as NRC

    ###########################################################################
    # GRAB OBSERVATIONS (KEY) LISTS AND THEN SET ROI BOUNDARIES  AND CM
    # SEQUENCE IF NEEDED TO TRACK A DRIFING TARGET (SYS I)
    ###########################################################################   
    ROI,obskeys,CM=NRC.NEDF_ROI_collections(collection=collection,writecatalog=writecatalog)
   
    ###########################################################################
    # SET MIN AND MAX VALUES FOR PLOT AXES AND COLOR BARS
    ###########################################################################
    if ctbls[0]=="jet":
        fNH3low=60
        fNH3high=160
        PCldlow=1200
        PCldhigh=2000
    elif ctbls[0]=="terrain_r":
        if coefs[0] != 0.0:
            fNH3low=60
            fNH3high=180
            PCldlow=1400
            PCldhigh=2400
        else:
            fNH3low=40
            fNH3high=160
            PCldlow=1000
            PCldhigh=2200          

    micronlow=0.5
    micronhigh=3.5

    ###########################################################################
    # LOOP OVER OBSERVATION KEYS AND CREATE INDIVIDUAL MAP SETS AND STATISTICS.
    # CLOSE ALL PLOTS IF REQUESTED TO REDUCE SCREEN CLUTTER
    #
    # !!!!NEED TO HAVE INPUTS FOR THE LIMB CORRECTION COEFFICIENC VALUES 
    # !!!!AND VET THEM
    ###########################################################################
    First=True

    for o in obskeys:
        dateobs,roilabel,mean1,stdv1,mean2,stdv2,meanamf=\
            L3JMP.L3_Jup_Map_Plot(obskey=o,imagetype='Map',target="Jupiter",
                        Smoothing=False,LatLims=LatLims,LonRng=LonRng,CMpref=CM[0],
                        LonSys=LonSys,showbands=False,coef=coefs,
                        subproj='NEZ',figxy=[8.0,4.0],FiveMicron=False,ROI=ROI,
                        ctbls=ctbls)
            
        print("*********",dateobs,roilabel,mean1,stdv1,mean2,stdv2)
        if First:
            datearr=[dateobs]
            roilabelarr=[roilabel]
            mean1arr=[mean1]
            stdv1arr=[stdv1]
            mean2arr=[mean2]
            stdv2arr=[stdv2]
            meanamfarr=[meanamf]
        else:
            datearr.append(dateobs)
            roilabelarr.append(roilabel)
            mean1arr.append(mean1)
            stdv1arr.append(stdv1)
            mean2arr.append(mean2)
            stdv2arr.append(stdv2)
            meanamfarr.append(meanamf)
        First=False
            
    if close:
        pl.close('all')

    ###########################################################################
    # CONVERT DATA TO NUMPY ARRAYS FOR PLOTTING
    ###########################################################################
    for k in range(0,len(datearr)):    
        datearr[k]=datetime.strptime(datearr[k], '%Y-%m-%dT%H:%M:%SZ')
    datearray=np.array(datearr)
    
    mean1array=np.array(mean1arr)
    mean2array=np.array(mean2arr)
    meanamfarray=np.array(meanamfarr)
    stdv1array=np.array(stdv1arr)
    stdv2array=np.array(stdv2arr)
    roilabelarray=np.array(roilabelarr)

    meanmean1=np.mean(mean1array,axis=0)
    stdvstdv1=np.std(mean1array,axis=0)
    meanmean2=np.mean(mean2array,axis=0)
    stdvstdv2=np.std(mean2array,axis=0)
    
    ###########################################################################
    # SET UP PLOTS FOR TIMELINE AND SCATTER
    ###########################################################################
    pathout='C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Studies/NEZ/'

    figtime,axstime=pl.subplots(1,2,figsize=(8,4), dpi=150, facecolor="white")
    figtime.suptitle(collection+" ROIs vs Time",x=0.5,ha='center',color='k')
    figscat,axsscat=pl.subplots(1,2,figsize=(8,4), dpi=150, facecolor="white",sharey=True)
    figscat.suptitle(collection+" ROI Cloud Top Pressure vs Ammonia Abundance",x=0.5,ha='center',color='k')
    figamf,axsamf=pl.subplots(1,2,figsize=(8,4), dpi=150, facecolor="white")
    figamf.suptitle(collection+" AMF",x=0.5,ha='center',color='k')

    ###########################################################################
    # PLOT DATA
    for i in range(0,mean1array.shape[1]):
        date_sort = datearray.argsort()
        amf_sort = meanamfarray[:,i].argsort()
        lbl_sort = roilabelarray[:,i].argsort()
        
        clr='C0'
        if "NH3" in roilabelarray[0,i]:
            clr='k'
        if "Plume" in roilabelarray[0,i]:
            clr='r'
        if "NEDF" in roilabelarray[0,i]:
            clr='y'
        if "NEB ref" in roilabelarray[0,i]:
            clr='brown'

        
        axstime[0].errorbar(datearray[date_sort],mean1array[date_sort,i],yerr=stdv1array[date_sort,i],
                            label=roilabelarray[0,i],
                        capsize=3.0,elinewidth=0.5,capthick=0.5,color=clr)
        axstime[1].errorbar(datearray[date_sort],mean2array[date_sort,i],yerr=stdv2array[date_sort,i],
                            label=roilabelarray[0,i],
                        capsize=3.0,elinewidth=0.5,capthick=0.5,color=clr)

        axsscat[0].errorbar(mean2array[:,i],mean1array[:,i],xerr=stdv2array[:,i],yerr=stdv1array[:,i],
                         label=roilabelarray[0,i],linewidth=0.0,
                         elinewidth=0.5,marker="o",markersize=3.0,color=clr)
        axsscat[1].errorbar(meanmean2[i],meanmean1[i],xerr=stdvstdv2[i],yerr=stdvstdv1[i],
                         label=roilabelarray[0,i],
                         elinewidth=0.5,marker="o",markersize=3.0,color=clr)
    
        axsamf[0].errorbar(meanamfarray[amf_sort,i],mean1array[amf_sort,i],yerr=stdv1array[amf_sort,i],
                         label=roilabelarray[0,i],linewidth=1.0,
                         elinewidth=0.5,marker="o",markersize=3.0,color=clr)
        axsamf[1].errorbar(meanamfarray[amf_sort,i],mean2array[amf_sort,i],yerr=stdv2array[amf_sort,i],
                         label=roilabelarray[0,i],linewidth=1.0,
                         elinewidth=0.5,marker="o",markersize=3.0,color=clr)
    
    ###########################################################################
    # FORMAT PLOTS
    for i in range(0,2):
        axstime[i].tick_params(axis='both', which='major', labelsize=8)
        axstime[i].grid(linewidth=0.2)
        axstime[i].xaxis.set_major_formatter(mdates.DateFormatter(dateformat))
        axstime[i].xaxis.set_minor_formatter(mdates.DateFormatter(dateformat))
        axstime[i].set_xlabel("Date-Time (UT)",fontsize=10)     

        axsscat[i].tick_params(axis='both', which='major', labelsize=8)
        axsscat[i].grid(linewidth=0.2)
        axsscat[i].set_ylim(PCldlow,PCldhigh)
        axsscat[i].set_xlim(fNH3low,fNH3high)
        axsscat[i].invert_yaxis()
       
        axsamf[i].tick_params(axis='both', which='major', labelsize=8)
        axsamf[i].set_xlim(1.0,2.0)
        axsamf[i].grid(linewidth=0.2)
        axsamf[i].set_xlabel("Air Mass Factor (one-way)",fontsize=10)     

    axstime[0].legend(ncols=2,fontsize=8)

    axstime[0].set_ylabel("Cloud Top Pressure (mb)",fontsize=10)  
    axstime[0].set_ylim(PCldlow,PCldhigh)
    axstime[0].invert_yaxis()
    axstime[0].set_title("Effective Cloud Top Pressure (mb)",fontsize=10)
            
    axstime[1].set_ylabel("Ammonia Mole Fraction (ppm)",fontsize=10)     
    axstime[1].set_ylim(fNH3low,fNH3high)   
    axstime[1].set_title("Ammonia Mole Fraction (ppm)",fontsize=10)
    
    figtime.subplots_adjust(left=0.08, right=0.97, top=0.85, bottom=0.12)
    figtime.savefig(pathout+collection+" ROI Time Series",dpi=300)

    axsscat[0].legend(ncols=2,fontsize=8)

    axsscat[0].set_ylabel("Cloud Top Pressure (mb)",fontsize=10)        
    axsscat[0].set_xlabel("Ammonia Mole Fraction (ppm)",fontsize=10) 
    axsscat[0].set_title("Individual ROI Sample Averages",fontsize=10)
    
    axsscat[1].set_xlabel("Ammonia Mole Fraction (ppm)",fontsize=10)    
    axsscat[1].set_title("ROIs Averaged over Time",fontsize=10)
    
    figscat.subplots_adjust(left=0.08, right=0.97, top=0.85, bottom=0.12)
    figscat.savefig(pathout+collection+" ROI Scatter",dpi=300)
    
    axsamf[0].legend(ncols=2,fontsize=8)

    axsamf[0].set_ylabel("Cloud Top Pressure (mb)",fontsize=10)  
    axsamf[0].set_ylim(PCldlow,PCldhigh)
    axsamf[0].invert_yaxis()
    axsamf[0].set_title("Effective Cloud Top Pressure (mb)",fontsize=10)
            
    axsamf[1].set_ylabel("Ammonia Mole Fraction (ppm)",fontsize=10)     
    axsamf[1].set_ylim(fNH3low,fNH3high)   
    axsamf[1].set_title("Ammonia Mole Fraction (ppm)",fontsize=10)

    figamf.subplots_adjust(left=0.08, right=0.97, top=0.85, bottom=0.12)
    figamf.savefig(pathout+collection+" AMF Scatter",dpi=300)
    

    """ #hold for mods if we want to do 5 micron stuff
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

    ###############################################################################
    # MAKE BLENDED MAP
    ###############################################################################
    print()
    print("##################")
    print("collection,[int(CM[0]-45),int(CM[0]+45)]=",collection,[int(CM[0]-45),int(CM[0]+45)])
    print("##################")
    print()

    MCM.MakeContiguousMap(False,False,False,obskeys,collection=collection,LonSys=LonSys,
                          FiveMicron=False,lats=LatLims,LonLims=[int(CM[0]-LonRng),int(CM[0]+LonRng)],
                          figsz=[3.0,6.0],ROI=ROI,variance=False,localmax=localmax,
                          proj="NEZ",ctbls=['terrain_r','Blues'])
