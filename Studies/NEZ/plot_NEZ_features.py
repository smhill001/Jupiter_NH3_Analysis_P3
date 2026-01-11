def plot_NEZ_features(feature_type='NH3',axsnp=False):
    """
    Plots data and trends from feature-regions (NH3, cloud plumes, or NEDFs)
    that were identified by find_blob.py and supporting functions, usually 
    executed as part of L4_Jup_Map_Plot.py run in a batch mode for a given
    apparition. The inputs for this routine are read from curated
    versions of the find_blob.py outputs that have been manually inspected
    to tag persistent features over time as they drift in longitude.
    
    The current curated data sets are for the 2024-24 appariation and are
    HARD-CODED below.

    Parameters
    ----------
    feature_type : String, optional
        This is a string tag indicating what type of feature is to be analyzed,
        e.g., NH3, Plume, or NEDF
    axsnp : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    import matplotlib.pyplot as pl
    import numpy as np
    from datetime import datetime
    import matplotlib.dates as mdates
    import wind_speed_from_longitude_pair as ws
    import csv
    
    ###########################################################################
    # Set path for data files and dictionary for curated data files
    ###########################################################################
    pth="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Studies/NEZ/"
    config={"NH3":{"file":"2024-2025 Mean Sys1 0-360 15N-15S blobs fNH3 Tracker.csv",
                    "title":"Enhanced Ammonia Regions Time Series",
                    "subtitles":[['fNH3 Lat. W. Centroid','fNH3 Mean','Mean Cloud Pressure'],
                                     ['fNH3 W. Lon. Centroid','fNH3 Area','Minimum Cloud Pressure']],
                    "ylabels":[['PG Latitude (deg)','fNH3 (ppm)','Pressure (mb)'],
                                      ['Longitude (deg)','Area (sq. deg.)','Pressure (mb)']],
                    "columns":[[8,10,13],[9,4,14]]},
             "Plume":{"file":"2024-2025 Mean Sys1 0-360 15N-15S blobs Plume Tracker.csv",
                      "title":"Cloud Plume Time Series",
                    "subtitles":[['Cloud Lat. W. Centroid','Cloud Mean','Mean fNH3'],
                                     ['Cloud Lon. W. Centroid','Cloud Area','Minimum fNH3']],
                    "ylabels":[['PG Latitude (deg)','Pressure (mb)','fNH3 (ppm)'],
                                      ['Longitude (deg)','Area (sq. deg.)','fNH3 (ppm)']],
                    "columns":[[8,10,13],[9,4,14]]},   
             "NEDF":{"file":"2024-2025 Mean Sys1 0-360 15N-15S blobs NEDF Tracker.csv",
                      "title":"NEDF Time Series",
                    "subtitles":[['NEDF Lat. W. Centroid','NEDF Mean','Mean fNH3'],
                                     ['NEDF Lon. W. Centroid','NEDF Area','Minimum fNH3']],
                    "ylabels":[['PG Latitude (deg)','Pressure (mb)','fNH3 (ppm)'],
                                      ['Longitude (deg)','Area (sq. deg.)','fNH3 (ppm)']],
                    "columns":[[8,10,13],[9,4,14]]}}
    
    datetimestart=datetime.strptime("2024-10-01", "%Y-%m-%d")
    datetimeend=datetime.strptime("2025-03-15", "%Y-%m-%d")
    
    infile=config[feature_type]["file"]
    title=config[feature_type]["title"]
    ylabels=np.array(config[feature_type]["ylabels"])
    subtitles=np.array(config[feature_type]["subtitles"])
    columns=np.array(config[feature_type]["columns"])
    
    fig,axs=pl.subplots(2,3,figsize=(12,6), dpi=150, facecolor="white",sharex=True)
    fig.suptitle(title, fontsize=14)
    
    data_array = np.genfromtxt(pth+infile, delimiter=',',dtype=None)
    
    date_form = mdates.DateFormatter("%y-%m-%d")
    
    pathmapplots="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Studies/NEZ/"
    lon_speed_by_label = {}
    lat_speed_by_label = {}
    area_by_label = {}
    strength_by_label = {}

    avglat=[]
    avglon=[]
    avgarea=[]
    avgstrength=[]
    avgfNH3=[]
    avgPCld=[]
    stdfNH3=[]
    stdPCld=[]
    labels=[]
    
    for pl_col in range(0,2):  #!!! I have row and column nomenclature swapped in the indices pl_col and pl_row
        for pl_row in range(0,3):
            col=columns[pl_col,pl_row]
            subtitle=subtitles[pl_col,pl_row]
            for i in range(1,11):  
                filtered_data = data_array[data_array[:,0] == str(i).encode('ascii')]
                if len(filtered_data)>0:
                    tstr=[]
                    t=[]
                    tjd=[]
                    for j in range(0,len(filtered_data)):
                        tstr.append(filtered_data[j,3].decode("ascii")[:19])
                        t.append(datetime.strptime(filtered_data[j,3].decode("ascii")[:19], "%Y-%m-%dT%H:%M:%S"))
                        tjd.append(filtered_data[j,2])
                    #print(tjd)
                    y=np.array(filtered_data[:,col],dtype=float)
                    tjd=np.array(tjd,dtype=float)
                    #print(tjd)
                    axs[pl_col,pl_row].scatter(t,y,s=5,label="Reg. "+str(i))
                    
                    #!!! Dont really need the if statements here - could just hardwire to data columns
                    #if pl_row in [0,1,2]: 
                    if pl_row == 0 and pl_col==0:  #Latitude Plot
                        avglat.append(np.mean(y))
                        labels.append(i)
                        xnp=np.array(filtered_data[:,10],dtype=float)
                        ynp=np.array(filtered_data[:,13],dtype=float)
                        latcoefficients = np.polyfit(tjd, y, 1)
                        p = np.poly1d(latcoefficients)
                        axs[pl_col,pl_row].plot(t,p(tjd),linewidth=0.5)

                        if len(t)>1:

                            lat_speed_by_label[i] = {
                                'label':i,
                                'avglat': avglat[labels.index(i)],
                                'latdrift': latcoefficients[0],
                                'latcoef0':latcoefficients[0],
                                'latcoef1':latcoefficients[1],
                                'avgfNH3':np.mean(avgfNH3),
                                'avgPCld':np.mean(avgPCld),
                                'stdfNH3':np.std(avgfNH3),
                                'stdPCld':np.std(avgPCld)

                            }

                        if feature_type=='NH3': #Order is reversed in source files
                            axsnp.scatter(xnp,ynp,s=5)
                            avgfNH3.append(np.mean(xnp))
                            avgPCld.append(np.mean(ynp))
                            stdfNH3.append(np.std(xnp))
                            stdPCld.append(np.std(ynp))
                        else:
                            axsnp.scatter(ynp,xnp,s=5)
                            avgfNH3.append(np.mean(ynp))
                            avgPCld.append(np.mean(xnp))
                            stdfNH3.append(np.std(ynp))
                            stdPCld.append(np.std(xnp))

                            
                    if pl_row == 0 and pl_col==1:  #Longitude Plot
                        avglon.append(np.mean(y))
                        loncoefficients = np.polyfit(tjd, y, 1)
                        p = np.poly1d(loncoefficients)
                        axs[pl_col,pl_row].plot(t,p(tjd),linewidth=0.5)
                        
                        if len(t)>1:
                            
                            mps,dpd=ws.wind_speed_from_long_pair(avglat[labels.index(i)],p(tjd.min()),
                                                                 p(tjd.max()),
                                                                 tstr[0],
                                                                 tstr[-1],
                                                                 system='I')
                            lon_speed_by_label[i] = {
                                'label':i,
                                'avglat': avglat[labels.index(i)],
                                'avglon': avglat[labels.index(i)],
                                'sys3speed': mps,
                                'sys1drift': dpd,
                                'sys1coef0':loncoefficients[0],
                                'sys1coef1':loncoefficients[1],
                                'avgfNH3':np.mean(avgfNH3),
                                'avgPCld':np.mean(avgPCld),
                                'stdfNH3':np.std(avgfNH3),
                                'stdPCld':np.std(avgPCld)

                            }

                            print("############## i,avglat,mps,dpd=",i,avglat[labels.index(i)], mps,dpd)
                            
                    if pl_row == 1 and pl_col==1:  #Area
                        avgarea.append(np.mean(y))
                        areacoefficients = np.polyfit(tjd, y, 1)
                        p = np.poly1d(areacoefficients)
                        axs[pl_col,pl_row].plot(t,p(tjd),linewidth=0.5)
                        
                        if len(t)>1:
                            
                            area_by_label[i] = {
                                'label':i,
                                'avglat': avglat[labels.index(i)],
                                'avglon': avglat[labels.index(i)],
                                'avgarea':avgarea[labels.index(i)],
                                'areacoef0':areacoefficients[0],
                                'areacoef1':areacoefficients[1],
                                'avgfNH3':np.mean(avgfNH3),
                                'avgPCld':np.mean(avgPCld),
                                'stdfNH3':np.std(avgfNH3),
                                'stdPCld':np.std(avgPCld)
                            }
                            
                    if pl_row == 1 and pl_col==0:  #Mean Strength
                        avgstrength.append(np.mean(y))
                        avgstrengthcoefficients = np.polyfit(tjd, y, 1)
                        p = np.poly1d(avgstrengthcoefficients)
                        axs[pl_col,pl_row].plot(t,p(tjd),linewidth=0.5)
                        
                        if len(t)>1:
                            
                            strength_by_label[i] = {
                                'label':i,
                                'avglat': avglat[labels.index(i)],
                                'avglon': avglat[labels.index(i)],
                                'avgfNH3':np.mean(avgfNH3),
                                'avgstr':avgstrength[labels.index(i)],
                                'strcoef0':avgstrengthcoefficients[0],
                                'strcoef1':avgstrengthcoefficients[1],
                                'avgPCld':np.mean(avgPCld),
                                'stdfNH3':np.std(avgfNH3),
                                'stdPCld':np.std(avgPCld)

                            }
                
            axs[pl_col,pl_row].set_title(subtitle)
            axs[pl_col,pl_row].set_ylabel(ylabels[pl_col,pl_row],fontsize=10)
            axs[pl_col,pl_row].tick_params('x', labelsize=8)
            axs[pl_col,pl_row].xaxis.set_major_formatter(date_form)
            axs[pl_col,pl_row].xaxis.set_major_locator(mdates.MonthLocator(interval=1))
            axs[pl_col,pl_row].set_xlim(datetimestart,datetimeend)
            axs[pl_col,pl_row].grid(linewidth=0.2)

            #print(ylabels[pl_col,pl_row])
            if "mb" in ylabels[pl_col,pl_row]:
                axs[pl_col,pl_row].invert_yaxis()

    axs[0,0].legend(fontsize=8,ncols=2,loc='lower left')
    axs[1,1].legend(fontsize=8,ncols=2,loc='upper left')



    #axs[1,1].legend(bbox_to_anchor=(0.5,-0.15),loc='center',borderaxespad=0,ncols=11,fontsize=8)
    
    fig.subplots_adjust(left=0.05, bottom=0.08, right=0.98, top=0.90,
                wspace=0.25, hspace=0.15)     
    fig.savefig(pathmapplots+"2024-25 "+feature_type+" time series.png",dpi=300)

    lonfields = [
        'label',
        'avglat',
        'avglon',
        'sys3speed',
        'sys1drift',
        'sys1coef0',
        'sys1coef1',
        'avgfNH3',
        'avgPCld',
        'stdfNH3',
        'stdPCld'
    ]

    with open(pathmapplots+feature_type+"lon_speed_by_label.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=lonfields)
        writer.writeheader()
        for item in lon_speed_by_label.values():
            writer.writerow(item)
            
    latfields = [
        'label',
        'avglat',
        'latdrift',
        'latcoef0',
        'latcoef1',
        'avglon',
        'avgfNH3',
        'avgPCld',
        'stdfNH3',
        'stdPCld'
    ]

    with open(pathmapplots+feature_type+"lat_speed_by_label.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=latfields)
        writer.writeheader()
        for item in lat_speed_by_label.values():
            writer.writerow(item)
            
    areafields = [
        'label',
        'avglat',
        'avglon',
        'avgarea',
        'areacoef0',
        'areacoef1',
        'avgfNH3',
        'avgPCld',
        'stdfNH3',
        'stdPCld'
    ]

    with open(pathmapplots+feature_type+"area_by_label.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=areafields)
        writer.writeheader()
        for item in area_by_label.values():
            writer.writerow(item)

    strengthfields = [
        'label',
        'avglat',
        'avglon',
        'avgfNH3',
        'avgstr',
        'strcoef0',
        'strcoef1',
        'avgPCld',
        'stdfNH3',
        'stdPCld'
    ]

    with open(pathmapplots+feature_type+"strength_by_label.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=strengthfields)
        writer.writeheader()
        for item in strength_by_label.values():
            writer.writerow(item)
    
    return np.array([avgfNH3, avgPCld,stdfNH3,stdPCld])
            
def plot_ALL_NEZ_features():
    """
    Plots data and trends from ALL feature-regions (NH3, cloud plumes, and NEDFs)
    that were identified by find_blob.py and supporting functions, usually 
    executed as part of L4_Jup_Map_Plot.py run in a batch mode for a given
    apparition. The inputs for this routine are read from curated
    versions of the find_blob.py outputs that have been manually inspected
    to tag persistent features over time as they drift in longitude.
       
    Calls
        plot_NEZ_features
        ZonalWinds
    Returns
    -------
    None.

    """
    import sys
    import matplotlib.pyplot as pl
    sys.path.append('../../Winds/')
    import ZonalHSTWinds as ZHSTW
    pathmapplots="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Studies/NEZ/"  

    ###########################################################################
    # Set up plot for ammonia-pressure plane for individual feature points
    ###########################################################################
    fignp,axsnp=pl.subplots(1,figsize=(6,6), dpi=150, facecolor="white",sharex=True)
    fignp.suptitle("Ammonia - Pressure Plane", fontsize=14)
    
    ###########################################################################
    # Run plot_NEZ_features three times, once each for ammonia, clouds, and NEDFS
    ###########################################################################
    avgsNH3=plot_NEZ_features(feature_type='NH3',axsnp=axsnp)
    avgsPlume=plot_NEZ_features(feature_type='Plume',axsnp=axsnp)
    avgsNEDF=plot_NEZ_features(feature_type='NEDF',axsnp=axsnp)
    axsnp.set_ylim(1500.,2200.)
    axsnp.set_xlim(60.,180.)
    axsnp.grid(linewidth=0.2)
    axsnp.set_title("Individual Features 2024-25",fontsize=14)
    axsnp.set_xlabel("Ammonia Mole Fraction (ppm)",fontsize=14)
    axsnp.set_ylabel("Cloud Pressure (mb)",fontsize=14)
    axsnp.legend()
    axsnp.invert_yaxis()
    fignp.savefig(pathmapplots+"fNH3-PCld Plane Individual Obs.png",dpi=300)


    ###########################################################################
    # Set up and plot average values on the ammonia-pressure plane for features
    ###########################################################################
    fignpavg,axsnpavg=pl.subplots(1,figsize=(6,6), dpi=150, facecolor="white",sharex=True)
    fignp.suptitle("Ammonia - Pressure Plane", fontsize=14)

    axsnpavg.errorbar(avgsNH3[0,:],avgsNH3[1,:],xerr=avgsNH3[2,:],yerr=avgsNH3[3,:],
                      color='C2',label="NH3",linewidth=0.0,
                      elinewidth=1.0,marker="o",markersize=3)
    axsnpavg.errorbar(avgsPlume[0,:],avgsPlume[1,:],xerr=avgsPlume[2,:],yerr=avgsPlume[3,:],
                     color='C0',label="Plume",linewidth=0.0,
                     elinewidth=1.0,marker="o",markersize=3)
    axsnpavg.errorbar(avgsNEDF[0,:],avgsNEDF[1,:],xerr=avgsNEDF[2,:],yerr=avgsNEDF[3,:],
                     color='r',label="NEDF",linewidth=0.0,
                     elinewidth=1.0,marker="o",markersize=3)
    
    labels,NH3_indices,Plume_indices,NEDF_indices=ZHSTW.NEZ_feature_quiver_plot()
    
    print(labels)
    axsnpavg.set_ylim(1500.,2200.)
    axsnpavg.set_xlim(60.,180.)
    axsnpavg.scatter(avgsNH3[0,NH3_indices],avgsNH3[1,NH3_indices],s=50,color='C2')
    axsnpavg.scatter(avgsPlume[0,Plume_indices],avgsPlume[1,Plume_indices],s=50,color='C0')
    axsnpavg.scatter(avgsNEDF[0,NEDF_indices],avgsNEDF[1,NEDF_indices],s=50,color='r')
    axsnpavg.invert_yaxis() 
    axsnpavg.grid(linewidth=0.2)
    axsnpavg.set_title("Feature Averages 2024-25",fontsize=14)
    axsnpavg.set_xlabel("Ammonia Mole Fraction (ppm)",fontsize=14)
    axsnpavg.set_ylabel("Cloud Pressure (mb)",fontsize=14)
    axsnpavg.legend()

    fignpavg.savefig(pathmapplots+"fNH3-PCld Plane Average Obs.png",dpi=300)


    figw,axsw=ZHSTW.ZonalHSTWinds("2024c_f631n","2024d_f631n")
    ZHSTW.overplot_NEZ_feature_speeds(axsw)
    axsw.legend()
    pathmapplots="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Studies/NEZ/"  
    
    figw.savefig(pathmapplots+"WindProfile.png",dpi=300)
