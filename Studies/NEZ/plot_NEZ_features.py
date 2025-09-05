def plot_NEZ_features(feature_type='NH3',axsnp=False):
    import matplotlib.pyplot as pl
    import numpy as np
    from datetime import datetime
    import matplotlib.dates as mdates
    import wind_speed_from_longitude_pair as ws
    import csv

    
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
    speed_by_label = {}

    avglat=[]
    avglon=[]
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
                        coefficients = np.polyfit(tjd, y, 1)
                        p = np.poly1d(coefficients)
                        axs[pl_col,pl_row].plot(t,p(tjd),linewidth=0.5)
                        
                        if len(t)>1:
                            
                            mps,dpd=ws.wind_speed_from_long_pair(avglat[labels.index(i)],p(tjd.min()),
                                                                 p(tjd.max()),
                                                                 tstr[0],
                                                                 tstr[-1],
                                                                 system='I')
                            speed_by_label[i] = {
                                'label':i,
                                'avglat': avglat[labels.index(i)],
                                'avglon': avglon[labels.index(i)],
                                'sys3speed': mps,
                                'drift': dpd,
                                'coef0':coefficients[0],
                                'coef1':coefficients[1],
                                'avgfNH3':np.mean(avgfNH3),
                                'avgPCld':np.mean(avgPCld),
                                'stdfNH3':np.std(avgfNH3),
                                'stdPCld':np.std(avgPCld)

                            }

                            print("############## i,avglat,mps,dpd=",i,avglat[labels.index(i)], mps,dpd)
                
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
    
    axs[1,1].legend(bbox_to_anchor=(0.5,-0.15),loc='center',borderaxespad=0,ncols=11,fontsize=8)
    fig.subplots_adjust(left=0.05, bottom=0.08, right=0.98, top=0.90,
                wspace=0.25, hspace=0.15)     
    fig.savefig(pathmapplots+"2024-25 "+feature_type+" time series.png",dpi=300)

    fields = [
        'label',
        'avglat',
        'avglon',
        'sys3speed',
        'drift',
        'coef0',
        'coef1',
        'avgfNH3',
        'avgPCld',
        'stdfNH3',
        'stdPCld'
    ]

    with open(pathmapplots+feature_type+"speed_by_label.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for item in speed_by_label.values():
            writer.writerow(item)
    
    return np.array([avgfNH3, avgPCld,stdfNH3,stdPCld])
            
def plot_ALL_NEZ_features():
    
    import matplotlib.pyplot as pl
    import ZonalWinds as ZW
    fignp,axsnp=pl.subplots(1,figsize=(6,6), dpi=150, facecolor="white",sharex=True)
    fignp.suptitle("N-P Plane", fontsize=14)
    
    avgsNH3=plot_NEZ_features(feature_type='NH3',axsnp=axsnp)
    avgsPlume=plot_NEZ_features(feature_type='Plume',axsnp=axsnp)
    avgsNEDF=plot_NEZ_features(feature_type='NEDF',axsnp=axsnp)
    
    axsnp.invert_yaxis()

    fignpavg,axsnpavg=pl.subplots(1,figsize=(6,6), dpi=150, facecolor="white",sharex=True)
    #axsnpavg.set_title("NH3-PCloud Plane", fontsize=14)
    axsnpavg.errorbar(avgsNH3[0,:],avgsNH3[1,:],xerr=avgsNH3[2,:],yerr=avgsNH3[3,:],
                      color='C2',label="NH3",linewidth=0.0,
                      elinewidth=1.0,marker="o",markersize=3)
    axsnpavg.errorbar(avgsPlume[0,:],avgsPlume[1,:],xerr=avgsPlume[2,:],yerr=avgsPlume[3,:],
                     color='C0',label="Plume",linewidth=0.0,
                     elinewidth=1.0,marker="o",markersize=3)
    axsnpavg.errorbar(avgsNEDF[0,:],avgsNEDF[1,:],xerr=avgsNEDF[2,:],yerr=avgsNEDF[3,:],
                     color='r',label="NEDF",linewidth=0.0,
                     elinewidth=1.0,marker="o",markersize=3)
    
    labels,NH3_indices,Plume_indices,NEDF_indices=ZW.NEZ_feature_quiver_plot()
    print(labels)
    axsnpavg.set_ylim(1500.,2200.)
    axsnpavg.set_xlim(60.,180.)
    axsnpavg.scatter(avgsNH3[0,NH3_indices],avgsNH3[1,NH3_indices],s=50,color='C2')
    axsnpavg.scatter(avgsPlume[0,Plume_indices],avgsPlume[1,Plume_indices],s=50,color='C0')
    axsnpavg.scatter(avgsNEDF[0,NEDF_indices],avgsNEDF[1,NEDF_indices],s=50,color='r')
    axsnpavg.invert_yaxis() 
    axsnpavg.grid(linewidth=0.2)
    axsnpavg.set_xlabel("Ammonia Mole Fraction (ppm)",fontsize=14)
    axsnpavg.set_ylabel("Cloud Pressure (mb)",fontsize=14)

    axsnpavg.legend()

    pathmapplots="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Studies/NEZ/"  
    fignpavg.savefig(pathmapplots+"ThermoCycle.png",dpi=300)

