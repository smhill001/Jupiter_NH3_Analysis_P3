def nez_analysis():
    """
    Created on Sat Dec 30 21:50:01 2023
    
    @author: smhil
    """
    
    import csv
    import numpy as np
    import pylab as pl
    
    path="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Studies/NEZ/"
    
    ###########################################################################
    # READ IN HILL MEASUREMENTS (DATE,LON,LAT) OF NH3 ENHANCEMENTS AND
    # CH4 DEEP CLOUDS.
    ###########################################################################
    fn="Hill_ammonia_mea.CSV"
    #CH4 deep cloud locations
    lsCH4obs=[]
    lsCH4lon=[]
    lsCH4lat=[]
    #NH3 enhancement locations
    lsNH3obs=[]
    lsNH3lon=[]
    lsNH3lat=[]
    
    with open(path+fn) as csv_file:
        csv_reader = csv.DictReader(csv_file)
        line_count = 0
        for row in csv_reader:
            if line_count == 0:
                #print(f'Column names are {", ".join(row)}')
                line_count += 1
            #print(f'\t{row["REC_NO"]} works in the {row["OBJECT"]} department, and was born in {row["REGION"]}.')
            if row["OBJECT"]=="DC1_SPOT":
                lsCH4obs.append(row["DATE"])
                lsCH4lon.append(float(row["LONG_1"]))
                lsCH4lat.append(float(row["LAT"]))
            if row["OBJECT"]=="WC1_SPOT":
                lsNH3obs.append(row["DATE"])
                lsNH3lon.append(float(row["LONG_1"]))
                lsNH3lat.append(float(row["LAT"]))
            line_count += 1
    print(f'Processed {line_count} lines.')
    
    ###########################################################################
    # READ IN HILL MEASUREMENTS (DATE,LON,LAT) OF NH3 ENHANCEMENTS AND
    # CH4 DEEP CLOUDS.
    ###########################################################################
    
    #WINJupos measurer PROJection locations
    lsPRJobs=[]
    lsPRJlon=[]
    lsPRJlat=[]
    lsPRJdlon=[]
    lsPRJlonsel=[]
    lsPRJdlonsel=[]
    lsPRJlatsel=[]
    
    fn="DC_PROJ_2022-23_sel.CSV"
    with open(path+fn) as csv_file:
        csv_reader = csv.DictReader(csv_file)
        line_count = 0
        for row in csv_reader:
            if line_count == 0:
                #print(f'Column names are {", ".join(row)}')
                line_count += 1
            #print(f'\t{row["REC_NO"]} works in the {row["OBJECT"]} department, and was born in {row["REGION"]}.')
            if row["OBJECT"]=="DC1_PROJ":
                lsPRJobs.append(row["DATE"])
                lsPRJlon.append(float(row["LONG_1"]))
                lsPRJlat.append(float(row["LAT"]))
            line_count += 1
    print(f'Processed {line_count} lines.')
    #print(lsPRJobs)
    
    ###########################################################################
    # LOOP OVER NH3 ENHANCEMENT MEASUREMENTS TO FIND *ALL* CH4 DEEP CLOUD
    # MEASUREMENTS MADE ON THE SAME UT DATE. OF THOSE DEEP CLOUD MEASUREMENTS,
    # FIND THE *ONE* THAT IS THE CLOSEST IN LONGITUDE TO THE SELECTED
    # NH3 ENHANCEMENT.
    ###########################################################################
    lsNH3dobs=[]    #list of NH3 obs dates with corresponding CH4 obs
    lsNH3dlon=[]    #list of corresponding delta longidutes (NH3lon-CH4lon)
    lsNH3dlat=[]    #list of corresponding delta latitudes (NH3lat-CH4lat)
    lsCH4latSel=[]  #list of CH4 latitudes that have corresponding NH3
                    #enhancements (for determining mean CH4 latitude)
                    
    for j in range(0,len(lsNH3obs)):
        #Get set of indices for the CH4 records with the same UT dates 
        #  (not times) as the current (jth) NH3 enhancement. With multiple
        #  NH3 enhancements on a given date, the set of 'matching' CH4
        #  records can be repeated multiple times.
        indices=[i for i, e in enumerate(lsCH4obs) if e == lsNH3obs[j]]
        #print(lsCH4obs[j],indices[:])
        #compute the closest CH4 DR in longitude and the index of that value 
        #closest=np.min(np.abs(lsNH3lon[j]-np.take(lsCH4lon,indices)))
        closesti=np.argmin(np.abs(lsNH3lon[j]-np.take(lsCH4lon,indices)))
        print(closesti,np.take(lsCH4lon,indices[closesti]))
        closest=lsNH3lon[j]-np.take(lsCH4lon,indices[closesti])
        #Filter for CH4 DFs within 10 deg lon of the NH3 enhancement and compute
        #proper sign for dlons that are over 180 deg
        if np.abs(closest)<10:
            if closest<180:
                lsNH3dlon.append(closest)
            if closest>=180:
                lsNH3dlon.append(closest-360.)
            lsNH3dobs.append(lsCH4obs[j])
            lsNH3dlat.append(lsNH3lat[j]-np.take(lsCH4lat,indices[closesti]))
            #lsNH3dlat.append(lsNH3lat[j])#
            #Load the latitude for the selected CH4 DFs
            lsCH4latSel.append(np.take(lsCH4lat,indices[closesti]))
            #print(lsNH3lon[j],np.take(lsCH4lon,indices),closest,lsNH3dlon)
            
        #REPEAT FOR WINJUPOS PROJs
        indicesPRJ=[i for i, e in enumerate(lsPRJobs) if e == lsNH3obs[j]]
        #print(j,lsNH3obs[j])
        #if isinstance(indicesPRJ, list): 
            #print(lsPRJobs)
        #compute the closest CH4 DR in longitude and the index of that value 
        #closest=np.min(np.abs(lsNH3lon[j]-np.take(lsCH4lon,indices)))
                            
        lsPRJlonsel.append(np.take(lsPRJlon,indicesPRJ))
        lsPRJdlonsel.append(np.take(lsPRJlon,indicesPRJ)-\
                            np.take(lsCH4lon,indices[closesti]))
            
        lsPRJdlon.append(np.array(lsPRJlon)-\
                            np.take(lsCH4lon,indices[closesti]))
        lsPRJlatsel.append(np.take(lsPRJlat,indicesPRJ))
    
    
    
    lsPRJlatsel=np.array(np.concatenate(lsPRJlatsel, axis=0 ))
    lsPRJlonsel=np.array(np.concatenate(lsPRJlonsel, axis=0 ))
    lsPRJdlonsel=np.array(np.concatenate(lsPRJdlonsel, axis=0 ))
    
    #lsPRJdlon=lsPRJdlon(np.where(np.abs(lsPRJdlon)<=10.))
    #lsPRJlatsel=lsPRJdlon(np.where(np.abs(lsPRJdlon)<=10.))
    
    #lsPRJdlon=np.concatenate( lsPRJdlon, axis=0 )
    print(np.mean(lsNH3dlon),np.std(lsNH3dlon))
    print(np.mean(lsNH3dlat),np.std(lsNH3dlat))
    print(np.mean(lsPRJdlonsel),np.std(lsPRJdlonsel))
    print(np.mean(lsPRJlatsel),np.std(lsPRJlatsel))
    ###############################################################################
    # CREATE FIGURES
    ###############################################################################        
    fig1,axs1=pl.subplots(2,1,figsize=(6.0,6.0), dpi=150, facecolor="white")
    fig1.suptitle("Relative Positions of NH3 Enhancement and 'Hot Spots'")
    axs1[0].set_title("From TBD to TBD")
    axs1[0].scatter(np.array(lsCH4lon),np.array(lsCH4lat),color='C0',s=5)
    axs1[0].scatter(np.array(lsNH3lon),np.array(lsNH3lat),color='C1',s=5)
    print(np.array(lsPRJlon).shape,np.array(lsPRJlat).shape)
    axs1[0].scatter(np.array(lsPRJlon),np.array(lsPRJlat),color='0.5',s=0.5)
    axs1[0].scatter(np.array(lsPRJlonsel),np.array(lsPRJlatsel),color='k',s=2)
    axs1[0].set_xlim(0,360)
    axs1[0].set_ylim(0,10)
    axs1[0].set_xlabel("System 1 Longitude (deg)")
    axs1[0].set_ylabel("PG Latitude (deg)")
    #print(np.array(lsNH3dlon).shape,np.array(lsNH3dlat).shape)
    
    axs1[1].scatter(np.array([0.0]),np.mean(np.array(lsCH4latSel)),color='C0',s=50,
                    alpha=0.5)
    axs1[1].scatter(np.mean(np.array(lsNH3dlon)),np.mean(np.array(lsNH3lat)),color='C1',s=50,
                    alpha=0.5)
    axs1[1].scatter(np.array(lsNH3dlon),np.array(lsNH3dlat)+np.mean(np.array(lsCH4latSel)),color='C1',s=5)
    #axs1[1].scatter(np.array(lsPRJdlonsel),np.array(lsPRJlatsel),color='k',s=2)
    print(np.array(lsPRJdlon).shape,np.array(lsPRJlat).shape)
    
    lsPRJdlonsela=lsPRJdlonsel[np.where(np.abs(lsPRJdlonsel)<=10.)]
    lsPRJlatsela=lsPRJlatsel[np.where(np.abs(lsPRJdlonsel)<=10.)]
    #axs1[1].scatter(np.array(lsPRJdlon),np.array(lsPRJlat),color='k',s=0.5)
    axs1[1].scatter(np.array(lsPRJdlonsela),np.array(lsPRJlatsela),color='k',s=2)
    axs1[1].scatter(np.mean(lsPRJdlonsela),np.mean(lsPRJlatsela),color='k',s=50,
                    alpha=0.5)
    axs1[1].set_xlim(-15,15)
    axs1[1].invert_xaxis()
    axs1[1].set_ylim(0,10)
    axs1[1].set_xlabel("System 1 Longitude Offset (deg)")
    axs1[1].set_ylabel("PG Latitude (deg)")
    
    fig1.savefig(path+'DF_Analysis.png',dpi=300)
