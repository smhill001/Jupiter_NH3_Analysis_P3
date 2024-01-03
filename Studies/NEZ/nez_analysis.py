def nez_analysis():
    """
    Created on Sat Dec 30 21:50:01 2023
    @author: smhil
    
    PURPOSE: TO ANALYZE THE RELATIVE POSITIONS IN THE NEZ OF OBSERVED
             NH3 ENHANCEMENTS AND CH4 DEEP CLOUDS. IN MY DATA, IT APPEARS
             THAT THESE ARE OFTEN (ALWAYS?) CLOSELY COUPLED IN AND AROUND
             CLASSIC VISIBLE LIGHT DARK FORMATIONS (DARK PROJECTIONS IN
             WINJUPOS). SO A WORKING ASSUMPTION IS THAT THERE IS A ONE-TO-ONE
             CORRESPONDENCE.
    """
    
    import csv
    import numpy as np
    import pylab as pl
    
    path="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Studies/NEZ/"
    
    pairs={} #original bland dictionary
    # NH3 rec, NH3 date, CH4 rec, CH4 date ...
    ###########################################################################
    # READ IN HILL MEASUREMENTS (DATE,LON,LAT) OF NH3 ENHANCEMENTS AND
    # CH4 DEEP CLOUDS.
    ###########################################################################
    fn="Hill_ammonia_mea.CSV"
    #CH4 deep cloud locations
    lsCH4rec=[]
    lsCH4obs=[]
    lsCH4lon=[]
    lsCH4lat=[]
    #NH3 enhancement locations
    lsNH3rec=[]
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
                lsCH4rec.append(row["REC_NO"])
                lsCH4obs.append(row["DATE"])
                lsCH4lon.append(float(row["LONG_1"]))
                lsCH4lat.append(float(row["LAT"]))
            if row["OBJECT"]=="WC1_SPOT":
                lsNH3rec.append(row["REC_NO"])
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
    lsPRJrec=[]
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
                print(f'Column names are {", ".join(row)}')
                line_count += 1
            if "PROJ" in row["OBJECT"]: #=="DC1_PROJ":
                lsPRJrec.append(row["REC_NO"])
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
        #print("j,lsNH3obs[j]=",j,lsNH3obs[j],lsNH3lon[j])
        indices=[i for i, e in enumerate(lsCH4obs) if e == lsNH3obs[j]]
        #print(lsNH33obs[j],indices[:],np.take(lsCH4lon,indices),)
        #compute the closest CH4 DR in longitude and the index of that value 
        closesti=np.argmin(np.abs(lsNH3lon[j]-np.take(lsCH4lon,indices)))
        #print(closesti,np.take(lsCH4lon,indices[closesti]))
        CH4darklon=np.take(lsCH4lon,indices[closesti])
        dlonmin=lsNH3lon[j]-np.take(lsCH4lon,indices[closesti])
        #print(dlonmin)
        #Filter for CH4 DFs within 10 deg lon of the NH3 enhancement and compute
        #proper sign for dlons that are over 180 deg
        if np.abs(dlonmin)<15:
            if dlonmin<180:
                lsNH3dlon.append(dlonmin)
            if dlonmin>=180:
                lsNH3dlon.append(dlonmin-360.)
            #lsNH3dobs.append(lsNH3obs[j])
            #Save the differential NH3-CH4 latitude along with the latitude
            # of the CH4 deep cloud spot (for later averaging)
            lsNH3dlat.append(lsNH3lat[j]-np.take(lsCH4lat,indices[closesti]))
            lsCH4latSel.append(np.take(lsCH4lat,indices[closesti]))

            #print("match")
            #print("NH3obs,NH3lon,CH4lon,dlon=",lsNH3obs[j],lsNH3lon[j],np.take(lsCH4lon,indices[closesti]),dlonmin)
            
        #print()
#REPEAT FOR WINJUPOS PROJs
        indicesPRJ=[i for i, e in enumerate(lsPRJobs) if e == lsNH3obs[j]]
        #print(j,lsNH3obs[j])
        #if isinstance(indicesPRJ, list): 
            #print(lsPRJobs)
        #compute the closest CH4 DR in longitude and the index of that value 
        #closest=np.min(np.abs(lsNH3lon[j]-np.take(lsCH4lon,indices)))
        #print("j,lsNH3obs[j]=",j,lsNH3obs[j],lsNH3lon[j])
        #print(indicesPRJ[:],np.take(lsPRJobs,indicesPRJ),
        #      np.take(lsPRJlon,indicesPRJ))
        # Load all Sys 1 longitude for PROJs on the given NH3obs date
        lsPRJlonsel.append(np.take(lsPRJlon,indicesPRJ))
        # Load all differential longitudes for PROJs on the given NH3obs date
        lsPRJdlonsel.append(np.take(lsPRJlon,indicesPRJ)-CH4darklon)
                     # Load all PG latitudes for PROJs on the given NH3obs day
        lsPRJlatsel.append(np.take(lsPRJlat,indicesPRJ))
    
    lsPRJlatsel=np.array(np.concatenate(lsPRJlatsel, axis=0 ))
    lsPRJlonsel=np.array(np.concatenate(lsPRJlonsel, axis=0 ))
    lsPRJdlonsel=np.array(np.concatenate(lsPRJdlonsel, axis=0 ))
    
    #lsPRJdlon=lsPRJdlon(np.where(np.abs(lsPRJdlon)<=10.))
    #lsPRJlatsel=lsPRJdlon(np.where(np.abs(lsPRJdlon)<=10.))
    
    #lsPRJdlon=np.concatenate( lsPRJdlon, axis=0 )
    print()
    print(np.mean(lsNH3dlon),np.std(lsNH3dlon))
    print(np.mean(lsNH3dlat),np.std(lsNH3dlat))
    #print(np.mean(lsPRJdlonsel),np.std(lsPRJdlonsel))
    #print(np.mean(lsPRJlatsel),np.std(lsPRJlatsel))
    ###############################################################################
    # CREATE Plot 0: PG Latitude vs System 1 Longitude
    ###############################################################################        
    fig1,axs1=pl.subplots(2,1,figsize=(6.0,6.0), dpi=150, facecolor="white")
    fig1.suptitle("Relative Positions of NH3 Enhancement and 'Hot Spots'")
    axs1[0].set_title('From '+str(lsNH3obs[0])+' to '+str(lsNH3obs[len(lsNH3obs)-1]))
    axs1[0].scatter(np.array(lsCH4lon),np.array(lsCH4lat),color='C0',s=5)
    axs1[0].scatter(np.array(lsNH3lon),np.array(lsNH3lat),color='C1',s=5)
    #print(np.array(lsPRJlon).shape,np.array(lsPRJlat).shape)
    #axs1[0].scatter(np.array(lsPRJlon),np.array(lsPRJlat),color='0.5',s=0.5)
    axs1[0].scatter(np.array(lsPRJlonsel),np.array(lsPRJlatsel),color='k',s=2)
    axs1[0].set_xlim(0,360)
    axs1[0].set_ylim(0,10)
    axs1[0].set_xlabel("System 1 Longitude (deg)")
    axs1[0].set_ylabel("PG Latitude (deg)")
    #print(np.array(lsNH3dlon).shape,np.array(lsNH3dlat).shape)
    
    ###############################################################################
    # CREATE Plot 1: PG Latitude vs Differential Longitude
    ###############################################################################        
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
    
    fig1.savefig(path+'nez_analysis.png',dpi=300)
