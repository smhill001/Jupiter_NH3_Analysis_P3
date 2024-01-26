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
    from scipy.stats import gaussian_kde
    
    path="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Studies/NEZ/"
    
    pairs={} #original bland dictionary
    # NH3 rec, NH3 date, CH4 rec, CH4 date ...
    ###########################################################################
    # READ IN HILL MEASUREMENTS (DATE,LON,LAT) OF NH3 ENHANCEMENTS AND
    # CH4 DEEP CLOUDS.
    ###########################################################################
    fn="Hill_ammonia_mea.CSV"
    #Set up empty lists
    #CH4 deep cloud locations
    CH4rec=[] #Record number
    CH4obs=[] #Observation date (UT)
    CH4lon=[] #System 1 Longitude
    CH4lat=[] #Planetographic Latitude
    #NH3 enhancement locations
    NH3rec=[] #Record number
    NH3obs=[] #Observation date (UT)
    NH3lon=[] #System 1 Longitude
    NH3lat=[] #Planetographic Latitude
    
    with open(path+fn) as csv_file:
        csv_reader = csv.DictReader(csv_file)
        line_count = 0
        for row in csv_reader:
            if line_count == 0:
                #print(f'Column names are {", ".join(row)}')
                line_count += 1
            #print(f'\t{row["REC_NO"]} works in the {row["OBJECT"]} department, and was born in {row["REGION"]}.')
            if row["OBJECT"]=="DC1_SPOT":
                CH4rec.append(row["REC_NO"])
                CH4obs.append(row["DATE"])
                CH4lon.append(float(row["LONG_1"]))
                CH4lat.append(float(row["LAT"]))
            if row["OBJECT"]=="WC1_SPOT":
                NH3rec.append(row["REC_NO"])
                NH3obs.append(row["DATE"])
                NH3lon.append(float(row["LONG_1"]))
                NH3lat.append(float(row["LAT"]))
            line_count += 1
    print(f'Processed {line_count} lines.')
    
    ###########################################################################
    # READ IN WINJupos MEASUREMENTS (DATE,LON,LAT) OF DARK PROJECTIONS 
    # (AND FESTOONS?)
    ###########################################################################
    
    #WINJupos measurer PROJections - all locations
    PRJrec=[] #Record number
    PRJobs=[] #Observation date (UT)
    PRJlon=[] #System 1 Longitude
    PRJlat=[] #Planetographic Latitude
    
    fn="DC_PROJ_2022-23_sel.CSV"
    with open(path+fn) as csv_file:
        csv_reader = csv.DictReader(csv_file)
        line_count = 0
        for row in csv_reader:
            if line_count == 0:
                print(f'Column names are {", ".join(row)}')
                line_count += 1
            if "PROJ" in row["OBJECT"]: #=="DC1_PROJ":
                PRJrec.append(row["REC_NO"])
                PRJobs.append(row["DATE"])
                PRJlon.append(float(row["LONG_1"]))
                PRJlat.append(float(row["LAT"]))
            line_count += 1

    print(f'Processed {line_count} lines.')
    
    #print(PRJobs)
    
    ###########################################################################
    # LOOP OVER NH3 ENHANCEMENT MEASUREMENTS TO FIND *ALL* CH4 DEEP CLOUD
    # MEASUREMENTS MADE ON THE SAME UT DATE. OF THOSE DEEP CLOUD MEASUREMENTS,
    # FIND THE *ONE* THAT IS THE CLOSEST IN LONGITUDE TO THE SELECTED
    # NH3 ENHANCEMENT.
    ###########################################################################
    NH3dobs=[]    #list of NH3 obs dates with corresponding CH4 obs
    NH3dlon=[]    #list of corresponding delta longidutes (NH3lon-CH4lon)
    NH3dlat=[]    #list of corresponding delta latitudes (NH3lat-CH4lat)
    CH4latSel=[]  #list of CH4 latitudes that have corresponding NH3
                    #enhancements (for determining mean CH4 latitude)

    PRJdlon=[]    #Relative Long to mean CH4 Deep Cloud
    PRJlonsel=[]  #Long of Dark Projections on same UT dates as CH4 Deep Cloud
    PRJdlonsel=[] #Relative Long to mean CH4 Deep Cloud on same UT dates as CH4  
    PRJlatsel=[]  #Lat of Dark Projections on same UT dates as CH4 Deep Cloud
    PRJlatClose=[]                    
    for j in range(0,len(NH3obs)):
    #for j in range(0,20):
        #Get set of indices for the CH4 records with the same UT dates 
        #  (not times) as the current (jth) NH3 enhancement. With multiple
        #  NH3 enhancements on a given date, the set of 'matching' CH4
        #  records can be repeated multiple times.
        #print("j,NH3obs[j]=",j,NH3obs[j],NH3lon[j])
        indices=[i for i, e in enumerate(CH4obs) if e == NH3obs[j]]
        #print(NH33obs[j],indices[:],np.take(CH4lon,indices),)
        #compute the index of CH4 deep cloud closest in longitude to the
        # jth NH3 enhancement
        closesti=np.argmin(np.abs(NH3lon[j]-np.take(CH4lon,indices)))
        #print(closesti,np.take(CH4lon,indices[closesti]))
        CH4darklon=np.take(CH4lon,indices[closesti])
        CH4darklat=np.take(CH4lat,indices[closesti])
        #dlonmin=NH3lon[j]-np.take(CH4lon,indices[closesti])
        dlonmin=NH3lon[j]-CH4darklon
        #print(dlonmin)
        
        #Filter for CH4 DFs within 45 deg lon of the NH3 enhancement and compute
        #proper sign for dlons that are over 180 deg
        if np.abs(dlonmin)<45:
            if dlonmin<180:
                NH3dlon.append(dlonmin)
            if dlonmin>=180:
                NH3dlon.append(dlonmin-360.)
            #NH3dobs.append(NH3obs[j])
            #Save the differential NH3-CH4 latitude along with the latitude
            # of the CH4 deep cloud spot (for later averaging)
            NH3dlat.append(NH3lat[j]-CH4darklat)
            CH4latSel.append(CH4darklat)

            #print("match")
            #print("NH3obs,NH3lon,CH4lon,dlon=",NH3obs[j],NH3lon[j],np.take(CH4lon,indices[closesti]),dlonmin)
            
        #print()
        
        #Get set of indices for the CH4 records with the same UT dates 
        #  (not times) as the current (jth) NH3 enhancement. With multiple
        #  NH3 enhancements on a given date, the set of 'matching' CH4
        #  records can be repeated multiple times.
        indicesPRJ=[i for i, e in enumerate(PRJobs) if e == NH3obs[j]]
        #if isinstance(indicesPRJ, list): 
            #print(PRJobs)
        #compute the closest CH4 DR in longitude and the index of that value 
        #closest=np.min(np.abs(NH3lon[j]-np.take(CH4lon,indices)))
        #print("j,NH3obs[j]=",j,NH3obs[j],NH3lon[j])
        #print(indicesPRJ[:],np.take(PRJobs,indicesPRJ),
        #      np.take(PRJlon,indicesPRJ))
        # Load all Sys 1 longitude for PROJs on the given NH3obs date
        
        if len(indicesPRJ)>>0:
            closesti=np.argmin(np.abs(NH3lon[j]-np.take(PRJlon,indicesPRJ)))
            print()
            print("j,NH3obs,PRJobs=",j,NH3obs[j],np.take(PRJobs,indicesPRJ),
                  np.take(PRJlon,indicesPRJ))
            print("*")
            print(closesti,np.take(PRJobs,indicesPRJ[closesti]))
            PRJdarklon=np.take(PRJlon,indicesPRJ[closesti])
            PRJdarklat=np.take(PRJlat,indicesPRJ[closesti])
            #dlonmin=NH3lon[j]-np.take(CH4lon,indices[closesti])
            dlonmin=NH3lon[j]-PRJdarklon
            print("NH3lon[j],indices,np.take(PRJlon,indices)=",
                  NH3lon[j],indices,np.take(PRJlon,indices))
            print("j,dlonmin=",j,dlonmin)
            
            if np.abs(dlonmin)<45:
                if dlonmin<180:
                    PRJdlon.append(dlonmin)
                if dlonmin>=180:
                    PRJdlon.append(dlonmin-360.)
     
                #NH3dlat.append(NH3lat[j]-CH4darklat)
                PRJlatClose.append(PRJdarklat)
    
            PRJlonsel.append(np.take(PRJlon,indicesPRJ))
            # Load all differential longitudes for PROJs on the given NH3obs date
            PRJdlonsel.append(np.take(PRJlon,indicesPRJ)-CH4darklon)
                         # Load all PG latitudes for PROJs on the given NH3obs day
            PRJlatsel.append(np.take(PRJlat,indicesPRJ))
            print(j,NH3obs[j],len(indicesPRJ),)

        #!!! I THINK I MAY NEED TO DO A SEPARATE LOOP HERE SINCE I NEED
        #!!! THE PROJECTIONS LOCATIONS RELATIVE TO THE CH4 HOLES, NOT
        #!!! RELATIVE TO THE NH3 ENHANCEMENTS. EITHER THAT OR I NEED OT ADD
        #!!! AN OFFSET FOR THE DIFFERENCE IN THE MEAN NH3 POSITION VERSUS.
        #!!! I ALSO DON'T SEEM TO ACCOUNT FOR THE MODULO 360 AT ALL HERE
        #!!! WHEREAS I DO FOR THE NH3-CH4 DIFFERENCE.
    
    CH4lon=np.array(CH4lon)
    CH4lat=np.array(CH4lat)
    NH3lon=np.array(NH3lon)
    NH3lat=np.array(NH3lat)

    PRJlatsel=np.array(np.concatenate(PRJlatsel, axis=0 ))
    PRJlonsel=np.array(np.concatenate(PRJlonsel, axis=0 ))
    PRJdlonsel=np.array(np.concatenate(PRJdlonsel, axis=0 ))
    print("PRJlonsel[0:2],PRJdlonsel[0:2]=",PRJlonsel[0:2],PRJdlonsel[0:2])
    PRJdlonsel[np.where(np.abs(PRJdlonsel)>180.)]=360.-PRJdlonsel[np.where(np.abs(PRJdlonsel)>180.)]
    #PRJdlon=PRJdlon(np.where(np.abs(PRJdlon)<=10.))
    #PRJlatsel=PRJdlon(np.where(np.abs(PRJdlon)<=10.))
    
    #PRJdlon=np.concatenate( PRJdlon, axis=0 )
    print()
    print("Mean and STD NH3 delta long=",np.mean(NH3dlon),np.std(NH3dlon))
    print("Mean and STD NH3 delta lat=",np.mean(NH3dlat),np.std(NH3dlat))
    #print(np.mean(PRJdlonsel),np.std(PRJdlonsel))
    #print(np.mean(PRJlatsel),np.std(PRJlatsel))
    ###############################################################################
    # CREATE Plot 0: PG Latitude vs System 1 Longitude
    ###############################################################################        
    fig1,axs1=pl.subplots(2,1,figsize=(6.0,6.0), dpi=150, facecolor="white")
    fig1.suptitle("NH3, CH4, & Dark Projections from "+str(NH3obs[0])+' to '+str(NH3obs[len(NH3obs)-1]))
    
    axs1[0].scatter(CH4lon,CH4lat,color='C0',s=5,label="CH4 Deep Cloud")
    axs1[0].scatter(NH3lon,NH3lat,color='C1',s=5,label="NH3 Enhancement")
    #print(np.array(PRJlon).shape,np.array(PRJlat).shape)
    #axs1[0].scatter(np.array(PRJlon),np.array(PRJlat),color='0.5',s=0.5)
    axs1[0].scatter(PRJlonsel,PRJlatsel,color='k',s=2,label="Proj/Fest")
    axs1[0].set_xlim(0,360)
    axs1[0].set_ylim(0,10)
    axs1[0].set_xlabel("System 1 Longitude (deg)")
    axs1[0].set_ylabel("PG Latitude (deg)")
    axs1[0].legend(ncol=3,fontsize=6,loc=8)
    #print(np.array(NH3dlon).shape,np.array(NH3dlat).shape)
    
    ###############################################################################
    # CREATE Plot 1: PG Latitude vs Differential Longitude
    ###############################################################################
    
    ############
    x=np.array(NH3dlon)
    y=np.array(NH3dlat)+np.mean(np.array(CH4latSel))
    heatmap, xedges, yedges = np.histogram2d(x,y,bins=(20,10))
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    #axs1[1].imshow(heatmap.T, extent=extent, origin='lower',cmap='gist_heat_r')
    ############
    nbins = 20
    k = gaussian_kde([x,y])
    xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))
    #axs1[1].set_title('Calculate Gaussian KDE')
    #axs1[1].pcolormesh(xi, yi, zi.reshape(xi.shape), cmap='gist_heat_r')
    axs1[1].pcolormesh(xi, yi, zi.reshape(xi.shape), shading='gouraud', 
                       cmap='gist_heat_r',label="Density")
    #axs1[1].contour(xi, yi, zi.reshape(xi.shape) )

    axs1[1].scatter(np.array([0.0]),np.mean(np.array(CH4latSel)),color='C0',s=50,
                    alpha=0.8,label="<CH4 Deep Cloud>")
    axs1[1].scatter(np.mean(np.array(NH3dlon)),np.mean(np.array(NH3lat)),color='C1',s=50,
                    alpha=0.8,label="<NH3 Enhancement>")
    axs1[1].scatter(np.array(NH3dlon),np.array(NH3dlat)+np.mean(np.array(CH4latSel)),
                    color='C1',s=5)
    #axs1[1].scatter(np.array(PRJdlonsel),np.array(PRJlatsel),color='k',s=2)
    
    PRJdlonsela=PRJdlonsel[np.where(np.abs(PRJdlonsel)<=45.)]
    PRJlatsela=PRJlatsel[np.where(np.abs(PRJdlonsel)<=45.)]
    #axs1[1].scatter(np.array(PRJdlon),np.array(PRJlat),color='k',s=0.5)
    #axs1[1].scatter(np.array(PRJdlonsela),np.array(PRJlatsela),color='k',s=2)
    axs1[1].scatter(np.array(PRJdlon),np.array(PRJlatClose),color='k',s=2)
    #axs1[1].scatter(np.array(PRJlonsel),np.array(PRJlatsel),color='r',s=2)
    axs1[1].scatter(np.mean(PRJdlonsela),np.mean(PRJlatsela),color='k',s=50,
                    alpha=0.8,label="<Proj/Fest>")
    axs1[1].set_xlim(-45,45)
    axs1[1].invert_xaxis()
    axs1[1].set_ylim(0,10)
    #axs1[1].set_title("Longitude Relative to Closest CH4 Deep Cloud (deg)")
    axs1[1].set_xlabel("Longitude Relative to Closest CH4 Deep Cloud (deg)")
    axs1[1].set_ylabel("PG Latitude (deg)")
    axs1[1].legend(ncol=1,fontsize=6)

    print()
    print("np.array(PRJdlon).shape,np.array(PRJlat).shape",
          np.array(PRJdlon).shape,np.array(PRJlat).shape)
    print()
    print("np.array(PRJdlonsel).shape,np.array(PRJlatsel).shape",
          np.array(PRJdlonsel).shape,np.array(PRJlatsel).shape)
    print()
    print("np.array(PRJdlonsela).shape,np.array(PRJlatsela).shape",
          np.array(PRJdlonsela).shape,np.array(PRJlatsela).shape)

    fig1.subplots_adjust(left=0.10, bottom=0.08, right=0.98, top=0.93,
                wspace=0.25, hspace=0.25)     

    fig1.savefig(path+'nez_analysis.png',dpi=300)
    
    
