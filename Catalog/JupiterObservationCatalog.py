"""
Created on Thu Oct 21 08:33:27 2021

PURPOSE: Plot datetime versus CM (I, II, III) of Jupiter observations. 
         (will be generalizable to all planetary observations eventually)


@author: Steven Hill
"""

def JupiterObservationCatalog():
    import sys
    drive='c:'
    sys.path.append(drive+'/Astronomy/Python Play')
    sys.path.append(drive+'/Astronomy/Python Play/Util_P3')
    sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy')
    sys.path.append(drive+'/Astronomy/Python Play/SPLibraries_P3')
    import os
    import ConfigFiles as CF
    from datetime import datetime, timedelta
    from astropy.table import QTable, Table, Column
    import numpy as np
    from astropy.io import fits, ascii

    date_list=[]
    map_list=[]
    cam_list=[]
    Filter_List=['380NUV','450BLU','467HeII','501OIII','550GRN','550OPN',
                 '620CH4','632OI', '647', '650RED', '656HIA','658NII',
                 '672SII','685NIR','730OII', '742NIR', '807NIR','889CH4',
                 '940NIR','1000NIR']
    path='c:/Astronomy/Projects/Planets/Jupiter/Imaging Data/'
    DIRarray=os.listdir(path)
    mappath='c:/Astronomy/Projects/Planets/Jupiter/Imaging Data/Mapping/NH3 Map Plots/'
    maparraytemp=os.listdir(mappath)

    #print(DIRarray)
    temp1=[k for k in DIRarray if os.path.isdir(path+k)]
    dateUTarray=[k for k in temp1 if "UT" in k]
    pathout='/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/'
    t = Table(names=('Date', 'Camera','Map',
                     '380NUV','450BLU','467HeII','501OIII','550GRN','550OPN',
                     '620CH4','632OI', '647', '650RED', '656HIA','658NII',
                     '672SII','685NIR','730OII', '742NIR', '807NIR','889CH4',
                     '940NIR','1000NIR'),
              dtype=('S10', 'S4','S30',
                     'B','B','B','B','B','B',
                     'B','B','B','B','B','B',
                     'B','B','B','B','B','B',
                     'B','B'))
    #print(dateUTarray)
    for i in dateUTarray:  #Session Date
        df=i[0:4]+"-"+i[4:6]+"-"+i[6:8]
        date_list.append(df)
        temp2=[k for k in maparraytemp if df in k]
        temp3=[k for k in temp2 if "Map.png" in k]
        map_list=map_list+temp3
        if len(temp3)==0:
            MapName=""
        else:
            MapName=str(temp3[0])

        wv_bool=np.zeros((20), dtype=bool)
        Camera=''
        sessionlist=os.listdir(path+i)
        for j in sessionlist: #j is a file in a session
            if "csv" in j:    #Here we are finding if there's  csv (video)
                if "Jupiter_VideoMetaData" in j: #And confirming video data
                    Camera="CMOS"
                    metadata=CF.video_metadata_list(path+i+"/"+j) #read the metadata
                    metadata.load_records()
                    for k in range(0,metadata.nrecords-1):
                        print(metadata.VideoFile[k])
                        strdate=metadata.VideoFile[k][0:15]
                        for ii in range(0,len(wv_bool)):
                            wv_bool[ii]=(Filter_List[ii][0:4] in metadata.VideoFile[k] or wv_bool[ii])
                            #print(ii, waveln,Filter_List[ii][0:3],wv_bool[ii])
                        #print strdate
                        #date=datetime.strptime(strdate,"%Y-%m-%d-%H%M")
                        #print date
                        #date_list.append(date)
                        #cam_list.append('CMOS')

            elif "fit" in j:
                if "Aligned" in j:
                    #print j
                    strdate=j[0:15]
                    Camera="CCD"
                    for k in range(0,metadata.nrecords-1):
                        #print(metadata.VideoFile[k])
                        #strdate=metadata.VideoFile[k][0:15]
                        for ii in range(0,len(wv_bool)):
                            wv_bool[ii]=(Filter_List[ii][0:4] in j or wv_bool[ii])
                    #print strdate
                    #date=datetime.strptime(strdate,"%Y-%m-%d-%H%M")
                    #print date
                    #date_list.append(date)
                    #cam_list.append('CCD')

        t.add_row((i,Camera,MapName,
                   wv_bool[0],wv_bool[1],wv_bool[2],wv_bool[3],wv_bool[4],wv_bool[5],
                   wv_bool[6],wv_bool[7],wv_bool[8],wv_bool[9],wv_bool[10],wv_bool[11],
                   wv_bool[12],wv_bool[13],wv_bool[14],wv_bool[15],wv_bool[16],wv_bool[17],
                   wv_bool[18],wv_bool[19]))
        print(i,
                   wv_bool[1],wv_bool[2],wv_bool[3])
    print()
    print(date_list)
    #print()
    #print(map_list)
    #print()
    #print(t)

    ascii.write(t,pathout+'JupiterObservationCatalog.csv',format='basic',
                overwrite=True,delimiter=',')

###############################################################################
def PlotJupiterObservations(CM=1):
    import matplotlib.pyplot as pl
    from datetime import datetime
    fig,ax=pl.subplots(nrows=4,ncols=1,figsize=(6.0,6.0), dpi=150, 
                       sharex=True,facecolor="white")
    fig.suptitle("Jupiter NH3 Observations",x=0.5,ha='center',color='k',
                 fontsize=16)

    startdate=datetime.strptime("2023-08-01","%Y-%m-%d")
    enddate=datetime.strptime("2024-01-01","%Y-%m-%d")
    subPlotJupiterObservations(ax[0],startdate,enddate,CM=CM)

    startdate=datetime.strptime("2022-07-19","%Y-%m-%d")
    enddate=datetime.strptime("2023-01-18","%Y-%m-%d")
    subPlotJupiterObservations(ax[1],startdate,enddate,CM=CM)
    
    startdate=datetime.strptime("2021-06-21","%Y-%m-%d")
    enddate=datetime.strptime("2021-12-07","%Y-%m-%d")
    subPlotJupiterObservations(ax[2],startdate,enddate,CM=CM)
    
    startdate=datetime.strptime("2020-06-21","%Y-%m-%d")
    enddate=datetime.strptime("2020-12-07","%Y-%m-%d")
    subPlotJupiterObservations(ax[3],startdate,enddate,CM=CM,xtitle=True)
    pathout='c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Catalog/'
    pl.savefig(pathout+'Catalog.png',dpi=320)


def subPlotJupiterObservations(ax,startdate,enddate,CM=1,xtitle=False):
    import sys
    drive='c:'
    sys.path.append(drive+'/Astronomy/Python Play')
    sys.path.append(drive+'/Astronomy/Python Play/Util_P3')
    sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy')
    sys.path.append(drive+'/Astronomy/Python Play/SPLibraries_P3')
    import ephem
    import os
    import numpy as np
    import ConfigFiles as CF
    from datetime import datetime, timedelta
    import EWLibV006_P3 as EWL
    import matplotlib.pyplot as pl
    import matplotlib.dates as mdates

    observer = ephem.Observer()
    #Location from Google Maps at 483 S Oneida Way, Denver 80224
    observer.lon = ephem.degrees('-104.907985')
    observer.lat = ephem.degrees('39.708200')

    date_list=[]
    cam_list=[]
    pathin='c:/Astronomy/Projects/Planets/Jupiter/Imaging Data/'
    path="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/"
    dateUTarray=os.listdir(pathin)
    for i in dateUTarray:
        Video=False
        if os.path.isdir(pathin+i) and str(i)[0:2]=="20":  #Only look in date directories
            print(i)
            sessionlist=os.listdir(pathin+i)  #List of files for a given date
            for j in sessionlist:
                if "Jupiter_VideoMetaData" in j and "csv" in j:
                    Video=True
            for j in sessionlist:
                print(j)
                if "csv" in j:
                    if "Jupiter_VideoMetaData" in j:
                        #print j
                        metadata=CF.video_metadata_list(pathin+i+"/"+j)
                        metadata.load_records()
                        for k in range(0,metadata.nrecords-1):
                            print(metadata.VideoFile[k])
                            strdate=metadata.VideoFile[k][0:15]
                            #print strdate
                            date=datetime.strptime(strdate,"%Y-%m-%d-%H%M")
                            #print date
                            date_list.append(date)
                            cam_list.append('CMOS')
                elif ("fit" in j or "FIT"  in j) and \
                     ("2020" or "2021" or "2022" in i[0:4]) and not(Video):
                    if "Aligned" in j:
                        print("* ",j)
                        strdate=j[0:15]
                        print("** ",strdate)
                        date=datetime.strptime(strdate,"%Y-%m-%d-%H%M")
                        print("*** ",date)
                        date_list.append(date)
                        cam_list.append('CCD')
                        
    date_list_datetime,elev_list,airmass_list,CMI_list,CMII_list,\
        Io_vis_list,Europa_vis_list, \
        Ganymede_vis_list,Callisto_vis_list= \
        EWL.JupiterEphemLists(date_list,observer)
    print(len(date_list_datetime))
    #fig,ax=pl.subplots(nrows=1,ncols=1,figsize=(6.0,4.0), dpi=150, facecolor="white")

    print(cam_list)
    CMOS_indices = [k for k, x in enumerate(cam_list) if x == 'CMOS']
    CCD_indices = [k for k, x in enumerate(cam_list) if x == 'CCD']
    if CM==1:
        ax.scatter(np.array(CMI_list)[CMOS_indices]*180./np.pi,np.array(date_list_datetime)[CMOS_indices],
                   s=2.0,color='C0',label="CMOS Obs")                            
        ax.scatter(np.array(CMI_list)[CCD_indices]*180./np.pi,np.array(date_list_datetime)[CCD_indices],
                   s=2.0,color='C1',label="CCD Obs")      
    elif CM==2:                      
        ax.scatter(np.array(CMII_list)[CMOS_indices]*180./np.pi,np.array(date_list_datetime)[CMOS_indices],
                   s=2.0,color='C0',label="CMOS Obs")                            
        ax.scatter(np.array(CMII_list)[CCD_indices]*180./np.pi,np.array(date_list_datetime)[CCD_indices],
                   s=2.0,color='C1',label="CCD Obs")      

    perijoveTime=['2020-07-25-0625','2020-09-16-0220',
                  '2021-07-21-0814','2021-09-02-2242','2021-10-16-1713',
                  '2022-08-17-1446','2022-09-29-1711','2022-11-06-2138','2022-12-15-0323']
    periJdatetimelist=[]
    for pjt in  range (0,len(perijoveTime)):
        periJdatetimelist.append(datetime.strptime(perijoveTime[pjt],"%Y-%m-%d-%H%M"))
    juno_list_datetime,Juno_elev_list,Juno_airmass_list,Juno_CMI_list,Juno_CMII_list,\
        Juno_Io_vis_list,Juno_Europa_vis_list, \
        Juno_Ganymede_vis_list,Juno_Callisto_vis_list= \
        EWL.JupiterEphemLists(periJdatetimelist,observer)

    JunoCMIEqX=[277.3,216.4,
                107.,225.,132.,
                319.,357.,37.,130.]

    for k in range(0,metadata.nrecords-1):
        #print metadata.VideoFile[k]
        strdate=metadata.VideoFile[k][0:15]
        #print strdate
        date=datetime.strptime(strdate,"%Y-%m-%d-%H%M")


    if CM==1:    
        ax.scatter(JunoCMIEqX,periJdatetimelist,color='k',s=10,
                   label="PJXX Sys 1 EqX")      
        ax.scatter(np.array(Juno_CMI_list)*180./np.pi,juno_list_datetime,color='r',
                   s=10,label="PJXX CMI")      
    elif CM==2:
        delta=np.array(Juno_CMII_list)-np.array(Juno_CMI_list)
        ax.scatter(JunoCMIEqX+delta,periJdatetimelist,color='k',s=10,
                   label="PJXX Sys 2 EqX")      
        ax.scatter(np.array(Juno_CMII_list)*180./np.pi,juno_list_datetime,color='r',
                   s=10,label="PJXX CMII")      
        
    ax.autoscale(enable=False)
    ax.yaxis_date()
    ax.set_adjustable('box')
    ax.tick_params(axis='both',labelsize=8)
    ax.set_xticks(np.linspace(0.,360.,13, endpoint=True))
    print([startdate,enddate])

    ytks = np.arange(startdate, enddate, timedelta(days=28)).astype(datetime)    
    ax.set(xlim=(0.0, 360.0), ylim=(startdate, enddate))
    ax.set_yticks(ytks)
    ax.yaxis.set_major_formatter(mdates.DateFormatter('%b-%d'))
    ax.set_ylabel(str(startdate)[0:4],fontsize=12)

    if xtitle:
        ax.set_xlabel("Sys "+str(CM)+" Longitude (deg)",fontsize=12)
        ax.legend(ncol=4,fontsize=7)
    ax.grid(linewidth=0.2)
    pl.subplots_adjust(left=0.15, bottom=0.12, right=0.95, top=0.92)

