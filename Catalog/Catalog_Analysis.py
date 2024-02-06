def Catalog_Analysis(writecsv=False,writejson=False,plothist=False):
    """
    Created on Fri Dec 15 10:50:45 2023
    
    PURPOSE: Create a catalog of observations as CSV and/or JSON and to plot
             statistical analyses of the observation catelog. Filtering
             on various parameters will be a future addition.
    @author: smhil
    """


    import sys
    drive='c:'
    sys.path.append(drive+'/Astronomy/Python Play')
    sys.path.append(drive+'/Astronomy/Python Play/Util_P3')
    sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy')
    sys.path.append('./Services')

    import numpy as np
    import json
    from astropy.table import  Table
    from astropy.io import ascii
    sys.path.append('./Services')
    import get_obs_list as getlist
    import get_WINJupos_ephem as WJ_ephem
    import matplotlib.pyplot as pl


    sourcefiles=getlist.get_obs_list()

    if writejson:
        filename="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Catalog/Catalog.json"
        with open(filename, "w") as fp:
            json.dump(sourcefiles , fp) 
    print("TEST0")

    if writecsv or plothist:
        filename="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Catalog/Catalog.csv"
    
        t = Table(names=('ObsID', 'NH3 Time','CM1','CM2','CM3',
                         'Telescope','FL','Camera','Seeing','Transparency',
                         'CH4 Abs / PCloud','CH4 Qual','NH3 Abs / fc(NH3)','NH3Qual','RGB Context'),
                  dtype=('S10','S10','S10','S10','S10',
                         'S10','S10','S10','S10','S10',
                         'S10','B','S10','B','S10'))
    
    
        for ObsID in sourcefiles:
            if int(ObsID[0:4])>2021:
                print(ObsID,ObsID[0:4])
                #for item in sourcefiles[ObsID]:
                #    print(item)
                NH3sec=str(int(str(sourcefiles[ObsID]["NH3file"][16:17]))*6)
                NH3time=(sourcefiles[ObsID]["NH3file"][0:10]+"_"+
                         sourcefiles[ObsID]["NH3file"][11:13]+":"+
                         sourcefiles[ObsID]["NH3file"][13:15]+":"+NH3sec.zfill(2))
                eph=WJ_ephem.get_WINJupos_ephem(NH3time)
                NH3_CM1=float(eph[0].strip())
                NH3_CM2=float(eph[1].strip())
                NH3_CM3=float(eph[2].strip())
    
                print(ObsID,sourcefiles[ObsID]["Telescope"],
                           sourcefiles[ObsID]["FL"],sourcefiles[ObsID]["Camera"],
                           sourcefiles[ObsID]["Seeing"],sourcefiles[ObsID]["Transparency"],
                           sourcefiles[ObsID]["CH4file"],True,#sourcefiles[ObsID]["CH4Qual"],
                           sourcefiles[ObsID]["NH3file"],True,#sourcefiles[ObsID]["NH3Qual"],
                           sourcefiles[ObsID]["RGBfile"])
                t.add_row((ObsID,NH3time,str(NH3_CM1),str(NH3_CM2),str(NH3_CM3),
                           sourcefiles[ObsID]["Telescope"],
                           sourcefiles[ObsID]["FL"],sourcefiles[ObsID]["Camera"],
                           "'"+sourcefiles[ObsID]["Seeing"]+"'","'"+sourcefiles[ObsID]["Transparency"]+"'",
                           sourcefiles[ObsID]["CH4file"],True,#sourcefiles[ObsID]["CH4Qual"],
                           sourcefiles[ObsID]["NH3file"],True,#sourcefiles[ObsID]["NH3Qual"],
                           sourcefiles[ObsID]["RGBfile"]))
        if writecsv:
            ascii.write(t,filename,format='basic',overwrite=True,delimiter=',')

    if plothist:
        CM1array=np.array(t['CM1']).astype(float)
        CM2array=np.array(t['CM2']).astype(float)
        CM3array=np.array(t['CM3']).astype(float)
        #for i in range(0,tarray.size):
        #    print(tarray[i])
        #print("** ",tarray.size)
        histCM1,binedgesCM1=np.histogram(CM1array,bins=np.linspace(0,360,num=36,
                                                           endpoint=True))
        histCM2,binedgesCM2=np.histogram(CM2array,bins=np.linspace(0,360,num=36,
                                                           endpoint=True))
        histCM3,binedgesCM3=np.histogram(CM3array,bins=np.linspace(0,360,num=36,
                                                           endpoint=True))
        
        bincenters = [(binedgesCM1[i]+binedgesCM1[i+1])/2. for i in range(len(binedgesCM1)-1)]
        fig,axs=pl.subplots(3,1,figsize=(6.0,6.0), dpi=150, facecolor="white",
                            sharey=True,sharex=True)

        axs[0].step(bincenters, histCM1)
        axs[0].set_xlim(0,360)
        axs[0].set_xticks(np.linspace(0,360,num=13))
        axs[0].set_ylim(0,10)
        axs[0].set_yticks(np.linspace(0,10,num=6))
        axs[0].set_title("System 1 Longitude")
        axs[1].step(bincenters, histCM2)
        axs[1].set_title("System 2 Longitude")
        axs[2].step(bincenters, histCM3)
        axs[2].set_title("System 3 Longitude")
        axs[2].set_xlabel("Longitude (deg)")       
        fig.suptitle("Distribution of Observations 2022-23")
        
        filename="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Catalog/Catalog_Hist.png"
        fig.savefig(filename,dpi=300)
        
    return(0)