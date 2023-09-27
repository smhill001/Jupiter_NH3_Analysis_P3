# -*- coding: utf-8 -*-
"""
Created on Sun Dec 05 08:48:16 2021

@author: Steven Hill
"""

def plot_TEXES_Groups(ax,clr="C2",prs=0.43798,mult=1.0):
    """
    PURPOSE:    This code reads and plots the NH3 mole fraction data from
                Fletcher et al., 2016 at 440mb pressure. It loops over the
                seven groups of observations and computes an average and
                standard deviation for plotting.
    """
    import matplotlib.pyplot as pl
    import scipy
    import numpy as np
    import matplotlib.pyplot as pl
    import ComputeNetRateJupiter_P3 as CNRJ

    print('***************** mult=',mult)
    pth="c:/Astronomy/Projects/SAS 2021 Ammonia/GeminiTEXES2017/ZonalResults/"
    pressure = scipy.fromfile(file=pth+"zmean_g1_retnh3_pressure.txt", dtype=float, count=-1, sep=" ")
    data=np.zeros((7,181))
    print("data=",data.shape)

    ind=np.where(np.abs(pressure-prs)/pressure<0.01)    #Pressure index
    PL=np.ndarray.flatten(np.array(ind))[0]             #Pressure
    #print('&&&&&&&&&&&&&&&&&&&&&&&&',ind)
    #print('&&&',np.ndarray.flatten(np.array(ind))[0],'***PL=',PL)
    
    PL=np.ndarray.flatten(np.array(ind))[0]
    for i in range(1,8):
        latc = scipy.fromfile(file=pth+"zmean_g"+str(i)+"_retnh3_lat.txt", dtype=float, count=-1, sep=" ")
        latg=Centric_to_Graphic(latc)
        latsize=len(latg)
        #print("latsize=========",latsize)
        Start=latsize*PL
        End=latsize*(PL+1)
        """if i==4:
            Start=59*PL
            End=59*(PL+1)
            latc = scipy.fromfile(file=pth+"zmean_g4_retnh3_lat.txt", dtype=float, count=-1, sep=" ")
            latg=Centric_to_Graphic(latc)
        else:
            Start=61*PL
            End=61*(PL+1)
            latc = scipy.fromfile(file=pth+"zmean_g1_retnh3_lat.txt", dtype=float, count=-1, sep=" ")
            latg=Centric_to_Graphic(latc)"""
        #print("Start,End=",Start,End)
        print(i)
        tmp = scipy.fromfile(file=pth+"zmean_g"+str(i)+"_retnh3_data.txt", dtype=float, count=-1, sep=" ")
        dat=tmp[Start:End]
        latgrid,tmpsig=CNRJ.uniform_lat_grid(latg,dat,Fine=True)

        #print("dat=",tmpsig.shape, tmpsig)
        data[i-1,:]=tmpsig
    #print(Start,End,pressure[PL])
    scaled_data_mean=np.mean(data,axis=0)*mult#*8.0e4
    scaled_data_std=np.std(data,axis=0)*mult#*8.0e4
    ax.plot(latgrid,scaled_data_mean,linewidth=1.0,
            label='Fletcher etal, 2020 ('+str(int(prs*1000.))+'mb)',color=clr)
    ax.fill_between(latgrid, scaled_data_mean-scaled_data_std, scaled_data_mean+scaled_data_std,
                    color=clr,alpha=0.08)
    #ax.plot(latgrid,data[3,:]*mult)
    #for i in range(1,8):
    #    ax.plot(latgrid,data[i-1,:]*mult)

    
def plot_Teifel(ax,clr='C0',width=1.5):
    """
    PURPOSE:    This code reads and plots the zonally averaged NH3 absorption
                EW at 645nm from data scanned from Teifel et al., 2018 
                figure 7.
    """
    import numpy as np
    from numpy import genfromtxt
    Lats=[-24.0,-12.5,0.0,12.5,24.0]
    EWs=np.array([5.92,6.78,6.75,6.35,5.38])*0.1
    #ax.scatter(Lats,EWs)
    pth="c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/"
    Teifel = np.array(genfromtxt(pth+"Teifel2018-Fig7.txt", delimiter=','))
    ax.plot(Teifel[4:24,0],Teifel[4:24,1]*0.1,label='Teifel etal, 2018',
            linewidth=width,color=clr)

def plot_Historical(ax,reference,clr='C0'):
    """
    PURPOSE:    This code reads and plots the zonally averaged NH3 absorption
                EW at 645nm from data scanned from Teifel et al., 2018 
                figure 7.
    """
    import matplotlib.pyplot as pl
    import numpy as np
    data={"Vdov2021":{"Region":["SPR","STB","STrZ","SEB","EZ","NEB","NTrZ","NTB","NPR","GRS"],
                      "Center_pgLat":[-45.00,-29.75,-23.40,-13.45,-0.15,15.60,20.80,27.8,45.00,-23.00],
                      "645EW":[4.2,5.1,5.7,6.0,6.2,5.5,4.4,4.7,4.6,4.7]},
          "Teif2018":{"Region":["STrZ","SEB","EZ","NEB","NTrZ"],
                            "Center_pgLat":[-23.40,-13.45,-0.15,15.60,20.8],
                            "645EW":[5.92,6.78,6.75,6.35,5.38]},
          "More1991":{"Region":["SPR","STB","STrZ","SEB","EZ","NEB","NTB","NPR"],
                            "Center_pgLat":[-45.00,-29.75,-23.40,-13.45,-0.15,20.80,27.8,45],
                            "645EW":[5.6,5.7,7.8,9.8,7.7,4.9,5.7,7.2]},
          "L&O1980":{"Region":["STrZ","SEB","EZ","NEB","GRS"],
                            "Center_pgLat":[-23.40,-13.45,-0.15,15.60,-23.0],
                            "645EW":[11.00,12.90,11.7,8.3,9.3]}}
    ax.scatter(data[reference]["Center_pgLat"],np.array(data[reference]["645EW"])*0.1,label=reference,color=clr)


def plot_profile_L2(ax,reference,LatLims=[45,135],LonRng=45.,band="CH4",
                    profile="Meridional",clr='C0',
                    width=1.0,style='solid',smooth=False):
    #Single module for aggregating profiles (meridional for now)
    #Also, the disk integrated transmissions
    #     are hardcoded!!! I need to have a master reference file for these 
    #     transmissions!
    # This is really profile aggregation- should I rename the program?
    
    import numpy as np
    from astropy.convolution import convolve, Box1DKernel
    import sys
    sys.path.append('./Services')
    import read_master_calibration
    import extract_profile as EP
    import get_L2_abs_data as GAOD
    import pylab as pl
    import get_batch_lists as GBL

    calibration,K_eff=read_master_calibration.read_master_calibration()
    sourcefiles=GAOD.get_L2_abs_data()
    DataSets=GBL.get_batch_lists()

    pth="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Analysis Data/L2 FITS/"

    Trans2EW={"VLT":{"EW_slope":{"NH3":-12.15343491,"CH4":-12.74115968},
                     "EW_const":{"NH3":12.15195684,"CH4":12.73323535}},
              "SCT":{"EW_slope":{"NH3":-12.87087479,"CH4":-13.52007859},
                     "EW_const":{"NH3":12.86940675,"CH4":13.51315382}}}

    if profile=="Meridional":
        AvgSum=np.zeros(180)
        StdSum=np.zeros(180)
        xlabel="Planetographic Latitude (deg)"
    elif profile=="Zonal":
        AvgSum=np.zeros(2*LonRng)
        StdSum=np.zeros(2*LonRng)
        xlabel="Longitude from CM (deg)"
    
    suffix={"CH4":"_Jupiter_620CH4AbsMap",
            "NH3":"_Jupiter_647NH3AbsMap"}
    
    figspaghetti,axsspaghetti=pl.subplots(1,1,figsize=(6.0,6.0), dpi=150, 
                                          facecolor="white")
    First=True
    Count=0
    print(reference)
    print(DataSets[reference])
    for ID in DataSets[reference]:
        print("*******ID=",ID)
        sourceindex=ID[0:4]+ID
        if len(ID)==11:
            version=ID[10]
            dataset=ID[0:10].replace('-','')+'UT'+version+"_Map"
        else:
            dataset=ID.replace('-','')+'UT'+"_Map"
        print("******dataset=",dataset)
        print(sourcefiles[dataset])
        file=sourcefiles[dataset][band+"file"]+suffix[band]
        try:
            file=file+sourcefiles[dataset]['Metadata']['Variation']
        except:
            print('No Variation')
        print("file=",file)
        
        
        Lats,AvgProf,StdProf=EP.extract_profile(pth,file+".fits",LonRng=LonRng,
                                                  profile=profile)
        Avg_EW=Trans2EW[reference[0:3]]["EW_slope"][band]*AvgProf+Trans2EW[reference[0:3]]["EW_const"][band]
        Std_EW=Trans2EW[reference[0:3]]["EW_slope"][band]*AvgProf+Trans2EW[reference[0:3]]["EW_const"][band]
        
        try:
            lbl=file[5:17]+" "+sourcefiles[dataset]['Metadata']['Variation']
        except:
            lbl=file[5:17]
        axsspaghetti.plot(Lats,Avg_EW,linewidth=0.5,
                          label=lbl)#+"_"+file[5:20])
        if First:
            Sum=Avg_EW
            Count=1
            First=False
        else:
            Sum=Sum+Avg_EW
            Count=Count+1
            print("No File")
    AvgSum=Sum/Count
    #Profile=np.loadtxt(data[reference]["path"]+data[reference]["fn"],usecols=range(2),delimiter=",")
    axsspaghetti.set_xlim(-45,45)
    axsspaghetti.set_ylim(0,2)
    if band=="CH4":
        axsspaghetti.set_ylim(0,2)
        axsspaghetti.set_title(reference+" Methane")

    elif band=="NH3":
        axsspaghetti.set_ylim(0,1)
        axsspaghetti.set_title(reference+" Ammonia")

    axsspaghetti.set_ylabel("Equivalent Width EW (nm)",fontsize=10)

    axsspaghetti.legend(fontsize=8,ncol=3)
    axsspaghetti.set_xlabel(xlabel,fontsize=10)

    if smooth:
        ax.plot(Lats,convolve(AvgSum, Box1DKernel(3)),color=clr,
                linewidth=width,linestyle=style,label=reference)        
        #ax.plot(Profile[:,0]-45.,convolve(EW, Box1DKernel(7)),color=clr,linewidth=width,linestyle=style,label=reference)        ax.plot(Lats,AvgMeridEW,color=clr,linewidth=width,linestyle=style,label=reference)
    else:
        #ax.plot(Profile[:,0]-45.,EW,color=clr,linewidth=width,linestyle=style,label=reference)
        ax.plot(Lats,AvgSum,color=clr,linewidth=width,
                linestyle=style,label=reference)
        
    path="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/"

    figspaghetti.savefig(path+"Analysis Data/Profiles/Profile_"+reference+"_"+band+"_"+profile+"_Absorption.png",dpi=300)
    
    return(0)

def plot_profiles_L3(ax,reference,LatLims=[45,135],LonRng=45.,profile="Meridional",
                         clr='C0',width=1.0,
                         style='solid',smooth=False,param='PCloud'):
    #     Also, the disk integrated transmissions
    #     are hardcoded!!! I need to have a master reference file for these 
    #     transmissions!
    #     !!!!And...what are we hardcoded to specific manually generated profile
    #     files?
    
    import numpy as np
    from astropy.convolution import convolve, Box1DKernel
    from astropy.io import fits
    import RetrievalLibrary as RL
    import sys
    sys.path.append('./Services')
    #import read_master_calibration
    import RetrievalLibrary as RL
    import extract_profile as EP
    import get_L2_abs_data as GAOD
    import pylab as pl

    sourcefiles=GAOD.get_L2_abs_data()

    data={"SCT 2022":{"path":"C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Analysis Data/L3 FITS/",
                              #"fn":"Profile of Jupiter GRS AVG NH3 Transmission 8 files-Celestron11.csv",
                              "fn":[],
                                    #why is 10/21 missing??? I've used it before!!!
                                    #"2022-10-21-0358_6-Jupiter-fNH3_Sys2.fits"],
                              "EW_slope":-12.82831873,
                              "EW_const":12.82685019},
          "VLTMUSE 2022":{"path":"C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Analysis Data/L3 FITS/",
                              "fn":["2022-09-19-0352_3"],
                              "EW_slope":-12.15343491,
                              "EW_const":12.15195684},
          "SCT 2023":{"path":"C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Analysis Data/L3 FITS/",
                              "fn":[],
                              "EW_slope":-12.82831873,
                              "EW_const":12.82685019}}
    if param=="PCloud":
        band="CH4"
        suffix="-Jupiter_PCloud_Sys2.fits"
        factor=1.0
    elif param=="fNH3":
        band="NH3"
        suffix="-Jupiter_fNH3_Sys2.fits"
        factor=1.0e6
    
    figspaghetti,axsspaghetti=pl.subplots(1,1,figsize=(6.0,6.0), dpi=150, facecolor="white")

    for ID in sourcefiles:
        if "Map" in ID:
            print("*******ID=",ID)
            if sourcefiles[ID][band+'Qual']:
                if 202207<int(ID[0:6])<202302 and sourcefiles[ID]["Metadata"]["Telescope"]=='C11':
                    print(int(ID[0:6]))
                    print("*******sourcefiles[ID]['CH4file'][0:4]=",sourcefiles[ID][band+"file"])
                    data["SCT 2022"]["fn"].append(sourcefiles[ID][band+"file"])
                elif 202307<int(ID[0:6])<202403:
                    print(int(ID[0:6]))
                    print("*******sourcefiles[ID]['CH4file'][0:4]=",sourcefiles[ID][band+"file"])
                    data["SCT 2023"]["fn"].append(sourcefiles[ID][band+"file"])

    pth=data[reference]["path"]
    if profile=="Meridional":
        AvgSum=np.zeros(180)
        StdSum=np.zeros(180)
    elif profile=="Zonal":
        AvgSum=np.zeros(2*LonRng)
        StdSum=np.zeros(2*LonRng)

    First=True
    for L3file in data[reference]["fn"]:                             
        fn=L3file+suffix
        print("$$$$$$$$$$$$$$fn=",fn)
        Lats,AvgMerid,StdMerid=EP.extract_profile(pth,fn,LonRng=LonRng,profile=profile)       
        #Convert from fraction to ppm by multiplying by 10^6
        if First and profile=="Meridional":
            AvgArr=AvgMerid
            AvgArr=np.reshape(AvgArr,(1,180))
            First=False
        if First and profile=="Zonal":
            AvgArr=AvgMerid
            AvgArr=np.reshape(AvgArr,(1,90))
            First=False
        else:
            AvgArr=np.vstack((AvgArr,AvgMerid))
            
        axsspaghetti.plot(Lats,AvgMerid*factor,linewidth=0.5,label=L3file[5:20])

    print("AvgArr.shape=",AvgArr.shape)
    AvgPro=np.mean(AvgArr[:,:],axis=0)*factor
    AvgStd=np.std(AvgArr[:,:],axis=0)*factor
    print("AvgPro.shape=",AvgPro.shape)
    print("AvgPro=",AvgPro)
    
    #Profile=np.loadtxt(data[reference]["path"]+data[reference]["fn"],usecols=range(2),delimiter=",")
    axsspaghetti.set_xlim(-45,45)
    if param=="PCloud":
        axsspaghetti.set_ylim(400,1000)
        axsspaghetti.set_ylabel("Effective Cloud-Top Pressure (mb)",fontsize=10)
    elif param=="fNH3":
        axsspaghetti.set_ylim(0,300)
        axsspaghetti.set_ylabel("Column-Average Ammonia Abundance (ppm)",fontsize=10)
    axsspaghetti.legend(fontsize=6,ncol=3)
    axsspaghetti.set_title(reference)
    axsspaghetti.set_xlabel("Planetographic Latitude (deg)",fontsize=10)

    if smooth:
        ax.plot(Lats,convolve(AvgPro, Box1DKernel(3)),color=clr,
                linewidth=width,linestyle=style,label=reference)  
        ax.fill_between(Lats, AvgPro-AvgStd, AvgPro+AvgStd,color=clr,alpha=.1)

        #ax.plot(Profile[:,0]-45.,convolve(EW, Box1DKernel(7)),color=clr,linewidth=width,linestyle=style,label=reference)        ax.plot(Lats,AvgMeridEW,color=clr,linewidth=width,linestyle=style,label=reference)
    else:
        #ax.plot(Profile[:,0]-45.,EW,color=clr,linewidth=width,linestyle=style,label=reference)
        ax.plot(Lats,AvgPro,color=clr,linewidth=width,
                linestyle=style,label=reference)
        ax.fill_between(Lats, AvgPro-AvgStd, AvgPro+AvgStd,color=clr,alpha=.1)

    path="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/"

    figspaghetti.savefig(path+"Analysis Data/Profiles/Profile_"+reference+"_"+param+"_"+profile+".png",dpi=300)

        
def Centric_to_Graphic(Latc):
    #Formula used is from Simon and Beebe, 1996
    import numpy as np
    Req=71492.0
    Rp=66845.0
    Latg=Latc
    for i in range(len(Latc)):
        Latg[i]=np.arctan(((Req/Rp)**2)*np.tan(Latc[i]*np.pi/180.))*180.0/np.pi
    #print(Latg)
    return Latg