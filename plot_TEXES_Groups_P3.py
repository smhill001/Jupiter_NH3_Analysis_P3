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

def plot_VLTMUSEandC11_EW_profiles(ax,reference,clr='C0',width=1.0,style='solid',smooth=False):
    #!!!This looks like it plots ABSORPTION and can do so as and EW for 
    #     comparision to prior work. However, there is only one linear fit
    #     for converting transmission to EW when separate ones are needed
    #     for SCT and VLT filters. Also, the disk integrated transmissions
    #     are hardcoded!!! I need to have a master reference file for these 
    #     transmissions!
    #     !!!!And...what are we hardcoded to specific manually generated profile
    #     files?
    
    import numpy as np
    from astropy.convolution import convolve, Box1DKernel
    
    data={"C11 2022":{"path":"/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/",
                              #"fn":"Profile of Jupiter GRS AVG NH3 Transmission 8 files-Celestron11.csv",
                              "fn":"Profile of Jupiter GRS AVG NH3 Transmission 5 files-Celestron11 +-20deg.csv",
                              "Disk_Int_Trans":0.972,
                              "EW_slope":-12.82831873,
                              "EW_const":12.82685019},
          "VLTMUSE 2022":{"path":"C:/Astronomy/Projects/Planets/Jupiter/Imaging Data/20220919UT/",
                              "fn":"Profile of 2022-09-19-0352_3-Jupiter-647NH3AbsMap-final1.csv",                           
                              "Disk_Int_Trans":0.962,
                              "EW_slope":-12.15343491,
                              "EW_const":12.15195684}}
                             
    Profile=np.loadtxt(data[reference]["path"]+data[reference]["fn"],usecols=range(2),delimiter=",")
    print("Profile.shape=",Profile.shape)
    Transmission=(Profile[:,1])*data[reference]["Disk_Int_Trans"]
    EW=data[reference]["EW_slope"]*Transmission+data[reference]["EW_const"]
    print("EW.shape=",EW.shape)
    if smooth:
        ax.plot(Profile[:,0]-45.,convolve(EW, Box1DKernel(7)),color=clr,linewidth=width,linestyle=style,label=reference)
    else:
        ax.plot(Profile[:,0]-45.,EW,color=clr,linewidth=width,linestyle=style,label=reference)


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