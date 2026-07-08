import socket
hostname = socket.gethostname()
import sys
sys.path.append('C:/Astronomy/Projects/SAS 2021 Ammonia/Visualization-and-Analysis/')

from config_VA import Fletcher_Profile,Plot_TEXES_code,Fletcher_Profile_Out

def Giles2017(dataset='4b'):
    import numpy as np
    fn="Giles2000 Figure "+dataset+".csv"
    path="C:/Astronomy/Projects/SAS 2021 Ammonia/Visualization-and-Analysis/"
    tmp = np.loadtxt(path+fn, delimiter=',')
    return(tmp)

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

def get_GEMINI_TEXES_Fletcher2020(prs=0.43798,mult=1000000):
    """
    This code gets fNH3 from March 12-14 2017 observations of Jupiter with the 
    GEMINI/TEXES instrument. The data is available in the supplementary information
    of Fletcher et al. (2020). This routine gets meridional profile data for
    a single pressure level for each of the seven observation groups in the
    paper.    

    Parameters
    ----------
    prs : float, optional
        DESCRIPTION. Pressure in bar. The default is 0.43798 bar
    mult : float, optional
        DESCRIPTION. Multiplicative factor to convert fNH3 from fraction to ppm.
        The default is 1.0E6.

    Returns
    -------
    latgrid
    scaled_data_mean
    scaled_data_std

    """
    import sys
    import numpy as np
    sys.path.append('./Photometry/code')
    sys.path.append(Plot_TEXES_code[hostname]+'/')
    import ComputeNetRateJupiter_P3 as CNRJ

    pth=Fletcher_Profile[hostname]+'/'
    pressure = np.fromfile(file=pth+"zmean_g1_retnh3_pressure.txt", dtype=float, count=-1, sep=" ")
    data=np.zeros((7,181))

    ind=np.where(np.abs(pressure-prs)/pressure<0.01)    #Pressure index
    PL=np.ndarray.flatten(np.array(ind))[0]             #Pressure
    #print('&&&&&&&&&&&&&&&&&&&&&&&&',ind)
    #print('&&&',np.ndarray.flatten(np.array(ind))[0],'***PL=',PL)
    for i in range(1,8):
        latc = np.fromfile(file=pth+"zmean_g"+str(i)+"_retnh3_lat.txt", dtype=float, count=-1, sep=" ")
        latg=Centric_to_Graphic(latc)
        latsize=len(latg)
        print("latsize=========",latsize)
        Start=latsize*PL
        End=latsize*(PL+1)
        tmp = np.fromfile(file=pth+"zmean_g"+str(i)+"_retnh3_data.txt", dtype=float, count=-1, sep=" ")
        dat=tmp[Start:End]
        latgrid,tmpsig=CNRJ.uniform_lat_grid(latg,dat,Fine=True)
        #return
        #print("dat=",tmpsig.shape, tmpsig)
        data[i-1,:]=tmpsig
    #print(Start,End,pressure[PL])
    scaled_data_mean=np.mean(data,axis=0)*mult
    scaled_data_std=np.std(data,axis=0)*mult
    
    return(latgrid,scaled_data_mean,scaled_data_std)

def get_all_GEMINI_TEXES_Fletcher2020(plot=True):
    """
    Created on Fri Nov 24 11:55:45 2023

    @author: smhil
    """
    #import sys
    import socket
    hostname = socket.gethostname()
    #drive='C:'
    #ys.path.append(drive+'/Astronomy/Python Play')
    #sys.path.append(drive+'/Astronomy/Python Play/Util_P3')
    #sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy_P3')
    #sys.path.append(drive+'/Astronomy/Python Play/SPLibraries_P3')
    import pylab as pl
    import numpy as np
    import matplotlib.ticker as ticker
   
    clr='C2'
    ###########################################################################
    # GET FLETCHER DATA
    ###########################################################################
    pth=Fletcher_Profile[hostname]+'/'

    pressure = np.fromfile(file=pth+"zmean_g1_retnh3_pressure.txt", dtype=float, count=-1, sep=" ")
    pressure=np.array(pressure)
    size=pressure.size
    data=np.zeros((181,size))
    std=np.zeros((181,size))
    dataavg=np.zeros(size)
    datastd=np.zeros(size)
    CH4=np.ones(size)*1810.

    for i in np.arange(0,size):
        plevel=pressure[i]
        latgrid,scaled_data_mean,scaled_data_std=get_GEMINI_TEXES_Fletcher2020(prs=plevel,
                                                               mult=1000000.)

        data[:,i]=scaled_data_mean
        std[:,i]=scaled_data_std
        dataavg[i]=np.nanmean(scaled_data_mean)
        #dataavg[i]=np.nanmean(scaled_data_mean[90:100])
        datastd[i]=np.nanmean(scaled_data_std)
        #print("pres, fNH3= ",plevel,dataavg[i])

    return latgrid,pressure,plevel,size,data,std,dataavg,datastd,CH4

def Profile_Vertical_Fletcher(plot=True):
    """
    Created on Fri Nov 24 11:55:45 2023

    @author: smhil
    """
    #import sys
    import socket
    hostname = socket.gethostname()
    import pylab as pl
    import numpy as np
    import matplotlib.ticker as ticker
   
    clr='C2'
    ###########################################################################
    # GET FLETCHER DATA
    ###########################################################################
    latgrid,pressure,plevel,size,data,std,dataavg,datastd,CH4=get_all_GEMINI_TEXES_Fletcher2020()
    print(latgrid)
    ###########################################################################
    # PLOT MERIDIONAL PROFILES (ALL)
    ###########################################################################
    if plot:
        figavgprof,axsavgprof=pl.subplots(1,1,figsize=(6.0,6.0), dpi=150, facecolor="white")

        for i in np.arange(0,size):
            axsavgprof.plot(latgrid,data[:,i],linewidth=1.0,
                    label='Fletcher etal, 2020 ('+str(int(plevel*1000.))+'mb)',color=clr)
            axsavgprof.fill_between(latgrid, data[:,i]-std[:,i], data[:,i]+std[:,i],
                            color=clr,alpha=0.08)
        axsavgprof.set_yscale('log')
        axsavgprof.set_ylim([0.001,100.])
        axsavgprof.invert_yaxis()
        #axsavgprof.legend()
    
    ###########################################################################
    # PLOT VERTICAL PROFILE
    ###########################################################################
    if plot:
        figvertprof,axsvertprof=pl.subplots(1,1,figsize=(6.0,6.0), dpi=150, facecolor="white")
        axsvertprof.plot(dataavg,pressure,label='Ammonia - Fletcher etal. [2020]')
        axsvertprof.plot(np.mean(data[90:100,:],axis=0),pressure,label='Ammonia - Fletcher etal. [2020]')
        #print(np.nanmean(data,axis=0))
        #print(dataavg)
        axsvertprof.fill_betweenx(pressure,dataavg-datastd,dataavg+datastd,alpha=0.1)
        #print("#################################")
        #print(dataavg,dataavg-datastd,dataavg+datastd)

        axsvertprof.plot(CH4,pressure,label='Methane')

        #Upper Haze
        #axsvertprof.fill_between([1,10000],[0.6,0.6],[0.2,0.2],color='0.9',
        #                         linewidth=[0,0],label='Haze')
        #Sheet Cloud
        #axsvertprof.fill_between([1,10000],[0.65,0.65],[0.67,0.67],color='0.5',
        #                         linewidth=[0,0],label='Sheet Cloud')
        
        axsvertprof.set_xscale('log')
        axsvertprof.set_xlim([1,10000.])
        axsvertprof.set_yscale('log')
        axsvertprof.set_ylim([0.01,10.])
        axsvertprof.tick_params(axis='both', which='both', labelsize=8)
        axsvertprof.yaxis.set_major_formatter(ticker.FormatStrFormatter('%1.2f'))
        #axsvertprof.yaxis.set_minor_formatter(ticker.FormatStrFormatter('%1.1f'))
        axsvertprof.invert_yaxis()
        
        axsvertprof.set_title('Abundances and Clouds')
        axsvertprof.set_xlabel('Abundance (ppm)')
        axsvertprof.set_ylabel('Pressure (bar)')
        axsvertprof.legend()
        
        path=Fletcher_Profile_Out[hostname]+'/'       
        figvertprof.savefig(path+"Profile Vertical Fletcher.png",dpi=300)
        
        tmp=np.array(Giles2017())
        print("#################",tmp[:,0]*1e6,tmp[:,1])
        axsvertprof.plot(tmp[:,0]*1e6,tmp[:,1],label='Giles++ 2017')


    return(pressure,dataavg,datastd)