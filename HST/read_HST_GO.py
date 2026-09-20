import sys
sys.path.append('C:/Astronomy/Projects/SAS 2021 Ammonia/Visualization-and-Analysis/')
sys.path.append('C:/Astronomy/Projects/SAS 2021 Ammonia/Data-Management-and-Access/processes/')

import matplotlib.pyplot as pl
import numpy as np
#import MakeContiguousMap as MCM
import make_patch as mp
import plot_patch as pp
import flatten_patch as fp
#from matplotlib.pyplot import imread
from astropy.io import fits
from scipy.ndimage import gaussian_filter
from matplotlib import cm, colors
from scipy.ndimage import zoom
import L4_Jup_Map_Plot_V2 as L4MP

#set constants and set gravity as a function of latitude (Jupiter is so oblate that gravity
#  varies significantly from the equator to the poles. 

amagat = 2.69e24 #Lodschmits number. (cm-2) This is really km-amagat - an arcane atmospheric science unit
mean_mol_wt = 3.85e-24 #gm/molecule, which is 2.22 gm/mole
fCH4 = 2.04e-3 #mole fraction of ammonia (old value from Galileo was 1.81e-3)
STP = 1.01e6  #dyne/cm^2 [(g-cm/s^2)/cm^2]`- standard pressure for one atmosphere
#K_eff_CH4620 = 0.427 #SCT & VLT value
#K_eff_NH3647 = 2.955 #SCT & VLT value
K_eff_CH4619 = 0.548 #HST value
K_eff_NH3645 = 3.090 #HST value

def gravity(planet="Jupiter"):
#This function calculates gravity (cm/s2) at a given planet latitude
    import numpy as np
    
    G=6.673*10**(-11)
    lats=90.0-np.linspace(0.5,179.5,180)
    garr=np.zeros((180,1))
    
    if planet=='Saturn':
        M=568*10**24
        Re=60268000
        epsilon=0.0980
        T=38362
    elif planet=='Jupiter':
        M=1901*10**24 #(1.898 × 10^27 kg)
        Re=71541000
        epsilon=0.06492
        T=35730
    
    latindx=0
    for lat in lats:
        print(latindx,lat)
        garr[latindx,0]=(100*G*M)/(Re**2*(1-epsilon*np.sin(np.deg2rad(lat))**2))-4*(np.pi**2)*Re*(1-epsilon*np.sin(np.deg2rad(lat))**2)*np.cos(np.deg2rad(lat))/T**2
        latindx=latindx+1
    #pl.plot(lats,garr[:,0])
    
    gmap=garr
    #print(gmap)
    for i in range(0,359):
        gmap=np.hstack((gmap,garr))
        
    #print(gmap.shape)
    #pl.imshow(gmap)
    return(gmap)


def computeCloudPressure(CH4data,LatLims,LonLims,dataversion=2):
    import make_patch as mp
    grav = gravity() #see file dropped in processes
    #compute ammonia optical depth according to Beer-Lambert law
    CH4_tau = -np.log(CH4data)
    #compute column abundance of molecules along the line-of-sight where 1000 corrects for units
    #  and K_eff are the effective absorption cross sections of each molecule, 0.427 and 2.955 for 
    #  methane and ammonia respectively
    CH4_Ncol = 1000*CH4_tau/K_eff_CH4619 
    #compute cloud pressure - CH4_Cloud_Press is computed in units of mb
    #!! zoom added to scale gravity grid to HST map file dimensions - SMH 1/15/2025
    gravity_patch=mp.make_patch(grav,LatLims,LonLims,180,180,pad=True)
    if dataversion==2:
        scale=1
    elif dataversion=='H':
        scale=20

    CH4_Cloud_Press = CH4_Ncol*amagat*zoom(gravity_patch,scale)*mean_mol_wt/(fCH4*STP)
   
    return CH4_Cloud_Press/2

def computeAmmoniaMoleFraction(CH4data, NH3data):
    #compute ammonia optical depth according to Beer-Lambert law
    CH4_tau = -np.log(CH4data)
    NH3_tau = -np.log(NH3data) 
  

    #compute column abundance of molecules along the line-of-sight where 1000 corrects for units
    #  and K_eff are the effective absorption cross sections of each molecule, 0.427 and 2.955 for 
    #  methane and ammonia respectively
    CH4_Ncol = 1000*CH4_tau/K_eff_CH4619 
    NH3_Ncol = 1000*NH3_tau/K_eff_NH3645
    fNH3=fCH4*(NH3_Ncol/CH4_Ncol)*1.0e6
    
    return fNH3

def get_HSTGO_filenames(obskeyHST):
    import socket
    hostname = socket.gethostname()
    from config_VA import Host_path
    pathHST=Host_path[hostname]+'/HST GO 18055/'+obskeyHST[:-1]+'/'+obskeyHST+'/unprocessed_L1/'

    if obskeyHST=='20251016UTa':
        fn619={'fn':'251015_619_2347_reg.fits','offset':[0,0]} #-9,0
        fn631={'fn':'251015_631_2352_reg.fits','offset':[0,0]}
        fn645={'fn':'251015_645_2354_reg.fits','offset':[0,0]} #-14,0
        fn275={'fn':'251015_275_2356_reg.fits','offset':[0,0]} #14,0
        fn395={'fn':'251016_395_0000_reg.fits','offset':[0,0]}
        fn502={'fn':'251016_502_0003_reg.fits','offset':[0,0]}
        fn673={'fn':'251015_673_2359_reg.fits','offset':[0,0]} #-13,0
        fn727={'fn':'251016_727_0010_reg.fits','offset':[0,0]} #7,0
        fn889={'fn':'251015_889_2339_reg.fits','offset':[0,0]}
    if obskeyHST=='20251016UTb':
        fn619={'fn':'251016_619_0015_reg.fits','offset':[0,0]} #-9,0
        fn631={'fn':'251016_631_0019_reg.fits','offset':[0,0]}
        fn645={'fn':'251016_645_0021_reg.fits','offset':[0,0]} #-14,0
        fn275={'fn':'251016_275_0005_reg.fits','offset':[0,0]} #14,0
        fn395={'fn':'251016_395_0000_reg.fits','offset':[0,0]}
        fn502={'fn':'251016_502_0003_reg.fits','offset':[0,0]}
        fn673={'fn':'251016_673_0007_reg.fits','offset':[0,0]} #-13,0
        fn727={'fn':'251016_727_0010_reg.fits','offset':[0,0]} #7,0
        fn889={'fn':'251015_889_2339_reg.fits','offset':[0,0]}
    if obskeyHST=='20251016UTc':
        fn619={'fn':'251016_619_0256_reg.fits','offset':[0,0]}
        fn631={'fn':'251016_631_0300_reg.fits','offset':[0,0]}
        fn645={'fn':'251016_645_0302_reg.fits','offset':[0,0]}
        fn275={'fn':'251016_275_0305_reg.fits','offset':[0,0]}#
        fn395={'fn':'251016_395_0309_reg.fits','offset':[0,0]}
        fn502={'fn':'251016_502_0311_reg.fits','offset':[0,0]}
        fn673={'fn':'251016_673_0307_reg.fits','offset':[0,0]} #4,0
        fn727={'fn':'251016_727_0318_reg.fits','offset':[0,0]} #5,0
        fn889={'fn':'251016_889_0251_reg.fits','offset':[0,0]} #0,-5
    if obskeyHST=='20251016UTd':
        fn619={'fn':'251016_619_0323_reg.fits','offset':[0,0]}
        fn631={'fn':'251016_631_0328_reg.fits','offset':[0,0]}
        fn645={'fn':'251016_645_0330_reg.fits','offset':[0,0]}
        fn275={'fn':'251016_275_0313_reg.fits','offset':[0,0]}#
        fn395={'fn':'251016_395_0309_reg.fits','offset':[0,0]}
        fn502={'fn':'251016_502_0311_reg.fits','offset':[0,0]}
        fn673={'fn':'251016_673_0316_reg.fits','offset':[0,0]} #4,0
        fn727={'fn':'251016_727_0318_reg.fits','offset':[0,0]} #5,0
        fn889={'fn':'251016_889_0251_reg.fits','offset':[0,0]} #0,-5
    if obskeyHST=='20251016UTe':
        fn619={'fn':'251016_619_0430_reg.fits','offset':[0,0]}
        fn631={'fn':'251016_631_0435_reg.fits','offset':[0,0]}
        fn645={'fn':'251016_645_0437_reg.fits','offset':[0,0]} ##8,0
        fn275={'fn':'251016_275_0439_reg.fits','offset':[0,0]}#
        fn395={'fn':'251016_395_0444_reg.fits','offset':[0,0]}
        fn502={'fn':'251016_502_0446_reg.fits','offset':[0,0]}
        fn673={'fn':'251016_673_0442_reg.fits','offset':[0,0]} #4,0
        fn727={'fn':'251016_727_0453_reg.fits','offset':[0,0]} #5,0
        fn889={'fn':'251016_889_0426_reg.fits','offset':[0,0]} #0,-5
    if obskeyHST=='20251016UTf':
        fn619={'fn':'251016_619_0458_reg.fits','offset':[0,0]}
        fn631={'fn':'251016_631_0502_reg.fits','offset':[0,0]}
        fn645={'fn':'251016_645_0504_reg.fits','offset':[0,0]} ##8,0
        fn275={'fn':'251016_275_0448_reg.fits','offset':[0,0]}#
        fn395={'fn':'251016_395_0444_reg.fits','offset':[0,0]}
        fn502={'fn':'251016_502_0446_reg.fits','offset':[0,0]}
        fn673={'fn':'251016_673_0450_reg.fits','offset':[0,0]} #4,0
        fn727={'fn':'251016_727_0453_reg.fits','offset':[0,0]} #5,0
        fn889={'fn':'251016_889_0426_reg.fits','offset':[0,0]} #0,-5
        
    if obskeyHST=='20251120UTa':
        fn619={'fn':'251120_619_0932_reg.fits','offset':[0,0]}
        fn631={'fn':'251120_631_0948_reg.fits','offset':[0,0]}
        fn645={'fn':'251120_645_0951_reg.fits','offset':[0,0]} ##8,0
        fn275={'fn':'251120_275_0944_reg.fits','offset':[0,0]}#
        fn395={'fn':'251120_395_0920_reg.fits','offset':[0,0]}
        fn502={'fn':'251120_502_0924_reg.fits','offset':[0,0]}
        fn673={'fn':'251120_673_0917_reg.fits','offset':[0,0]} #4,0
        fn727={'fn':'251120_727_0912_reg.fits','offset':[0,0]} #5,0
        fn889={'fn':'251120_889_0928_reg.fits','offset':[0,0]} #0,-5
    if obskeyHST=='20251120UTb':
        fn619={'fn':'251120_619_1414_reg.fits','offset':[0,0]}
        fn631={'fn':'251120_631_1430_reg.fits','offset':[0,0]}
        fn645={'fn':'251120_645_1434_reg.fits','offset':[0,0]} ##8,0
        fn275={'fn':'251120_275_1426_reg.fits','offset':[0,0]}#
        fn395={'fn':'251120_395_1402_reg.fits','offset':[0,0]}
        fn502={'fn':'251120_502_1406_reg.fits','offset':[0,0]}
        fn673={'fn':'251120_673_1359_reg.fits','offset':[0,0]} #4,0
        fn727={'fn':'251120_727_1354_reg.fits','offset':[0,0]} #5,0
        fn889={'fn':'251120_889_1410_reg.fits','offset':[0,0]} #0,-5
    if obskeyHST=='20251120UTc':
        fn619={'fn':'251120_619_1723_reg.fits','offset':[0,0]}
        fn631={'fn':'251120_631_1738_reg.fits','offset':[0,0]}
        fn645={'fn':'251120_645_1742_reg.fits','offset':[0,0]} ##8,0
        fn275={'fn':'251120_275_1735_reg.fits','offset':[0,0]}#
        fn395={'fn':'251120_395_1711_reg.fits','offset':[0,0]}
        fn502={'fn':'251120_502_1715_reg.fits','offset':[0,0]}
        fn673={'fn':'251120_673_1707_reg.fits','offset':[0,0]} #4,0
        fn727={'fn':'251120_727_1703_reg.fits','offset':[0,0]} #5,0
        fn889={'fn':'251120_889_1719_reg.fits','offset':[0,0]} #0,-5
    if obskeyHST=='20251221UTa':
        fn619={'fn':'251221_619_1850_reg.fits','offset':[0,0]}
        fn631={'fn':'251221_631_1905_reg.fits','offset':[0,0]}
        fn645={'fn':'251221_645_1909_reg.fits','offset':[0,0]} ##8,0
        fn275={'fn':'251221_275_1902_reg.fits','offset':[0,0]}#
        fn395={'fn':'251221_395_1838_reg.fits','offset':[0,0]}
        fn502={'fn':'251221_502_1842_reg.fits','offset':[0,0]}
        fn673={'fn':'251221_673_1834_reg.fits','offset':[0,0]} #4,0
        fn727={'fn':'251221_727_1830_reg.fits','offset':[0,0]} #5,0
        fn889={'fn':'251221_889_1846_reg.fits','offset':[0,0]} #0,-5
    if obskeyHST=='20251221UTb':
        fn619={'fn':'251222_619_1518_reg.fits','offset':[0,0]}
        fn631={'fn':'251222_631_1534_reg.fits','offset':[0,0]}
        fn645={'fn':'251222_645_1537_reg.fits','offset':[0,0]} ##8,0
        fn275={'fn':'251222_275_1530_reg.fits','offset':[0,0]}#
        fn395={'fn':'251222_395_1506_reg.fits','offset':[0,0]}
        fn502={'fn':'251222_502_1510_reg.fits','offset':[0,0]}
        fn673={'fn':'251222_673_1503_reg.fits','offset':[0,0]} #4,0
        fn727={'fn':'251222_727_1458_reg.fits','offset':[0,0]} #5,0
        fn889={'fn':'251222_889_1514_reg.fits','offset':[0,0]} #0,-5
    if obskeyHST=='20251221UTc':
        fn619={'fn':'251223_619_1622_reg.fits','offset':[0,0]}
        fn631={'fn':'251223_631_1638_reg.fits','offset':[0,0]}
        fn645={'fn':'251223_645_1641_reg.fits','offset':[0,0]} ##8,0
        fn275={'fn':'251223_275_1634_reg.fits','offset':[0,0]}#
        fn395={'fn':'251223_395_1610_reg.fits','offset':[0,0]}
        fn502={'fn':'251223_502_1614_reg.fits','offset':[0,0]}
        fn673={'fn':'251223_673_1607_reg.fits','offset':[0,0]} #4,0
        fn727={'fn':'251223_727_1602_reg.fits','offset':[0,0]} #5,0
        fn889={'fn':'251223_889_1618_reg.fits','offset':[0,0]} #0,-5
        
    return pathHST,fn619,fn631,fn645,fn275,fn395,fn502,fn673,fn727,fn889
        
def read_HSTGO_fits(path,filename,LonSys,plot=True,dataunit=0):
    import sys
    sys.path.append('..\Maps')
    import convert_system3_to_I_II_spice as clong
    
    hdulist=fits.open(path+filename['fn'])
    hdulist.info()
    hdr=hdulist[0].header
    if 'PHOTIF' in hdr:
        data=hdulist[dataunit].data*hdr['PHOTIF']
    else:
        data=hdulist[dataunit].data*4.70214815131e-05
    emi=hdulist[1].data
    inc=hdulist[1].data
    hdulist.close()
    if '_275_' in filename['fn']:
        #k=0.65
        k=1.0
    else:
        k=1.0
    print("##################################################################")
    print("filename=",filename)
    print("k= ",k)
    print("##################################################################")
    minnaert=(np.cos(np.deg2rad(inc))**k)*(np.cos(np.deg2rad(emi))**(k-1.0))
    data=data/minnaert
    #import matplotlib.pyplot as pl
    #pl.imshow(minnaert)
    
    """
    print("##########data.shape=",data.shape)
    print("##########hdr['LAT_TOP']=",hdr['LAT_TOP'])

    print("##########hdr['LAT_STEP']=",hdr['LAT_STEP'])
    print("##########hdr['LAT_BOT']=",hdr['LAT_BOT'])
    print("##########hdr['NAXIS1']=",hdr['NAXIS1'])
    print("##########hdr['NAXIS2']=",hdr['NAXIS2'])
    print("##########hdr['LON_LEFT']=",hdr['LON_LEFT'])
    print("##########hdr['LON_STEP']=",hdr['LON_STEP'])
    print("##########hdr['LON_RIGH']=",hdr['LON_RIGH'])
    """
    HST=np.zeros((3600,7200))
    INC=np.zeros((3600,7200))
    EMI=np.zeros((3600,7200))
    leftindx=int((360-hdr['LON_LEFT'])/(-hdr['LON_STEP'])) % 7200
    loncount=hdr['NAXIS1']
    rightindx=(leftindx+loncount) % 7200
    #topindx=int((90-hdr['LAT_TOP'])/hdr['LAT_STEP'])   
    topindx=int((90+hdr['LAT_BOT'])/hdr['LAT_STEP'])   
    latcount=hdr['NAXIS2']
    """
    print()
    print("----------topindx=",topindx)
    print("----------latcount=",latcount)
    print("----------leftindx=",leftindx)
    print("----------loncount=",loncount)
    """
    if rightindx<leftindx:
        HST[topindx:topindx+latcount,leftindx:7200]=data[:,0:loncount-rightindx]
        HST[topindx:topindx+latcount,0:rightindx]=data[:,loncount-rightindx:loncount]
        INC[topindx:topindx+latcount,leftindx:7200]=inc[:,0:loncount-rightindx]
        INC[topindx:topindx+latcount,0:rightindx]=inc[:,loncount-rightindx:loncount]
        EMI[topindx:topindx+latcount,leftindx:7200]=emi[:,0:loncount-rightindx]
        EMI[topindx:topindx+latcount,0:rightindx]=emi[:,loncount-rightindx:loncount]
    else:
        HST[topindx:topindx+latcount,leftindx:rightindx]=data
        INC[topindx:topindx+latcount,leftindx:rightindx]=inc
        EMI[topindx:topindx+latcount,leftindx:rightindx]=emi
        
    HST=np.flipud(HST)
    INC=np.flipud(INC)
    EMI=np.flipud(EMI)
    #HSTy=np.roll(HST,filename['offset'][0],axis=0)
    #HSTx=np.roll(HSTy,filename['offset'][1],axis=1)
    
    #Roll to LonSys
    CM3=hdr['TRG_LON']
    dateobs=hdr["date-obs"]+"T"+hdr["time-obs"]
    longitudes=clong.convert_system3_to_I_II_spice(dateobs, CM3)
    #print("########### ROLL ##########")
    #print(hdr["date-obs"])
    #print("LonSys,CM3,Sys I Long=",LonSys,CM3,longitudes["System I"])

    if LonSys=='1':
        roll=CM3-longitudes["System I"]
    elif LonSys=='2':
        roll=CM3-longitudes["System II"]
    else:
        roll=0

    print(filename,roll)
    print("########### ROLL ##########")

    HST_rolled=np.roll(HST,int(roll)*20,axis=1)
    INC_rolled=np.roll(INC,int(roll)*20,axis=1)
    EMI_rolled=np.roll(EMI,int(roll)*20,axis=1)
    
    if plot:
        fig,ax=pl.subplots(1,figsize=(8,6), dpi=150, facecolor="white")
        ax.imshow(HST_rolled,cmap='gray')
        ax.set_title(filename['fn'][-24:])
    
    return HST_rolled,hdr,INC_rolled,EMI_rolled

#def make_L2_HSTGO_abs_data(pathHST,fn619,fn631,fn645,LonSys,plot=True):
    
def make_L2_HSTGO_abs_data(pathHST,hdr619,HST619,hdr631,HST631,
                           hdr645,HST645,LonSys,LatLims,LonLims,
                           dataversion=2,plot=False):
    """
    Get Minnaert-corrected reflectivity data from the raw radiance files
    provided by Mike Wong. Compute the the ammonia and methane absorption.
    Then compute the cloud pressure and ammonia mole fraction.


    Parameters
    ----------
    pathHST : TYPE
        DESCRIPTION.
    hdr619 : TYPE
        DESCRIPTION.
    HST619 : TYPE
        DESCRIPTION.
    hdr631 : TYPE
        DESCRIPTION.
    HST631 : TYPE
        DESCRIPTION.
    hdr645 : TYPE
        DESCRIPTION.
    HST645 : TYPE
        DESCRIPTION.
    LonSys : TYPE
        DESCRIPTION.
    LatLims : TYPE
        DESCRIPTION.
    LonLims : TYPE
        DESCRIPTION.
    dataversion : TYPE, optional
        DESCRIPTION. The default is 2.
    plot : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    PCld : TYPE
        DESCRIPTION.
    fNH3 : TYPE
        DESCRIPTION.
    PCldhdr : TYPE
        DESCRIPTION.
    fNH3hdr : TYPE
        DESCRIPTION.
    CH4abs : TYPE
        DESCRIPTION.
    NH3abs : TYPE
        DESCRIPTION.
    CH4abshdr : TYPE
        DESCRIPTION.
    NH3abshdr : TYPE
        DESCRIPTION.

    """

    import copy
    import process_L1Y_helpers as hp
    
    ###########################################################################
    # Compute methane and ammonia absorption, CH4abs and NH3abs, respectively.
    # and update the header information
    ###########################################################################
    CH4abs=HST619/HST631
    
    CH4abshdr = copy.deepcopy(hdr619)
    CH4abshdr["DATE-OBS"] = hp.averageDates(hdr619["DATE-OBS"]+'T'+hdr619["TIME-OBS"],
                                          hdr631["DATE-OBS"]+'T'+hdr631["TIME-OBS"], 
                                          "%Y-%m-%dT%H:%M:%S.%f")
    hp.averageHdrNum(CH4abshdr, hdr619, hdr631, "TRG_LAT")
    hp.averageHdrNum(CH4abshdr, hdr619, hdr631, "TRG_LON")
    hp.averageHdrNum(CH4abshdr, hdr619, hdr631, "SUN_LAT")
    hp.averageHdrNum(CH4abshdr, hdr619, hdr631, "SUN_LON")
    #hp.averageHdrNum(PCldhdr, hdr619, hdr631, "HIERARCH PLANMAP LIGHT-TIME")
    hp.averageHdrNum(CH4abshdr, hdr619, hdr631, "TRG_D")
    CH4abshdr["BUNIT"]='619nm Methane Transmission'
    #del CH4abshdr["MISSVAL"]

    NH3abs=HST645/HST631
    NH3abshdr = copy.deepcopy(hdr645)
    NH3abshdr["DATE-OBS"] = hp.averageDates(hdr645["DATE-OBS"]+'T'+hdr645["TIME-OBS"],
                                          hdr631["DATE-OBS"]+'T'+hdr631["TIME-OBS"], "%Y-%m-%dT%H:%M:%S.%f")
    hp.averageHdrNum(NH3abshdr, hdr645, hdr631, "TRG_LAT")
    hp.averageHdrNum(NH3abshdr, hdr645, hdr631, "TRG_LON")
    hp.averageHdrNum(NH3abshdr, hdr645, hdr631, "SUN_LAT")
    hp.averageHdrNum(NH3abshdr, hdr645, hdr631, "SUN_LON")
    #hp.averageHdrNum(fNH3hdr, hdr619, hdr631, "HIERARCH PLANMAP LIGHT-TIME")
    hp.averageHdrNum(NH3abshdr, hdr645, hdr631, "TRG_D")
    NH3abshdr["BUNIT"]='Mole Fraction'
    #del NH3abshdr["MISSVAL"]

    ###########################################################################
    # Compute cloud pressure, PCld, and ammonia mole fraction, fNH3

    PCld=computeCloudPressure(CH4abs,LatLims,LonLims,dataversion=dataversion)
    PCldhdr = copy.deepcopy(CH4abshdr)
    PCldhdr["BUNIT"]='Cloud-top Press'
    #del PCldhdr["MISSVAL"]

    fNH3=computeAmmoniaMoleFraction(CH4abs, NH3abs)
    fNH3hdr = copy.deepcopy(hdr645)
    fNH3hdr["DATE-OBS"] = hp.averageDates(hdr645["DATE-OBS"]+'T'+hdr645["TIME-OBS"],
                                          hdr631["DATE-OBS"]+'T'+hdr631["TIME-OBS"], "%Y-%m-%dT%H:%M:%S.%f")
    hp.averageHdrNum(fNH3hdr, hdr645, hdr631, "TRG_LAT")
    hp.averageHdrNum(fNH3hdr, hdr645, hdr631, "TRG_LON")
    hp.averageHdrNum(fNH3hdr, hdr645, hdr631, "SUN_LAT")
    hp.averageHdrNum(fNH3hdr, hdr645, hdr631, "SUN_LON")
    #hp.averageHdrNum(fNH3hdr, hdr619, hdr631, "HIERARCH PLANMAP LIGHT-TIME")
    hp.averageHdrNum(fNH3hdr, hdr645, hdr631, "TRG_D")
    fNH3hdr["BUNIT"]='Mole Fraction'
    #del fNH3hdr["MISSVAL"]

    return PCld,fNH3,PCldhdr,fNH3hdr,CH4abs,NH3abs,CH4abshdr,NH3abshdr

def make_L2_HSTGO_AOI_CI(pathHST,hdr275,HST275,hdr889,HST889,
                         hdr395,HST395,hdr631,HST631,LonSys,plot=True):
    """
    

    Parameters
    ----------
    pathHST : TYPE
        DESCRIPTION.
    hdr275 : TYPE
        DESCRIPTION.
    HST275 : TYPE
        DESCRIPTION.
    hdr889 : TYPE
        DESCRIPTION.
    HST889 : TYPE
        DESCRIPTION.
    hdr395 : TYPE
        DESCRIPTION.
    HST395 : TYPE
        DESCRIPTION.
    hdr631 : TYPE
        DESCRIPTION.
    HST631 : TYPE
        DESCRIPTION.
    LonSys : TYPE
        DESCRIPTION.
    plot : TYPE, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    AOI : TYPE
        DESCRIPTION.
    CI : TYPE
        DESCRIPTION.

    """
    import copy
    import process_L1Y_helpers as hp

    #!!!! Need to update headers here including adding input headers to call
    AOI=HST889/HST275
    AOIhdr = copy.deepcopy(hdr889)
    
    AOIhdr["DATE-OBS"] = hp.averageDates(hdr889["DATE-OBS"]+'T'+hdr889["TIME-OBS"],
                                          hdr275["DATE-OBS"]+'T'+hdr275["TIME-OBS"], 
                                          "%Y-%m-%dT%H:%M:%S.%f")
    hp.averageHdrNum(AOIhdr, hdr889, hdr275, "TRG_LAT")
    hp.averageHdrNum(AOIhdr, hdr889, hdr275, "TRG_LON")
    hp.averageHdrNum(AOIhdr, hdr889, hdr275, "SUN_LAT")
    hp.averageHdrNum(AOIhdr, hdr889, hdr275, "SUN_LON")
    #hp.averageHdrNum(PCldhdr, hdr619, hdr631, "HIERARCH PLANMAP LIGHT-TIME")
    hp.averageHdrNum(AOIhdr, hdr889, hdr275, "TRG_D")
    AOIhdr["BUNIT"]='Altitude Opacity Index'
    #del CH4abshdr["MISSVAL"]
    
    CI=HST395/HST631
    CIhdr = copy.deepcopy(hdr631)
    CIhdr["DATE-OBS"] = hp.averageDates(hdr631["DATE-OBS"]+'T'+hdr631["TIME-OBS"],
                                          hdr395["DATE-OBS"]+'T'+hdr395["TIME-OBS"], 
                                          "%Y-%m-%dT%H:%M:%S.%f")
    hp.averageHdrNum(CIhdr, hdr631, hdr395, "TRG_LAT")
    hp.averageHdrNum(CIhdr, hdr631, hdr395, "TRG_LON")
    hp.averageHdrNum(CIhdr, hdr631, hdr395, "SUN_LAT")
    hp.averageHdrNum(CIhdr, hdr631, hdr395, "SUN_LON")
    #hp.averageHdrNum(PCldhdr, hdr619, hdr631, "HIERARCH PLANMAP LIGHT-TIME")
    hp.averageHdrNum(CIhdr, hdr631, hdr395, "TRG_D")
    AOIhdr["BUNIT"]='Color Index'
    #del CH4abshdr["MISSVAL"]
    if plot:
        fig,ax=pl.subplots(1,figsize=(8,6), dpi=150, facecolor="white")
        ax.imshow(AOI)

    return AOI,CI

def normalizeBrightness(radianceArr, emissionArr):
    #From Leah Tiktin L1Y helpers code
   
    #mask = (0 < emissionArr) & (emissionArr < 75)
    mask = emissionArr < 60
    print("normalizeBrightness: radianceArr.shape, emissionArr.shape, mask.shape=",
          radianceArr.shape, emissionArr.shape, mask.shape)
    radianceMax = radianceArr[mask].max()
    radianceMean=radianceArr[mask].mean()
    #normRadiance = radianceArr / radianceMax
    normRadiance = radianceArr / radianceMean
  
    return np.array(normRadiance)

def make_L2_HSTGO_RGB_data(pathHST,dataR,dataG,dataB,emiR,emiG,emiB,LonSys,plot=True):
    
    #HSTRnorm=normalizeBrightness(dataR,emiR)
    #HSTGnorm=normalizeBrightness(dataG,emiG)
    #HSTBnorm=normalizeBrightness(dataB,emiB)
    #RGB=np.stack((HSTRnorm,HSTGnorm,HSTBnorm),axis=2)
    RGB=np.stack((dataR,dataG,dataB),axis=2)
    if plot:
        fig,ax=pl.subplots(1,figsize=(8,6), dpi=150, facecolor="white")
        ax.imshow(RGB)
        #ax.set_title(filename['fn'][-24:])

    return RGB

def plot_HST_global_maps(CH4abs,NH3abs):
    figCH4,axsCH4=pl.subplots(1,figsize=(8,6), dpi=150, facecolor="white")
    axsCH4.imshow(CH4abs,cmap='Blues_r',vmin=1400,vmax=3600)
    axsCH4.set_title("CH4 Absorption")
    
    figNH3,axsNH3=pl.subplots(1,figsize=(8,6), dpi=150, facecolor="white")
    axsNH3.imshow(CH4abs/NH3abs,cmap='terrain_r',vmin=0,vmax=200)
    #axsNH3.imshow(NH3abs,cmap='Blues_r',vmin=1.1,vmax=1.6)
    axsNH3.set_title("NH3 Absorption")

def plot_HSTGO_abs_patches(obskeyHST,LatLims,LonLims,CH4abs,NH3abs,axsCH4abspatch,
                           axsNH3patch,minmax1=[1400,3600],minmax2=[0,200],
                           cbar_reverse1=True,cbar_reverse2=False,
                           titles=['fNH3 (ppm)','PCloud (mbar)'],
                           ctbls=['Blues','terrain_r']):

    ###########################################################################
    # MAKE PATCH
    ###########################################################################
    CH4patch=mp.make_patch(CH4abs,LatLims,LonLims,180,180,pad=True)
    CH4patchflat=fp.flatten_patch(CH4patch)
    NH3patch=mp.make_patch(NH3abs,LatLims,LonLims,180,180,pad=True)
    NH3patchflat=fp.flatten_patch(NH3patch)

    ###########################################################################
    # PLOT SECTION
    ###########################################################################   
    cbttl="Mean="+str(np.mean(CH4patchflat))[:4]+" $\pm$ "+str(np.std(CH4patchflat))[:3]
    pp.plot_patch(CH4patchflat,LatLims,LonLims,180,180,ctbls[0],axsCH4abspatch,
                   cbarplot=True,cbar_title=cbttl,cbar_reverse=cbar_reverse1,
                   vn=minmax1[0],vx=minmax1[1])  
    axsCH4abspatch.set_title(titles[1],fontsize=10)
        
    cbttl="Mean="+str(np.mean(NH3patchflat))[:3]+" $\pm$ "+str(np.std(NH3patchflat))[:2]
    pp.plot_patch(NH3patchflat,LatLims,LonLims,180,180,ctbls[1],axsNH3patch,
                   cbarplot=True,cbar_title=cbttl,cbar_reverse=cbar_reverse2,
                   vn=minmax2[0],vx=minmax2[1])
    axsNH3patch.set_title(titles[0],fontsize=10)
    
    return CH4patchflat,NH3patchflat
    
def plot_HSTGO_RGB_patches(obskeyHST,LatLims,LonLims,RGBpatch,axsRGBpatch,
                           title='RGB (673/502/395)'):
    RGB4Display=np.power(np.array(RGBpatch).astype(float),1.3)
    #RGB4Display=RGB4Display/RGB4Display.max()
    show=axsRGBpatch.imshow(RGB4Display,
               extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],
                       90-LatLims[0]],
                       aspect="equal")
    im_ratio = RGBpatch.shape[0]/RGBpatch.shape[1]
    cbar = pl.colorbar(show, 
               orientation='vertical',cmap='gist_heat',
               ax=axsRGBpatch,fraction=0.046*im_ratio, pad=0.04)
    #cbar.ax.set_yticklabels(np.around(tx,3))
    cbar.ax.tick_params(labelsize=6,color="k")#if iSession >1:
    #cbar.ax.set_ylabel(cbar_title,size=8)#,loc="top")
    #cbar.ax.yaxis.set_label_coords(-0.7, 0.5)
    cbar.ax.set_visible(False)

    axsRGBpatch.set_title(title,fontsize=10)
    
    
    
def plot_HSTGO_abs_3dpatches(obskeyHST,LatLims,LonLims,CH4abspatchflat,NH3patchflat):
    fig3dCld, ax3dCld = pl.subplots(1,figsize=(12,8),subplot_kw={"projection": "3d"})
    fig3dboth, ax3dboth = pl.subplots(1,figsize=(12,8),subplot_kw={"projection": "3d"})
    #fig3dNH3, ax3dNH3 = pl.subplots(1,figsize=(12,8),subplot_kw={"projection": "3d"})
    # Make data.
    Lons = np.arange((360-LonLims[0])/0.05,(360-LonLims[1])/0.05,-1.0)
    Lats = np.arange((90-LatLims[0])/0.05,(90-LatLims[1])/0.05, -1.0)
    print(Lats)
    X, Y = np.meshgrid(Lons, Lats)
    Z = CH4abspatchflat
    W = NH3patchflat
    
    norm = colors.Normalize(vmin=0.60, vmax=0.70)
    cmap = cm.get_cmap('terrain_r') 
    facecolors = cmap(norm(W))
    
    ###############################################################################
    # Plot Pressure on Pressure surface
    surfPCld = ax3dCld.plot_surface(X, Y, Z, cmap="Blues_r",linewidth=0, 
                                    antialiased=False,
                                    vmin=0.78,vmax=0.92)
    ax3dCld.set_zlim(0.5,1.2)
    ax3dCld.invert_xaxis()
    ax3dCld.set_box_aspect((np.ptp(Lons), np.ptp(Lats), 30/0.05))
    ax3dCld.view_init(45, -60, 0) 
    ax3dCld.set_xlabel(' Longitude deg')
    ax3dCld.set_ylabel('PG Latitude (deg)')
    
    ###############################################################################
    # Plot Ammonia on Pressure surface
    surfboth = ax3dboth.plot_surface(X, Y, Z, facecolors=facecolors, 
                                     rstride=1, cstride=1, antialiased=False)
    ax3dboth.set_zlim(0.5,1.2)
    ax3dboth.invert_xaxis()
    ax3dboth.set_box_aspect((np.ptp(Lons), np.ptp(Lats), 30/0.05))
    ax3dboth.view_init(45, -60, 0) 
    ax3dboth.set_xlabel(' Longitude deg')
    ax3dboth.set_ylabel('PG Latitude (deg)')

def write_HST_fits_patch(obskeyHST,LonSys,patchflat,hdr,pathout,patchtype,
                         LatLims=False,LonLims=False,dmin=0.0,dmax=1.0):
    """
    If LatLims and LonLims are provided, then the full map array is saved,
    but only the patch area has data. This is done as a 16 bit int.
    Otherwise, the FITS file is of the patch only and is a 64-bit float.

    Parameters
    ----------
    obskeyHST : TYPE
        DESCRIPTION.
    LonSys : TYPE
        DESCRIPTION.
    patchflat : TYPE
        DESCRIPTION.
    pathout : TYPE
        DESCRIPTION.
    patchtype : TYPE
        DESCRIPTION.
    LatLims : TYPE, optional
        DESCRIPTION. The default is False.
    LonLims : TYPE, optional
        DESCRIPTION. The default is False.
    dmin : TYPE, optional
        DESCRIPTION. The default is 0.0.
    dmax : TYPE, optional
        DESCRIPTION. The default is 1.0.

    Returns
    -------
    None.

    """
    import os
    print("DEBUG DEBUG DEBUG LATLIMS LONLIMS")
    print(LatLims, LonLims)
    if LatLims and LonLims:
        fullmap=np.zeros((3600,7200))
        ylims=np.array(LatLims)*20
        xlims=np.array(LonLims)*20
        print(xlims,ylims)
        print(patchflat.shape)
        fullmap[ylims[0]:ylims[1],xlims[0]:xlims[1]]=patchflat
        hdudata = fits.PrimaryHDU(data=fullmap,header=hdr)
        bscale = (dmax - dmin) / 65535.0
        bzero  = dmin + 32768 * bscale   
        hdudata.scale('int16', bscale=bscale, bzero=bzero)
    else:
        hdudata = fits.PrimaryHDU(patchflat)
      
    hdul = fits.HDUList([hdudata])
    print(hdr)

    fnout=pathout+'/'+obskeyHST+' HST Sys'+ LonSys +' '+patchtype+'.fits'
    try:
        os.remove(fnout)
    except: 
        print("file doesn't exist")
    hdul.writeto(fnout)
    hdul.close()

###############################################################################
def HSTGO_process_and_plot(obskeyHST,LatLims,LonLimsInput,LonSys='3',
                           cont=False,waveplot=False):
    
    import os
    import sys
    from scipy.ndimage import gaussian_filter
    sys.path.append('C:/Astronomy/Projects/SAS 2021 Ammonia/Visualization-and-Analysis/')

    import L4_Jup_Map_Plot_V2 as L4MP 
    #import figsize_and_aspect,set_up_figure
    from astropy.io import fits
    
    LonLims=[(360-LonLimsInput[1])%360,(360-LonLimsInput[0])%360]
    if LonLimsInput[0]==0:
        LonLims[1]=360
    print(LonLims)
    
    # Get file names and navigation adjustments
    pathHST,fn619,fn631,fn645,fn275,fn395,fn502,fn673,fn727,fn889=\
        get_HSTGO_filenames(obskeyHST)

    # Make raw data patches
    data275,hdr275,inc275,emi275=read_HSTGO_fits(pathHST,fn275,LonSys,plot=False,dataunit=0)
    data275patch=mp.make_patch(data275,LatLims,LonLims,180,180,pad=True)
    emi275patch=mp.make_patch(emi275,LatLims,LonLims,180,180,pad=True)
    #data275patchflat=fp.flatten_patch(data275patch)
    del hdr275["MISSVAL"]
    if fn395 !=' ':
        data395,hdr395,inc395,emi395=read_HSTGO_fits(pathHST,fn395,LonSys,plot=False,dataunit=0)
        data395patch=mp.make_patch(data395,LatLims,LonLims,180,180,pad=True)
        emi395patch=mp.make_patch(emi395,LatLims,LonLims,180,180,pad=True)
        #data395patchflat=fp.flatten_patch(data395patch)
        del hdr395["MISSVAL"]
    if fn502!=' ':
        data502,hdr502,inc502,emi502=read_HSTGO_fits(pathHST,fn502,LonSys,plot=False,dataunit=0)
        data502patch=mp.make_patch(data502,LatLims,LonLims,180,180,pad=True)
        emi502patch=mp.make_patch(emi502,LatLims,LonLims,180,180,pad=True)
        #data502patchflat=fp.flatten_patch(data502patch)
        del hdr502["MISSVAL"]
    data619,hdr619,inc619,emi619=read_HSTGO_fits(pathHST,fn619,LonSys,plot=False,dataunit=0)
    data619patch=mp.make_patch(data619,LatLims,LonLims,180,180,pad=True)
    emi619patch=mp.make_patch(emi619,LatLims,LonLims,180,180,pad=True)
    data619patchflat=fp.flatten_patch(data619patch)
    del hdr619["MISSVAL"]
    data631,hdr631,inc631,emi631=read_HSTGO_fits(pathHST,fn631,LonSys,plot=False,dataunit=0)
    data631patch=mp.make_patch(data631,LatLims,LonLims,180,180,pad=True)
    emi631patch=mp.make_patch(emi631,LatLims,LonLims,180,180,pad=True)
    data631patchflat=fp.flatten_patch(data631patch)
    del hdr631["MISSVAL"]
    data645,hdr645,inc645,emi645=read_HSTGO_fits(pathHST,fn645,LonSys,plot=False,dataunit=0)
    data645patch=mp.make_patch(data645,LatLims,LonLims,180,180,pad=True)
    emi645patch=mp.make_patch(emi645,LatLims,LonLims,180,180,pad=True)
    data645patchflat=fp.flatten_patch(data645patch)
    del hdr645["MISSVAL"]
    data673,hdr673,inc673,emi673=read_HSTGO_fits(pathHST,fn673,LonSys,plot=False,dataunit=0)
    data673patch=mp.make_patch(data673,LatLims,LonLims,180,180,pad=True)
    emi673patch=mp.make_patch(emi673,LatLims,LonLims,180,180,pad=True)
    #data673patchflat=fp.flatten_patch(data673patch)
    del hdr673["MISSVAL"]
    if fn727!=' ':
        data727,hdr727,inc727,emi727=read_HSTGO_fits(pathHST,fn727,LonSys,plot=False,dataunit=0)
        data727patch=mp.make_patch(data727,LatLims,LonLims,180,180,pad=True)
        emi727patch=mp.make_patch(emi727,LatLims,LonLims,180,180,pad=True)
        #data727patchflat=fp.flatten_patch(data727patch)
        del hdr727["MISSVAL"]
    if fn889!=' ':
        data889,hdr889,inc889,emi889=read_HSTGO_fits(pathHST,fn889,LonSys,plot=False,dataunit=0)
        data889patch=mp.make_patch(data889,LatLims,LonLims,180,180,pad=True)
        emi889patch=mp.make_patch(emi889,LatLims,LonLims,180,180,pad=True)
        #data889patchflat=fp.flatten_patch(data889patch)
        del hdr889["MISSVAL"]
    
    ###########################################################################
    # Compute PCld and fNH3 from Mike Wongs radiance data with a Minnaert 
    # correction applied.
    #PCld,fNH3,PCldhdr,fNH3hdr,CH4abs,NH3abs,CH4abshdr,NH3abshdr = \
    #    make_L2_HSTGO_abs_data(pathHST,fn619,fn631,fn645,LonSys,plot=False)
    PCldpatch,fNH3patch,PCldhdr,fNH3hdr,CH4abs,NH3abs,CH4abshdr,NH3abshdr = \
        make_L2_HSTGO_abs_data(pathHST,hdr619,data619patch,hdr631,data631patch,
                               hdr645,data645patch,LonSys,LatLims,LonLims,
                               dataversion='H',plot=False)
    #plot_HST_global_maps(CH4abs,NH3abs)
    #CH4patch=mp.make_patch(CH4abs,LatLims,LonLims,180,180,pad=True)
    PCldpatchflat=fp.flatten_patch(PCldpatch)
    #NH3patch=mp.make_patch(NH3abs,LatLims,LonLims,180,180,pad=True)
    fNH3patchflat=fp.flatten_patch(fNH3patch)

    ###########################################################################
    # Compute AOI and CI from Mike Wongs radiance data with a Minnaert 
    # correction applied.
    AOIpatch,CIpatch=make_L2_HSTGO_AOI_CI(pathHST,hdr275,data275patch,hdr889,data889patch,
                                          hdr395,data395patch,hdr631,data631patch,
                                          LonSys,plot=False)
    AOIpatchflat=fp.flatten_patch(AOIpatch)
    CIpatchflat=fp.flatten_patch(CIpatch)

    # Make standard RGB and false color 'methane' RGB
    RGB=make_L2_HSTGO_RGB_data(pathHST,data673patch,data502patch,data395patch,
                               emi673patch,emi502patch,emi395patch,LonSys,plot=False)
    #RGBpatch=mp.make_patch(RGB,LatLims,LonLims,180,180)
    RGBpatch=RGB

    RGBMeth=make_L2_HSTGO_RGB_data(pathHST,data673patch,data727patch,data889patch,
                                   emi673patch,emi727patch,emi889patch,LonSys,plot=False)
    #RGBMethpatch=mp.make_patch(RGBMeth,LatLims,LonLims,180,180)
    RGBMethpatch=RGBMeth

    ###########################################################################
    # Write FITS Patches
    ########################################################################### 
    import socket
    hostname = socket.gethostname()
    from config_VA import Host_path
    #pathout="C:/Astronomy/Projects/SAS 2021 Ammonia/Data/HST GO 18055/"+obskeyHST[:-1]+"/"+obskeyHST
    pathout=Host_path[hostname]+'/HST GO 18055/'+obskeyHST[:-1]+'/'+obskeyHST

    if not os.path.exists(pathout+"/L1"):
        os.makedirs(pathout+"/L1")
    if not os.path.exists(pathout+"/L3"):
        os.makedirs(pathout+"/L3")

    # Write L3 logitudinally flattened rnvironmental parameters and indices 
    write_HST_fits_patch(obskeyHST,LonSys,PCldpatchflat,PCldhdr,pathout+"/L3",'PCld',
                         LatLims=LatLims, LonLims=LonLims,dmin=0.0,dmax=4000.0)
    write_HST_fits_patch(obskeyHST,LonSys,fNH3patchflat,fNH3hdr,pathout+"/L3",'fNH3',
                         LatLims=LatLims, LonLims=LonLims,dmin=-100.0,dmax=600.0)
    write_HST_fits_patch(obskeyHST,LonSys,CIpatchflat,fNH3hdr,pathout+"/L3",'CI',
                         LatLims=LatLims, LonLims=LonLims,dmin=0.0,dmax=1.0)
    write_HST_fits_patch(obskeyHST,LonSys,AOIpatchflat,fNH3hdr,pathout+"/L3",'AOI',
                         LatLims=LatLims, LonLims=LonLims,dmin=0.0,dmax=1.0)


    #Write L1 Minnaert-corrected reflectances
    write_HST_fits_patch(obskeyHST,LonSys,data275patch,hdr275,pathout+"/L1",'275 Refl',
                         LatLims=LatLims, LonLims=LonLims,dmin=0.0,dmax=10.0)

    write_HST_fits_patch(obskeyHST,LonSys,data395patch,hdr395,pathout+"/L1",'395 Refl',
                         LatLims=LatLims, LonLims=LonLims,dmin=0.0,dmax=10.0)

    write_HST_fits_patch(obskeyHST,LonSys,data502patch,hdr502,pathout+"/L1",'502 Refl',
                         LatLims=LatLims, LonLims=LonLims,dmin=0.0,dmax=10.0)

    write_HST_fits_patch(obskeyHST,LonSys,data619patch,hdr619,pathout+"/L1",'619 Refl',
                         LatLims=LatLims, LonLims=LonLims,dmin=0.0,dmax=10.0)

    write_HST_fits_patch(obskeyHST,LonSys,data631patch,hdr631,pathout+"/L1",'631 Refl',
                         LatLims=LatLims, LonLims=LonLims,dmin=0.0,dmax=10.0)

    write_HST_fits_patch(obskeyHST,LonSys,data645patch,hdr645,pathout+"/L1",'645 Refl',
                         LatLims=LatLims, LonLims=LonLims,dmin=0.0,dmax=10.0)

    write_HST_fits_patch(obskeyHST,LonSys,data673patch,hdr673,pathout+"/L1",'673 Refl',
                         LatLims=LatLims, LonLims=LonLims,dmin=0.0,dmax=10.0)

    write_HST_fits_patch(obskeyHST,LonSys,data727patch,hdr727,pathout+"/L1",'727 Refl',
                         LatLims=LatLims, LonLims=LonLims,dmin=0.0,dmax=10.0)

    write_HST_fits_patch(obskeyHST,LonSys,data889patch,hdr889,pathout+"/L1",'889 Refl',
                         LatLims=LatLims, LonLims=LonLims,dmin=0.0,dmax=10.0)

    ###########################################################################
    # PLOT SECTION - this is like an L4 plot format
    ###########################################################################   
    figsz,aspect,plot_adjust=L4MP.figsize_and_aspect(LatLims,LonLims)
    # SETUP FIRST FIGURE
    titles=['fNH3 (ppm)','PCloud (mbar)']
    RGBtitle='RGB (673/502/395)'
    fig1,axs1=L4MP.set_up_figure(figsz,obskeyHST,LonSys,RGBaxs=2)
    # PLOT fNH3
    cbttlNH3="Mean="+str(np.mean(fNH3patchflat))[:4]+" $\pm$ "+str(np.std(fNH3patchflat))[:3]
    fNH3factor=1.5
    pp.plot_patch(fNH3patchflat*fNH3factor,LatLims,LonLims,180,180,'terrain_r',axs1[0],
                   cbarplot=True,cbar_title=cbttlNH3,cbar_reverse=False,
                   vn=0,vx=300)  
    axs1[1].set_title('fNH3 (ppm)',fontsize=10)
    # PLOT PCld    
    cbttlCH4="Mean="+str(np.mean(PCldpatchflat))[:4]+" $\pm$ "+str(np.std(PCldpatchflat))[:3]
    pp.plot_patch(PCldpatchflat,LatLims,LonLims,180,180,'Blues',axs1[1],
                   cbarplot=True,cbar_title=cbttlCH4,cbar_reverse=True,
                   vn=1000,vx=3000)  
    axs1[1].set_title('PCloud (mbar)',fontsize=10)
    # PLOT RGB
    print("################",RGB.shape)
    plot_HSTGO_RGB_patches(obskeyHST,LatLims,LonLims,RGBpatch,axs1[2],title=RGBtitle)
    fig1.subplots_adjust(**plot_adjust)    
    
    # SETUP SECOND FIGURE
    titles=['Color Index (CI)','Altitude Opacity Index (AOI)']
    RGBtitle="'Methane' RGB (673/727/889)"
    fig2,axs2=L4MP.set_up_figure(figsz,obskeyHST,LonSys,RGBaxs=2)
    # PLOT CI    
    cbttlCI="Mean="+str(np.mean(CIpatchflat))[:4]+" $\pm$ "+str(np.std(CIpatchflat))[:3]
    pp.plot_patch(CIpatchflat,LatLims,LonLims,180,180,'Spectral',axs2[0],
                   cbarplot=True,cbar_title=cbttlCI,cbar_reverse=False,
                   vn=0.35,vx=0.75)  
    axs2[0].set_title('Color Index (CI)',fontsize=10)    
    # PLOT AOI
    cbttlAOI="Mean="+str(np.mean(AOIpatchflat))[:4]+" $\pm$ "+str(np.std(AOIpatchflat))[:3]
    pp.plot_patch(AOIpatchflat,LatLims,LonLims,180,180,'Greys_r',axs2[1],
                   cbarplot=True,cbar_title=cbttlAOI,cbar_reverse=False,
                   vn=0.1,vx=0.4)  
    axs2[1].set_title('Altitude Opacity Index (AOI)',fontsize=10)
    # PLOT METHANE RGB
    #print("RGBMethpatch.shape,[RGBMethpatch[:,:,0],RGBMethpatch[:,:,1],RGBMethpatch[:,:,2]*4].shape",
    #      RGBMethpatch.shape,np.stack((RGBMethpatch[:,:,0],RGBMethpatch[:,:,1],RGBMethpatch[:,:,2]*4),axis=2).shape)
    plot_HSTGO_RGB_patches(obskeyHST,LatLims,LonLims,
                           np.stack((RGBMethpatch[:,:,0]*.8,RGBMethpatch[:,:,1]*1.5,RGBMethpatch[:,:,2]*7),axis=2),
                           axs2[2],title=RGBtitle)
    fig2.subplots_adjust(**plot_adjust)     

    if cont:
        tx_fNH3=[100,150,200,250]
        tx_PCld=[1600,2000,2400,2800]
        smoothed_fNH3 = gaussian_filter(fNH3patchflat, sigma=5)
        smoothed_PCld = gaussian_filter(PCldpatchflat, sigma=5)
        L4MP.ApplyContours(axs1,2,smoothed_fNH3,tx_fNH3,smoothed_PCld,tx_PCld,
                          LatLims,LonLims,IRTFcollection=False,IRTFaxs='',
                          CH4889collection=False,CH4889plot=False,CH4889axs='',
                          HST=True)
        L4MP.ApplyContours(axs2,2,smoothed_fNH3,tx_fNH3,smoothed_PCld,tx_PCld,
                          LatLims,LonLims,IRTFcollection=False,IRTFaxs='',
                          CH4889collection=False,CH4889plot=False,CH4889axs='',
                          HST=True)

    path="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Studies/"
    if waveplot:
        L4MP.RossbyWavePlot(obskeyHST,LonLims,fNH3patchflat,PCldpatchflat,
                       [figsz[0],figsz[1]],path,LonSys,HST=True)

    print(pathout)
    fig1.savefig(pathout+'/'+obskeyHST+' HST Sys'+ LonSys +' fNH3+PCld.png',dpi=150)
    fig2.savefig(pathout+'/'+obskeyHST+' HST Sys'+ LonSys +' CI+AOI.png',dpi=150)


