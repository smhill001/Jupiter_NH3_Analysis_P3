def get_JALPO_collectionlists(collectionlists):
    
    from os import listdir
    
    pathRGB="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/5micron/JALPO JPEGS/"
      
    if collectionlists=="5um":
        RGBfiles=["j22mapsL3_08.jpg","j22mapsL3_09a.jpg","j22mapsL3_09b.jpg","j22mapsL3_11.jpg","j22mapsL3_13.jpg"]
        
        yranges={"j22mapsL3_08.jpg":[1350,1580],
                 "j22mapsL3_09a.jpg":[2360,2590],
                 "j22mapsL3_09b.jpg":[1603,1833],
                 "j22mapsL3_11.jpg":[1099,1329],
                 "j22mapsL3_13.jpg":[ 340, 570]}
        
        collections={"j22mapsL3_08.jpg":"01",
                 "j22mapsL3_09a.jpg":"02",
                 "j22mapsL3_09b.jpg":"03",
                 "j22mapsL3_11.jpg":"04",
                 "j22mapsL3_13.jpg":"05"}
        
        outnames={"j22mapsL3_08.jpg":"2022-08-02-0000_0-Jupiter_5micron_RGB_CM2_L360_MAP-BARE",
                 "j22mapsL3_09a.jpg":"2022-08-16-1200_0-Jupiter_5micron_RGB_CM2_L360_MAP-BARE",
                 "j22mapsL3_09b.jpg":"2022-08-18-0000_0-Jupiter_5micron_RGB_CM2_L360_MAP-BARE",
                 "j22mapsL3_11.jpg":"2022-09-03-1200_0-Jupiter_5micron_RGB_CM2_L360_MAP-BARE",
                 "j22mapsL3_13.jpg":"2022-09-29-1200_0-Jupiter_5micron_RGB_CM2_L360_MAP-BARE"}
        
    elif collectionlists=="889CH4":
        
        RGBfiles=["j22mapsL3_08.jpg","j22mapsL3_09a.jpg","j22mapsL3_13.jpg",
                       "j23-24mapsL3-13a.jpg","j23-24mapsL3-13b.jpg","j23-24mapsL3-20.jpg",
                       "j25-26mapsL3_04.jpg","j25-26mapsL3-06.jpg"]
        
        yranges={"j22mapsL3_08.jpg":[1856,2086],
                      "j22mapsL3_09a.jpg":[1856,2086],
                      "j22mapsL3_13.jpg":[ 88, 318],
                      "j23-24mapsL3-13a.jpg":[ 851, 1081],
                      "j23-24mapsL3-13b.jpg":[ 1356, 1586],             
                      "j23-24mapsL3-20.jpg":[ 94, 324],             
                      "j25-26mapsL3_04.jpg":[598,828],
                      "j25-26mapsL3-06.jpg":[852,1082]}
        
        collections={"j22mapsL3_08.jpg":"20220802-20220803",
                 "j22mapsL3_09a.jpg":"20220816-20220816",
                 "j22mapsL3_13.jpg":"20220929-20220929",
                 "j23-24mapsL3-13a.jpg":"20231017-20231017",
                 "j23-24mapsL3-13b.jpg":"20231018-20231018",
                 "j23-24mapsL3-20.jpg":"20240131-20240201",             
                 "j25-26mapsL3_04.jpg":"20250919-20250920",
                 "j25-26mapsL3-06.jpg":"20251017-20251018"}
        
        outnames={"j22mapsL3_08.jpg":"2022-08-02-0000_0-Jupiter_889CH4_RGB_CM2_L360_MAP-BARE",
                 "j22mapsL3_09a.jpg":"2022-08-16-1200_0-Jupiter_889CH4_RGB_CM2_L360_MAP-BARE",
                 "j22mapsL3_13.jpg":"2022-09-29-1200_0-Jupiter_889CH4_RGB_CM2_L360_MAP-BARE",
                 "j23-24mapsL3-13a.jpg":"2023-10-17-1200_0-Jupiter_889CH4_RGB_CM2_L360_MAP-BARE",
                 "j23-24mapsL3-13b.jpg":"2023-10-18-1200_0-Jupiter_889CH4_RGB_CM2_L360_MAP-BARE",
                 "j23-24mapsL3-20.jpg":"2024-02-01-0000_0-Jupiter_889CH4_RGB_CM2_L360_MAP-BARE",
                 "j25-26mapsL3_04.jpg":"2025-09-20-0000_0-Jupiter_889CH4_RGB_CM2_L360_MAP-BARE",
                 "j25-26mapsL3-06.jpg":"2025-10-17-0000_0-Jupiter_889CH4_RGB_CM2_L360_MAP-BARE"}
        
    return RGBfiles,yranges,collections,outnames,pathRGB

def write_fits_JALPO_map(data,env_data_type,mean_time,
                              LonSys,lats,LonLims,collection,
                              proj='../../Data/L4 FITS (cont maps)/',
                              target="Jupiter"):
    
    import numpy as np
    from astropy.io import fits
    from time import gmtime, strftime
    import os
    import get_spice_ephem as sp_ephem
    import make_L4_cont_maps as M4M

    pathout="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Studies/"
    hdudata = fits.PrimaryHDU(data)
    hdul = fits.HDUList([hdudata])   

    fn,epochjd,dt_string,dt_string_fits=M4M.createL4FileName(mean_time, env_data_type,
                                                         LonSys,lats,LonLims,collection)
    
    ephem=sp_ephem.get_spice_ephem(dt_string_fits,planet=target)

    hdul[0].header['BUNIT']=env_data_type    
    hdul[0].header['BITPIX']=-32
    hdul[0].header['DATE-OBS']=dt_string_fits
    hdul[0].header['AUTHOR']='Hill, S. M.'
    hdul[0].header['FILENAME']=fn

    hdul[0].header['OBJECT']='Jupiter'
    hdul[0].header['TIME']=strftime("%Y-%m-%d %H:%M:%S", gmtime())
    
    hdul[0].header['TELESCOP']='JALPO'
    hdul[0].header['INSTRUME']='JALPO'
    hdul[0].header['SEEING']='N/A'
    hdul[0].header['TRANSPAR']='N/A'
    #hdul[0].header['VERSION']=('1.2','Lat. Dep. Gravity; eta=2 (not 4)')
    hdul[0].header['CTYPE1']=('Sys. '+LonSys+' Longitude','deg')
    hdul[0].header['CTYPE2']=('PG Latitude','deg')
    hdul[0].header['CM1']=float(ephem[0].strip())
    hdul[0].header['CM2']=float(ephem[1].strip())
    hdul[0].header['CM3']=float(ephem[2].strip())
        
    fnout=pathout+proj+"/"+fn
    
    try:
        os.remove(fnout)
    except: 
        print("file doesn't exist")
    hdul.writeto(fnout)
    hdul.close()
        
    return
from matplotlib.pyplot import imread
import pylab as pl
import numpy as np
from os import listdir
import copy
import scipy.ndimage as ndimage
from imageio import imwrite
import get_WINJupos_ephem as WJ_ephem


def make_JALPO_maps(collectionlists,env_data_type):
    """
    Takes L3 cylindrical maps from Shinji Mizumoto at JALPO and extracts 
    60S-60N latitudes and recenters them with zero longitude at the right and
    360 longitude at the left.
    
    https://alpo-j.sakura.ne.jp/Latest/j_Cylindrical_Maps/j_Cylindrical_Maps.htm
    
    

    Parameters
    ----------
    RGBfiles : TYPE
        DESCRIPTION.
    yranges : TYPE
        DESCRIPTION.
    outnames : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    from astropy.time import Time
    
    RGBfiles,yranges,collections,outnames,pathRGB=get_JALPO_collectionlists(collectionlists)

    figysize=len(RGBfiles)
    fig1,axs1=pl.subplots(figysize,1,figsize=(6.0,6.0), dpi=150, facecolor="white",
                          sharex=True,sharey=True)
    fig2,axs2=pl.subplots(figysize,1,figsize=(4.0,6.0), dpi=150, facecolor="white",
                          sharex=True,sharey=True)
    #lon_back_plane=np.zeros(230,939,dtype=float)
    colat=np.arange(30,150,1)
    rows=np.array((1.0-(np.cos(colat*np.pi/180.)/np.cos(30*np.pi/180.)))*115,dtype=int)
    mymap=np.zeros((180,360,3),dtype=float)
    #print(colat)
    #print()
    #print(rows)
    subplot=0
    for RGBfile in RGBfiles:
        print("## ",RGBfile)
        if '.jpg' in RGBfile:
            print()
            print(RGBfile)
            print()
            RGB=imread(pathRGB+RGBfile)
            print(RGB.shape)
            if "j22" in RGBfile:
                roll=120
            elif "j23" in RGBfile:
                roll=240
            elif "j25" in RGBfile:
                roll=160
            
            outname=outnames[RGBfile]
            RGBtime=(outname[0:10]+"_"+outname[11:13]+":"+outname[13:15]+":00")
            eph=WJ_ephem.get_WINJupos_ephem(RGBtime,planet="Jupiter")
            #time.sleep(5)
            RGB_CM1=float(eph[0].strip())
            RGB_CM2=float(eph[1].strip())
            RGB_CM3=float(eph[2].strip())
    
            y0=yranges[RGBfile][0]
            y1=yranges[RGBfile][1]
            RGBpatch=copy.deepcopy(RGB[y0:y1,294:1231,:]) #-> row, col, depth
            axs1[subplot].imshow(RGBpatch,origin='upper')
    
            NH3LonLims=[360.,0.]
            axs1[subplot].xlim=[360-NH3LonLims[0],360-NH3LonLims[1]]
            axs1[subplot].set_xticks(np.linspace(450,0,31), minor=False)
            axs1[4].xticklabels=np.array(np.mod(np.linspace(450,0,31),360))
    
    
            longscaledmap=ndimage.zoom(RGBpatch, (1, 0.3842, 1))
            print(longscaledmap.shape)
            colatindx=30
            for row in rows:
                mymap[colatindx,:,:]=longscaledmap[row,:,:]
                #print(colatindx,row,np.max(longscaledmap[row,:,:]),np.max(mymap[colatindx,:,:]))
    
                colatindx=colatindx+1
            
            print("################ roll=",roll)
            mymap=np.roll(mymap,-roll,axis=1)
            #RGBroll=RGB_CM2-RGB_CM3
            #print("########## RGBroll=",RGBroll)
            #mymap=np.roll(mymap,int(-RGBroll),axis=1)   
            ###################################################################
            # Write output PNG file
            imwrite(pathRGB+outnames[RGBfile]+'.png', mymap)#.astype(np.uint16))
            
            ###################################################################
            # Write FITS files
            print(RGBtime,RGBtime.replace("_"," "))
            temptimestr=RGBtime.replace("_"," ")
            mean_time=Time(temptimestr,format='iso',scale='utc')
            write_fits_JALPO_map(np.mean(mymap,axis=2),env_data_type,mean_time,
                                          '3',[0,180],[0,360],collections[RGBfile],
                                          target="Jupiter")
    
            ###################################################################
            # Plot for diagnostic purposes
            axs2[subplot].imshow(np.array(mymap,dtype=int),origin='upper')
            #axs2[subplot].xlim=[360-NH3LonLims[0],360-NH3LonLims[1]]
            axs2[subplot].set_xticks(np.linspace(360,0,25), minor=False)
            axs2[4].xticklabels=np.array(np.mod(np.linspace(360,0,25),360))
            axs2[subplot].tick_params(axis='both',which='major',labelsize=7)
            
            subplot=subplot+1
    
    
    pathmapplots=pathRGB
    
    fig1.subplots_adjust(left=0.05, bottom=0.05, right=0.97, top=0.93,
                wspace=0.10, hspace=0.03)     
    fig1.savefig(pathmapplots+" "+env_data_type+" maps.png",dpi=300)
    
    fig2.subplots_adjust(left=0.05, bottom=0.05, right=0.97, top=0.93,
                wspace=0.10, hspace=0.03)     
    fig2.savefig(pathmapplots+" "+env_data_type+" maps remapped.png",dpi=300)
