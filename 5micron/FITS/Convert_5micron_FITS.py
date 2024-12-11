def Convert_5micron_FITS():
    """
    Created on Mon Nov 11 08:06:50 2024
    
    @author: smhil
    """
    import sys
    drive='c:'
    sys.path.append(drive+'/Astronomy/Python Play')
    sys.path.append(drive+'/Astronomy/Python Play/Util_P3')
    sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy')
    sys.path.append('./Services')

    from astropy.io import fits
    sys.path.append('./Services')
    import pylab as pl
    import numpy as np
    import os
    import scipy
    import get_WINJupos_ephem as WJ_ephem
    import planetmapper as pm
    from scipy import interpolate
    import copy
    
    pathFITS="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/5micron/FITS/"
    
    for fname in os.listdir(pathFITS):
        #name, cmap, ext = os.path.splitext(fname)
        #print("name, cmap, ext",name, cmap, ext)
        print("fname",fname)
        print("fname[:-9]",fname[-9:])
        tail=fname[-9:]
        print("pathFITS+fname",pathFITS+fname)
        if tail=="cmap.fits":
            print("fname=",fname)
            Sourcehdulist=fits.open(pathFITS+fname)
            Sourcehdulist.info()
            Sourcehdr=Sourcehdulist[0].header
            Sourcedata=Sourcehdulist[0].data
            Sourcehdulist.close()
            print("DATE_OBS",Sourcehdr["DATE_OBS"])
            print("TIME_OBS",Sourcehdr["TIME_OBS"])
            
            secs=float(Sourcehdr["TIME_OBS"][6:])
            minfrac=str(round(secs/6.))
            secsint=str(round(secs))
            #print("secs",secs)
            #print("minfrac",minfrac)
        
            fnskeleton="-Jupiter_5micron.fits"
            fnout=Sourcehdr["DATE_OBS"]+"-"+Sourcehdr["TIME_OBS"][0:2]+\
                  Sourcehdr["TIME_OBS"][3:5]+"_"+minfrac+fnskeleton
            #print("fnout",fnout)
        
            dateobs=Sourcehdr["DATE_OBS"]+"T"+Sourcehdr["TIME_OBS"][0:6]+secsint.zfill(2)+"Z"
            print("dateobs",dateobs)
        
            sec=str(int(str(fnout[16:17]))*6) #COMPUTE FROM FRACTIONAL WINJUPOS MINUTE
            time=(fnout[0:10]+"_"+fnout[11:13]+":"+fnout[13:15]+":"
                     +sec.zfill(2))
            eph=WJ_ephem.get_WINJupos_ephem(time,planet="Jupiter")
            CM1=float(eph[0].strip())
            CM2=float(eph[1].strip())
            CM3=float(eph[2].strip())
        
            LCMIII=Sourcehdr["LCMIII"]
            
            
            if Sourcehdr["NAXIS1"]==720:
                #data=np.resize(Sourcedata,[180,360])
                #data = Sourcedata.resize((180,360), PIL.Image.ANTIALIAS)
                data=scipy.ndimage.zoom(Sourcedata, 0.5, order=0)
            else:
                data=Sourcedata
        
            flipdata=np.flipud(data)
            Finaldata=copy.deepcopy(flipdata)
            latpc=np.arange(89.5,-89.6,-1)
            j=pm.Body("Jupiter")
            latpg=j.centric2graphic_lonlat(latpc,latpc)[1]
            #print("latpg",latpg)
            #print Wavelength.size,Signal.size
            for long in range(0,360,1):
                Interp=interpolate.interp1d(latpg,flipdata[:,long],kind='linear', 
                                            copy=True,bounds_error=False, 
                                            fill_value=np.NaN,axis=0)  
                SignalonGrid=Interp(latpc)
                Finaldata[:,long]=SignalonGrid
            
            rollangle=LCMIII-CM3
            if np.abs(rollangle)>2.0:
                print("################# rollangle ####",dateobs,rollangle)
            Finaldata=np.roll(Finaldata,int(rollangle),axis=1)

            hdu = fits.PrimaryHDU(Finaldata)
            hdu.header=Sourcehdr
            hdul = fits.HDUList([hdu])
            hdul[0].header['FILENAME']=fnout
            hdul[0].header['DATE-OBS']=dateobs
            hdul[0].header['AUTHOR']='Hill, S. M.'
            hdul[0].header['CM1']=(CM1,'Sys. 1 Long. Central Meridian')
            hdul[0].header['CM2']=(CM2,'Sys. 2 Long. Central Meridian')
            hdul[0].header['CM3']=(CM3,'Sys. 3 Long. Central Meridian')
            hdul[0].header['CTYPE1']=('Sys. 3 Longitude','deg')
            hdul[0].header['CTYPE2']=('PG Latitude','deg')
        
            out=pathFITS+"output/"+fnout
            
            try:
                os.remove(out)
            except: 
                print("file doesn't exist")
            hdul.writeto(out)
            hdul.close()
        
        
            pl.imshow(Finaldata)
