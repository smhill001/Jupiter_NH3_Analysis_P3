def blendstack(stackdata,stackweights):
    """
    Blend individual map patches into a larger map using using weights based
    on center-to-limb distance.

    Parameters
    ----------
    stackdata : TYPE
        DESCRIPTION.
    stackweights : TYPE
        DESCRIPTION.

    Returns
    -------
    blendweightdata : TYPE
        DESCRIPTION
    stdvdata : TYPE
        DESCRIPTION
    fracdata : TYPE
        DESCRIPTION

    """
    import numpy as np
    indzf=np.where(stackdata==0)
    stackdata[indzf]=np.nan
    blenddata=np.nanmean(stackdata,axis=2)
    
    stackdatamasked = np.ma.MaskedArray(stackdata, mask=np.isnan(stackdata))
    print("$$$$$$$$")
    print(stackdatamasked.shape,stackweights.shape)
    blendweightdata=np.ma.average(stackdatamasked, axis=2, weights=stackweights) 
    
    stdvdata=np.nanstd(stackdata,axis=2)
    fracdata=stdvdata/blenddata
    
    return(blendweightdata,stdvdata,fracdata)

def make_map_frame_lims(hdr,lonhalfwidth,boxcar,LonSys):
    """
    Create the masking longitudes for data and weights for blended for the
    input map frames for contiguous maps

    Parameters
    ----------
    hdr : TYPE
        DESCRIPTION.
    lonhalfwidth : TYPE
        DESCRIPTION.
    boxcar : int
        DESCRIPTION.
    LonSys : int
        DESCRIPTION.

    Returns
    -------
    ll_0 : int
        DESCRIPTION.
    ll_1 : int
        DESCRIPTION.
    wl_0 : int
        DESCRIPTION.
    wl_1 : int
        DESCRIPTION.

    """
        # Map LonSys to appropriate central meridian key
    
    cm_key = {
        '1': 'CM1',
        '2': 'CM2',
        '3': 'CM3'#,
        #'subobs':180
    }.get(LonSys)

    if cm_key:
        cm = hdr[cm_key]
        ll_0 = int(360 - (cm - lonhalfwidth))
        ll_1 = int(360 - (cm + lonhalfwidth))
        wl_0 = int(360 - (cm - lonhalfwidth + boxcar))
        wl_1 = int(360 - (cm + lonhalfwidth - boxcar))
    
    return ll_0,ll_1,wl_0,wl_1

def make_5u_cont_map(IRTFcollection,LonSys,boxcar,lats):
    """
    

    Parameters
    ----------
    IRTFdataset : TYPE
        DESCRIPTION.
    LonSys : TYPE
        DESCRIPTION.
    boxcar : TYPE
        DESCRIPTION.
    lats : TYPE
        DESCRIPTION.

    Returns
    -------
    blendweightIRTF : TYPE
        DESCRIPTION.
    stdvIRTF : TYPE
        DESCRIPTION.
    fracIRTF : TYPE
        DESCRIPTION.

    """
    import sys
    drive='c:'
    sys.path.append(drive+'/Astronomy/Python Play')
    sys.path.append('C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/5micron')

    import numpy as np
    from astropy.io import fits
    import scipy.ndimage as ndi
    from astropy.time import Time
    import get_5um_collection as g5u
    

    ###########################################################################
    # Loop over 5um FITS files to create maps
    ###########################################################################
    pathIRTF_FITS="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/5micron/FITS/output/"
    stackIRTF,stackweightsIRTF,stackTimeIRTF=\
        [np.zeros((180, 360)) for _ in range(3)]

    IRTFdataset=g5u.get_5um_collection(IRTFcollection)

    print("IRTFdataset=",IRTFdataset)
    FirstIRTF=True
    for IRTFfile in IRTFdataset:
        #Read IRTF data file
        print("*******obsdate=",IRTFfile)
        
        IRTFhdulist=fits.open(pathIRTF_FITS+IRTFfile)
        IRTFhdulist.info()
        IRTFhdr=IRTFhdulist[0].header
        IRTFdata=IRTFhdulist[0].data
        IRTFhdulist.close()
        
        #Roll to LonSys other than 3 if necessary
        if LonSys=='1':
            IRTFroll=IRTFhdr["CM3"]-IRTFhdr["CM1"]
            IRTFdatar=np.roll(IRTFdata,int(IRTFroll),axis=1)
        if LonSys=='2':
            IRTFroll=IRTFhdr["CM3"]-IRTFhdr["CM2"]
            IRTFdatar=np.roll(IRTFdata,int(IRTFroll),axis=1)
        else:
            IRTFdatar=IRTFdata
                       
        ll_0,ll_1,wl_0,wl_1=make_map_frame_lims(IRTFhdr,70,boxcar,LonSys)
        
        #######################################################################
        # RESET EMPTY FRAMES FOR NEWEST LIMB-CORRECTED MAP REGION
        #######################################################################
        outputweightsIRTF,outputmask,datamask=\
            [np.zeros((180, 360)) for _ in range(3)]

        epoch=Time(IRTFhdr["Date-Obs"], format='isot', scale='utc')

        if wl_1>0 and wl_0<360:
            outputmask[lats[0]:lats[1],wl_1:wl_0]= 1.0
            print("$$$$$$$$$$-------- wl_1,wl_0,=",wl_1,wl_0)
        elif wl_0>360:
            outputmask[lats[0]:lats[1],wl_1:360]= 1.0
            outputmask[lats[0]:lats[1],0:np.mod(wl_0,360)]= 1.0
            print("$$$$$$$$$$-------- wl_1:360,0:np.mod(wl_0,360)=",wl_1,360,0,np.mod(wl_0,360))
        elif wl_1<0:
            outputmask[lats[0]:lats[1],np.mod(wl_1,360):360]= 1.0
            outputmask[lats[0]:lats[1],0:wl_0]= 1.0
            print("$$$$$$$$$$-------- wl_1:360,0:np.mod(wl_0,360)=",wl_1,360,0,np.mod(wl_0,360))

        outputweightsIRTF=ndi.uniform_filter1d(outputmask,boxcar,1)
        datamask[outputweightsIRTF>0.0]=1.0

        if FirstIRTF:
            stackIRTF=IRTFdatar*datamask
            stackIRTF=np.reshape(stackIRTF,(180,360,1))
            stackweightsIRTF=outputweightsIRTF
            stackweightsIRTF=np.reshape(stackweightsIRTF,(180,360,1))
            stackTimeIRTF=epoch.jd*datamask
            stackTimeIRTF=np.reshape(epoch.jd*datamask,(180,360,1))
        else:
            stackIRTF=np.dstack((stackIRTF,IRTFdatar*datamask))
            tempweightsIRTF=np.reshape(outputweightsIRTF,(180,360,1))
            stackweightsIRTF=np.dstack((stackweightsIRTF,tempweightsIRTF))
            tempTimeIRTF=np.reshape(epoch.jd*datamask,(180,360,1))
            stackTimeIRTF=np.dstack((stackTimeIRTF,tempTimeIRTF))
        FirstIRTF=False
            
    #######################################################################
    # Flatten cubes by taking the weighted mean and ignoring zero values
    #######################################################################
    indzf=np.where(stackIRTF==0)
    stackIRTF[indzf]=np.nan
    blendweightIRTF,stdvIRTF,fracIRTF=blendstack(stackIRTF,stackweightsIRTF)
    stackTimeIRTF=np.array(stackTimeIRTF)
    blendweightTimeIRTF,stdvTimeIRTF,fracTimeIRTF=blendstack(stackTimeIRTF,stackweightsIRTF)

    return blendweightIRTF,stdvIRTF,fracIRTF,blendweightTimeIRTF,stdvTimeIRTF,fracTimeIRTF


def createL4FileName(mean_time, env_data_type,LonSys,lats,LonLims,collection):
    from astropy.time import Time
    import numpy as np
    import make_lat_lon_str as MLLS
    
    mean_epoch_jd=Time(np.mean(mean_time),format='jd',scale='utc')
    dt_object = mean_epoch_jd.to_datetime()
    dt_string = dt_object.strftime('%Y%m%dT%H%M%S')
    dt_string_fits = dt_object.strftime('%Y-%m-%dT%H:%M:%S')

    latstr,lonstr=MLLS.make_lat_lon_str(lats,LonLims)
    print(dt_string+'_Jupiter_')
    fnout=dt_string+'_Jupiter_'+collection+'_ContMap_L4'+env_data_type+' Sys'+LonSys+' '+lonstr+' '+latstr+'.fits'

    return fnout,mean_epoch_jd,dt_string,dt_string_fits

def write_fits_contiguous_map(data,variance,fract,env_data_type,mean_time,
                              LonSys,lats,LonLims,collection,
                              first_fits_hdr,proj='../../Data/L4 FITS (cont maps)/',
                              target="Jupiter"):
    
    import numpy as np
    from astropy.io import fits
    from time import gmtime, strftime
    import os
    import get_spice_ephem as sp_ephem

    pathout="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Studies/"
    hdudata = fits.PrimaryHDU(data)
    variancedata=fits.ImageHDU(variance)
    fractdata=fits.ImageHDU(fract)
    timedata=fits.ImageHDU(mean_time)
    hdul = fits.HDUList([hdudata,variancedata,fractdata,timedata])

    fn,epochjd,dt_string,dt_string_fits=createL4FileName(mean_time, env_data_type,
                                                         LonSys,lats,LonLims,collection)
    
    ephem=sp_ephem.get_spice_ephem(dt_string_fits,planet=target)

    hdul[0].header['BUNIT']=env_data_type    
    hdul[0].header['BITPIX']=-32
    hdul[0].header['DATE-OBS']=dt_string_fits
    hdul[0].header['AUTHOR']='Hill, S. M.'
    hdul[0].header['FILENAME']=fn

    hdul[0].header['OBJECT']='Jupiter'
    hdul[0].header['TIME']=strftime("%Y-%m-%d %H:%M:%S", gmtime())
    
    hdul[0].header['TELESCOP']=first_fits_hdr['TELESCOP']
    hdul[0].header['INSTRUME']=first_fits_hdr['INSTRUME']
    hdul[0].header['SEEING']=first_fits_hdr['SEEING']
    hdul[0].header['TRANSPAR']=first_fits_hdr['TRANSPAR']
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

def make_L4_cont_maps(collection="20251016-20251017",obskeys=False,LonSys='3',
                      FiveMicron='fits',Five_obskey='',IRTFcollection='20240205-20240205',
                      lats=[0,180],LonLims=[0,360],variance=True,proj='../../Data/L4 FITS (cont maps)/',
                      bare_maps=False,cb=False,LimbCorrection=True,
                      lonhalfwidth=45,boxcar=9):
    """
    Rewriten Nov 2025 to only produce global, cylindrical FITS files in
        system 3 longitude and pg latitude. All other functions of the original
        routine - plotting/subplotting, extrema, etc, have been removed and
        decoupled.
    """
    import sys
    drive='c:'
    sys.path.append(drive+'/Astronomy/Python Play')
    sys.path.append(drive+'/Astronomy/Python Play/Util_P3')
    sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy')
    sys.path.append('./Services')

    import numpy as np
    import scipy.ndimage as ndi
    sys.path.append('./Maps')
    from astropy.time import Time

    import read_fits_map_L2_L3 as RFM
    import get_map_collection as gmc
    
    if not(obskeys):
        obskeys,dummy=gmc.get_map_collection(collection)
        print("################################",obskeys)

    ###########################################################################
    # Establish empty arrays for stacked contiguous maps
    ###########################################################################
    stackfNH3,stackPCloud,stackweights,stackR,stackG,stackB,stackTime=\
        [np.zeros((180, 360)) for _ in range(7)]

    ###########################################################################
    # Loop over observations in each data set and create 3D cubes
    ###########################################################################
    print("obskeys=",obskeys)
    First=True
    for obskey in obskeys:
        print("*******obsdate in MakeContiguousMap=",obskey)
        PCloudhdr,PClouddata,fNH3hdr,fNH3data,sza,eza,RGB,RGB_CM2,RGBtime= \
                        RFM.read_fits_map_L2_L3(obskey=obskey,LonSys=LonSys,
                                                imagetype="Map",Level="L3",
                                                LimbCorrection=LimbCorrection)
                        
        #Set longitude limits for data (ll_x) and weights (wl_x)
        ll_0,ll_1,wl_0,wl_1=make_map_frame_lims(fNH3hdr,lonhalfwidth,boxcar,LonSys)
        
        #######################################################################
        # RESET EMPTY FRAMES FOR NEWEST LIMB-CORRECTED MAP REGION
        #######################################################################
        outputmask,outputweights,datamask=\
            [np.zeros((180, 360)) for _ in range(3)]
        epoch=Time(fNH3hdr["Date-Obs"], format='isot', scale='utc')

        #Calculate mask and weights arrays
        # Link for uniform_filter1d
        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.uniform_filter1d.html#scipy.ndimage.uniform_filter1d       
        # 2025-11-11: UPDATED WEIGHTING SCHEME
        if wl_1>0 and wl_0<360:
            outputmask[lats[0]:lats[1],wl_1:wl_0]= 1.0
            print("$$$$$$$$$$-------- wl_1,wl_0,=",wl_1,wl_0)
        elif wl_0>360:
            outputmask[lats[0]:lats[1],wl_1:360]= 1.0
            outputmask[lats[0]:lats[1],0:np.mod(wl_0,360)]= 1.0
            print("$$$$$$$$$$-------- wl_1:360,0:np.mod(wl_0,360)=",wl_1,360,0,np.mod(wl_0,360))
        elif wl_1<0:
            outputmask[lats[0]:lats[1],np.mod(wl_1,360):360]= 1.0
            outputmask[lats[0]:lats[1],0:wl_0]= 1.0
            print("$$$$$$$$$$-------- wl_1:360,0:np.mod(wl_0,360)=",wl_1,360,0,np.mod(wl_0,360))

        outputweights=ndi.uniform_filter1d(outputmask,boxcar,1)
        datamask[outputweights>0.0]=1.0
        
        if First:
            stackfNH3=fNH3data*datamask
            stackfNH3=np.reshape(stackfNH3,(180,360,1))
            stackPCloud=PClouddata*datamask
            stackPCloud=np.reshape(stackPCloud,(180,360,1))
            stackweights=outputweights
            stackweights=np.reshape(stackweights,(180,360,1))
            stackR=RGB[:,:,0]*datamask
            stackR=np.reshape(stackR,(180,360,1))
            stackG=RGB[:,:,1]*datamask
            stackG=np.reshape(stackG,(180,360,1))
            stackB=RGB[:,:,2]*datamask
            stackB=np.reshape(stackB,(180,360,1))
            stackTime=epoch.jd*datamask
            stackTime=np.reshape(stackTime,(180,360,1))
        else:
            stackfNH3=np.dstack((stackfNH3,fNH3data*datamask))
            stackPCloud=np.dstack((stackPCloud,PClouddata*datamask))
            tempweights=np.reshape(outputweights,(180,360,1))
            stackweights=np.dstack((stackweights,tempweights))
            tempR=np.reshape(RGB[:,:,0]*datamask,(180,360,1))
            stackR=np.dstack((stackR,tempR))
            tempG=np.reshape(RGB[:,:,1]*datamask,(180,360,1))
            stackG=np.dstack((stackG,tempG))
            tempB=np.reshape(RGB[:,:,2]*datamask,(180,360,1))
            stackB=np.dstack((stackB,tempB))
            tempTime=np.reshape(epoch.jd*datamask,(180,360,1))
            stackTime=np.dstack((stackTime,tempTime))
        First=False
        
    ###########################################################################
    # Flatten cubes by taking the mean and ignoring zero values
    #!!!  WRITE FITS FILES OF MEAN AND STANDARD DEVIATIONS OF FULL MAPS!!!!
    #
    # Link for taking weighted averages of masked arrays
    # https://numpy.org/doc/2.2/reference/generated/numpy.ma.average.html
    # https://stackoverflow.com/questions/21113384/python-numpy-weighted-average-with-nans
    ###########################################################################
    #Blend science data
    print(stackfNH3.shape)
    blendweightfNH3,stdvfNH3,fracfNH3=blendstack(stackfNH3,stackweights)
    blendweightPCloud,stdvPCloud,fracPCloud=blendstack(stackPCloud,stackweights)
  
    ###########################################################################
    #Blend RGB
    blendweightR,stdvR,fracR=blendstack(stackR,stackweights)
    blendweightG,stdvG,fracG=blendstack(stackG,stackweights)
    blendweightB,stdvB,fracB=blendstack(stackB,stackweights)

    blendRGBweight=np.zeros([180,360,3])
    blendRGBweight[:,:,0]=blendweightR
    blendRGBweight[:,:,1]=blendweightG
    blendRGBweight[:,:,2]=blendweightB
    
    ###########################################################################
    #Blend DateTime
    stackTime=np.array(stackTime)
    blendweightTime,stdvTime,fracTime=blendstack(stackTime,stackweights)

    ###########################################################################

    write_fits_contiguous_map(blendweightfNH3,stdvfNH3,fracfNH3,'fNH3',blendweightTime,
                                  LonSys,lats,LonLims,collection,
                                  fNH3hdr,proj=proj)
    
    write_fits_contiguous_map(blendweightPCloud,stdvPCloud,fracPCloud,'PCld',blendweightTime,
                                  LonSys,lats,LonLims,collection,
                                  fNH3hdr,proj=proj)
    
    write_fits_contiguous_map(blendRGBweight,blendRGBweight,blendRGBweight,'RGB',blendweightTime,
                                  LonSys,lats,LonLims,collection,
                                  fNH3hdr,proj=proj)
    
    if FiveMicron=='fits':
        blendweightIRTF,stdvIRTF,fracIRTF,blendweightTimeIRTF,stdvTimeIRTF,fracTimeIRTF=\
            make_5u_cont_map(IRTFcollection,LonSys,boxcar,lats)
        write_fits_contiguous_map(blendweightIRTF,stdvIRTF,fracIRTF,'IRTF',blendweightTimeIRTF,
                                      LonSys,lats,LonLims,IRTFcollection,
                                      fNH3hdr,proj=proj)

    return
