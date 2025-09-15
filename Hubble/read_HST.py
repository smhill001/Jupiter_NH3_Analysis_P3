import matplotlib.pyplot as pl
import numpy as np
import find_blob as fb
import MakeContiguousMap as MCM
import make_patch as mp
import make_patch_RGB as mpr
import plot_patch as pp
from imageio import imread


def getFilename(files, keyword):
    """
    filters observations by substring in filename

    Parameters:
    data(Object): keys with values being list of filenames
    keyword (string): substring to filter filenames by

    Returns:
    Object: keys with filename arrays filtered by keyword
    """
    filteredData=[]
    for file in files:
        if keyword in file:
            print("############",file)
            filteredData=file
    return filteredData

def read_HST(obskey="2024c_f631",LonSys='1'):
    """
    Reads an HST OPAL map (FITS or TIF) from 2024 and returns a map in
    system I, II, or III coordinates.

    Parameters
    ----------
    obskey : TYPE, optional
        DESCRIPTION. The default is "2024c_f631".
    LonSys : TYPE, optional
        DESCRIPTION. The default is '1'.

    Returns
    -------
    datar : TYPE
        DESCRIPTION.

    """
    import sys
    drive='c:'
    sys.path.append(drive+'/Astronomy/Python Play')
    sys.path.append(drive+'/Astronomy/Python Play/Util_P3')
    sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy')
    sys.path.append('./Services')

    #from matplotlib.pyplot import imread
    from imageio import imread
    from PIL import Image

    from astropy.io import fits
    sys.path.append('./Services')
    import os
    import numpy as np
    import convert_system3_to_I_II_spice as clong
    
    pathHST='C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Hubble/Data/'
    filenames=os.listdir(pathHST)
    
    filename=getFilename(filenames,obskey)
    print("filename=====",filename)
    if "fits" in filename:
        
        hdulist=fits.open(pathHST+filename)
        hdulist.info()
        hdr=hdulist[0].header
        data=hdulist[0].data
        hdulist.close()
    
        dateobs=hdr["date-obs"]
        longitudes=clong.convert_system3_to_I_II_spice(dateobs, 180.)           

    if "tif" in filename:
        data=imread(pathHST+filename)
        temp=Image.open(pathHST+filename).convert('RGB')
        data=np.array(temp)

        if "2024c" in filename:
            hdulist=fits.open(pathHST+"hlsp_opal_hst_wfc3-uvis_jupiter-2024c_f275w_v1_globalmap.fits")
        if "2024d" in filename:
            hdulist=fits.open(pathHST+"hlsp_opal_hst_wfc3-uvis_jupiter-2024d_f275w_v1_globalmap.fits")

        hdulist.info()
        hdr=hdulist[0].header
        hdulist.close()
    
        dateobs=hdr["date-obs"]
        longitudes=clong.convert_system3_to_I_II_spice(dateobs, 180.)
     
    if LonSys=='1':
        roll=180.-longitudes["System I"]
        print("LonSys=",LonSys,longitudes["System I"])
    elif LonSys=='2':
        roll=180.-longitudes["System II"]
    else:
        roll=0

    if "fits" in filename:
        datar=np.roll(data,int(roll)*10,axis=1)
    if "tif" or "png" in filename:
        datar=np.roll(data,int(roll)*10,axis=1)

    print(filename)
    #pl.imshow(data)
    return datar,dateobs

###############################################################################
def HST_feature_overlay(HSTkey="2024d_f395n-f502n",LonSys='1'):
    """
    Runs MakeContinguousMap on 2024-11-18 data with the feature overlay.
    Then reads an HST OPAL TIF image (2024d) and plots the feature overlay.
    This is all in System I longitude and covers NEZ features.

    Parameters
    ----------
    LonSys : TYPE, optional
        DESCRIPTION. The default is '1'.

    Returns
    -------
    None.

    """
    LatLims=[75,95]
    CM=270
    LonRng=90
    LatLims=np.array(LatLims)
    LonLims=np.array([360-int(CM+LonRng),360-int(CM-LonRng)])
    
    print("########",LatLims,LonLims,CM,LonRng)
    latstr,lonstr=MCM.make_lat_and_lon_str(LatLims,LonLims)

    lats,blendweightPCloud,blendweightfNH3,blendRGBweight,labeled_fNH3, props_fNH3,labeled_Plum, props_Plum,labeled_NEDF, props_NEDF=\
        MCM.MakeContiguousMap(collection="20241118-20241118",obskeys=False,LonSys='1',
                         FiveMicron=False,Five_obskey='',IRTFdataset='', 
                         lats=LatLims,LonLims=[180,360],figsz=[6.0,6.0],
                         ROI=False, variance=False,localmax=False,segment=True,proj='NEZ',
                         ctbls=['terrain_r','Blues'], cont=True,bare_maps=False,cb=True, 
                         axNH3=False,axCH4=False,axRGB=False,LimbCorrection=True,lonhalfwidth=45,boxcar=9)
    
    #mapHST=read_HST(obskey="2024d_f631",LonSys='1')
    mapHST,dateHST=read_HST(obskey="2024d_f395n-f502n",LonSys='1')
    print(mapHST.shape)
    patch=mp.make_patch(mapHST,LatLims,LonLims,CM,LonRng,pad=True)
    
    figH,axsH=pl.subplots(1,figsize=(10,4.5), dpi=150, facecolor="white")
    axsH.grid(linewidth=0.2)
    axsH.ylim=[-90.,90.]
    axsH.xlim=[0.,360.]
    axsH.set_xticks(np.linspace(450.,0.,31), minor=False)
    xticklabels=np.array(np.mod(np.linspace(450,0,31),360))
    axsH.set_xticklabels(xticklabels.astype(int))
    axsH.set_yticks(np.linspace(-90,90,13), minor=False)
    yticklabels=np.array(np.linspace(-90,90,13))
    axsH.set_yticklabels(yticklabels.astype(int))
    axsH.tick_params(axis='both', which='major', labelsize=7)
    axsH.set_ylabel("PG Latitude (deg)")
    axsH.set_adjustable('box') 
    
    pp.plot_patch(patch,LatLims,LonLims,CM,LonRng,"gray",axsH,cbarplot=False,
                   cbar_title="Test",cbar_reverse=False,vn=100,vx=300,n=6)
    """
    fb.plot_regions_on_axis(axsH, labeled_fNH3, props_fNH3,lon_lims=np.array(LonLims)+180,lats=lats,
                 plot_contours=False, plot_masks=True,plot_labels=False,contour_color='C0')
    fb.plot_regions_on_axis(axsH, labeled_Plum, props_Plum,lon_lims=np.array(LonLims)+180,lats=lats,
                 plot_contours=False, plot_masks=True,plot_labels=False, contour_color='white')
    fb.plot_regions_on_axis(axsH, labeled_NEDF, props_NEDF,lon_lims=np.array(LonLims)+180,lats=lats,
                 plot_contours=False, plot_masks=True,plot_labels=False, contour_color='black')
    """
    pathmapplots="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Studies/NEZ/"
    figH.savefig(pathmapplots+"HST "+HSTkey+" Mean Sys"+LonSys+" "+lonstr+" "+latstr+" map.png",dpi=300)
