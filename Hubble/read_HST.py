import matplotlib.pyplot as pl
import numpy as np
import find_blob as fb
import MakeContiguousMap as MCM


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
            filteredData=file
    return filteredData

def read_HST(obskey="2024c_f631",LonSys='1'):
    import sys
    drive='c:'
    sys.path.append(drive+'/Astronomy/Python Play')
    sys.path.append(drive+'/Astronomy/Python Play/Util_P3')
    sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy')
    sys.path.append('./Services')

    from matplotlib.pyplot import imread
    from astropy.io import fits
    sys.path.append('./Services')
    import os
    import numpy as np
    import convert_system3_to_I_II_spice as clong
    
    pathHST='C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Hubble/Data/'
    filenames=os.listdir(pathHST)
    
    filename=getFilename(filenames,obskey)
    print(filename)
    if "fits" in filename:
        
        hdulist=fits.open(pathHST+filename)
        hdulist.info()
        hdr=hdulist[0].header
        data=hdulist[0].data
        hdulist.close()
    
        print(hdr["date-obs"])
        print("LonSys=",LonSys)
        print()
        longitudes=clong.convert_system3_to_I_II_spice(hdr["date-obs"], 180.)
        print("LonSys=",LonSys)
        print(longitudes)
            

    if "tif" in filename:
        data=imread(pathHST+filename)
        if "2024c" in filename:
            hdulist=fits.open(pathHST+"hlsp_opal_hst_wfc3-uvis_jupiter-2024c_f275w_v1_globalmap.fits")
        if "2024d" in filename:
            hdulist=fits.open(pathHST+"hlsp_opal_hst_wfc3-uvis_jupiter-2024d_f275w_v1_globalmap.fits")

        print("TEST")
        hdulist.info()
        hdr=hdulist[0].header
        hdulist.close()
    
        print(hdr["date-obs"])
        print("LonSys=",LonSys)
        print()
        longitudes=clong.convert_system3_to_I_II_spice(hdr["date-obs"], 180.)
     
    if LonSys=='1':
        roll=180.-longitudes["System I"]
        print("LonSys=",LonSys,longitudes["System I"])
    elif LonSys=='2':
        roll=180.-longitudes["System II"]
    else:
        roll=0

    if "fits" in filename:
        datar=np.roll(data,int(roll)*10,axis=1)
    if "tif" in filename:
        datar=np.roll(data,int(roll)*10,axis=1)

    print(filename)
    #pl.imshow(data)
    return datar
    
def make_patch_HST(Map,LatLims,LonLims,CM2deg,LonRng,pad=True):
    """
    Purpose: 
        Make a map patch and handle the case where the data overlap
        the map edges. This is designed for a map with Jovian longitude
        conventions that with the left boundary at 360 ascending from
        the right boundary at 0. In WinJUPOS, the actual map setting
        shows the left boundary at zero, which is of course, also 360.
    
    Parameters
    ----------
    Map : NUMPY ARRAY [180,360]
        DESCRIPTION.
    LatLims : NUMPY ARRAY [2]
        DESCRIPTION. Colatitudes of patch boundary. 
        !!! Need details of convention.
    LonLims : NUMPY ARRAY [2]
        DESCRIPTION. Initial and final longitudes of patch boundary. 
        !!! Need details of convention. (colongitudes?)
    CM2deg : TYPE
        DESCRIPTION. Central Meridian to center patch on
    LonRng : TYPE
        DESCRIPTION.
    pad : INTEGER, optional
        DESCRIPTION. The default is True. Doesn't seem like I've used this in
        ages and it's commented out. Appears to deal with array wrapping.

    Returns
    -------
    patch : TYPE
        DESCRIPTION.

    """
    import numpy as np
    print("**************")
    print(LatLims[0],LatLims[1],LonLims[0],LonLims[1])
    patch=np.copy(Map[LatLims[0]:LatLims[1],LonLims[0]:LonLims[1]])
    if CM2deg<LonRng:
        print("******************  CM2deg<LonRng")
        patch=np.concatenate((np.copy(Map[LatLims[0]:LatLims[1],LonLims[0]-1:3600]),
                              np.copy(Map[LatLims[0]:LatLims[1],0:LonLims[1]-3600])),axis=1)
    if CM2deg>3600-LonRng:
        patch=np.concatenate((np.copy(Map[LatLims[0]:LatLims[1],3600+LonLims[0]:3600]),
                              np.copy(Map[LatLims[0]:LatLims[1],0:LonLims[1]])),axis=1)

    return patch    

def plot_patch_HST(patch,LatLims,LonLims,CM2,LonRng,colorscale,axis,
               cbarplot=True,cbar_title="Test",cbar_reverse=False,vn=0.10,vx=0.20,n=6):
    """
    Created on Thu Apr 11 18:54:35 2024
    
    @author: smhil
    """    
    import numpy as np
    import pylab as pl
    #import make_patch as MP

    print("@@@@@@@@@@@@ LatLims, LonLims, CM2, LonRng",LatLims, LonLims, CM2, LonRng)
    #patch=MP.make_patch(fullmap,LatLims,LonLims,CM2,LonRng)
    np.nan_to_num(patch, copy=False, nan=0.0, posinf=0.0, neginf=0.0)
    tx=np.linspace(vn,vx,n,endpoint=True)

    print("******** extent=",[360-LonLims[0],360-LonLims[1],90-LatLims[1],
            90-LatLims[0]])
    show=axis.imshow(patch, colorscale, origin='upper',vmin=vn,vmax=vx,  
               extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],
                       90-LatLims[0]],
                       aspect="equal",
                       interpolation='none')        # No smoothing/interpolation
                

    im_ratio = patch.shape[0]/patch.shape[1]
    if cbarplot:
        cbar = pl.colorbar(show, ticks=tx, 
                   orientation='vertical',cmap='gist_heat',
                   ax=axis,fraction=0.046*im_ratio, pad=0.04)
        cbar.ax.set_yticklabels(np.around(tx,3))
        cbar.ax.tick_params(labelsize=6,color="k")#if iSession >1:
        cbar.ax.set_ylabel(cbar_title,size=8)#,loc="top")
        cbar.ax.yaxis.set_label_coords(-0.7, 0.5)
        if cbar_reverse:
            cbar.ax.invert_yaxis()
    #if colorscale=="Greys":
    #    cbar.set_label('Cloud Top Pressure (mb)',fontsize=7)

    return patch,vn,vx,tx

LatLims=[74,94]
CM=270
LonRng=90
LatLims=np.array(LatLims)
LonLims=np.array([360-int(CM+LonRng),360-int(CM-LonRng)])

LatLimsx=LatLims*10
LonRngx=LonRng*10
CMx=CM*10
LonLimsx=LonLims*10
print("########",LatLimsx,LonLimsx,CMx,LonRngx)
#LonLims=np.array(LonLims)*10

#mapHST=read_HST(obskey="2024d_f631",LonSys='1')
mapHST=read_HST(obskey="2024d_f395n-f502n",LonSys='1')
patch=make_patch_HST(mapHST,LatLimsx,LonLimsx,CMx,LonRngx,pad=True)

#fig,axs=pl.subplots(1,figsize=(8,3), dpi=150, facecolor="white")
#axs.imshow(
#    patch,
#    extent=[360, 180, -5, 15],  # [x_min, x_max, y_min, y_max]
#    aspect='equal',              # Prevent distortion
#    interpolation='none',cmap='gray')        # No smoothing/interpolation
    
figH,axsH=pl.subplots(1,figsize=(8,3), dpi=150, facecolor="white")

#axsH.grid(linewidth=0.2)
#axsH.ylim=[-5.,15.]
#axsH.xlim=[180.,360.]
#axsH.set_xticks(np.linspace(450.,0.,31), minor=False)
#xticklabels=np.array(np.mod(np.linspace(450,0,31),360))
#axsH.set_xticklabels(xticklabels.astype(int))
#axsH.set_yticks(np.linspace(-90,90,13), minor=False)
#yticklabels=np.array(np.linspace(-90,90,13))
#axsH.set_yticklabels(yticklabels.astype(int))
#axsH.tick_params(axis='both', which='major', labelsize=7)
#axsH.set_ylabel("PG Latitude (deg)")
#axsH.set_adjustable('box') 

#!!!!!!Need to reconcile differences in LonLims below versus that computed above for Hubble    
lats,blendweightPCloud,blendweightfNH3,blendRGBweight,labeled_fNH3, props_fNH3,labeled_Plum, props_Plum,labeled_NEDF, props_NEDF=\
    MCM.MakeContiguousMap(collection="20241118-20241118",obskeys=False,LonSys='1',
                     FiveMicron=False,Five_obskey='',IRTFdataset='', 
                     lats=LatLims,LonLims=[180,360],figsz=[6.0,6.0],
                     ROI=False, variance=False,localmax=False,segment=True,proj='maps',
                     ctbls=['terrain_r','Blues'], cont=True,bare_maps=False,cb=True, 
                     axNH3=False,axCH4=False,axRGB=False,LimbCorrection=True,lonhalfwidth=45,boxcar=9)

CM=270
LonRng=90
LatLims=np.array(LatLims)
CM=CM
LonRng=LonRng
LonLims=[360-int(CM+LonRng),360-int(CM-LonRng)]
plot_patch_HST(patch,LatLims,LonLims,CM,LonRng,"gray",axsH,cbarplot=True,
               cbar_title="Test",cbar_reverse=False,vn=100,vx=300,n=6)
"""axs1.grid(linewidth=0.2)
#axs1.ylim=[-450.,450.]
axs1.xlim=[3600-LonLims[0],3600-LonLims[1]]
axs1.set_xticks(np.linspace(4500,0,31), minor=False)
xticklabels=np.array(np.mod(np.linspace(4500,0,31),3600))
axs1.set_xticklabels(xticklabels.astype(int))
axs1.set_yticks(np.linspace(-450,450,7), minor=False)
axs1.tick_params(axis='both', which='major', labelsize=9)
"""



fb.plot_regions_on_axis(axsH, labeled_fNH3, props_fNH3,lon_lims=np.array(LonLims)+180,lats=lats,
             plot_contours=False, plot_masks=True,plot_labels=False,contour_color='C0')
fb.plot_regions_on_axis(axsH, labeled_Plum, props_Plum,lon_lims=np.array(LonLims)+180,lats=lats,
             plot_contours=False, plot_masks=True,plot_labels=False, contour_color='white')
fb.plot_regions_on_axis(axsH, labeled_NEDF, props_NEDF,lon_lims=np.array(LonLims)+180,lats=lats,
             plot_contours=False, plot_masks=True,plot_labels=False, contour_color='black')
#figz,axsz=pl.subplots(1,figsize=(8,3), dpi=150, facecolor="white")
#axsz.imshow(labeled_fNH3)
