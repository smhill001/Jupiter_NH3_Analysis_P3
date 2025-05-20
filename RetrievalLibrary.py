# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 14:43:04 2023

@author: smhil
"""
"""
def make_patch(Map,LatLims,LonLims,CM2deg,LonRng,pad=True):
    """
    Purpose: Make a map patch and handle the case where the data overlap
             the map edges. This is designed for a map with Jovian longitude
             conventions that with the left boundary at 360 ascending from
             the right boundary at 0. In WinJUPOS, the actual map setting
             shows the left boundary at zero, which is of course, also 360.
    """
    import numpy as np
    print("**************")
    print(LatLims[0],LatLims[1],LonLims[0],LonLims[1])
    patch=np.copy(Map[LatLims[0]:LatLims[1],LonLims[0]:LonLims[1]])
    if CM2deg<LonRng:
        print("******************  CM2deg<LonRng")
        patch=np.concatenate((np.copy(Map[LatLims[0]:LatLims[1],LonLims[0]-1:360]),
                              np.copy(Map[LatLims[0]:LatLims[1],0:LonLims[1]-360])),axis=1)
    if CM2deg>360-LonRng:
        patch=np.concatenate((np.copy(Map[LatLims[0]:LatLims[1],360+LonLims[0]:360]),
                              np.copy(Map[LatLims[0]:LatLims[1],0:LonLims[1]])),axis=1)
    #if pad:
    #    patch_pad=np.pad(patch,5,mode='reflect')
    #print("####################### Patch shape",patch.shape)

    return patch
"""

"""
def make_contours_CH4_patch(ax,CH4Abs_conv,LatLims,LonLims,lvls=[0.71,0.73,0.75,0.77,0.79],frmt='%3.1e',clr='w'):
    """
    PURPOSE: Overlay countours of NH3 absorption data on Jovian maps.
             Specifically designed for equivalent widths with mean values of
             ~0.55nm
    """
    clrs = []
    for i in range(6):
        clrs.append(clr)
    cs=ax.contour(CH4Abs_conv,origin='upper', 
                  extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],90-LatLims[0]],
                  colors=clrs, alpha=0.8,levels=lvls,
                  #linewidths=[0.5,0.5,0.5,0.5,0.5,0.5],
                  linewidths=[1.0,1.0,1.0,1.0,1.0,1.0],
                  linestyles=['dashed','dashed','dashed','dashed','dashed'])
    #ax.clabel(cs,[19.0,19.5,20.0,20.5,21.0],inline=True,fmt='%2.1f',fontsize=8)
    #print(lvls)
    ax.clabel(cs,lvls,inline=True,fmt=frmt,fontsize=9)
    
"""
def make_contours_CH4(ax,CH4Abs_conv,lvls=[0.71,0.73,0.75,0.77,0.79],frmt='%3.1e'):
    """
    PURPOSE: Overlay countours of NH3 absorption data on Jovian maps.
             Specifically designed for equivalent widths with mean values of
             ~0.55nm
    """
    cs=ax.contour(CH4Abs_conv,origin='upper', 
                  colors=['w','w','w','w','w'], alpha=0.5,levels=lvls,
                  linewidths=[0.5,0.5,1.0,0.5,0.5],
                  linestyles=['dashed','dashed','solid','dashed','dashed'])
    #ax.clabel(cs,[19.0,19.5,20.0,20.5,21.0],inline=True,fmt='%2.1f',fontsize=8)
    ax.clabel(cs,lvls,inline=True,fmt=frmt,fontsize=8)

def sharpen(image):
    from astropy.convolution import Gaussian2DKernel
    from astropy.convolution import convolve
    kernel = Gaussian2DKernel(x_stddev=5)
    blurred=convolve(image,kernel)
    tst=image+0.99*(image-blurred) #Need to revalidate that this is correct
    return tst

def patchstats(patch):
    import numpy as np
    from numpy import nanmean,nanstd,nanmin,nanmax

    mean=np.nanmean(patch)
    stdev=np.nanstd(patch)
    minimum=np.nanmin(patch)
    maximum=np.nanmax(patch)
    
    return [mean,stdev,minimum,maximum]