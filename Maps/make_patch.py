# -*- coding: utf-8 -*-
"""
Created on Tue May 13 08:24:10 2025

@author: smhil
"""

def make_patch(Map,LatLims,LonLims,CM2deg,LonRng,pad=True):
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
    #print("**************")
    #print(LatLims[0],LatLims[1],LonLims[0],LonLims[1])
    patch=np.copy(Map[LatLims[0]:LatLims[1],LonLims[0]:LonLims[1]])
    if CM2deg<LonRng:
        #print("******************  CM2deg<LonRng")
        patch=np.concatenate((np.copy(Map[LatLims[0]:LatLims[1],LonLims[0]-1:360]),
                              np.copy(Map[LatLims[0]:LatLims[1],0:LonLims[1]-360])),axis=1)
    if CM2deg>360-LonRng:
        patch=np.concatenate((np.copy(Map[LatLims[0]:LatLims[1],360+LonLims[0]:360]),
                              np.copy(Map[LatLims[0]:LatLims[1],0:LonLims[1]])),axis=1)
    #if pad:
    #    patch_pad=np.pad(patch,5,mode='reflect')
    #print("####################### Patch shape",patch.shape)

    return patch