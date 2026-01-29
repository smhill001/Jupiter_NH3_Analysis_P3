def make_patch(Map,LatLims,LonLims,CM,LonRng,pad=True):
    """
    Purpose: 
        Make a map patch and handle the case where the data overlap
        the map edges. This is designed for a map with Jovian longitude
        conventions that with the left boundary at 360 ascending from
        the right boundary at 0. In WinJUPOS, the actual map setting
        shows the left boundary at zero, which is of course, also 360.
        
    Updates:
        2025-08-13: Accommodated full map sizes other than 180x360 to
                    facilitate use of Hubble data (1800x3600). Tested 
                    also on [1800,3600,3] HST *.TIF data and worked fine.
                    Can probably replace "make_patch_RGB.py"
    
    Parameters
    ----------
    Map : NUMPY ARRAY [nx180,nx360,<1>]
        DESCRIPTION.
    LatLims : NUMPY ARRAY [2]
        DESCRIPTION. Colatitudes of patch boundary in degrees. 
    LonLims : NUMPY ARRAY [2]
        DESCRIPTION. Initial and final longitudes of patch boundary. 
        !!! Need details of convention. (colongitudes?)
    CM : TYPE
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
    import copy
    
    scale=int(Map.shape[0]/180)
    print("######## scale=",scale)
    lon_max=360*scale
    LatLims=np.array(LatLims)*scale
    LonRng=LonRng*scale
    CM=CM*scale
    LonLims=np.array(LonLims)*scale
    
    print(lon_max,LatLims,LonRng,CM,LonLims)
    
    patch=np.copy(Map[LatLims[0]:LatLims[1],LonLims[0]:LonLims[1]])
    if CM<LonRng:
        print("******************  CM2deg<LonRng")
        patch=np.concatenate((np.copy(Map[LatLims[0]:LatLims[1],LonLims[0]-1:lon_max]),
                              np.copy(Map[LatLims[0]:LatLims[1],0:LonLims[1]-lon_max])),axis=1)
    if CM>lon_max-LonRng:
        print("******************  CM2deg>LonRng")

        patch=np.concatenate((np.copy(Map[LatLims[0]:LatLims[1],lon_max+LonLims[0]:lon_max]),
                              np.copy(Map[LatLims[0]:LatLims[1],0:LonLims[1]])),axis=1)
        print("lon_max+LonLims[0]:lon_max,0:LonLims[1]=",lon_max+LonLims[0],lon_max,0,LonLims[1])
    return patch    
