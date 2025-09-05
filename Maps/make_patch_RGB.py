def make_patch_RGB(Map,LatLims,LonLims,CM2deg,LonRng,pad=True):
    """
    Purpose: Make a map patch and handle the case where the data overlap
             the map edges. This is designed for a map with Jovian longitude
             conventions that with the left boundary at 360 ascending from
             the right boundary at 0. In WinJUPOS, the actual map setting
             shows the left boundary at zero, which is of course, also 360.
    """
    import numpy as np
    
    patch=np.copy(Map[LatLims[0]:LatLims[1],LonLims[0]:LonLims[1],:])
    if CM2deg<LonRng:
        patch=np.concatenate((np.copy(Map[LatLims[0]:LatLims[1],LonLims[0]-1:360,:]),
                              np.copy(Map[LatLims[0]:LatLims[1],0:LonLims[1]-360,:])),axis=1)
    if CM2deg>360-LonRng:
        patch=np.concatenate((np.copy(Map[LatLims[0]:LatLims[1],360+LonLims[0]:360,:]),
                              np.copy(Map[LatLims[0]:LatLims[1],0:LonLims[1],:])),axis=1)
    #if pad:
    #    patch_pad=np.pad(patch,5,mode='reflect')

    return patch