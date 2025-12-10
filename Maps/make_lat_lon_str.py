def make_lat_lon_str(lats,LonLims):
    """
    Creates strings for insertion into file names for the latitude and
    longitude range of a give map.

    Parameters
    ----------
    lats : numpy array
        DESCRIPTION.
    LonLims : numpy array
        DESCRIPTION.

    Returns
    -------
    latstr : string
        Latitude range description string for use in filenames.
    lonstr : string
        Longitude range description string for use in filenames.

    """
    if int(lats[0])<90:
        latstr=str(90-lats[0])+"N"
    if int(lats[0])==90:
        latstr=str(90-lats[0])
    if int(lats[0])>90:
        latstr=str(lats[0]-90)+"S"
        
    if int(lats[1])<90:
        latstr=latstr+"-"+str(90-lats[1])+"N"
    if int(lats[1])==90:
        latstr=latstr+"-"+str(90-lats[1])
    if int(lats[1])>90:
        latstr=latstr+"-"+str(lats[1]-90)+"S"
        
    lonstr=str(LonLims[0])+"-"+str(LonLims[1])
    return latstr,lonstr