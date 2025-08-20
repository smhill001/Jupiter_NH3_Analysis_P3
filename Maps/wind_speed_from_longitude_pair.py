import numpy as np
from datetime import datetime

def wind_speed_from_long_pair(lat_deg, lon1_deg, lon2_deg,
                      t1_str, t2_str,
                      system='III',
                      radius_equatorial_km=71492.0,
                      flattening=0.06487):
    """
    Compute zonal wind speed at a given planetographic latitude
    from two longitudes and times.

    Parameters
    ----------
    lat_deg : float
        Planetographic latitude in degrees (+N).
    lon1_deg, lon2_deg : float
        Longitudes at time1 and time2, in the chosen system (I, II, or III).
        Input convention: positive WESTWARD (Jovian standard).
    t1, t2 : datetime
        Two observation times.
    system : str
        'I', 'II', or 'III'. Input longitude system.
    radius_equatorial_km : float
        Equatorial radius of Jupiter [km].
    flattening : float
        1 - (polar_radius/equatorial_radius).

    Returns
    -------
    speed_mps : float
        Zonal wind speed in m/s, positive eastward relative to System III.
    """

    # --- Parse times
    try:
        t1 = datetime.strptime(t1_str, "%Y-%m-%d %H:%M:%S")
        t2 = datetime.strptime(t2_str, "%Y-%m-%d %H:%M:%S")
    except ValueError:
        raise ValueError("Times must be in format 'YYYY-MM-DD hh:mm:ss'")

    # --- Time delta in seconds
    dt = (t2 - t1).total_seconds()


    # --- Rotation periods (IAU 2000 values)
    P_sysI   = 9*3600 + 50*60 + 30.0     # sec
    P_sysII  = 9*3600 + 55*60 + 40.6     # sec
    P_sysIII = 9*3600 + 55*60 + 29.71    # sec

    omega_I   = 2*np.pi / P_sysI
    omega_II  = 2*np.pi / P_sysII
    omega_III = 2*np.pi / P_sysIII

    # --- Raw longitude displacement (positive westward input, flip sign)
    dlon = lon2_deg - lon1_deg
    dlon = (dlon + 180) % 360 - 180  # normalize to -180..180
    dlon_east = -dlon  # convert to positive eastward convention

    # --- Planetographic latitude to planetocentric latitude
    lat_rad = np.radians(lat_deg)
    e2 = 2*flattening - flattening**2
    lat_centric = np.arctan((1 - e2) * np.tan(lat_rad))

    # --- Local radius of rotation at given latitude
    R_eq = radius_equatorial_km * 1000.0
    R_lat = R_eq * np.cos(lat_centric)

    # --- Linear speed (m/s) relative to chosen system
    angular_rate = np.radians(dlon_east) / dt
    speed_rel_system = angular_rate * R_lat

    # --- Add/subtract base rotation rate difference vs System III
    if system.upper() == 'I':
        delta_omega = omega_I - omega_III
    elif system.upper() == 'II':
        delta_omega = omega_II - omega_III
    elif system.upper() == 'III':
        delta_omega = 0.0
    else:
        raise ValueError("system must be 'I', 'II', or 'III'")

    speed_mps = speed_rel_system + delta_omega * R_lat

    return speed_mps,dlon/dt*86400
