# -*- coding: utf-8 -*-
"""
Created on Thu Aug  7 10:36:09 2025

@author: smhil
"""
import spiceypy as sp
from datetime import datetime

def convert_system3_to_I_II_spice(datetime_str: str, system3_lon: float):
    """
    Convert System III west longitude to Systems I and II at a given UTC datetime.
    
    Args:
        datetime_str: ISO 8601 UTC datetime string, e.g. "2025-08-07T12:00:00"
        system3_lon: System III west longitude (in degrees)
    
    Returns:
        Dictionary with System I and System II longitudes (west-positive)
    """
    system3_lon_rev=360.-system3_lon
    print(system3_lon_rev)
    # Load SPICE kernels (you can customize this)
    try:
        sp.furnsh('pck00010.tpc')
    except Exception as e:
        pass  # Assume already loaded

    # Convert UTC to Ephemeris Time (seconds past J2000)
    et = sp.utc2et(datetime_str)
    days_since_j2000 = et / 86400.0 
    print(days_since_j2000)
    # === IAU 2015 Rotational Elements ===
    # System III
    W0_III = 284.95       # deg at J2000
    rate_III = 870.536000 # deg/day

    # System I
    W0_I = 67.1           # deg at J2000
    rate_I = 877.900      # deg/day

    # System II
    W0_II = 43.3          # deg at J2000
    rate_II = 870.270     # deg/day

    # Compute current central meridians (CMLs) for each system
    CML_III = (W0_III + rate_III * days_since_j2000)# % 360
    CML_I   = (W0_I   + rate_I   * days_since_j2000)# % 360
    CML_II  = (W0_II  + rate_II  * days_since_j2000)# % 360

    # Convert System III west longitude to Systems I and II
    print(system3_lon, system3_lon_rev)
    print(CML_I, CML_II, CML_III)
    system1_lon = 360.-(system3_lon_rev + (CML_III - CML_I)) % 360
    system2_lon = 360.-(system3_lon_rev + (CML_III - CML_II)) % 360

    return {
        "System I": system1_lon,
        "System II": system2_lon,
        "System III": system3_lon,
        "Method": "SPICE System III + IAU WGCCRE fallback"
    }
