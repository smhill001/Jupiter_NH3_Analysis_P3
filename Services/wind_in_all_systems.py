# -*- coding: utf-8 -*-
"""
Created on Wed Apr  9 08:17:26 2025

@author: smhil

From ChatGPT 4/11/2025
You can call it “Jovian Atmospheric Extrema Tracker” — or just "Jovian Extrema" 
for short. I’ve saved it, so whenever you return, just mention that name and 
we’ll pick up right where we left off!
"""


import numpy as np

# Jupiter's 1-bar level radii (meters)
a = 71492e3  # Equatorial
b = 66854e3  # Polar

# IAU authoritative rotation periods (in seconds)
rotation_periods = {
    "System I": 35430.0,
    "System II": 35740.6,
    "System III": 35729.71
}

def jupiter_radius(phi_deg):
    """Compute radius at a given planetographic latitude."""
    phi_rad = np.radians(phi_deg)
    sin_phi = np.sin(phi_rad)
    cos_phi = np.cos(phi_rad)
    num = (a**2 * cos_phi)**2 + (b**2 * sin_phi)**2
    denom = (a * cos_phi)**2 + (b * sin_phi)**2
    r = np.sqrt(num / denom)
    return r

def compute_wind_speed(phi_deg, lon1_deg, lon2_deg, jd1, jd2):
    """Compute physical wind speed in m/s from longitudinal drift."""
    r = jupiter_radius(phi_deg)
    delta_lon_rad = np.radians(lon2_deg - lon1_deg)  # preserve sign
    distance = -r * delta_lon_rad  # negative: eastward drift is positive wind
    delta_t_sec = (jd2 - jd1) * 86400
    return distance / delta_t_sec

def system_velocity(phi_deg, system):
    """Compute zonal speed of a rotation system at given latitude."""
    T = rotation_periods[system]
    r = jupiter_radius(phi_deg)
    return (2 * np.pi / T) * r * np.cos(np.radians(phi_deg))  # m/s

def convert_to_deg_per_day(v_mps, phi_deg):
    """Convert zonal speed from m/s to degrees/day."""
    r = jupiter_radius(phi_deg)
    phi_rad = np.radians(phi_deg)
    if np.isclose(np.cos(phi_rad), 0):
        return np.nan  # near poles
    omega_radps = v_mps / (r * np.cos(phi_rad))
    return omega_radps * (180 / np.pi) * 86400  # deg/day

def wind_in_all_systems(phi_deg, lon1_deg, lon2_deg, jd1, jd2, reference_system):
    """
    Compute wind speed in the reference system and translate to all three systems.
    Returns: dict with wind speeds relative to System I, II, III in m/s and deg/day.
    """
    if reference_system not in rotation_periods:
        raise ValueError("Invalid reference system. Choose from: 'System I', 'System II', 'System III'.")

    physical_wind = compute_wind_speed(phi_deg, lon1_deg, lon2_deg, jd1, jd2)
    v_ref = system_velocity(phi_deg, reference_system)
    absolute_motion = v_ref + physical_wind  # inertial speed

    results = {}
    for system in rotation_periods:
        v_sys = system_velocity(phi_deg, system)
        wind_mps = absolute_motion - v_sys
        wind_deg_day = convert_to_deg_per_day(wind_mps, phi_deg)
        results[system] = {
            "m/s": wind_mps,
            "deg/day": wind_deg_day
        }

    return results

import csv

def export_wind_results_to_csv(filename, phi_deg, lon1_deg, lon2_deg, jd1, jd2, reference_system, wind_results):
    """
    Export wind results to a CSV file.
    """
    with open(filename, mode='w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([
            "Planetographic Latitude (deg)",
            "Longitude 1 (deg)",
            "Longitude 2 (deg)",
            "JD1",
            "JD2",
            "Reference System",
            "Relative System",
            "Wind Speed (m/s)",
            "Wind Speed (deg/day)"
        ])
        for system, vals in wind_results.items():
            writer.writerow([
                phi_deg,
                lon1_deg,
                lon2_deg,
                jd1,
                jd2,
                reference_system,
                system,
                f"{vals['m/s']:.4f}",
                f"{vals['deg/day']:.4f}"
            ])
    print(f"✅ Wind results saved to {filename}")
    
#EXAMPLE
wind_results=wind_in_all_systems(6, 308, 275, 2460642.77565138, 2460647.7011267,"System III")

path="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Studies/maps/"
export_wind_results_to_csv(
    filename=path+"jupiter_wind_speeds.csv",
    phi_deg=6,
    lon1_deg=308,
    lon2_deg=275,
    jd1=2460642.77565138,
    jd2=2460647.7011267,
    reference_system="System III",
    wind_results=wind_results
)