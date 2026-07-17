# -*- coding: utf-8 -*-
"""
Created on Fri Dec 15 09:07:32 2023

@author: smhil
"""

import planetmapper as pm
import spiceypy as spice
import matplotlib.pyplot as pl

import sys
sys.path.append('../../Visualization-and-Analysis/') #!!!!!! New change

import socket
hostname = socket.gethostname()
from config_VA import spice_path

drive='c:'
sys.path.append(drive+'/Astronomy/Python Play')
#sys.path.append("C:/Users/smhil/spice_kernels/")
#sys.path.append("C:/Users/smhil/spice_kernels/naif/generic_kernels/lsk/")
#path="C:/Astronomy/Projects/SAS 2021 Ammonia/VLT MUSE/"

#fn="2022-09-19-0440_Jupiter_620nm.fits"
#pm.set_kernel_path("C:/Users/smhil/spice_kernels")
#jtest=pm.observation(path+fn)

pm.base.prevent_kernel_loading()
spice.furnsh(spice_path[hostname]+"/naif/generic_kernels/lsk/naif0012.tls.pc")
#spice.furnsh(spice_path[hostname]+"/naif/generic_kernels/pck/earth_200101_990825_predict.bpc")
spice.furnsh(spice_path[hostname]+"/naif/generic_kernels/pck/pck00011.tpc")
spice.furnsh(spice_path[hostname]+"/naif/generic_kernels/pck/pck00010.tpc")
spice.furnsh(spice_path[hostname]+"/naif/generic_kernels/pck/gm_de440.tpc")
spice.furnsh(spice_path[hostname]+"/naif/generic_kernels/pck/gm_de431.tpc")

spice.furnsh(spice_path[hostname]+"naif/generic_kernels/spk/planets/de430.bsp")
spice.furnsh(spice_path[hostname]+"/naif/generic_kernels/spk/satellites/jup365.bsp")
spice.furnsh(spice_path[hostname]+"/naif/generic_kernels/spk/satellites/sat459.bsp")

spice.furnsh(spice_path[hostname]+"/naif/JUNO/kernels/spk/spk_ref_251119_281001_260129.bsp")
spice.furnsh(spice_path[hostname]+"/naif/JUNO/kernels/spk/spk_ref_250724_281001_250724.bsp")
spice.furnsh(spice_path[hostname]+"/naif/JUNO/kernels/spk/spk_ref_250526_281001_250505.bsp")
spice.furnsh(spice_path[hostname]+"/naif/JUNO/kernels/spk/spk_ref_231110_251016_231110.bsp")

body = pm.BodyXY('Jupiter', '2026-03-02 04:59:00', sz=500,observer='EARTH')
y=body.set_disc_params(x0=250, y0=250, r0=200)
x=body.get_backplane_map('INCIDENCE')

print("body.subpoint_lon=",body.subpoint_lon)
print("body.subpoint_lat=",body.subpoint_lat)
print("body.subpoint_distance=",body.subpoint_distance)

import numpy as np
times = np.arange(
    np.datetime64('2025-10-17T14:30:00'),
    np.datetime64('2025-10-17T14:50:00'),
    np.timedelta64(1, 'm')
)
time_strings = np.datetime_as_string(times)
for time in time_strings:
    bodyJUNO = pm.BodyXY('Jupiter', time, sz=500,observer='JUNO')
    print("body.subpoint_lon=",bodyJUNO.subpoint_lon)
    print("body.subpoint_lat=",bodyJUNO.subpoint_lat)
    print("body.subpoint_distance=",bodyJUNO.subpoint_distance)
bodyJUNO = pm.BodyXY('Jupiter', ['2025-10-17T14:40:00'][0], sz=500,observer='JUNO')
print("body.subpoint_lon=",bodyJUNO.subpoint_lon)
print("body.subpoint_lat=",bodyJUNO.subpoint_lat)
print("body.subpoint_distance=",bodyJUNO.subpoint_distance)


#fig,axs=pl.subplots(figsize=(8.0,4.0), dpi=150, facecolor="white")
#axs.imshow(x)
#fn=path="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Maps/PlanetmapperTest.png"
#fig.savefig(fn,dpi=300)

#gui = pm.gui.GUI(J)
"""
import glob
import planetmapper

for path in sorted(glob.glob('data/*.fits')):
    # Running from Python allows you to customise SPICE settings like the aberration correction
    observation = planetmapper.Observation(path, aberration_correction='CN+S')

    # Run some custom setup
    observation.add_other_bodies_of_interest('Io', 'Europa', 'Ganymede', 'Callisto')
    observation.set_plate_scale_arcsec(42) # set a custom plate scale
    observation.rotation_from_wcs() # get the disc rotation from the header's WCS info

    # Run the GUI to fit the observation interactively
    # This will open a GUI window every loop
    coords = observation.run_gui()

    # More custom code can go here to use the fitted observation...
    # for example, we can print some values for the last click location
    if coords:
        x, y = coords[-1]
        print(observation.xy2lonlat(x, y))
        
        
import planetmapper

body = planetmapper.Body('Jupiter', '2020-01-01')
body.plot_wireframe_radec(show=True)     
"""   