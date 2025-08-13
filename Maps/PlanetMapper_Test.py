# -*- coding: utf-8 -*-
"""
Created on Fri Dec 15 09:07:32 2023

@author: smhil
"""

import planetmapper as pm
import spiceypy as spice
import matplotlib.pyplot as pl

import sys
drive='c:'
sys.path.append(drive+'/Astronomy/Python Play')
#sys.path.append("C:/Users/smhil/spice_kernels/")
#sys.path.append("C:/Users/smhil/spice_kernels/naif/generic_kernels/lsk/")
#path="C:/Astronomy/Projects/SAS 2021 Ammonia/VLT MUSE/"

#fn="2022-09-19-0440_Jupiter_620nm.fits"
#pm.set_kernel_path("C:/Users/smhil/spice_kernels")
#jtest=pm.observation(path+fn)

pm.base.prevent_kernel_loading()
spice.furnsh("C:/Astronomy/Python Play/spice_kernels/naif/generic_kernels/lsk/naif0012.tls.pc")
spice.furnsh("C:/Astronomy/Python Play/spice_kernels/naif/generic_kernels/pck/earth_200101_990825_predict.bpc")
spice.furnsh("C:/Astronomy/Python Play/spice_kernels/naif/generic_kernels/pck/pck00011.tpc")
spice.furnsh("C:/Astronomy/Python Play/spice_kernels/naif/generic_kernels/pck/pck00010.tpc")
spice.furnsh("C:/Astronomy/Python Play/spice_kernels/naif/generic_kernels/pck/gm_de440.tpc")
spice.furnsh("C:/Astronomy/Python Play/spice_kernels/naif/generic_kernels/pck/gm_de431.tpc")

spice.furnsh("C:/Astronomy/Python Play/spice_kernels/naif/generic_kernels/spk/planets/de430.bsp")
spice.furnsh("C:/Astronomy/Python Play/spice_kernels/naif/generic_kernels/spk/satellites/jup365.bsp")
spice.furnsh("C:/Astronomy/Python Play/spice_kernels/naif/generic_kernels/spk/satellites/sat453.bsp")


body = pm.BodyXY('Jupiter', '2025-03-02 04:59:00', sz=500)
y=body.set_disc_params(x0=250, y0=250, r0=200)
x=body.get_backplane_map('INCIDENCE')

print("body.subpoint_lon=",body.subpoint_lon)
print("body.subpoint_lat=",body.subpoint_lat)

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