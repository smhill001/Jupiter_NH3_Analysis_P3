def make_sza_eza_planes(dateobs='2023-10-17 04:04:42',First=True):
    """
    Created on Sun Dec 17 07:08:46 2023
    
    PURPOSE: Make mu and mu_0 backplanes for L3 FITS files.  These backplanes 
             should probably be applied further upstream in the processsing
             pipeline, but for the time being, they can serve as part of
             a semiempirical correction when plotting maps of fNH3 and PCloud.
             
             
    @author: smhil
    """
    
    import planetmapper as pm
    import spiceypy as spice
    import matplotlib.pyplot as pl
    import numpy as np
    
    import sys
    drive='c:'
    sys.path.append(drive+'/Astronomy/Python Play')
    #sys.path.append("C:/Users/smhil/spice_kernels/")
    #sys.path.append("C:/Users/smhil/spice_kernels/naif/generic_kernels/lsk/")
    #path="C:/Astronomy/Projects/SAS 2021 Ammonia/VLT MUSE/"
    
    #fn="2022-09-19-0440_Jupiter_620nm.fits"
    #pm.set_kernel_path("C:/Users/smhil/spice_kernels")
    #jtest=pm.observation(path+fn)
    if First:
        pm.base.prevent_kernel_loading()
        spice.furnsh("C:/Astronomy/Python Play/spice_kernels/naif/generic_kernels/lsk/naif0012.tls.pc")
        spice.furnsh("C:/Astronomy/Python Play/spice_kernels/naif/generic_kernels/pck/earth_200101_990825_predict.bpc")
        spice.furnsh("C:/Astronomy/Python Play/spice_kernels/naif/generic_kernels/pck/pck00011.tpc")
        spice.furnsh("C:/Astronomy/Python Play/spice_kernels/naif/generic_kernels/pck/gm_de440.tpc")
        spice.furnsh("C:/Astronomy/Python Play/spice_kernels/naif/generic_kernels/pck/gm_de431.tpc")
        
        spice.furnsh("C:/Astronomy/Python Play/spice_kernels/naif/generic_kernels/spk/planets/de430.bsp")
        spice.furnsh("C:/Astronomy/Python Play/spice_kernels/naif/generic_kernels/spk/satellites/jup365.bsp")
    
    
    body = pm.BodyXY('Jupiter', dateobs, sz=500)
    #y=body.set_disc_params(x0=250, y0=250, r0=200)
    eza=np.flipud(np.cos(body.get_backplane_map('EMISSION')*np.pi/180.))
    sza=np.flipud(np.cos(body.get_backplane_map('INCIDENCE')*np.pi/180.))
    #fig,axs=pl.subplots(3,1,figsize=(8.0,8.0), dpi=150, facecolor="white")

    #axs[0].imshow(1.0/eza,vmin=-5.,vmax=5.)  
    #axs[1].imshow(1.0/sza,vmin=-5.,vmax=5.)
    #amfdata=(1.0/sza+1.0/eza)/2.0
    #axs[2].imshow(amfdata,vmin=-5.,vmax=5.)
    #fn=path="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Maps/PlanetmapperTest.png"
    #fig.savefig(fn,dpi=300)
    #print("################")
    #print(dateobs)
    
    return(eza,sza)