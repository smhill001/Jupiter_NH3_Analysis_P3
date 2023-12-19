def sza_eza_planes():
    """
    Created on Sun Dec 17 07:08:46 2023
    
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
    spice.furnsh("C:/Astronomy/Python Play/spice_kernels/naif/generic_kernels/pck/gm_de440.tpc")
    spice.furnsh("C:/Astronomy/Python Play/spice_kernels/naif/generic_kernels/pck/gm_de431.tpc")
    
    spice.furnsh("C:/Astronomy/Python Play/spice_kernels/naif/generic_kernels/spk/planets/de430.bsp")
    spice.furnsh("C:/Astronomy/Python Play/spice_kernels/naif/generic_kernels/spk/satellites/jup365.bsp")
    
    
    body = pm.BodyXY('Jupiter', '2023-10-17 04:04:42', sz=500)
    y=body.set_disc_params(x0=250, y0=250, r0=200)
    eza=body.get_backplane_map('EMISSION')
    sza=body.get_backplane_map('INCIDENCE')
    fig,axs=pl.subplots(2,1,figsize=(8.0,8.0), dpi=150, facecolor="white")
    fig1,axs1=pl.subplots(figsize=(8.0,4.0), dpi=150, facecolor="white")

    axs[0].imshow(eza)  
    axs[1].imshow(sza)
    axs1.imshow(eza/sza)
    fn=path="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Maps/PlanetmapperTest.png"
    fig.savefig(fn,dpi=300)
    return(eza,sza)