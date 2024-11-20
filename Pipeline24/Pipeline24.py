def Pipeline24(obskey,planet="Saturn"):

    import sys
    drive='c:'
    sys.path.append(drive+'/Astronomy/Python Play')
    sys.path.append(drive+'/Astronomy/Python Play/Util_P3')
    sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy')
    sys.path.append(drive+'/Astronomy/Python Play/SPLibraries_P3')
    import os
    import ConfigFiles as CF
    import get_L0A_png_list as GL
    import gui_fitter as GF
    
    dataset=GL.get_L0A_png_list(obskey,planet=planet)
    GF.gui_fitter(dataset,planet=planet)