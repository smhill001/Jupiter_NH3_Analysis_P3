# -*- coding: utf-8 -*-

def make_HST_fits_script(obskeyHST,LonSys):
    
    import sys
    import socket
    hostname = socket.gethostname()
    from config_VA import Host_path
    from HST_Study_Maps import Study_Maps

    #sys.path.append('C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/HST/')
    #sys.path.append('C:/Astronomy/Projects/SAS 2021 Ammonia/Visualization-and-Analysis/')
    #sys.path.append()
    import read_HST_GO as HGO

    CoLatLims=Study_Maps[obskeyHST][LonSys]['CoLatLims']
    LonRng=Study_Maps[obskeyHST][LonSys]['LonRng']
    CMpref=Study_Maps[obskeyHST][LonSys]['CMpref']
     

    LonLimsWest=[CMpref-LonRng,CMpref+LonRng]
    HGO.HSTGO_process_and_plot(obskeyHST,CoLatLims,LonLimsWest,LonSys=LonSys)
        
