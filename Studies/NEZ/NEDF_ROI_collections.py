def NEDF_ROI_collections(collection="20241202-20241202 NEDF 340"):

    import sys
    drive='c:'
    sys.path.append(drive+'/Astronomy/Python Play')
    sys.path.append(drive+'/Astronomy/Python Play/Util_P3')
    sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy')
    sys.path.append('./Services')
    
    import os
    import pylab as pl
    import numpy as np
    import matplotlib.dates as mdates
    from datetime import datetime

    import L3_Jup_Map_Plot as L3JMP
    import get_map_collection as gmc
    import MakeContiguousMap as MCM

    ###########################################################################
    # GRAB OBSERVATIONS (KEY) LISTS AND THEN SET ROI BOUNDARIES  AND CM
    # SEQUENCE IF NEEDED TO TRACK A DRIFING TARGET (SYS I)
    ###########################################################################
    obslist,dummy=gmc.get_map_collection(collection)
    
    if collection=="20230827-20240301 NEDF":
        ROI={"NEDF1":[80,83,207.0,2.0],
             "Fest1":[85,88,215.0,2.0],
             "Plume1":[85,88,207.0,2.0],
             "+NH3":[84,89,210,10.0],
             "-NH3":[84,89,185,10.0],
             "STrZ Ref":[91,93,200.0,3.0]}     
        obslist=["20230827UTa","20230905UTa","20230922UTa","20230929UTa",
                 "20231006UTa"]
        CM=[200]

    if collection=="20241115-20241115":
        ROI={"CORE":[91,93,200.0,2.0],
             "NE Plume":[91,93,200.0,2.0],
             "SW Plume":[91,93,200.0,2.0],
             "SEB Ref1":[91,93,200,2.0],
             "SEB Ref2":[91,93,200,2.0],
             "STrZ Ref":[91,93,200.0,3.0]}     
        obslist=["20241115UTa","20241115UTb","20241115UTc","20241115UTd","20241115UTe"]
        CM=[200]

    if collection=="20241129-20241129 NEDF 235":
        ROI={"NEDF1":[80,82,237.0,2.0],
             #"Fest1":[82,84,243.0,2.0],
             #"Fest2":[92,94,243.0,2.0],
             "Plume1":[83,85,229.0,2.0],
             "+NH3":[85,88,236,2.0],
             "NEB Ref":[75,77,232,8.0]}
        obslist=["20241129UTc","20241129UTd","20241129UTe","20241129UTf",
                 "20241129UTg"]
        CM=[235]
        
    if collection=="20241129-20241129 NEDF 280":
        ROI={"NEDF1":[80,83,282.0,2.0],
             "Fest1":[84,87,288.0,2.0],
             "Fest2":[92,94,288.0,2.0],
             "Plume1":[83,85,274.0,2.0],
             "+NH3":[84,87,279,2.0],
             "NEB Ref":[75,77,277,8.0]}
        obslist=["20241129UTg","20241129UTh","20241129UTi","20241129UTj","20241129UTk",
                 "20241129UTl","20241129UTm"]
        CM=[280]

    if collection=="20241129-20241129 NEDF 322":
        ROI={"NEDF1":[80,83,327.0,2.0],
             "Fest1":[84,87,333.0,2.0],
             "Fest2":[92,94,333.0,2.0],
             "Plume1":[83,85,319.0,2.0],
             "+NH3":[84,87,324,2.0],
             "NEB Ref":[75,77,322,8.0]}
        obslist=["20241129UTj","20241129UTk","20241129UTl","20241129UTm","20241129UTn"]
        CM=[322]

    if collection=="20241129-20241129 NEDF 340":
        ROI={"NEDF1":[80,83,353.0,2.0],
             "Fest1":[84,87,357.0,2.0],
             "Fest2":[92,94,353.0,2.0],
             "Plume1":[82,84,341.0,2.0],
             "+NH3":[84,87,351,2.0],
             "NEB Ref":[75,77,347,8.0]}
        obslist=["20241129UTj","20241129UTk","20241129UTl","20241129UTm","20241129UTn"]
        CM=[340]
        
    if collection=="20241202-20241202 NEDF 235":
        ROI={"NEDF1":[80,82,237.0,2.0],
             "Fest1":[82,84,243.0,2.0],
             "Fest2":[92,94,243.0,2.0],
             "Plume1":[83,85,229.0,2.0],
             "+NH3":[85,88,236,2.0],
             "NEB Ref":[75,77,232,8.0]}
        obslist=["20241202UTa","20241202UTb","20241202UTc"]
        CM=[235]
            
    if collection=="20241202-20241202 NEDF 280":
        ROI={"NEDF1":[80,83,282.0,2.0],
             "Fest1":[84,87,288.0,2.0],
             "Fest2":[92,94,288.0,2.0],
             "Plume1":[83,85,274.0,2.0],
             "+NH3":[84,87,279,2.0],
             "NEB Ref":[75,77,277,8.0]}
        obslist=["20241202UTa","20241202UTb","20241202UTc","20241202UTd",
                 "20241202UTe","20241202UTf"]
        CM=[280]

    if collection=="20241202-20241202 NEDF 322":
        ROI={"NEDF1":[80,83,327.0,2.0],
             "Fest1":[84,87,333.0,2.0],
             "Fest2":[92,94,333.0,2.0],
             "Plume1":[83,85,319.0,2.0],
             "+NH3":[84,87,324,2.0],
             "NEB Ref":[75,77,322,8.0]}
        obslist=["20241202UTd","20241202UTe",
                 "20241202UTf","20241202UTg"]
        CM=[322]

    if collection=="20241202-20241202 NEDF 340":
        ROI={"NEDF1":[80,83,353.0,2.0],
             "Fest1":[84,87,357.0,2.0],
             "Fest2":[92,94,353.0,2.0],
             "Plume1":[82,84,341.0,2.0],
             "+NH3":[84,87,351,2.0],
             "NEB Ref":[75,77,347,8.0]}
        obslist=["20241202UTe",
                 "20241202UTf","20241202UTg",
                 "20241202UTh"]
        CM=[340]
        
    if collection=="20241105-20241105 NEDF A":
        ROI={"NEDF1":[80,83,25.0,2.0],
             "Plume1":[82,84,15.0,2.0],
             "+NH3":[84,87,25,7.0],
             "NEB Ref":[75,77,27,8.0]}
        obslist=['20241105UTa',
                 '20241105UTb',
                 '20241105UTc',
                 '20241105UTd',
                 '20241105UTe'#,
                 #'20241105UTf'#,
                 #'20241105UTg',
                 #'20241105UTh',
                 #'20241105UTi',
                 #'20241105UTj',               
                 ]
        CM=[25]
        
    if collection=="20241118-20241118 NEDF A":
        ROI={"NEDF1":[79,82,35.0,2.0],
             "Plume1":[83,85,21.0,2.0],
             "+NH3":[83,86,29,3.0],
             "NEB Ref":[75,77,29,8.0]}
        obslist=[#'20241118UTa',
                 #'20241118UTb',
                 #'20241118UTc',
                 #'20241118UTd',
                 #'20241118UTe',
                 #'20241118UTf',
                 #'20241118UTg',
                 #'20241118UTh',
                 #'20241118UTi',
                 #'20241118UTj',
                 #'20241118UTk',
                 '20241118UTl',
                 '20241118UTm',
                 '20241118UTn',
                 '20241118UTo',
                 '20241118UTp'               
                 ]
        CM=[30]
        
    if collection=="20241128-20241128 NEDF A":
        ROI={"NEDF1":[78,81,43.0,2.0],
             "Plume1":[85,87,29.0,2.0],
             "+NH3":[82,85,36,3.0],
             "NEB Ref":[75,77,35,8.0]}
        obslist=['20241128UTa',
                 '20241128UTb',
                 '20241128UTc',
                 '20241128UTd',
                 '20241128UTe',
                 '20241128UTf'#,
                 #'20241128UTg',
                 #'20241128UTh',
                 #'20241128UTi',
                 #'20241128UTj',
                 #'20241128UTk',
                 #'20241128UTl',
                 #'20241128UTm',
                 #'20241128UTn',
                 #'20241128UTo'               
                 ]
        CM=[35]

    if collection=="20241203-20241203 NEDF A":
        ROI={"NEDF1":[79,82,41.0,2.0],
             "Plume1":[82,84,21.0,2.0],
             "+NH3":[82,85,32,3.0],
             "NEB Ref":[75,77,35,8.0]}
        obslist=['20241203UTa',
                 '20241203UTb',
                 '20241203UTc',
                 '20241203UTd'#,
                 #'20241203UTe',
                 #'20241203UTf',
                 #'20241203UTg',
                 #'20241203UTh',
                 #'20241203UTi',
                 #'20241203UTj',
                 #'20241203UTk',
                 #'20241203UTl',
                 #'20241203UTm',
                 #'20241203UTn',
                 #'20241203UTo',
                 #'20241203UTp'               
                 ]
        CM=[30]
        
    if collection=="20241205-20241205 NEDF A":
        ROI={"NEDF1":[80,83,42.0,2.0],
             "Plume1":[81,84,16.0,2.0],
             "+NH3":[83,86,32,3.0],
             "NEB Ref":[75,77,35,8.0]}
        obslist=[#'20241205UTa',
                 #'20241205UTb',
                 '20241205UTc',
                 '20241205UTd',
                 '20241205UTe',
                 '20241205UTf',
                 '20241205UTg'#,
                 #'20241205UTh'#,
                 #'20241205UTi'#,
                 #'20241205UTj',
                 #'20241205UTk',
                 #'20241205UTl',
                 #'20241205UTm',
                 #'20241205UTn',
                 #'20241205UTo',
                 #'20241205UTp'               
                 ]
        CM=[30]
        
    if collection=="20250106-20250106 NEDF CD":
        ROI={"NEDF C":[80,83,145.0,2.0],
             "Plume C":[83,86,135.0,4.0],
             "+NH3 C":[85,88,144,3.0],
             
             "NEDF D":[80,84,176.0,2.0],
             "Plume D":[82,86,165.0,3.0],
             "+NH3 D":[84,88,173,3.0],
             
             "NEB Ref":[75,77,155,10.0]}
        obslist=[#'20250106UTa', #Gradient issue
                 #'20250106UTb', #some gradient, and maybe an artifact over GRS (at least in IGB)
                 #'20250106UTc',
                 #'20250106UTd',
                 #'20250106UTe', #FF artifact?
                 #'20250106UTf', 
                 #'20250106UTg', #some gradient
                 #'20250106UTh',
                 #'20250106UTi', #a little offset time with the NH3 images, but seems okay
                 #'20250106UTj', #FF artifacts (2x)
                 #'20250106UTk',
                 #'20250106UTl', #FF artifacts (4x?)
                 #'20250106UTm',
                 #'20250106UTn', #Only a single Ha image
                 #'20250106UTo',
                 #'20250106UTp', #Gradient issues
                 '20250106UTq', #Very sharp! - great views of SEB outbreak
                 '20250106UTr', #Very sharp! - great views of SEB outbreak
                 '20250106UTs', #Very sharp - great views of SEB outbreak; FF artifacts in NEB
                 #'20250106UTt'            
                 ]
        CM=[155]
        
    ###########################################################################
    # SEB Collections
    ###########################################################################
        
    if collection=="20241115-20241115 SEB OB":
        ROI={"CORE":[103,106,324.0,2.0],
             "NE Plume":[99,102,320.0,2.0],
             "SW Plume":[105,108,329.0,2.0],
             "SEB Ref1":[102,104,334,2.0],
             "SEB Ref2":[100,108,311,2.0],
             "STrZ Ref":[111,115,323.0,3.0]}     
        obslist=["20241115UTa","20241115UTb","20241115UTc","20241115UTd","20241115UTe"]
        CM=[322]

    if collection=="20241129-20241129 SEB OB":
        ROI={"CORE":[102,105,323.0,2.0],
             "NE Plume":[98,101,320.0,2.0],
             "SW Plume":[104,107,329.0,2.0],
             "SEB Ref1":[102,104,334,2.0],
             "SEB Ref2":[105,108,313,3.0],
             "STrZ Ref":[111,113,323.0,5.0]}
        obslist=["20241129UTg","20241129UTh","20241129UTi",
                 "20241129UTj","20241129UTk","20241129UTl","20241129UTm"]
        CM=[322]
    
    if collection=="20241202-20241202 SEB OB":
        ROI={"CORE":[102,105,326.0,2.0],
             "NE Plume":[101,104,318.0,2.0],
             "SW Plume":[106,109,333.0,2.0],
             "SEB Ref1":[102,104,336,2.0],
             "SEB Ref2":[105,108,313,3.0],
             "STrZ Ref":[111,113,326.0,5.0]}
        obslist=["20241202UTc","20241202UTd","20241202UTe","20241202UTf",
                 "20241202UTg","20241202UTh"]
        CM=[325]
   
    return(ROI,obslist,CM)