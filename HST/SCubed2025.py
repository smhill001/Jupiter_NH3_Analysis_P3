# -*- coding: utf-8 -*-
#!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!! TO BE DEPRECATED ON 7/1/2026
# REPLACED BY:
    # HST_Study_Maps.py
    # HST_Analysis_Script.py
    # make_HST_fits_script.py
    
def SCubed2025(obskeyHST,LonSys,ROI_ID=False,makefits=False,HST=True,SCT=False):
    
    import sys
    import socket
    hostname = socket.gethostname()
    from config_VA import Host_path
    #sys.path.append('C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/HST/')
    #sys.path.append('C:/Astronomy/Projects/SAS 2021 Ammonia/Visualization-and-Analysis/')
    #sys.path.append()
    import L3_Jup_Map_Plot_V2 as L3MP
    #import read_HST_GO as HGO

    ROI=False
    segment=False
    compare=False

    ###############################################################################
    ###############################################################################
    #
    # PJ 77 - HST data 2025-10-16
    #
    ###############################################################################
    ###############################################################################
    # 20251015UTa HST
    ###############################################################################  
    
    Study_Maps={'20251016UTa':{'collection':'20251016-20251016',
                                '3':{'plotoptions':["contours"],
                                     'CoLatLims':[45,135],
                                     'LonRng':45,
                                     'CMpref':255},
                                '1':{'plotoptions':["scatter","resid","correl"],
                                     'CoLatLims':[75,105],
                                     'LonRng':45,
                                     'CMpref':280,
                                     'NEZ-East':{'ROI':{"Hot Spot":[82,84,267.0,2.0],
                                                       "Gyre":[85,87,270.0,3.0],
                                                       "Cloud Plume":[82,84,262.0,2.0],
                                                       "Reference":[75,77,270.0,5.0]},
                                                'CoLatLims':[75,90],
                                                'LonRng':15,
                                                'CMpref':270
                                                },
                                     'NEZ-West':{'ROI':{"Hot Spot":[81,82,297.0,2.0],
                                                       "Gyre":[83,86,299.0,4.0],
                                                       "Cloud Plume":[81,83,293.0,2.0],
                                                       "Reference":[75,77,295.0,5.0]},
                                                'CoLatLims':[75,90],
                                                'LonRng':15,
                                                'CMpref':295
                                                }
                                     }

                                },
                '20251016UTc':{'collection':'20251017-20251017',
                               '3':{'plotoptions':["contours"],
                                    'CoLatLims':[45,135],
                                    'LonRng':30,
                                    'CMpref':30},
                               '1':{'plotoptions':["scatter","resid","correl"],
                                    'CoLatLims':[75,105],
                                    'LonRng':45,
                                    'CMpref':45,
                                    'NEZ-East':{'ROI':{"Hot Spot":[82,83,14.0,2.0],
                                                      "Gyre":[84,86,15.0,3.0],
                                                      "Cloud Plume":[82,84,5.0,3.0],
                                                      "Reference":[76,78,15,4.0]},
                                               'CoLatLims':[75,90],
                                               'LonRng':15,
                                               'CMpref':15
                                               },
                                    'NEZ-West':{'ROI':{"Hot Spot":[83,84,53.0,2.0],
                                                      "Gyre":[84,85,53.0,3.0],
                                                      "Cloud Plume":[84,86,34.0,3.0],
                                                      "Reference":[76,78,45,4.0]},
                                               'CoLatLims':[75,90],
                                               'LonRng':15,
                                               'CMpref':45
                                               },
                                    'SEZ-SED':{'ROI':{"Hot Spot":[98,100,61.0,2.0],
                                                      "Gyre":[95,97,63.0,4.0],
                                                      "Cloud Plume":[96,98,53.0,2.0],
                                                      "Reference":[104,105,60,8.0]},
                                               'CoLatLims':[90,105],
                                               'LonRng':15,
                                               'CMpref':60
                                               }
                                    }

                               },
                '20251016UTf':{'collection':'20251017-20251017',
                               '3':{'plotoptions':["contours"],
                                    'CoLatLims':[45,135],
                                    'LonRng':45,
                                    'CMpref':90},
                               '1':{'plotoptions':["scatter","resid","correl"],
                                    'CoLatLims':[75,105],
                                    'LonRng':45,
                                    'CMpref':105,
                                    'NEZ-East':{'ROI':{"Hot Spot":[83,85,96.0,3.0],
                                                      "Gyre":[85,87,98.0,4.0],
                                                      "Cloud Plume":[82,84,79.0,3.0],
                                                      "Reference":[76,78,90.0,4.0]},
                                               'CoLatLims':[75,90],
                                               'LonRng':15,
                                               'CMpref':90
                                               },
                                    'NEZ-West':{'ROI':{"Hot Spot":[82,84,125.0,3.0],
                                                      "Gyre":[84,86,126.0,3.0],
                                                      "Cloud Plume":[82,84,113.0,3.0],
                                                      "Reference":[76,78,120.0,4.0]},
                                               'CoLatLims':[75,90],
                                               'LonRng':15,
                                               'CMpref':120
                                               },
                                    'SEZ-SED':{'ROI':{"Hot Spot":[98,100,61.0,2.0],
                                                      "Gyre":[95,97,63.0,4.0],
                                                      "Cloud Plume":[96,98,53.0,2.0],
                                                      "Reference":[104,105,60,8.0]},
                                               'CoLatLims':[90,105],
                                               'LonRng':15,
                                               'CMpref':60
                                               }
                                    }

                               },
                '20251120UTa':{'collection':'20251116-20251116',
                               '3':{'plotoptions':["contours"],
                                    'CoLatLims':[45,135],
                                    'LonRng':30,
                                    'CMpref':135},
                               '1':{'plotoptions':["scatter","resid","correl"],
                                    'CoLatLims':[75,105],
                                    'LonRng':30,
                                    'CMpref':65
                                    }

                               },
                '20251120UTb':{'collection':'20251119-20251119',
                               '3':{'plotoptions':["contours"],
                                    'CoLatLims':[45,135],
                                    'LonRng':30,
                                    'CMpref':300},
                               '1':{'plotoptions':["scatter","resid","correl"],
                                    'CoLatLims':[75,105],
                                    'LonRng':30,
                                    'CMpref':230,
                                    'NEZ-East':{'ROI':{"Hot Spot":[81,82,215.0,3.0],
                                                      "Gyre":[83,85,215.0,4.0],
                                                      "Cloud Plume":[82,84,203.0,3.0],
                                                      "Reference":[76,78,215.0,4.0]},
                                               'CoLatLims':[75,90],
                                               'LonRng':15,
                                               'CMpref':215
                                               },

                                    }

                               },
                '20251120UTc':{'collection':'20251119-20251119',
                               '3':{'plotoptions':["contours"],
                                    'CoLatLims':[45,135],
                                    'LonRng':30,
                                    'CMpref':50},
                               '1':{'plotoptions':["scatter","resid","correl"],
                                    'CoLatLims':[75,105],
                                    'LonRng':30,
                                    'CMpref':330
                                    }

                               },
                '20251221UTa':{'collection':'20251222-20251222',
                               '3':{'plotoptions':["contours"],
                                    'CoLatLims':[45,135],
                                    'LonRng':30,
                                    'CMpref':90},
                               '1':{'plotoptions':["scatter","resid","correl"],
                                    'CoLatLims':[75,105],
                                    'LonRng':30,
                                    'CMpref':250
                                    }

                               },
                '20251221UTb':{'collection':'20251222-20251222',
                               '3':{'plotoptions':["contours"],
                                    'CoLatLims':[45,135],
                                    'LonRng':30,
                                    'CMpref':120},
                               '1':{'plotoptions':["scatter","resid","correl"],
                                    'CoLatLims':[75,105],
                                    'LonRng':30,
                                    'CMpref':285
                                    }

                               },
                '20251221UTc':{'collection':'20251222-20251222',
                               '3':{'plotoptions':["contours"],
                                    'CoLatLims':[45,135],
                                    'LonRng':30,
                                    'CMpref':310},
                               '1':{'plotoptions':["scatter","resid","correl"],
                                    'CoLatLims':[75,105],
                                    'LonRng':30,
                                    'CMpref':125
                                    }

                               }

                 }
    
    
  
    ###############################################################################
    ###############################################################################
    ###############################################################################
    
    print("@@@@@@@@@@@@@@@@@@@@@ ROI_ID=",ROI_ID)

    collection=Study_Maps[obskeyHST]['collection']
    plotoptions=Study_Maps[obskeyHST][LonSys]['plotoptions']
    if ROI_ID in Study_Maps[obskeyHST][LonSys]:
        ROI=Study_Maps[obskeyHST][LonSys][ROI_ID]['ROI']
        CoLatLims=Study_Maps[obskeyHST][LonSys][ROI_ID]['CoLatLims']
        LonRng=Study_Maps[obskeyHST][LonSys][ROI_ID]['LonRng']
        CMpref=Study_Maps[obskeyHST][LonSys][ROI_ID]['CMpref']
    else:
        CoLatLims=Study_Maps[obskeyHST][LonSys]['CoLatLims']
        LonRng=Study_Maps[obskeyHST][LonSys]['LonRng']
        CMpref=Study_Maps[obskeyHST][LonSys]['CMpref']
     

    if makefits:
        LonLimsWest=[CMpref-LonRng,CMpref+LonRng]
        HGO.HSTGO_process_and_plot(obskeyHST,CoLatLims,LonLimsWest,LonSys=LonSys)
        
    if HST:
        print("@@@@@@@@@@@@@@@@@@@@@ ROI_ID=",ROI_ID)
        L3MP.L3_Jup_Map_Plot_V2(obskey=obskeyHST, 
                                CoLatLims=CoLatLims,LonRng=LonRng,CMpref=CMpref,LonSys=LonSys, 
                                subproj='SCubed 2025/'+obskeyHST,plotoptions=plotoptions,
                                dataversion='H',smoothcont=5,segment=segment,ROI_ID=ROI_ID,
                                ROI=ROI,compare=compare)

    if SCT:
        L3MP.L3_Jup_Map_Plot_V2(obskey=collection, 
                                CoLatLims=CoLatLims,LonRng=LonRng,CMpref=CMpref,LonSys=LonSys, 
                                subproj='SCubed 2025/'+obskeyHST,plotoptions=plotoptions,
                                dataversion=2,smoothcont=0,segment=segment,ROI_ID=ROI_ID,
                                ROI=ROI,compare=compare)
