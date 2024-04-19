# -*- coding: utf-8 -*-
def read_master_calibration():
    """
    Created on Wed Jul 26 13:43:34 2023
    
    @author: smhil
    
    PURPOSE: To provide a single, configuration controlled source of 
             calibration and k_eff data. To be added will be the EW fits.
             
    UPDATES:
        4/15/2024: Updated VLT-Filter: 'CH4GlobalTrans': from 0.885 to 0.897
                   This reflects using the average of both 2022 MUSE observations
                   (7/30 & 9/19) instead of just the 9/19 MUSE observation
    """
    
    #!! ULTIMATELY THIS SHOULD COME FROM A MASTER CSV FILE WHERE THESE ARE 
    #   GENERATED OUTPUTS FROM OTHER VERSION CONTROLLED PROGRAMS.
    calibration={'AGU 2022':{'CH4GlobalTrans':0.858,'NH3GlobalTrans':0.928},
                 'Model 1':{'CH4GlobalTrans':0.880,'NH3GlobalTrans':0.960},
                 'Model 2':{'CH4GlobalTrans':0.878,'NH3GlobalTrans':0.940},
                 'Observed':{'CH4GlobalTrans':0.910,'NH3GlobalTrans':0.972},
                 'VLT-MUSE':{'CH4GlobalTrans':0.861,'NH3GlobalTrans':0.961},
                 'SCT-Obs-Final':{'CH4GlobalTrans':0.920,'NH3GlobalTrans':0.972},
                 'VLT-Obs-Final':{'CH4GlobalTrans':0.893,'NH3GlobalTrans':0.962},
                 'VLT-Filter':{'CH4GlobalTrans':0.897,'NH3GlobalTrans':0.964}}

    K_eff={'CH4_620':{'C11':0.427,'VLT':0.454},
           'NH3_647':{'C11':2.955,'VLT':3.129}}
    
    return(calibration, K_eff)
