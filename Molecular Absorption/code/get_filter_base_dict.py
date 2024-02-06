# -*- coding: utf-8 -*-

def get_filter_base_dict(Telescope,FilterList=['620','632','647','656'],
                         Inst=True):
    """
    Created on Wed Jul 19 14:00:10 2023

    A microservice to provide the base dictionary with filter information,
    including where to find filter transmission files.

    Used by:
        JupiterFilterPerformance.py
            ->NH3_Filter_Library_P3.py/compute_filter_Jupiter_transmissions
                ->get_filter_base_dict

    @author: smhil
    """
    import numpy as np
    import GeneralSpecUtils_P3 as GSU

    
    if Telescope=='SCT':
         path='c:/Astronomy/Projects/Techniques/InstrumentPerformance-P3/Filters/'
         filterdict={'620':{'transfile':'620CH4/620CH4_Transmission.txt',
                            'filtname':'620CH4','halfwdth':10.},
                      '632':{'transfile':'632OI/632OI_Transmission.txt',
                             'filtname':'632OI','halfwdth':10.},
                      '647':{'transfile':'647CNT/647CNT_Transmission.txt',
                             'filtname':'647NH3','halfwdth':10.},
                      '656':{'transfile':'656HIA/656HIA_Transmission.txt',
                             'filtname':'656HIA','halfwdth':10.},
                      '658':{'transfile':'658NII/658NII_Transmission.txt',
                             'filtname':'658NII','halfwdth':5.},
                      '672':{'transfile':'672SII/672SII_Transmission.txt',
                             'filtname':'672SII','halfwdth':10.},
                      '730':{'transfile':'730OII/730OII_Transmission.txt',
                             'filtname':'730OII','halfwdth':10.},
                      '889':{'transfile':'889CH4/889CH4_Transmission.txt',
                             'filtname':'889CH4','halfwdth':10.},
                      '940':{'transfile':'940NIR/940NIR_Transmission.txt',
                             'filtname':'940NIR','halfwdth':10.}}
    elif Telescope=='VLT':
         path='C:/Astronomy/Projects/SAS 2021 Ammonia/VLT MUSE/'
         filterdict={'620':{'transfile':'620CH4_MUSE_Transmission.txt',
                            'filtname':'620CH4','halfwdth':10.},
                      '632':{'transfile':'632OI_MUSE_Transmission.txt',
                             'filtname':'632OI','halfwdth':10.},
                      '647':{'transfile':'647NH3_MUSE_Transmission.txt',
                             'filtname':'647NH3','halfwdth':10.},
                      '656':{'transfile':'656HIA_MUSE_Transmission.txt',
                             'filtname':'656HIA','halfwdth':10.}}
    
    telepath='C:/Astronomy/Projects/Techniques/InstrumentPerformance-P3/'
    TelescopePerformance=np.loadtxt(telepath+'SystemResponseCLR-1260mm200lpm.txt',usecols=range(2))
    
    for filtr in FilterList:
        TempFiltTrans=np.loadtxt(path+filterdict[filtr]['transfile'],
                                 usecols=range(2))
        if Inst and Telescope=='SCT':
            filterdict[filtr]['FiltTrans']=GSU.SpectrumMath(TelescopePerformance,
                                                            TempFiltTrans,"Multiply")
        else:
            filterdict[filtr]['FiltTrans']=TempFiltTrans

    return(path,filterdict)