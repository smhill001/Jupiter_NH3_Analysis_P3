# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 09:54:10 2023

@author: smhil
"""
def get_keff(filterdata,FilterList,CH4,NH3):
    """
    PURPOSE: Loads effective absorption and center wavelength for each filter
             in the FilterList into the dictionary filterdata. Essentially 
             a wrapper for the k_eff computed in K_eff to be loaded into
             the filterdata dictionary.

    Parameters
    ----------
    filterdata : DICT
        Dictionary containing all metadata relevant to a given filter
    FilterList : LIST
        List of filters (strings) to be looped over for the k_eff computation
    CH4 : numpy float array(TBC?)
        Array of absorption cross section values versus wavelength for methane
    NH3 : numpy float array(TBC?)
        Array of absorption cross section values versus wavelength for ammonia.

    Returns
    -------
    filterdata : DICT
        Dictionary containing all metadata relevant to a given filter

    """
    import NH3_Filter_Library_P3 as NFL
    for filtr in FilterList:

        filterdata[filtr]['keff_CH4'],filterdata[filtr]['leff_CH4']= \
            NFL.K_eff(filterdata[filtr]['FiltTrans'],CH4,\
            filterdata[filtr]['halfwdth'],filtr)
        filterdata[filtr]['keff_NH3'],filterdata[filtr]['leff_NH3']= \
            NFL.K_eff(filterdata[filtr]['FiltTrans'],NH3,\
            filterdata[filtr]['halfwdth'],filtr)
                
    return(filterdata)
