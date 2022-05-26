# -*- coding: utf-8 -*-
"""
Created on Fri Feb 18 09:27:30 2022
    cont_absorption_calcs:
            PURPOSE: Compute the expected ratio between the an absorption
                     spectrum and a hypothetical continuum spectrum convolved
                     with a given filter passband
            INPUTS:  Convolution of the continuum and filter passband,
                     convolution of the absorption and filter passband,
                     beginning and ending wavelengths, and the filter name.
            OUTPUTS: Integral of the continuum convolved with the filter and
                     integral of the absorption convolved with the filter
                     (NOTE: currently (5/25/22) ratios are printed but not 
                      returned)

@author: smhil
"""

def cont_absorption_calcs(ContinuumProduct,AbsorptionProduct,wv1,wv2,filtername):
    import numpy as np
    StartIndex=np.where(ContinuumProduct[:,0]==wv1)
    EndIndex=np.where(ContinuumProduct[:,0]==wv2)
    ContinIntegral=sum(ContinuumProduct[StartIndex[0][0]:EndIndex[0][0],1])
    AbsorpIntegral=sum(AbsorptionProduct[StartIndex[0][0]:EndIndex[0][0],1])
    
    print("########### Jupiter "+filtername+" absorption/continuum")
    print("index=",StartIndex[0][0],EndIndex[0][0])
    print("Contin, Absorp=",ContinIntegral,AbsorpIntegral)
    print("Ratio, 1-Ratio=",AbsorpIntegral/ContinIntegral,1.0-AbsorpIntegral/ContinIntegral)
    print("1/(1-Ratio)=",1.0/(1.0-AbsorpIntegral/ContinIntegral))
    print(" ")

    return ContinIntegral,AbsorpIntegral

def K_eff(FilterTransmission,CH4_Crossection,wv1,wv2,filtername):
    import numpy as np
    import GeneralSpecUtils_P3 as GSU
    import copy

    keff_CH4Product=GSU.SpectrumMath(FilterTransmission,CH4_Crossection,"Multiply")
    lam=copy.deepcopy(FilterTransmission)
    lam[:,1]=FilterTransmission[:,0]
    lamprod=GSU.SpectrumMath(FilterTransmission,lam,"Multiply")
    StartIndex=np.where(keff_CH4Product[:,0]==wv1)
    EndIndex=np.where(keff_CH4Product[:,0]==wv2)
    #keff_CH4=sum(keff_CH4Product[StartIndex889[0][0]:EndIndex889[0][0],1])/(898.-880.)
    keff_CH4=sum(keff_CH4Product[StartIndex[0][0]:EndIndex[0][0],1])/\
                 sum(FilterTransmission[StartIndex[0][0]:EndIndex[0][0],1])
    leff_CH4=sum(lamprod[StartIndex[0][0]:EndIndex[0][0],1])/\
                 sum(FilterTransmission[StartIndex[0][0]:EndIndex[0][0],1])
    print("################################")
    #print(filtername,"K_eff=",keff_CH4)
    print(filtername,"K_eff : %4.3f, leff : %5.2f" % (keff_CH4, leff_CH4))
    print("")
    
