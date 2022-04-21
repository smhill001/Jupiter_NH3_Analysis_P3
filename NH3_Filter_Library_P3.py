# -*- coding: utf-8 -*-
"""
Created on Fri Feb 18 09:27:30 2022

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
