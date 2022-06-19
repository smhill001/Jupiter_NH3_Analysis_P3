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

def K_eff(P,FilterTransmission,CH4_Crossection,wv1,wv2,filtername,axis):
    import numpy as np
    import GeneralSpecUtils_P3 as GSU
    import copy
    import NH3_Filter_Library_P3 as NFL

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
    return(keff_CH4,leff_CH4)

def Weighting_Function(P,Keff,filtername,axis):
    import numpy as np

    #Note Hill mean_mol_wt is 2.31 where Mendikoa is 2.22
    #Results in a 0.96 ratio between the resulting optical depths

    amagat=2.69e24 #Lodschmits number?
    gravity=2228.0
    mean_mol_wt=3.85e-24
    fCH4=1.81e-3
    STP=1.01e6

    kmatm=(P/1.0e5)*STP*fCH4/(amagat*mean_mol_wt*gravity)
    tau_hill=kmatm*Keff
    tau_mend=22.4e4*(P/1.0e5)*fCH4*Keff/(2.22*gravity)
    print("tau_hill/tau_mend = ",tau_hill/tau_mend)
    trans=np.exp(-2.0*tau_hill)
    Dtrans=np.diff(trans)
    DP=np.diff(P)
    PDP=0.5*DP+P[:-1]

    return(tau_hill)
    #axis.plot(kmatm,P,label='km-atm')
    #axis.plot(trans,P,label='Trans')
    #axis.plot(-Dtrans/np.max(-Dtrans),PDP,linewidth=0.5,label=filtername+" CH4")

def Rayleigh_Function(P,leff,filtername,axis):
    import numpy as np

    amagat=2.69e24 #Lodschmits number?
    gravity=2228.0
    mean_mol_wt=3.85e-24
    fCH4=1.81e-3
    STP=1.01e6
    
    Na=6.022E23
    
    A=8.14*10E-29
    B=1.28*10E-30
    C=1.61*10E-32
    #Jupiter Data
    nH2=0.839
    nHe=0.156
    mu=2.22

    wvln=leff/1000.
    sigmaH2=A/wvln**4+B/wvln**6+C/wvln**8
    sigmaHe=0.05*sigmaH2
    
    print("nH2*sigmaH2+nHe*sigmaHe",nH2*sigmaH2+nHe*sigmaHe)
    
    tau_R=(1.0e6*(P/1e5)/(mean_mol_wt*gravity))*(nH2*sigmaH2+nHe*sigmaHe)
    tauR=(P/1e5)*Na*nH2*sigmaH2*10**6/(mu*gravity)+(P/1e5)*Na*nHe*sigmaHe*10**6/(mu*gravity)
    print("tau_R/tauR",tau_R/tauR)
    
    trans_R=np.exp(-2.0*tau_R)
    Dtrans_R=np.diff(trans_R)
    transR=np.exp(-2.0*tauR)
    DtransR=np.diff(transR)
    
    DP=np.diff(P)
    PDP=0.5*DP+P[:-1]
    #axis.plot(-Dtrans_R/np.max(-Dtrans_R),PDP,linewidth=0.5,label=filtername+" R")
    #axis.plot(-DtransR/np.max(-DtransR),PDP,linewidth=0.5,label=filtername+" R")

    #g=Gravity(planet,lat);
    return(tau_R)


def Compute_Transmission(P,tau_R,tau_CH4,filtername,axistrans,axisweight):
    import numpy as np
    
    trans_R=np.exp(-2.0*tau_R)
    Dtrans_R=np.diff(trans_R)
    trans_CH4=np.exp(-2.0*tau_CH4)
    Dtrans_CH4=np.diff(trans_CH4)
    trans=np.exp(-2.0*(tau_R+tau_CH4))
    Dtrans=np.diff(trans)
    
    DP=np.diff(P)
    PDP=0.5*DP+P[:-1]
    axistrans.plot(trans,P,linewidth=1.0,label=filtername)
    axisweight.plot(-Dtrans/np.max(-Dtrans),PDP,linewidth=1.0,label=filtername)
