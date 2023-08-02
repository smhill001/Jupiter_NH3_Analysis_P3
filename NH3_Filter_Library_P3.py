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
                     
     K_eff
        PURPOSE: Compute K_eff and l_eff for CH4 or NH3 by convolving the 
                 absorption corss section and filter transmission.
                 
     tau_gas_versus_P
        PURPOSE: Compute two-way optical depth, tau, due to gas absorption
                 as a function of pressure.
     tau_rayleigh_versus_P
     Compute_Transmission
     Get_Albedo_and_Absorption
     SpectralModeling
     MoonAlbedos

@author: smhil
"""

def cont_absorption_calcs(ContinuumProduct,AbsorptionProduct,wv1,wv2,filtername,prn=True):
    import numpy as np
    StartIndex=np.where(ContinuumProduct[:,0]==wv1)
    EndIndex=np.where(ContinuumProduct[:,0]==wv2)
    ContinIntegral=sum(ContinuumProduct[StartIndex[0][0]:EndIndex[0][0],1])
    AbsorpIntegral=sum(AbsorptionProduct[StartIndex[0][0]:EndIndex[0][0],1])
    TransIntegral=AbsorpIntegral/ContinIntegral
    
    if prn:
        print("########### Jupiter "+filtername+" absorption/continuum")
        print("index=",StartIndex[0][0],EndIndex[0][0])
        print("Contin, Absorp=",ContinIntegral,AbsorpIntegral)
        print("Ratio, 1-Ratio=",AbsorpIntegral/ContinIntegral,1.0-AbsorpIntegral/ContinIntegral)
        print("1/(1-Ratio)=",1.0/(1.0-AbsorpIntegral/ContinIntegral))
        print(" ")

    return ContinIntegral,AbsorpIntegral,TransIntegral

def K_eff(FilterTransmission,Abs_Crossection,halfwidth, filtr):
    """
    COMPUTE K_eff AND l_eff FOR AN ABSORBING GAS, CH4 OR NH3, GIVEN THE 
      CROSSECTION AND INPUT TRANSMISSION
    Calls:
        tbd
    Called by:
        -> get_keff 
    """
    import numpy as np
    import GeneralSpecUtils_P3 as GSU
    import copy
    import NH3_Filter_Library_P3 as NFL
    
    #Set index interval based on wavelength interval
    wv1,wv2=float(filtr)-float(halfwidth),float(filtr)+float(halfwidth)
    StartIndex=np.where(FilterTransmission[:,0]==wv1)
    EndIndex=np.where(FilterTransmission[:,0]==wv2)

    #Compute product of transmission and crossection and k_eff
    Abs_Product=GSU.SpectrumMath(FilterTransmission,Abs_Crossection,"Multiply")
    keff=sum(Abs_Product[StartIndex[0][0]:EndIndex[0][0],1])/ \
                 sum(FilterTransmission[StartIndex[0][0]:EndIndex[0][0],1])

    #Compute product of transmission and wavelength for center wavelength
    lam=copy.deepcopy(FilterTransmission)
    lam[:,1]=FilterTransmission[:,0]
    lamprod=GSU.SpectrumMath(FilterTransmission,lam,"Multiply")
    leff=sum(lamprod[StartIndex[0][0]:EndIndex[0][0],1])/ \
                 sum(FilterTransmission[StartIndex[0][0]:EndIndex[0][0],1])

    return(keff,leff)

def tau_gas_versus_P(P,Keff,filtername,axis,gas='CH4'):
    """
    ###########################################################################
    # COMPUTE TAU DUE TO GAS ABSORPTION AS A FUNCTION OF PRESSURE 
    #   (WEIGHTING FUNCTION)
    #   !!!! NEED TO ADD NH3 OPTION WITH SOME BASELINE ABUNDANCE OR PROFILE
    #   !!!! CAN COMPARE MENDIKOA TO HILL COMPUTATION
    #   !!!! HAS COMMENTED CODE FOR PLOTTING WHICH COULD BE USEFUL TO REACTIVATE
    #   !!!! SHOULD HAVE GRAVITY AS A FUNCTION OF LATITUDE?
    ###########################################################################
    """
    import numpy as np

    #Note Hill mean_mol_wt is 2.31 where Mendikoa is 2.22
    #Results in a 0.96 ratio between the resulting optical depths

    amagat=2.69e24 #Lodschmits number?
    gravity=2228.0
    mean_mol_wt=3.85e-24
    if gas=='CH4':
        fCH4=1.81e-3
    elif gas=="NH3":
        fCH4=1.5e-4
    STP=1.01e6

    kmatm=(P/1.0e5)*STP*fCH4/(amagat*mean_mol_wt*gravity)
    tau_gas=kmatm*Keff
    tau_mend=22.4e4*(P/1.0e5)*fCH4*Keff/(2.22*gravity)
    #print("tau_gas/tau_mend = ",tau_gas/tau_mend)
    trans=np.exp(-2.0*tau_gas)
    Dtrans=np.diff(trans)
    DP=np.diff(P)
    PDP=0.5*DP+P[:-1]

    return(tau_gas)
    #axis.plot(kmatm,P,label='km-atm')
    #axis.plot(trans,P,label='Trans')
    #axis.plot(-Dtrans/np.max(-Dtrans),PDP,linewidth=0.5,label=filtername+" CH4")

def tau_rayleigh_versus_P(P,leff,filtername,axis):
    ###########################################################################
    # COMPUTE TAU DUE TO RAYLEIGH SCATTERING AS A FUNCTION OF PRESSURE
    #   !!!! HAS COMMENTED CODE FOR PLOTTING WHICH COULD BE USEFUL TO REACTIVATE
    #   !!!! SHOULD HAVE GRAVITY AS A FUNCTION OF LATITUDE?
    ###########################################################################
    import numpy as np

    #amagat=2.69e24 #Lodschmits number?
    gravity=2228.0
    mean_mol_wt=3.85e-24
    #fCH4=1.81e-3
    #STP=1.01e6
    
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
    
    #print("nH2*sigmaH2+nHe*sigmaHe",nH2*sigmaH2+nHe*sigmaHe)
    
    tau_R=(1.0e6*(P/1e5)/(mean_mol_wt*gravity))*(nH2*sigmaH2+nHe*sigmaHe)
    tauR=(P/1e5)*Na*nH2*sigmaH2*10**6/(mu*gravity)+(P/1e5)*Na*nHe*sigmaHe*10**6/(mu*gravity)
    #print("tau_R/tauR",tau_R/tauR)
    
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


def Compute_Transmission(P,tau_R,tau_gas,filtername,axistrans,axisweight):
    import numpy as np
    
    trans_R=np.exp(-2.0*tau_R)
    Dtrans_R=np.diff(trans_R)
    trans_gas=np.exp(-2.0*tau_gas)
    Dtrans_gas=np.diff(trans_gas)
    trans=np.exp(-2.0*(tau_R+tau_gas))
    Dtrans=np.diff(trans)
    
    DP=np.diff(P)
    PDP=0.5*DP+P[:-1]
    axistrans.plot(trans,P,linewidth=1.0,label=filtername)
    axisweight.plot(-Dtrans/np.max(-Dtrans),PDP,linewidth=1.0,label=filtername)
    
    return(trans)
    
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 15:07:47 2022

Get Albedo and Molecular Data

@author: smhil
"""
        
def SpectralModeling(s_NH3=0.018,s_CH4=0.304,refl=0.53):
    """
    Model the spectrum of Jupiter using CH4 and NH3 column densities (km-atm)
    along with a constant reflectivity for the cloud tops. 
    !Currently does not include Rayleigh scattering and cloud top
    reflectivity is constant
    Also provides linear fit coefficients between transmission in a given
    filter and molecular band equivalent width
    """
    import sys
    sys.path.append('c:/Astronomy/Python Play')
    sys.path.append('c:/Astronomy/Python Play/Util_P3')
    sys.path.append('c:/Astronomy/Python Play/SPLibraries_P3')
    sys.path.append('c:/Astronomy/Python Play/SpectroPhotometry/Spectroscopy_P3')
    import matplotlib.pyplot as pl
    import numpy as np
    import GeneralSpecUtils_P3 as GSU
    import NH3_Filter_Library_P3 as NFL
    import copy
    sys.path.append('./Services')
    import get_albedo_continua_crossections as gACC

    path='c:/Astronomy/Projects/Techniques/InstrumentPerformance-P3/Filters/'
    filterdata={'620':{'transfile':'620CH4/620CH4_Transmission.txt',
                       'filtname':'620CH4','filtwdth':10.},
                '647':{'transfile':'647CNT/647CNT_Transmission.txt',
                    'filtname':'647NH3','filtwdth':10.},
                '620MUSE':{'transfile':'',
                          'filtname':'620MUSE','filtwdth':10},
                '647MUSE':{'transfile':'',
                          'filtname':'647MUSE','filtwdth':10}}

    #x0,x1,xtks=600.,1000.,9
    #y0,y1,ytks=0.0,0.7,8
    x0,x1,xtks=600.,680.,5
    y0,y1,ytks=0.4,0.6,5
    Albedo,Continua,CH4,NH3= \
        gACC.get_albedo_continua_crossections(x0,x1,xtks,y0,y1,ytks,
                                              Crossect=False)

    figtest,axtest=pl.subplots(1,1,figsize=(6.0,3.5), dpi=150, facecolor="white")

    axtest.plot(Albedo[:,0],Albedo[:,1],label='Jupiter Albedo (Karkoschka, 1994)',linewidth=1.0,color='C0')
    axtest.set_xlim(x0,x1)
    axtest.set_xticks(np.linspace(x0,x1,xtks, endpoint=True))
    axtest.set_ylim(y0,y1)
    axtest.set_yticks(np.linspace(y0,y1,ytks, endpoint=True))
    print(CH4.shape)
    CH4_trans=refl*(np.exp(-s_CH4*CH4[:,1]))
    NH3_trans=refl*(np.exp(-s_NH3*NH3[:,1]))
    gas_trans=refl*(np.exp(-s_CH4*CH4[:,1]-s_NH3*NH3[:,1]))
    print(CH4_trans.shape)
    axtest.plot(CH4[:,0],CH4_trans,label='CH4 '+str(s_CH4)[:5]+' km-atm',linewidth=0.5,color='C1')
    axtest.plot(NH3[:,0],NH3_trans,label='NH3 '+str(s_NH3)[:5]+' km-atm',linewidth=0.5,color='C2')
    axtest.plot(NH3[:,0],gas_trans,label='gas (CH4+NH3)',linewidth=1.0,color='C3')
    axtest.grid(linewidth=0.2)
    axtest.legend(fontsize=7, loc='best')
    axtest.tick_params(axis='both', which='major', labelsize=8)
    axtest.set_xlabel("Wavelength (nm)")
    axtest.set_title("Albedo Spectral Model")
    axtest.set_ylabel("Albedo")

    figtest.subplots_adjust(left=0.10, right=0.90, top=0.9, bottom=0.14)
    figtest.savefig('c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Albedo_Spectral_Model.png',dpi=320)

    ###########################################################################
    # Compute and plot transmission of the 647nm filter and the equivalent
    # width of the band as a function of the column abundance of NH3.
    ###########################################################################
    filterdata['647']['FiltTrans']=np.loadtxt(path+filterdata['647']['transfile'],usecols=range(2))

    F647_wvs=[642.,652.]
    F647_idx=[np.argmin(abs(F647_wvs[0]-filterdata['647']['FiltTrans'][:,0])),
              np.argmin(abs(F647_wvs[1]-filterdata['647']['FiltTrans'][:,0]))]
    filterdata['647MUSE']['FiltTrans']=copy.deepcopy(filterdata['647']['FiltTrans'])
    filterdata['647MUSE']['FiltTrans'][:,1]=np.zeros((len(filterdata['647MUSE']['FiltTrans'][:,1])))
    filterdata['647MUSE']['FiltTrans'][F647_idx[0]:F647_idx[1],1]=1.
    
    s_NH3_Arr=np.arange(0.0,0.031,0.002) #column abundance in km-atm
    NH3_trans_Arr=[]
    NH3_trans_Arr_MUSE=[]
    W_Arr=[]
    
    Model='Piecewise1'
    Tele='SCT'
    Continuum_Albedo=np.zeros((Continua[Model]['WaveGrid'].size,2))
    Continuum_Albedo[:,0]=Continua[Model]['WaveGrid']
    Continuum_Albedo[:,1]=Continua[Model]['Albedo']
    
    Continuum=copy.deepcopy(Continuum_Albedo)
    Continuum[:,1]=np.ones(Continuum_Albedo.shape[0])
    NH3_trans=copy.deepcopy(Continuum_Albedo)
    NH3_trans_MUSE=copy.deepcopy(Continuum_Albedo)
    print("C",Continuum_Albedo.shape)
    
    filterdata['647MUSE']['ContAlbedoProd']=GSU.SpectrumMath(filterdata['647MUSE']['FiltTrans'],Continuum_Albedo,"Multiply")
    filterdata['647MUSE']['AbsrProd']=GSU.SpectrumMath(filterdata['647MUSE']['FiltTrans'],Albedo,"Multiply")
    filterdata['647MUSE']['ContAlbedo_Int'],filterdata['647MUSE']['AbsrAlbedo_Int'],filterdata['647MUSE']['TransAlbedoInt']= \
        NFL.cont_absorption_calcs(filterdata['647MUSE']['ContAlbedoProd'],filterdata['647MUSE']['AbsrProd'], \
                                  float('647')-filterdata['647MUSE']['filtwdth'],\
                                  float('647')+filterdata['647MUSE']['filtwdth'], \
                                      filterdata['647MUSE']['filtname'])


    figslopeNH3,axslopeNH3=pl.subplots(1,1,figsize=(6.0,3.5), dpi=150, facecolor="white")
    for s in s_NH3_Arr:
        print("############ s=",s)
        temp=1.0*(np.exp(-s*NH3[:,1]))
        
        print('t',temp.shape)
        NH3_trans[:,1]=temp #NH3 transmission as a function of wavelength
        print(NH3_trans.shape)
        filterdata['647']['ContProd']=GSU.SpectrumMath(filterdata['647']['FiltTrans'],Continuum,"Multiply")
        filterdata['647']['AbsrProd']=GSU.SpectrumMath(filterdata['647']['FiltTrans'],NH3_trans,"Multiply")
        
        filterdata['647']['Cont_Int'],filterdata['647']['Absr_Int'],filterdata['647']['TransInt']= \
            NFL.cont_absorption_calcs(filterdata['647']['ContProd'],filterdata['647']['AbsrProd'], \
                                      float('647')-filterdata['647']['filtwdth'],\
                                      float('647')+filterdata['647']['filtwdth'], \
                                          filterdata['647']['filtname'])
                        
        filterdata['647']['Tau_Albedo']=-np.log(filterdata['647']['TransInt'])
        
        NH3_trans_Arr.append(filterdata['647']['TransInt'])

        ind1,ind2=np.argmin(abs(NH3_trans[:,0]-627.)),np.argmin(abs(NH3_trans[:,0]-667.))
        W=np.sum(1.0-NH3_trans[ind1:ind2,1])*0.5 #for 0.5 nm bins
        W_Arr.append(W)
        print("#### W=",W)
        
        ##########

        filterdata['647MUSE']['ContProd']=GSU.SpectrumMath(filterdata['647MUSE']['FiltTrans'],Continuum,"Multiply")
        filterdata['647MUSE']['AbsrProd']=GSU.SpectrumMath(filterdata['647MUSE']['FiltTrans'],NH3_trans,"Multiply")
        
        filterdata['647MUSE']['Cont_Int'],filterdata['647MUSE']['Absr_Int'],filterdata['647MUSE']['TransInt']= \
            NFL.cont_absorption_calcs(filterdata['647MUSE']['ContProd'],filterdata['647MUSE']['AbsrProd'], \
                                      float('647')-filterdata['647MUSE']['filtwdth'],\
                                      float('647')+filterdata['647MUSE']['filtwdth'], \
                                          filterdata['647MUSE']['filtname'])
                        
        filterdata['647MUSE']['Tau_Albedo']=-np.log(filterdata['647MUSE']['TransInt'])
        
        NH3_trans_Arr_MUSE.append(filterdata['647MUSE']['TransInt'])

    
    R_W=np.corrcoef(s_NH3_Arr,W_Arr)[0,1]
    R_trans=np.corrcoef(s_NH3_Arr,NH3_trans_Arr)[0,1]
    R_trans_MUSE=np.corrcoef(s_NH3_Arr,NH3_trans_Arr_MUSE)[0,1]
    NH3fit=np.polyfit(NH3_trans_Arr,W_Arr,1)
    NH3fitMUSE=np.polyfit(NH3_trans_Arr_MUSE,W_Arr,1)

    axslopeNH3.plot(s_NH3_Arr,NH3_trans_Arr,color='C0',label='Transmission, R='+str(R_trans)[:5])
    axslopeNH3.plot(s_NH3_Arr,NH3_trans_Arr_MUSE,color='C2',label='Transmission MUSE, R='+str(R_trans_MUSE)[:5])
    axWNH3=axslopeNH3.twinx()
    axWNH3.plot(s_NH3_Arr,W_Arr,color='C1',label='Eq. Width (nm), R='+str(R_W)[:4])
    axslopeNH3.grid(linewidth=0.2)
    axslopeNH3.set_xlabel("NH3 Column Abundance (km-atm)")
    axslopeNH3.set_ylabel("Transmission")
    axWNH3.set_ylabel("Equivalent Width (nm)")
    axslopeNH3.tick_params(axis='both', which='major', labelsize=8)
    axWNH3.tick_params(axis='both', which='major', labelsize=8)
    axslopeNH3.set_xlim(0.0,0.03)
    axslopeNH3.set_xticks(np.linspace(0.0,0.03,7, endpoint=True))
    axslopeNH3.set_ylim(0.9,1.02)
    axslopeNH3.set_yticks(np.linspace(0.9,1.02,7, endpoint=True))
    axWNH3.set_ylim(0.0,1.2)
    axWNH3.set_yticks(np.linspace(0.0,1.2,7, endpoint=True))
    axslopeNH3.legend(fontsize=7, loc=2)
    axWNH3.legend(fontsize=7, loc=1)
    axslopeNH3.set_title("Transmission and Equivalent Width for *Pure* NH3")

    figslopeNH3.subplots_adjust(left=0.10, right=0.90, top=0.9, bottom=0.14)
    figslopeNH3.savefig('c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/NH3_Transmission_and_EW.png',dpi=320)
    
    ###########################################################################
    # Compute and plot transmission of the 620nm filter and the equivalent
    # width of the band as a function of the column abundance of CH4.
    ###########################################################################
    filterdata['620']['FiltTrans']=np.loadtxt(path+filterdata['620']['transfile'],usecols=range(2))

    F620_wvs=[615.,625.]
    F620_idx=[np.argmin(abs(F620_wvs[0]-filterdata['620']['FiltTrans'][:,0])),
              np.argmin(abs(F620_wvs[1]-filterdata['620']['FiltTrans'][:,0]))]
    filterdata['620MUSE']['FiltTrans']=copy.deepcopy(filterdata['647']['FiltTrans'])
    filterdata['620MUSE']['FiltTrans'][:,1]=np.zeros((len(filterdata['620MUSE']['FiltTrans'][:,1])))
    filterdata['620MUSE']['FiltTrans'][F620_idx[0]:F620_idx[1],1]=1.

    s_CH4_Arr=np.arange(0.0,0.500,0.02) #column abundance in km-atm
    CH4_trans_Arr=[]
    CH4_trans_Arr_MUSE=[]
    W_Arr=[]
    Continuum=np.copy(Continuum_Albedo)
    Continuum[:,1]=np.ones(Continuum_Albedo.shape[0])
    CH4_trans=copy.deepcopy(Continuum_Albedo)
    CH4_trans_MUSE=copy.deepcopy(Continuum_Albedo)
    print("C",Continuum_Albedo.shape)
    
    filterdata['620MUSE']['ContAlbedoProd']=GSU.SpectrumMath(filterdata['620MUSE']['FiltTrans'],Continuum_Albedo,"Multiply")
    filterdata['620MUSE']['AbsrProd']=GSU.SpectrumMath(filterdata['620MUSE']['FiltTrans'],Albedo,"Multiply")
    filterdata['620MUSE']['ContAlbedo_Int'],filterdata['620MUSE']['AbsrAlbedo_Int'],filterdata['620MUSE']['TransAlbedoInt']= \
        NFL.cont_absorption_calcs(filterdata['620MUSE']['ContAlbedoProd'],filterdata['620MUSE']['AbsrProd'], \
                                  float('620')-filterdata['620MUSE']['filtwdth'],\
                                  float('620')+filterdata['620MUSE']['filtwdth'], \
                                      filterdata['620MUSE']['filtname'])

    figslopeCH4,axslopeCH4=pl.subplots(1,1,figsize=(6.0,3.5), dpi=150, facecolor="white")
    for s in s_CH4_Arr:
        print("############ s=",s)
        temp=1.0*(np.exp(-s*CH4[:,1]))
        
        print('t',temp.shape)
        CH4_trans[:,1]=temp
        print(CH4_trans.shape)
        filterdata['620']['ContProd']=GSU.SpectrumMath(filterdata['620']['FiltTrans'],Continuum,"Multiply")
        filterdata['620']['AbsrProd']=GSU.SpectrumMath(filterdata['620']['FiltTrans'],CH4_trans,"Multiply")
        
        filterdata['620']['Cont_Int'],filterdata['620']['Absr_Int'],filterdata['620']['TransInt']= \
            NFL.cont_absorption_calcs(filterdata['620']['ContProd'],filterdata['620']['AbsrProd'], \
                                      float('620')-filterdata['620']['filtwdth'],\
                                      float('620')+filterdata['620']['filtwdth'], \
                                          filterdata['620']['filtname'])
                        
        filterdata['620']['Tau_Albedo']=-np.log(filterdata['620']['TransInt'])
        
        CH4_trans_Arr.append(filterdata['620']['TransInt'])

        ind1,ind2=np.argmin(abs(CH4_trans[:,0]-600.)),np.argmin(abs(CH4_trans[:,0]-640.))
        W=np.sum(1.0-CH4_trans[ind1:ind2,1])*0.5 #for 0.5 nm bins
        W_Arr.append(W)
        print("#### W=",W)

        ########
        
        filterdata['620MUSE']['ContProd']=GSU.SpectrumMath(filterdata['620MUSE']['FiltTrans'],Continuum,"Multiply")
        filterdata['620MUSE']['AbsrProd']=GSU.SpectrumMath(filterdata['620MUSE']['FiltTrans'],CH4_trans,"Multiply")
        
        filterdata['620MUSE']['Cont_Int'],filterdata['620MUSE']['Absr_Int'],filterdata['620MUSE']['TransInt']= \
            NFL.cont_absorption_calcs(filterdata['620MUSE']['ContProd'],filterdata['620MUSE']['AbsrProd'], \
                                      float('620')-filterdata['620MUSE']['filtwdth'],\
                                      float('620')+filterdata['620MUSE']['filtwdth'], \
                                          filterdata['620MUSE']['filtname'])
                        
        filterdata['620MUSE']['Tau_Albedo']=-np.log(filterdata['620MUSE']['TransInt'])
        
        CH4_trans_Arr_MUSE.append(filterdata['620MUSE']['TransInt'])

    
    R_trans=np.corrcoef(s_CH4_Arr,CH4_trans_Arr)[0,1]
    R_W=np.corrcoef(s_CH4_Arr,W_Arr)[0,1]
    CH4fit=np.polyfit(CH4_trans_Arr,W_Arr,1)
    CH4fitMUSE=np.polyfit(CH4_trans_Arr_MUSE,W_Arr,1)

    axslopeCH4.plot(s_CH4_Arr,CH4_trans_Arr,color='C0',label='Transmission, R='+str(R_trans)[:5])
    axslopeCH4.plot(s_CH4_Arr,CH4_trans_Arr_MUSE,color='C2',label='Transmission MUSE, R='+str(R_trans)[:5])
    axWCH4=axslopeCH4.twinx()
    axWCH4.plot(s_CH4_Arr,W_Arr,color='C1',label='Eq. Width (nm), R='+str(R_W)[:4])
    axslopeCH4.grid(linewidth=0.2)
    axslopeCH4.set_xlabel("CH4 Column Abundance (km-atm)")
    axslopeCH4.set_ylabel("Transmission")
    axWCH4.set_ylabel("Equivalent Width (nm)")
    axslopeCH4.tick_params(axis='both', which='major', labelsize=8)
    axWCH4.tick_params(axis='both', which='major', labelsize=8)
    axslopeCH4.set_xlim(0.0,0.5)
    axslopeCH4.set_xticks(np.linspace(0.0,0.5,6, endpoint=True))
    axslopeCH4.set_ylim(0.9,1.1)
    axslopeCH4.set_yticks(np.linspace(0.8,1.1,7, endpoint=True))
    axWCH4.set_ylim(0.0,1.2)
    axWCH4.set_yticks(np.linspace(0.0,3.0,7, endpoint=True))
    axslopeCH4.legend(fontsize=7, loc=2)
    axWCH4.legend(fontsize=7, loc=1)
    axslopeCH4.set_title("Transmission and Equivalent Width for *Pure* CH4")

    figslopeCH4.subplots_adjust(left=0.10, right=0.90, top=0.9, bottom=0.14)
    figslopeCH4.savefig('c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/CH4_Transmission_and_EW.png',dpi=320)    
    
    print("filterdata['620MUSE']['TransAlbedoInt']=",filterdata['620MUSE']['TransAlbedoInt'])
    print("filterdata['647MUSE']['TransAlbedoInt']=",filterdata['647MUSE']['TransAlbedoInt'])
    print()
    print("NH3fit=",NH3fit)
    print("NH3fitMUSE=",NH3fitMUSE)
    print("CH4fit=",CH4fit)
    print("CH4fitMUSE=",CH4fitMUSE)
    print()
    
def MoonAlbedos(Moon):
    import numpy as np
    import GeneralSpecUtils_P3 as GSU
    
    #Io_leading1980 = np.fromfile(file="c:/Astronomy/Projects/Planets/JovianMoons/References/io.leading_no_header.txt", dtype=float, count=-1, sep=" ")    
    #Io_leading1980=np.reshape(Io_leading1980,[int(Io_leading1980.size/3),3])
    moonspath="c:/Astronomy/Projects/Planets/JovianMoons/References/"
    moonfiles={'Io':'io.trailing_no_header.txt','Europa':'europa_no_header.txt',
               'Ganymede':'ganymede_no_header.txt','Callisto':'callisto_no_header.txt'}
    Moon1980 = np.fromfile(file=moonspath+moonfiles[Moon], dtype=float, count=-1, sep=" ")    
    Moon1980=np.reshape(Moon1980,[int(Moon1980.size/3),3])

    WaveGrid,SignalonGrid=GSU.uniform_wave_grid(Moon1980[:,0]*1000.,Moon1980[:,1],
                                            Extend=False,Fine=False)
    MoonGrid=np.zeros((WaveGrid.size,2))
    MoonGrid[:,0]=WaveGrid
    MoonGrid[:,1]=SignalonGrid

    return MoonGrid

def vert_profile_quad_plot(SupTitle="SupTitle"):
    # Simply sets up quad plot for Transmission and Contribution functions
    import matplotlib.pyplot as pl
    fig_trans,axs_trans=pl.subplots(2,2,figsize=(6.0,6.0), dpi=150, facecolor="white",
                                  sharex=True,sharey=True)
    fig_trans.suptitle(SupTitle)
    for i in range(0,2):
        for j in range(0,2):
            axs_trans[i,j].set_ylim(100.,5000000.)
            axs_trans[i,j].set_yscale('log')
            axs_trans[i,j].set_ylim(axs_trans[i,j].get_ylim()[::-1]) #reverse y-axis
            axs_trans[i,j].set_xlim(0.,1.) #reverse y-axis
            axs_trans[i,j].set_xscale('linear')
            axs_trans[i,j].grid(which='both')
            if i==1:
                axs_trans[i,j].set_xlabel("Transmission")
            if j==0:
                axs_trans[i,j].set_ylabel("Pressue (Pa)")
    axs_trans[0,0].set_title("Gas+Rayleigh")
    axs_trans[0,1].set_title("CH4")
    axs_trans[1,0].set_title("Rayleigh")
    axs_trans[1,1].set_title("NH3")
    return fig_trans,axs_trans

def compute_vertical_transmission_profiles(Jupiterdata,FilterList,CH4,NH3,Pres,
                                           fout_sfx="Test"):
    """
    PURPOSE: Computes vertical profiles of the transmission and weighting 
             functions for gas absorption (methane and ammonia) along with 
             Rayleigh scattering then plots the results.

    Parameters
    ----------
    Jupiterdata : dict
        This dictionary will have the basic filter information...TBD
    FilterList : TYPE
        DESCRIPTION.
    CH4 : TYPE
        DESCRIPTION.
    NH3 : TYPE
        DESCRIPTION.
    P : TYPE
        DESCRIPTION.
    fnout : TYPE, optional
        DESCRIPTION. The default is 'filtereffectivedata.csv'.

    Returns
    -------
    Jupiterdata : TYPE
        DESCRIPTION.

    """
    import numpy as np
    import NH3_Filter_Library_P3 as NFL
    import sys
    sys.path.append('./Services')
    import get_keff
    
    fig_trans,axs_trans=NFL.vert_profile_quad_plot(SupTitle="Two-way Transmission for Absorbers and Scatterers")
    fig_Keff,axs_Keff=NFL.vert_profile_quad_plot(SupTitle="Weighting for Absorbers and Scatterers")
    
    # Write file header for filter data csv file
    #pth='c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/'
    #filtereffectivedata = open(pth+fnout, 'w')
    #tmp="Wavelength (nm),Filter Name,k_eff (NH3),l_eff (NH3),k_eff (CH4),l_eff (CH4),Trans,Tau,NH3 (m-atm),CH4 (m-atm)\n"
    #filtereffectivedata.write(tmp)
    # Compute effective absorption coefficients for each filter
    #Jupiterdata=get_keff.get_keff(Jupiterdata,FilterList,CH4,NH3)
    Jupiterdata=get_keff.get_keff(Jupiterdata,FilterList,CH4,NH3)

    for filtr in FilterList:
        
        #Compute optical depth (tau) independently for methane, ammonia, Rayleigh 
        #  scattering, and total gas (CH4+NH3) for each filter as a function of
        #  pressure level.
            
        Jupiterdata[filtr]['tau_CH4']= \
            NFL.tau_gas_versus_P(Pres,Jupiterdata[filtr]['keff_CH4'],
                                 Jupiterdata[filtr]['filtname'],axs_Keff,gas='CH4')
        Jupiterdata[filtr]['tau_NH3']= \
            NFL.tau_gas_versus_P(Pres,Jupiterdata[filtr]['keff_NH3'],
                                 Jupiterdata[filtr]['filtname'],axs_Keff,gas='NH3')
        Jupiterdata[filtr]['tau_R']= \
            NFL.tau_rayleigh_versus_P(Pres,Jupiterdata[filtr]['leff_CH4'],
                                      Jupiterdata[filtr]['filtname'],axs_Keff)
        tau_gas=Jupiterdata[filtr]['tau_CH4']+Jupiterdata[filtr]['tau_NH3']
        
        #Compute and plot transmission and weighting functions
        
        tmp=NFL.Compute_Transmission(Pres,Jupiterdata[filtr]['tau_R'],
                                     tau_gas,Jupiterdata[filtr]['filtname'],
                                     axs_trans[0,0],axs_Keff[0,0])
        tmp=NFL.Compute_Transmission(Pres,Jupiterdata[filtr]['tau_R']*0.0,
                                     Jupiterdata[filtr]['tau_CH4'],
                                     Jupiterdata[filtr]['filtname']+' CH4',
                                     axs_trans[0,1],axs_Keff[0,1])
        tmp=NFL.Compute_Transmission(Pres,Jupiterdata[filtr]['tau_R'],
                                     tau_gas*0.0,Jupiterdata[filtr]['filtname']+' Ray',
                                     axs_trans[1,0],axs_Keff[1,0])
        tmp=NFL.Compute_Transmission(Pres,Jupiterdata[filtr]['tau_R']*0.0,
                                     Jupiterdata[filtr]['tau_NH3'],
                                     Jupiterdata[filtr]['filtname']+' NH3',
                                     axs_trans[1,1],axs_Keff[1,1])
               
        """Jupiterdata[filtr]['NH3ColDens']=1000.*Jupiterdata[filtr]['Tau_Albedo']/Jupiterdata[filtr]['keff_NH3']
        Jupiterdata[filtr]['CH4ColDens']=1000.*Jupiterdata[filtr]['Tau_Albedo']/Jupiterdata[filtr]['keff_CH4']
        
        #Create the CSV record string and append it to the output file
        
        tmp=filtr+","+Jupiterdata[filtr]['filtname']+","+str(Jupiterdata[filtr]['keff_NH3'])+","\
                +str(Jupiterdata[filtr]['leff_NH3'])+","+str(Jupiterdata[filtr]['keff_CH4'])+","\
                +str(Jupiterdata[filtr]['leff_CH4'])+","+str(Jupiterdata[filtr]['TransInt'])+","\
                +str(Jupiterdata[filtr]['Tau_Albedo'])+","+str(Jupiterdata[filtr]['NH3ColDens'])+","\
                +str(Jupiterdata[filtr]['CH4ColDens'])+"\n"
        filtereffectivedata.write(tmp)"""
    
    #filtereffectivedata.close()
    
    axs_trans[0,0].legend(loc=1,ncol=3, borderaxespad=0.,prop={'size':6})
    fig_trans.subplots_adjust(left=0.12, right=0.96, top=0.90, bottom=0.09)
            
    axs_Keff[0,0].legend(loc=1,ncol=3, borderaxespad=0.,prop={'size':6})
    fig_Keff.subplots_adjust(left=0.12, right=0.96, top=0.90, bottom=0.09)
    pthout='c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Contribution Functions/'
    fig_trans.savefig(pthout+'VertTransProf'+fout_sfx+'.png',dpi=320,bbox_inches = 'tight')
    fig_Keff.savefig(pthout+'VertContFunct'+fout_sfx+'.png',dpi=320,bbox_inches = 'tight')
    
    return Jupiterdata

def compute_filter_spectrum(x0,x1,xtks,y0,y1,ytks,FilterList,
                                         Albedo,Continuum_Albedo,ContinuumModel,
                                         Telescope='SCT',scale=True):
    """
    PURPOSE: PLOT FILTER TRANSMISSIONS CONVOLVED WITH DISK-INTEGRATED ALBEDO 
             AND CONTINUUM
    Calls:
    Called by:
        -> JupiterFilterPerformance
        -> AmmoniaTest_P3.py       

    Parameters
    ----------
    x0 : Float
        DESCRIPTION.
    x1 : Float
        DESCRIPTION.
    xtks : Int
        DESCRIPTION.
    y0 : Float
        DESCRIPTION.
    y1 : Float
        DESCRIPTION.
    ytks : Int
        DESCRIPTION.
    filterwavelength : TYPE
        DESCRIPTION.
    Albedo : TYPE
        DESCRIPTION.
    Continuum_Albedo : TYPE
        DESCRIPTION.
    ContinuumModel : TYPE
        DESCRIPTION.
    Telescope : Str, optional
        DESCRIPTION. The default is 'SCT'.
    scale : Bool, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    Jupiterdata : Dict
        DESCRIPTION.
    axs1 : Plot axis object
        DESCRIPTION.

    """
    import numpy as np
    import NH3_Filter_Library_P3 as NFL
    import matplotlib.pyplot as pl
    import GeneralSpecUtils_P3 as GSU
    import sys
    sys.path.append('./Services')
    import get_filter_base_dict

    # GET JUPITER BASE DICTIONARY AND SET UP PLOT
    path,Jupiterdata=get_filter_base_dict.get_filter_base_dict(Telescope)

    fig1,axs1=pl.subplots(1,1,figsize=(6.0,3.5), dpi=150, facecolor="white",
                          sharex=True)
    axs1.set_xlim(x0,x1)
    axs1.set_xticks(np.linspace(x0,x1,xtks, endpoint=True))
    axs1.set_ylim(y0,y1)
    axs1.set_yticks(np.linspace(y0,y1,ytks, endpoint=True))
    axs1.grid(linewidth=0.2)
    axs1.tick_params(axis='both', which='major', labelsize=8)
    axs1.set_ylabel("Albedo x Transmission",color="black")
    axs1.set_xlabel("Wavelength (nm)")

    # PLOT ALBEDO AND CONTINUUM MODEL
    axs1.plot(Albedo[:,0],Albedo[:,1],label='Albedo')
    axs1.plot(Continuum_Albedo[:,0],Continuum_Albedo[:,1],label='Continuum')
    if Telescope=='VLT':
        MPath='C:/Astronomy/Projects/SAS 2021 Ammonia/VLT MUSE/'
        MUSEAlbedoSpec=np.loadtxt(MPath+'MUSE_AlbedoSpec.txt',usecols=range(2))
        axs1.plot(MUSEAlbedoSpec[:,0],MUSEAlbedoSpec[:,1]*764.,label="MUSE Albedo")
        #MUSESmoothedSpec=np.loadtxt(MPath+'MUSE_SmoothedSpec.txt',usecols=range(2))
        #axs1.plot(MUSESmoothedSpec[:,0],MUSESmoothedSpec[:,1]*764.,label="MUSE Smoothed")

    # PLOT INSTRUMENT+FILTER TRANSMISSION PROFILES WITH CONTINUUM AND ALBEDO
    for filtr in FilterList:
        # COMPUTE PRODUCT OF INSTRUMENT+FILTER TRANSMISSION WITH CONTINUUM AND 
        # ALBEDO 
        Jupiterdata[filtr]['ContProd']=\
            GSU.SpectrumMath(Jupiterdata[filtr]['FiltTrans'],
                             Continuum_Albedo,"Multiply")
        Jupiterdata[filtr]['AbsrProd']=\
            GSU.SpectrumMath(Jupiterdata[filtr]['FiltTrans'],
                             Albedo,"Multiply")
    
        # COMPUTE INTEGRAL TRANSMITTED SIGNAL FOR EACH FILTER WITH CONTINUUM 
        # AND ALBEDO INPUT
        Jupiterdata[filtr]['Cont_Int'],Jupiterdata[filtr]['Absr_Int'],Jupiterdata[filtr]['TransInt']= \
            NFL.cont_absorption_calcs(Jupiterdata[filtr]['ContProd'],Jupiterdata[filtr]['AbsrProd'], \
                                      float(filtr)-Jupiterdata[filtr]['halfwdth'],\
                                      float(filtr)+Jupiterdata[filtr]['halfwdth'], \
                                          Jupiterdata[filtr]['filtname'])
        # COMPUTE OPTICAL DEPTH IN A FILTER FROM THE TRANSMISSION
        Jupiterdata[filtr]['Tau_Albedo']=-np.log(Jupiterdata[filtr]['TransInt'])
        
        #! WOULD BE A GREAT PLACE TO ADD COMPUTATION OF EQUIVALENT WIDTH FOR
        #! CH4 AND NH3 BANDS?
    
        zeros=np.zeros(Jupiterdata[filtr]['FiltTrans'].shape[0])
        
        # SCALE FILTER PEAK TO CONTINUUM OR NOT
        if scale:
            wvidx=np.argmin(np.abs(Jupiterdata[filtr]['ContProd'][:,0]-np.int(filtr)))
            maxcont=np.nanmax(Jupiterdata[filtr]['ContProd'][wvidx-20:wvidx+20,1])
            temp=np.nan_to_num(Jupiterdata[filtr]['ContProd'][:,:], copy=True, nan=0.0)
            idx=int(np.argmin(np.abs(temp[:,1]-maxcont)))
            normscale=Continuum_Albedo[idx,1]/maxcont
        else:
            normscale=1.0
            
        # ONLY WRITE LEGEND ONE TIME
        if str(filtr)[0:3]=='620':
            axs1.fill_between(Jupiterdata[filtr]['AbsrProd'][:,0],zeros, Jupiterdata[filtr]['AbsrProd'][:,1]*normscale,label='Jupiter Signal',color='C0',alpha=0.2)
            axs1.fill_between(Jupiterdata[filtr]['AbsrProd'][:,0],Jupiterdata[filtr]['AbsrProd'][:,1]*normscale, Jupiterdata[filtr]['ContProd'][:,1]*normscale,label='Gas-free Signal',color='C1',alpha=0.2)
        else:
            axs1.fill_between(Jupiterdata[filtr]['AbsrProd'][:,0],zeros, Jupiterdata[filtr]['AbsrProd'][:,1]*normscale,color='C0',alpha=0.2)
            axs1.fill_between(Jupiterdata[filtr]['AbsrProd'][:,0],Jupiterdata[filtr]['AbsrProd'][:,1]*normscale, Jupiterdata[filtr]['ContProd'][:,1]*normscale,color='C1',alpha=0.2)
    
    axs1.legend(fontsize=8,loc=2,ncol=2)
    
    # LABELING FOR PUBLICATION FIGURE
    if Telescope=='SCT':
        Label="a) "
    if Telescope=='VLT':
        Label="b) "
    axs1.set_title(Label+Telescope+" Filter Performance (Continuum "+ContinuumModel+")")
    
    fig1.subplots_adjust(left=0.10, right=0.90, top=0.9, bottom=0.14)
    pth='c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Filter Plots/'
    fig1.savefig(pth+'JupiterFilterPerformance-'+Telescope+'-'\
                 +ContinuumModel+'.png',dpi=320)
    
    return Jupiterdata,axs1
