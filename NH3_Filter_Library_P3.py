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

def K_eff(P,FilterTransmission,Abs_Crossection,wv1,wv2,filtername,axis):
    ###########################################################################
    # COMPUTE K_eff AND l_eff FOR AN ABSORBING GAS, CH4 OR NH3, GIVEN THE 
    #   CROSSECTION AND FILTER TRANSMISSION
    #   !!!! AGAIN, NEED TO HAVE TELESCOPE PERFORMANCE HERE ALSO
    #   !!!! Is pressure, P, actually used? Remove from passed parameters.
    ###########################################################################
    import numpy as np
    import GeneralSpecUtils_P3 as GSU
    import copy
    import NH3_Filter_Library_P3 as NFL

    Abs_Product=GSU.SpectrumMath(FilterTransmission,Abs_Crossection,"Multiply")
    lam=copy.deepcopy(FilterTransmission)
    lam[:,1]=FilterTransmission[:,0]
    lamprod=GSU.SpectrumMath(FilterTransmission,lam,"Multiply")
    StartIndex=np.where(Abs_Product[:,0]==wv1)
    EndIndex=np.where(Abs_Product[:,0]==wv2)
    #print(wv1,wv2)
    keff=sum(Abs_Product[StartIndex[0][0]:EndIndex[0][0],1])/\
                 sum(FilterTransmission[StartIndex[0][0]:EndIndex[0][0],1])
    leff=sum(lamprod[StartIndex[0][0]:EndIndex[0][0],1])/\
                 sum(FilterTransmission[StartIndex[0][0]:EndIndex[0][0],1])
    #print("################################")
    #print(filtername,"K_eff=",keff_CH4)
    #print(filtername,"K_eff : %4.3f, leff : %5.2f" % (keff, leff))
    #print("")
    return(keff,leff)

def tau_gas_versus_P(P,Keff,filtername,axis,gas='CH4'):
    ###########################################################################
    # COMPUTE TAU DUE TO GAS ABSORPTION AS A FUNCTION OF PRESSURE 
    #   (WEIGHTING FUNCTION)
    #   !!!! NEED TO ADD NH3 OPTION WITH SOME BASELINE ABUNDANCE OR PROFILE
    #   !!!! CAN COMPARE MENDIKOA TO HILL COMPUTATION
    #   !!!! HAS COMMENTED CODE FOR PLOTTING WHICH COULD BE USEFUL TO REACTIVATE
    #   !!!! SHOULD HAVE GRAVITY AS A FUNCTION OF LATITUDE?
    ###########################################################################
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
def Get_Albedo_and_Absorption(x0,x1,xtks,y0,y1,ytks,Lutz=True,Crossect=True,
                              ContMod='Both'):
    
    import sys
    sys.path.append('c:/Astronomy/Python Play')
    sys.path.append('c:/Astronomy/Python Play/Util_P3')
    sys.path.append('c:/Astronomy/Python Play/SPLibraries_P3')
    sys.path.append('c:/Astronomy/Python Play/SpectroPhotometry/Spectroscopy_P3')
    import matplotlib.pyplot as pl
    import numpy as np
    from scipy import interpolate
    import GeneralSpecUtils_P3 as GSU
    from numpy import genfromtxt
    import NH3_Filter_Library_P3 as NFL
    
    ###############################################################################
    # LOAD JOVIAN DISK-INTEGRATEDALBEDO DATA FROM KARKOSCHKA, 1994 (DATA FROM 1993)
    ###############################################################################
    Jupiter_Karkoschka1993 = np.fromfile(file="c:/Astronomy/Projects/Planets/Saturn/Spectral Data/Karkoschka/1993.tab.txt", dtype=float, count=-1, sep=" ")
    kark1993nrows=int(Jupiter_Karkoschka1993.size/8)
    Jupiter_Karkoschka1993=np.reshape(Jupiter_Karkoschka1993,[kark1993nrows,8])
    
    Albedo_KarkRef1993=np.zeros((kark1993nrows,2))
    Albedo_KarkRef1993[:,0]=Jupiter_Karkoschka1993[:,0]
    Albedo_KarkRef1993[:,1]=Jupiter_Karkoschka1993[:,3] #Albedo
    WaveGrid,SignalonGrid=GSU.uniform_wave_grid(Albedo_KarkRef1993[:,0],Albedo_KarkRef1993[:,1],
                                            Extend=False,Fine=False)
    
    Albedo=np.zeros((WaveGrid.size,2))
    Albedo[:,0]=WaveGrid
    Albedo[:,1]=SignalonGrid
    ###############################################################################
    # COMPUTE SPLINE FIT TO ALBEDO, PICKING POINTS WITH MINIMUM GAS ABSORPTION.
    ###############################################################################
    SplineWV1= np.array([560.0, 580.0, 600.0,610.0, 635.0, 660.0, 675.0, 690., 714.0,
                         745.0, 830.0, 945.0,1050.0])
    SplineWV2 = np.array([560.0, 580.0, 600.0, 632.0, 678.0, 
                         745.0, 830.0, 945.0,1050.0])
    Continua={'Model1':{'WV':SplineWV1},'Model2':{'WV':SplineWV2}}
    Models=['Model1','Model2']
    #SplineWV=SplineWV2
    
    for Model in Models:
        Continua[Model]['Mag']=np.ones(Continua[Model]['WV'].size)
        for i in range(0,Continua[Model]['WV'].size):
            Start=Continua[Model]['WV'][i]-.0000001
            End=Continua[Model]['WV'][i]+.0000001
            SplineWVIndices=np.where((Albedo[:,0] >Start) & \
                 (Albedo[:,0] < End))
            #print("i= ",i,SplineWVIndices)
            Continua[Model]['Mag'][i]=np.log10(Albedo[SplineWVIndices[0],1])
        
        #print(SplineMag)
        Continua[Model]['Mag'][Continua[Model]['Mag'].size-1]=np.log10(0.42)
        #x = np.arange(0, 2*np.pi+np.pi/4, 2*np.pi/8)
        #y = np.sin(x)
        tck = interpolate.splrep(Continua[Model]['WV'], Continua[Model]['Mag'], s=0)
        #xnew = np.arange(0, 2*np.pi, np.pi/50)
        Temp = 10**interpolate.splev(WaveGrid, tck, der=0)
        
        Continua[Model]['WaveGrid']=WaveGrid
        Continua[Model]['Albedo']=Temp
        Continuum_Albedo=np.zeros((WaveGrid.size,2))
        Continuum_Albedo[:,0]=WaveGrid
    if ContMod=='1' or ContMod=='Both':
        Continuum_Albedo[:,1]=Continua['Model1']['Albedo']
    if ContMod=='2' or ContMod=='Both':
        Continuum_Albedo[:,1]=Continua['Model2']['Albedo']

    ###############################################################################
    # LOAD METHANE ABSORPTION DATA FROM KARKOSCHKA, 1994 (DATA FROM 1993)
    ###############################################################################
    CH4_KarkRef1993=np.zeros((kark1993nrows,2))
    CH4_KarkRef1993[:,0]=Jupiter_Karkoschka1993[:,0]
    CH4_KarkRef1993[:,1]=Jupiter_Karkoschka1993[:,2] #CH4 Coef
    WaveGrid,SignalonGrid=GSU.uniform_wave_grid(CH4_KarkRef1993[:,0],CH4_KarkRef1993[:,1],
                                            Extend=False,Fine=False)
    CH4=np.zeros((WaveGrid.size,2))
    CH4[:,0]=WaveGrid
    CH4[:,1]=SignalonGrid
    ###############################################################################
    # LOAD AMMONIA ABSORPTION DATA EITHER FROM LUTZ & OWEN 1980 OR FROM EXOMOL 
    #   (IRWIN, 2022 - PERSONAL COMMUNICATION)
    ###############################################################################
    fn='Lutz&Owen1980_Figure5_AmmoniaCrossSection.csv'
    pth="c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/"
    NH3_Lutz_Owen_1980 = np.array(genfromtxt(pth+fn, delimiter=','))
    
    fn='Exomol_NH3.csv'
    pth="c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/"
    NH3_Exomol = np.array(genfromtxt(pth+fn, delimiter=','))
    #print("***",NH3_Exomol.shape)
    wv,csec=GSU.uniform_wave_grid(NH3_Exomol[:,0]*1000.,NH3_Exomol[:,1],Extend=False,Fine=False)
    NH3_Exomol_regrid=np.transpose(np.array([wv,csec]))
    
    NH3=NH3_Exomol_regrid
    #print("***",NH3_Exomol_regrid.shape,NH3_Exomol_regrid.max(1))
    ###############################################################################
    
    ###### Get reference regional I/F reflectivities from Dahl, 2021?
    #
    #Dahl_NEB="c:/Astronomy/Projects/SAS 2021 Ammonia/Dahl Spectra-NEB.txt"
    #NEB = genfromtxt(Dahl_NEB, delimiter=',')
    #Dahl_SEB="c:/Astronomy/Projects/SAS 2021 Ammonia/Dahl Spectra-SEB.txt"
    #SEB = genfromtxt(Dahl_SEB, delimiter=',')
    #Dahl_EZ="c:/Astronomy/Projects/SAS 2021 Ammonia/Dahl Spectra-EZ.txt"
    #EZ = genfromtxt(Dahl_EZ, delimiter=',')
    
    ###### Plot disk-integrated reference albedo and simulated ammonia-free albedo
    fig_molecules,ax_molecules=pl.subplots(1,1,figsize=(6.0,3.5), dpi=150, facecolor="white",
                          sharex=True)
    
    #Plot Layout Configuration
    #x0,x1,xtks=600.,1000.,9
    #y0,y1,ytks=0.0,0.7,8
    
    # Set x limits
    ax_molecules.set_xlim(x0,x1)
    # Set x ticks
    ax_molecules.set_xticks(np.linspace(x0,x1,xtks, endpoint=True))
    # Set y limits
    ax_molecules.set_ylim(y0,y1)
    ax_molecules.set_yticks(np.linspace(y0,y1,ytks, endpoint=True))
    # Set y ticks
    ax_molecules.grid(linewidth=0.2)
    ax_molecules.tick_params(axis='both', which='major', labelsize=8)
    ax_molecules.set_ylabel("Albedo",color="black")
    
    ax_molecules.plot(Albedo[:,0],Albedo[:,1],label='Jupiter Albedo (Karkoschka, 1994)',linewidth=1.0,color='C0')
    if ContMod=='1' or ContMod=='Both':
        ax_molecules.plot(Continua['Model1']['WaveGrid'],Continua['Model1']['Albedo'],label='Continuum Model 1',
                     linewidth=1,linestyle='--',color='C0')
    if ContMod=='2' or ContMod=='Both':
        ax_molecules.plot(Continua['Model2']['WaveGrid'],Continua['Model2']['Albedo'],label='Continuum Model 2',
                     linewidth=1,linestyle='-.',color='C1')
        
    #ax_molecules.plot(Continuum_Albedo[:,0],Continuum_Albedo[:,1],label='Continuum Albedo',
    #             linewidth=1,linestyle='--',color='C0')
    ax_molecules.set_title("Albedo and Molecular Absorption")
    
    if Crossect:
        axs1b = ax_molecules.twinx()  # instantiate a second axes that shares the same x-axis
        axs1b.ticklabel_format(axis='y')
        axs1b.tick_params(axis='y', which='major', labelsize=8)
        axs1b.set_yscale('log')
        axs1b.set_ylim(1e-4,1e3)
        axs1b.set_ylabel("Absorption Coefficient 1/(km-atm)")#,color="green")
        
        axs1b.plot(CH4_KarkRef1993[:,0],CH4_KarkRef1993[:,1],label='CH4 Abs. Coef. (Karkoschka, 1994) ',linewidth=1.0,color='C2')
        axs1b.plot(NH3[:,0],NH3[:,1],label='NH3 Abs. Coef. (ExoMol) ',linewidth=1.0,color='C3')
        if Lutz:
            axs1b.plot(NH3_Lutz_Owen_1980[:,0],NH3_Lutz_Owen_1980[:,1],label='NH3 Abs. Coef. (Lutz & Owen, 1980) ',linewidth=0.5,color='C3')
        axs1b.legend(fontsize=7, loc=1)

    ax_molecules.set_xlabel("Wavelength (nm)")
    ax_molecules.legend(fontsize=7, loc=2)
    fig_molecules.subplots_adjust(left=0.10, right=0.90, top=0.9, bottom=0.14)
    
    fig_molecules.savefig('c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Albedo_Molecular_Absorption.png',dpi=320)
    
    return(Albedo,Continuum_Albedo,CH4,NH3)
    
    ########## END OF FIRST FUNCTION AND PLOT ##############
    
    
def SpectralModeling(s_NH3=0.018,s_CH4=0.304,refl=0.53):
###############################################################################
# Model the spectrum of Jupiter using CH4 and NH3 column densities (km-atm)
# along with a constant reflectivity for the cloud tops. 
# !!!! Currently does not include Rayleigh scattering.
###############################################################################
    import sys
    sys.path.append('c:/Astronomy/Python Play')
    sys.path.append('c:/Astronomy/Python Play/Util_P3')
    sys.path.append('c:/Astronomy/Python Play/SPLibraries_P3')
    sys.path.append('c:/Astronomy/Python Play/SpectroPhotometry/Spectroscopy_P3')
    import matplotlib.pyplot as pl
    import numpy as np
    import GeneralSpecUtils_P3 as GSU
    import NH3_Filter_Library_P3 as NFL
    path='c:/Astronomy/Projects/Techniques/InstrumentPerformance-P3/Filters/'
    filterdata={'620':{'transfile':'620CH4/620CH4_Transmission.txt',
                       'filtname':'620CH4','filtwdth':10.},
                '647':{'transfile':'647CNT/647CNT_Transmission.txt',
                    'filtname':'647NH3','filtwdth':10.}}

    #x0,x1,xtks=600.,1000.,9
    #y0,y1,ytks=0.0,0.7,8
    x0,x1,xtks=600.,680.,5
    y0,y1,ytks=0.4,0.6,5
    Albedo,Continuum_Albedo,CH4,NH3=NFL.Get_Albedo_and_Absorption(x0,x1,xtks,y0,y1,ytks,
                                                                  ContMod='Both',Crossect=False)

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

    s_NH3_Arr=np.arange(0.0,0.031,0.002) #column abundance in km-atm
    NH3_trans_Arr=[]
    W_Arr=[]
    Continuum=np.copy(Continuum_Albedo)
    Continuum[:,1]=np.ones(Continuum_Albedo.shape[0])
    NH3_trans=np.copy(Continuum_Albedo)
    print("C",Continuum_Albedo.shape)
    figslope,axslope=pl.subplots(1,1,figsize=(6.0,3.5), dpi=150, facecolor="white")
    for s in s_NH3_Arr:
        print("############ s=",s)
        temp=1.0*(np.exp(-s*NH3[:,1]))
        
        print('t',temp.shape)
        NH3_trans[:,1]=temp
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
    
    R_trans=np.corrcoef(s_NH3_Arr,NH3_trans_Arr)[0,1]
    R_W=np.corrcoef(s_NH3_Arr,W_Arr)[0,1]
    NH3fit=np.polyfit(NH3_trans_Arr,W_Arr,1)
    
    axslope.plot(s_NH3_Arr,NH3_trans_Arr,color='C0',label='Transmission, R='+str(R_trans)[:5])
    axW=axslope.twinx()
    axW.plot(s_NH3_Arr,W_Arr,color='C1',label='Eq. Width (nm), R='+str(R_W)[:4])
    axslope.grid(linewidth=0.2)
    axslope.set_xlabel("NH3 Column Abundance (km-atm)")
    axslope.set_ylabel("Transmission")
    axW.set_ylabel("Equivalent Width (nm)")
    axslope.tick_params(axis='both', which='major', labelsize=8)
    axW.tick_params(axis='both', which='major', labelsize=8)
    axslope.set_xlim(0.0,0.03)
    axslope.set_xticks(np.linspace(0.0,0.03,7, endpoint=True))
    axslope.set_ylim(0.9,1.02)
    axslope.set_yticks(np.linspace(0.9,1.02,7, endpoint=True))
    axW.set_ylim(0.0,1.2)
    axW.set_yticks(np.linspace(0.0,1.2,7, endpoint=True))
    axslope.legend(fontsize=7, loc=2)
    axW.legend(fontsize=7, loc=1)
    axslope.set_title("Transmission and Equivalent Width for *Pure* NH3")

    figslope.subplots_adjust(left=0.10, right=0.90, top=0.9, bottom=0.14)
    figslope.savefig('c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/NH3_Transmission_and_EW.png',dpi=320)
    
    ###########################################################################
    # Compute and plot transmission of the 620nm filter and the equivalent
    # width of the band as a function of the column abundance of CH4.
    ###########################################################################
    filterdata['620']['FiltTrans']=np.loadtxt(path+filterdata['620']['transfile'],usecols=range(2))

    s_CH4_Arr=np.arange(0.0,0.500,0.02) #column abundance in km-atm
    CH4_trans_Arr=[]
    W_Arr=[]
    Continuum=np.copy(Continuum_Albedo)
    Continuum[:,1]=np.ones(Continuum_Albedo.shape[0])
    CH4_trans=np.copy(Continuum_Albedo)
    print("C",Continuum_Albedo.shape)
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
    
    R_trans=np.corrcoef(s_CH4_Arr,CH4_trans_Arr)[0,1]
    R_W=np.corrcoef(s_CH4_Arr,W_Arr)[0,1]
    CH4fit=np.polyfit(CH4_trans_Arr,W_Arr,1)

    axslopeCH4.plot(s_CH4_Arr,CH4_trans_Arr,color='C0',label='Transmission, R='+str(R_trans)[:5])
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
    
    print("NH3fit=",NH3fit)
    print("CH4fit=",CH4fit)
    
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