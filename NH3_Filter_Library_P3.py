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

def K_eff(P,FilterTransmission,Abs_Crossection,wv1,wv2,filtername,axis):
    ###########################################################################
    # COMPUTE K_eff AND l_eff FOR AN ABSORBING GAS, CH4 OR NH3, GIVEN THE 
    #   CROSSECTION AND FILTER TRANSMISSION
    #   !!!! AGAIN, NEED TO HAVE TELESCOPE PERFORMANCE HERE ALSO
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
    keff=sum(Abs_Product[StartIndex[0][0]:EndIndex[0][0],1])/\
                 sum(FilterTransmission[StartIndex[0][0]:EndIndex[0][0],1])
    leff=sum(lamprod[StartIndex[0][0]:EndIndex[0][0],1])/\
                 sum(FilterTransmission[StartIndex[0][0]:EndIndex[0][0],1])
    print("################################")
    #print(filtername,"K_eff=",keff_CH4)
    print(filtername,"K_eff : %4.3f, leff : %5.2f" % (keff, leff))
    print("")
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
    
    print("nH2*sigmaH2+nHe*sigmaHe",nH2*sigmaH2+nHe*sigmaHe)
    
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
    
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 15:07:47 2022

Get Albedo and Molecular Data

@author: smhil
"""
def Get_Albedo_and_Absorption(x0,x1,xtks,y0,y1,ytks):
    
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
    SplineWV1= np.array([560.0, 580.0, 600.0, 635.0, 660.0, 675.0, 690., 714.0,
                         745.0, 830.0, 945.0,1050.0])
    SplineWV2 = np.array([560.0, 580.0, 600.0, 677.0, 690., 
                         745.0, 830.0, 945.0,1050.0])
    SplineWV=SplineWV1
    SplineMag=np.ones(SplineWV.size)
    for i in range(0,SplineWV.size):
        Start=SplineWV[i]-.0000001
        End=SplineWV[i]+.0000001
        SplineWVIndices=np.where((Albedo[:,0] >Start) & \
             (Albedo[:,0] < End))
        #print("i= ",i,SplineWVIndices)
        SplineMag[i]=np.log10(Albedo[SplineWVIndices[0],1])
    
    #print(SplineMag)
    SplineMag[SplineMag.size-1]=np.log10(0.42)
    x = np.arange(0, 2*np.pi+np.pi/4, 2*np.pi/8)
    y = np.sin(x)
    tck = interpolate.splrep(SplineWV, SplineMag, s=0)
    #xnew = np.arange(0, 2*np.pi, np.pi/50)
    Temp = 10**interpolate.splev(WaveGrid, tck, der=0)
    
    Continuum_Albedo=np.zeros((WaveGrid.size,2))
    Continuum_Albedo[:,0]=WaveGrid
    Continuum_Albedo[:,1]=Temp
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
    ax_molecules.plot(Continuum_Albedo[:,0],Continuum_Albedo[:,1],label='Continuum Albedo',
                 linewidth=1,linestyle='--',color='C0')
    ax_molecules.set_title("Albedo and Molecular Absorption")
    
    axs1b = ax_molecules.twinx()  # instantiate a second axes that shares the same x-axis
    axs1b.ticklabel_format(axis='y')
    axs1b.tick_params(axis='y', which='major', labelsize=8)
    axs1b.set_yscale('log')
    axs1b.set_ylim(1e-4,1e3)
    axs1b.set_ylabel("Absorption Coefficient 1/(km-atm)")#,color="green")
    
    axs1b.plot(CH4_KarkRef1993[:,0],CH4_KarkRef1993[:,1],label='CH4 Abs. Coef. (Karkoschka, 1994) ',linewidth=1.0,color='C2')
    axs1b.plot(NH3[:,0],NH3[:,1],label='NH3 Abs. Coef. (ExoMol) ',linewidth=1.0,color='C3')
    axs1b.plot(NH3_Lutz_Owen_1980[:,0],NH3_Lutz_Owen_1980[:,1],label='NH3 Abs. Coef. (Lutz & Owen, 1980) ',linewidth=0.5,color='C3')
    
    #axs1b.plot(CH4[:,0],CH4[:,1],label='CH4 Abs. Coef. ',linewidth=1,color='k')
    ax_molecules.set_xlabel("Wavelength (nm)")
    ax_molecules.legend(fontsize=7, loc=2)
    axs1b.legend(fontsize=7, loc=1)
    fig_molecules.subplots_adjust(left=0.10, right=0.90, top=0.9, bottom=0.14)
    
    fig_molecules.savefig('c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Albedo_Molecular_Absorption.png',dpi=320)
    
    return(Albedo,Continuum_Albedo,CH4,NH3)
    
    ########## END OF FIRST FUNCTION AND PLOT ##############
