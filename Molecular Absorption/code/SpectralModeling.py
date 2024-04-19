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
    import get_filter_base_dict as GFBD

    ###########################################################################
    # Get albedo, continua, and cross section data and plot
    ###########################################################################

    x0,x1,xtks=600.,680.,5
    y0,y1,ytks=0.4,0.6,5
    Albedo,Continua,CH4,NH3,NH3_LO1980= \
        gACC.get_albedo_continua_crossections(x0,x1,xtks,y0,y1,ytks,
                                              Crossect=False)
    #NH3=NH3_LO1980
    
    ###########################################################################
    # Get albedo, continua, and cross section data and plot
    ###########################################################################

    figtest,axtest=pl.subplots(1,1,figsize=(6.0,3.5), dpi=150, facecolor="white")

    axtest.plot(Albedo[:,0],Albedo[:,1],label='Jupiter Albedo (Karkoschka, 1994)',linewidth=1.0,color='C0')
    axtest.set_xlim(x0,x1)
    axtest.set_xticks(np.linspace(x0,x1,xtks, endpoint=True))
    axtest.set_ylim(y0,y1)
    axtest.set_yticks(np.linspace(y0,y1,ytks, endpoint=True))
    #print(CH4.shape)
    CH4_trans=refl*(np.exp(-s_CH4*CH4[:,1]))
    NH3_trans=refl*(np.exp(-s_NH3*NH3[:,1]))
    gas_trans=refl*(np.exp(-s_CH4*CH4[:,1]-s_NH3*NH3[:,1]))
    #print(CH4_trans.shape)
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

    FilterList=['620','647']
    SCTpath,SCTdict=GFBD.get_filter_base_dict('SCT',FilterList=FilterList,
                             Inst=True)
    VLTpath,VLTdict=GFBD.get_filter_base_dict('VLT',FilterList=FilterList,
                             Inst=True)
    Model='Piecewise1'

    for tele in ['SCT','VLT']:
        #print(tele)
        if tele=='SCT':
            filterdata=SCTdict
            path=SCTpath
            i=0
        elif tele=='VLT':
            filterdata=VLTdict
            path=VLTpath
            i=0
        
        figsct,axssct=pl.subplots(1,1,figsize=(6.0,6.0), dpi=150, facecolor="white")
        figsct.suptitle(tele)
        axssct.plot(filterdata['647']['FiltTrans'][:,0],filterdata['647']['FiltTrans'][:,1])
        axssct.set_xlim(600.,680.)
            
        fig_tau,axs_tau=NFL.Tau_EW_quad_plot(SupTitle="Transmission and Equivalent Width for "+tele)
        
        for filtr in FilterList:
            #print(filtr)
            if filtr=='647':
                wvs=[642.,652.] #Original
                #F647_wvs=[636.,656.]
                Band_idx=[np.argmin(abs(wvs[0]-filterdata[filtr]['FiltTrans'][:,0])),
                          np.argmin(abs(wvs[1]-filterdata[filtr]['FiltTrans'][:,0]))]
                s_Arr=np.arange(0.0,0.031,0.002) #column abundance in km-atm
                Cont_Wave=[637.,657.] #Original
                #NH3_Cont_Wave=[636.,656.]
                ind1,ind2=np.argmin(abs(NH3[:,0]-Cont_Wave[0])),np.argmin(abs(NH3[:,0]-Cont_Wave[1]))
                j=0
                gas='Ammonia'
                csect=NH3
            elif filtr=='620':
                wvs=[615.,625.] #Original
                #F647_wvs=[636.,656.]
                Band_idx=[np.argmin(abs(wvs[0]-filterdata[filtr]['FiltTrans'][:,0])),
                          np.argmin(abs(wvs[1]-filterdata[filtr]['FiltTrans'][:,0]))]
                s_Arr=np.arange(0.0,0.500,0.02) #column abundance in km-atm
                Cont_Wave=[610.,630.] #Original
                ind1,ind2=np.argmin(abs(NH3[:,0]-Cont_Wave[0])),np.argmin(abs(NH3[:,0]-Cont_Wave[1]))
                j=1
                gas='Methane'
                csect=CH4

            Trans_Arr=[]
            W_Arr=[]
            """
            filterdata[filtr]['FiltTrans']=copy.deepcopy(filterdata[filtr]['FiltTrans'])
            filterdata[filtr]['FiltTrans'][:,1]=np.zeros((len(filterdata[filtr]['FiltTrans'][:,1])))
            filterdata[filtr]['FiltTrans'][Band_idx[0]:Band_idx[1],1]=1.
            """
            Continuum_Albedo=np.zeros((Continua[Model]['WaveGrid'].size,2))
            Continuum_Albedo[:,0]=Continua[Model]['WaveGrid']
            Continuum_Albedo[:,1]=Continua[Model]['Albedo']
            
            Continuum=copy.deepcopy(Continuum_Albedo)
            Continuum[:,1]=np.ones(Continuum_Albedo.shape[0])
            Trans=copy.deepcopy(Continuum_Albedo)
            #print("C",Continuum_Albedo.shape)
            
            filterdata[filtr]['ContAlbedoProd']=GSU.SpectrumMath(filterdata[filtr]['FiltTrans'],Continuum_Albedo,"Multiply")
            filterdata[filtr]['AbsrProd']=GSU.SpectrumMath(filterdata[filtr]['FiltTrans'],Albedo,"Multiply")
            filterdata[filtr]['ContAlbedo_Int'],filterdata[filtr]['AbsrAlbedo_Int'],filterdata[filtr]['TransAlbedoInt']= \
                NFL.cont_absorption_calcs(filterdata[filtr]['ContAlbedoProd'],filterdata[filtr]['AbsrProd'], \
                                          float(filtr)-filterdata[filtr]['halfwdth'],\
                                          float(filtr)+filterdata[filtr]['halfwdth'], \
                                              filterdata[filtr]['filtname'],prn=False)
                
            #figslopeNH3,axslopeNH3=pl.subplots(1,1,figsize=(6.0,3.5), dpi=150, facecolor="white")
            for s in s_Arr:
                #print("############ s=",s)
                temp=1.0*(np.exp(-s*csect[:,1]))
                
                #print('t',temp.shape)
                Trans[:,1]=temp #NH3 transmission as a function of wavelength
                #print(NH3_trans.shape)
                filterdata[filtr]['ContProd']=GSU.SpectrumMath(filterdata[filtr]['FiltTrans'],Continuum,"Multiply")
                filterdata[filtr]['AbsrProd']=GSU.SpectrumMath(filterdata[filtr]['FiltTrans'],Trans,"Multiply")
                
                filterdata[filtr]['Cont_Int'],filterdata[filtr]['Absr_Int'],filterdata[filtr]['TransInt']= \
                    NFL.cont_absorption_calcs(filterdata[filtr]['ContProd'],filterdata[filtr]['AbsrProd'], \
                                              float(filtr)-filterdata[filtr]['halfwdth'],\
                                              float(filtr)+filterdata[filtr]['halfwdth'], \
                                                  filterdata[filtr]['filtname'],prn=False)
                                
                filterdata[filtr]['Tau_Albedo']=-np.log(filterdata[filtr]['TransInt'])
                
                Trans_Arr.append(filterdata[filtr]['TransInt'])
        
                W=np.sum(1.0-Trans[ind1:ind2,1])*0.5 #for 0.5 nm bins
                W_Arr.append(W)
                #print("#### W=",W)
                            
            TauSCTNH3=-np.log(Trans_Arr)
            R_W=np.corrcoef(s_Arr,W_Arr)[0,1]
            R_trans=np.corrcoef(s_Arr,TauSCTNH3)[0,1]
            
            NH3fit=np.polyfit(TauSCTNH3,W_Arr,1)
            NH3transfit=np.polyfit(Trans_Arr,W_Arr,1)
            
            #x=copy.deepcopy(TauSCTNH3)
            #x = x[:,np.newaxis]
            #a, _, _, _ = np.linalg.lstsq(x, W_Arr)
            #print()
            #print(gas+" "+tele+" fit=",NH3fit)
            print()
            print(gas+" "+tele+" Trans fit=",NH3transfit)

            #print("************** i,j=",i,j)
            axs_tau[0,j].plot(s_Arr,TauSCTNH3,color='C0',label='Tau, R='+str(R_trans)[:5])
            axs_tau[i,j].grid(linewidth=0.2)
            #axs_tau[i,j].set_xlabel("Column Abundance (km-atm)")
            #axs_tau[0,0].set_ylabel("Opacity")
            axs_tau[0,j].set_xlim(0.0,0.03)
            axs_tau[0,j].set_xticks(np.linspace(0.0,0.03,7, endpoint=True))
            axs_tau[0,0].set_ylim(0.0,0.3)
            axs_tau[0,0].set_yticks(np.linspace(0.0,0.3,4, endpoint=True))
            axs_tau[0,1].set_xlim(0.0,0.5)
            axs_tau[0,1].set_xticks(np.linspace(0.0,0.5,6, endpoint=True))
            axs_tau[0,1].set_ylim(0.0,0.3)
            axs_tau[0,1].set_yticks(np.linspace(0.0,0.3,4, endpoint=True))

            axWNH3=axs_tau[i,j].twinx()
            axWNH3.plot(s_Arr,W_Arr,color='C1',label='Eq. Width (nm), R=')#+str(R_W)[:4])
            axWNH3.set_ylim(0.0,3.0)
            axWNH3.set_yticks(np.linspace(0.0,3,6, endpoint=True))
            
            ###################################################################
            #Plot EW vs Tau and *fit* of EW vs Tau
            ###################################################################
            axs_tau[1,j].plot(TauSCTNH3,W_Arr,color='C0',label='Model')#+str(R_trans)[:5])
            fit=NH3fit[0]*np.array(TauSCTNH3)+NH3fit[1]
            #FitLin=a*x

            axs_tau[1,j].plot(TauSCTNH3,fit,color='C1',label='Linear Fit')#+str(R_trans)[:5])
            axErr=axs_tau[1,j].twinx()
            axErr.plot(TauSCTNH3,W_Arr/fit,color='C2',label='Ratio')
            axErr.set_ylim(0.9,1.1)
            axErr.set_yticks(np.linspace(0.9,1.1,5, endpoint=True))

            axs_tau[1,0].set_xlim(0.0,0.1)
            axs_tau[1,0].set_xticks(np.linspace(0.0,0.1,5, endpoint=True))
            axs_tau[1,0].set_ylim(0.0,1.2)
            axs_tau[1,0].set_yticks(np.linspace(0.0,1.2,7, endpoint=True))
            
            axs_tau[1,1].set_xlim(0.0,0.25)
            axs_tau[1,1].set_xticks(np.linspace(0.0,0.25,6, endpoint=True))
            axs_tau[1,1].set_ylim(0.0,1.2)
            axs_tau[1,1].set_yticks(np.linspace(0.0,2.5,6, endpoint=True))
            
            axs_tau[i,j].tick_params(axis='both', which='major', labelsize=8)

            if j==0:
                axWNH3.tick_params(labelright=False)    
                axErr.tick_params(labelright=False)    

            if j==1:
                axWNH3.set_ylabel("Equivalent Width (nm)")
                axWNH3.tick_params(axis='both', which='major', labelsize=8)

            axs_tau[0,j].legend(fontsize=7, loc=2)
            axs_tau[1,j].legend(fontsize=7, loc=2)
            axWNH3.legend(fontsize=7, loc=1)
            axErr.legend(fontsize=7, loc=2)
    
        fig_tau.subplots_adjust(left=0.10, right=0.90, top=0.9, bottom=0.14)
        fig_tau.savefig('c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Tau_vs_EW.png',dpi=320)

    