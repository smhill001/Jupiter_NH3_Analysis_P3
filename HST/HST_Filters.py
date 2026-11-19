import os
os.chdir("/mnt/data/git_repos/Jupiter_NH3_Analysis_P3/")

def HST_Filters(filters):
    import sys
    #sys.path.append('C:/Astronomy/Python Play/SpectroPhotometry/Spectroscopy_P3/')
    sys.path.append('./Molecular Absorption/code/')
    sys.path.append('/mnt/data/git_repos/Visualization-and-Analysis/')
    from config_VA import config_VA, Host_path, Pipeline_path
    import socket
    hostname = socket.gethostname()
    
    import numpy as np
    import matplotlib.pyplot as pl
    import GeneralSpecUtils_P3 as GSU
    import NH3_Filter_Library_P3 as NFL
    import get_ammonia_data as get_NH3
    import get_karkoschka_data as get_CH4
    
    
    #path="C:/Astronomy/Projects/SAS 2021 Ammonia/Data/HST GO 18055/Filter Throughput/"
    path=Host_path[hostname]+"/HST GO 18055/Filter Throughput/"
    
    filterfiles={619:'wfc3_uvis1_fq619n.txt',
                 631:'wfc3_uvis1_f631n.txt',
                 645:'wfc3_uvis1_f645n.txt'}
    
    
    NH3=get_NH3.get_ammonia_data(projpath=Pipeline_path[hostname],Source="Irwin")
    CH4=get_CH4.get_karkoschka_data(projpath=Pipeline_path[hostname],Type='CH4')
       
    fig1,axs1=pl.subplots(1,figsize=(6,4), dpi=150, facecolor="white")
    fig1.suptitle("Transmission")
    
    for filt in filters:
    
        fn=filterfiles[filt]
        data = np.loadtxt(path+fn, skiprows=1)
        WG,SG=GSU.uniform_wave_grid(data[:,0]/10.,data[:,1])
        keff_NH3,leff_NH3=NFL.K_eff(np.array([WG,SG]).T,NH3,15,filt)
        print(filt)
        print(keff_NH3,leff_NH3)
        keff_CH4,leff_CH4=NFL.K_eff(np.array([WG,SG]).T,CH4,15,filt)
        print(keff_CH4,leff_CH4)
        
        #axs1.plot(data[:,0]/10.,data[:,1],label=str(filt))
        axs1.plot(WG,SG,label=str(filt)+"N")
    
    axs1.set_ylim(0,0.25)
    axs1.set_xlim(600.,680)
    axs1.legend(loc=1)
    
    axsabs=axs1.twinx()
    #axsabs.set_ylim(0,10)
    
    #axsabs.set_yscale('log')
    #axsabs.set_ylim(5e-3,1e2)
    axsabs.set_ylim(0,5)

    axsabs.plot(NH3[:,0],NH3[:,1],label="NH3",linestyle='dashed')
    axsabs.plot(CH4[:,0],CH4[:,1]*5,'C3',label="CH4",linestyle="dashed")
    axsabs.legend(loc=2)



HST_Filters([619,631,645])
#print()
#print('631')
#HST_Filters(631)
#print()
#print('645')
#HST_Filters(645)