def HST_Filters(filt):
    import sys
    sys.path.append('C:/Astronomy/Python Play/SpectroPhotometry/Spectroscopy_P3/')
    sys.path.append('C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Molecular Absorption/code/')
    
    import numpy as np
    import matplotlib.pyplot as pl
    import GeneralSpecUtils_P3 as GSU
    import NH3_Filter_Library_P3 as NFL
    import get_ammonia_data as get_NH3
    import get_karkoschka_data as get_CH4
    
    
    path="C:/Astronomy/Projects/SAS 2021 Ammonia/Data/HST GO 18055/Filter Throughput/"
    
    filterfiles={619:'wfc3_uvis1_fq619n.txt',
                 631:'wfc3_uvis1_f631n.txt',
                 645:'wfc3_uvis1_f645n.txt'}
    
    fn=filterfiles[filt]
    data = np.loadtxt(path+fn, skiprows=1)
    WG,SG=GSU.uniform_wave_grid(data[:,0]/10.,data[:,1])
    
    NH3=get_NH3.get_ammonia_data(Source="Irwin")
    CH4=get_CH4.get_karkoschka_data(Type='CH4')
    
    keff_NH3,leff_NH3=NFL.K_eff(np.array([WG,SG]).T,NH3,15,filt)
    print(keff_NH3,leff_NH3)
    keff_CH4,leff_CH4=NFL.K_eff(np.array([WG,SG]).T,CH4,15,filt)
    print(keff_CH4,leff_CH4)
    
    fig1,axs1=pl.subplots(1,figsize=(6,8), dpi=150, facecolor="white")
    fig1.suptitle("Transmission")
    
    axs1.plot(data[:,0]/10.,data[:,1])
    axs1.plot(WG,SG)
    axs1.set_xlim(600.,680.)
    axsabs=axs1.twinx()
    axsabs.set_yscale('log')
    
    axsabs.plot(NH3[:,0],NH3[:,1])
    axsabs.plot(CH4[:,0],CH4[:,1],'C3')

print('619')
HST_Filters(619)
print()
print('631')
HST_Filters(631)
print()
print('645')
HST_Filters(645)