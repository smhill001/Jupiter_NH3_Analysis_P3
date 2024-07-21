# -*- coding: utf-8 -*-

def JupiterFilterPerformance(TelescopeList=["SCT","VLT"],
                             ModelList=['Spline1','Spline2','Piecewise1',
                                        'Piecewise2','Linear 2pt.'],
                             FilterList=['620','632','647','656']):
    """
    Purpose is to compute and plot crossection and albedo data then to 
        compute and plot filter performance for all combinations of 
        telescopes, continuum models, and filters
    Calls
        ->NH3_Filter_Library_P3.py/Get_Albedo_Continua_Crossections 
            for albedo and crossection data
        ->NH3_Filter_Library_P3.py/compute_filter_spectrum
            to compute and plot filter data for all allowed filter-telescope
            combinations. As of 7/19/2023, telescopes include SCT and VLT

    Returns
    -------
    None.

    """
    import sys
    sys.path.append('c:/Astronomy/Python Play')
    sys.path.append('c:/Astronomy/Python Play/Util_P3')
    sys.path.append('c:/Astronomy/Python Play/SPLibraries_P3')
    sys.path.append('c:/Astronomy/Python Play/SpectroPhotometry/Spectroscopy_P3')
    import numpy as np
    import NH3_Filter_Library_P3 as NFL
    sys.path.append('./Services')
    import get_albedo_continua_crossections as gACC
    import get_keff
    projpath="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/"

    ###########################################################################
    #!!!!NEED A MICRO SERVICE FOR THIS!!!!
    ###########################################################################
    gravity=22.280 #m/s^2  ######!!!!!!Needs to be fixed!!!
    mmolwt=3.85e-27 #kg/molecule
    fCH4=1.81e-3
    Lodschmit=2.687E+25 #m-3
    STP=1.01325e5  #Pa

    # Get molecular and Jovian data   
    x0,x1,xtks=600.,680.,9
    y0,y1,ytks=0.0,0.7,8
    Albedo,Continua,CH4,NH3,NH3_Lutz_Owen_1980,fig_molecules,ax_molecules= \
        gACC.get_albedo_continua_crossections(x0,x1,xtks,y0,y1,ytks,Crossect=True)

    # PLOT FILTER TRANSMISSIONS CONVOLVED WITH DISK-INTEGRATED ALBEDO AND CONTINUUM
    pathout='Molecular Absorption/data output/'
    filtereffectivedata = open(projpath+pathout+'JupiterFilterPerformance.csv', 'w')
    tmp="Wavelength (nm),Filter Name,k_eff (NH3),l_eff (NH3),k_eff (CH4),"\
        +"l_eff (CH4),Telescope, ContModel,Transmission,Tau,EW (nm),"\
        +"N[NH3],N[CH4],Pres (mb),f(NH3)\n"
    filtereffectivedata.write(tmp)
    
    for Tele in TelescopeList:
        for ModelName in ModelList:       
            Continuum_Albedo=np.zeros((Continua[ModelName]['WaveGrid'].size,2))
            Continuum_Albedo[:,0]=Continua[ModelName]['WaveGrid']
            Continuum_Albedo[:,1]=Continua[ModelName]['Albedo']
            
            filterdata,axsFilt= \
                NFL.compute_filter_spectrum(x0,x1,xtks,y0,y1,ytks,FilterList,
                                            Albedo,Continuum_Albedo,
                                            ModelName,Telescope=Tele)
                
            filterdata=get_keff.get_keff(filterdata,FilterList,CH4,NH3)



            for filtr in FilterList:
                filterdata[filtr]['NH3ColDens']=\
                    1000.*filterdata[filtr]['Tau_Albedo']/ \
                        filterdata[filtr]['keff_NH3']
                filterdata[filtr]['CH4ColDens']=\
                    1000.*filterdata[filtr]['Tau_Albedo']/ \
                        filterdata[filtr]['keff_CH4']
                        
                #!!!Need to fix this an the paper's formula for pressure!!
                filterdata[filtr]['Pres']=(filterdata[filtr]['CH4ColDens'])\
                    *gravity*mmolwt*Lodschmit/(4.0*fCH4*STP)#-> It works 7/20/2023
                    
                if filtr=='647':
                    filterdata[filtr]['fNH3']=fCH4*filterdata['647']['NH3ColDens']/\
                        filterdata['620']['CH4ColDens'] #-> It works 7/20/2023
                else:
                    filterdata[filtr]['fNH3']='-'
                    
                        
                print(filtr)
                tmp=filtr+","+filterdata[filtr]['filtname']+","\
                         +str(filterdata[filtr]['keff_NH3'])+","\
                         +str(filterdata[filtr]['leff_NH3'])+","\
                         +str(filterdata[filtr]['keff_CH4'])+","\
                         +str(filterdata[filtr]['leff_CH4'])+","\
                    +Tele+','\
                    +ModelName+','+str(filterdata[filtr]['TransInt'])+","\
                    +str(filterdata[filtr]['Tau_Albedo'])+","\
                    +"  ,"\
                    +str(filterdata[filtr]['NH3ColDens'])+","\
                    +str(filterdata[filtr]['CH4ColDens'])+","\
                    +str(filterdata[filtr]['Pres'])+","\
                    +str(filterdata[filtr]['fNH3'])+"\n"
                        
                filtereffectivedata.write(tmp)
    filtereffectivedata.close()
    """          
        
        tmp=filtr+","+Jupiterdata[filtr]['filtname']+","+str(Jupiterdata[filtr]['keff_NH3'])+","\
                +str(Jupiterdata[filtr]['leff_NH3'])+","+str(Jupiterdata[filtr]['keff_CH4'])+","\
                +str(Jupiterdata[filtr]['leff_CH4'])+","+str(Jupiterdata[filtr]['TransInt'])+","\
                +str(Jupiterdata[filtr]['Tau_Albedo'])+","+str(Jupiterdata[filtr]['NH3ColDens'])+","\
                +str(Jupiterdata[filtr]['CH4ColDens'])+"\n"
    """
