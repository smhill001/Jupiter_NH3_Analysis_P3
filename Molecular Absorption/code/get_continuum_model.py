# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 16:20:06 2023

@author: smhil
"""
def get_continuum_model(Albedo,ModelList=['Spline1','Spline2','Piecewise1',
                                          'Piecewise2','Linear 2pt.']):
    import numpy as np
    from scipy import interpolate
    import GeneralSpecUtils_P3 as GSU

    WaveGrid=Albedo[:,0]
    SplineWV1= np.array([560.0, 580.0, 600.0,610.0, 632.0, 657.0, 677.0, 690., 714.0,
                         745.0, 830.0, 945.0,1050.0])
    SplineWV2 = np.array([560.0, 580.0, 600.0, 632.0, 678.0, 
                         745.0, 830.0, 945.0,1050.0])
    ExtrapWV = np.array([632.,657.])
    Continua={'Spline1':{'WV':SplineWV1,'MType':'Spline'},
              'Spline2':{'WV':SplineWV2,'MType':'Spline'},
              'Piecewise1':{'WV':SplineWV1,'MType':'Piecewise'},
              'Piecewise2':{'WV':SplineWV2,'MType':'Piecewise'},
              'Linear 2pt.':{'WV':ExtrapWV,'MType':'Linear'}}
    
    for Model in ModelList:
        print(Model)
        Continua[Model]['Mag']=np.ones(Continua[Model]['WV'].size)
        #Loop through spline wavelengths and look up albedo magnitude
        for i in range(0,Continua[Model]['WV'].size):
            Start=Continua[Model]['WV'][i]-.0000001
            End=Continua[Model]['WV'][i]+.0000001
            SplineWVIndices=np.where((Albedo[:,0] >Start) & \
                 (Albedo[:,0] < End))
            #print("i= ",i,SplineWVIndices)
            if Continua[Model]['MType']=='Spline':
                Continua[Model]['Mag'][i]=np.log10(Albedo[SplineWVIndices[0],1])
            elif Continua[Model]['MType']=='Piecewise' or Continua[Model]['MType']=='Linear':
                Continua[Model]['Mag'][i]=Albedo[SplineWVIndices[0],1]
        
        if Continua[Model]['MType']=='Spline':
            Continua[Model]['Mag'][Continua[Model]['Mag'].size-1]=np.log10(0.42) #what is this hard code?!!
            #Compute the knots and derivatives of the spline fit
            tck = interpolate.splrep(Continua[Model]['WV'], Continua[Model]['Mag'], s=0)
            #Evaluate the spline log magnitudes across the standard wave grid
            Temp = 10**interpolate.splev(WaveGrid, tck, der=0)
        elif Continua[Model]['MType']=='Piecewise':
            WaveGrid,SignalonGrid=GSU.uniform_wave_grid(Continua[Model]['WV'],
                                                        Continua[Model]['Mag'],
                                                        Extend=False,Fine=False)
            Temp=SignalonGrid
        elif Continua[Model]['MType']=='Linear':
            print("*********************Linear",Model)
            print(Continua[Model]['WV'])
            print(Continua[Model]['Mag'])
            Slope=(Continua[Model]['Mag'][1]-Continua[Model]['Mag'][0])/ \
                (Continua[Model]['WV'][1]-Continua[Model]['WV'][0])
            print("Slope=",Slope)
            Temp=Slope*(WaveGrid-Continua[Model]['WV'][0])+Continua[Model]['Mag'][0]

        Continua[Model]['WaveGrid']=WaveGrid
        Continua[Model]['Albedo']=Temp  #Continuum Albedo, that is

    return(Continua)