# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 13:21:32 2025

@author: smhil
"""



def gravity(planet="Jupiter"):
#This function calculates gravity (cm/s2) at a given planet latitude
    import numpy as np
    import matplotlib.pyplot as pl
    
    G=6.673*10**(-11)
    lats=90.0-np.linspace(0.5,179.5,180)
    garr=np.zeros((180,1))
    
    if planet=='Saturn':
        M=568*10**24
        Re=60268000
        epsilon=0.0980
        T=38362
    elif planet=='Jupiter':
        M=1901*10**24 #(1.898 × 10^27 kg)
        Re=71541000
        epsilon=0.06492
        T=35730
    
    latindx=0
    for lat in lats:
        print(latindx,lat)
        garr[latindx,0]=(100*G*M)/(Re**2*(1-epsilon*np.sin(np.deg2rad(lat))**2))-4*(np.pi**2)*Re*(1-epsilon*np.sin(np.deg2rad(lat))**2)*np.cos(np.deg2rad(lat))/T**2
        latindx=latindx+1
    #pl.plot(lats,garr[:,0])
    
    gmap=garr
    #print(gmap)
    for i in range(0,359):
        gmap=np.hstack((gmap,garr))
        
    #print(gmap.shape)
    #pl.imshow(gmap)
    return(gmap)