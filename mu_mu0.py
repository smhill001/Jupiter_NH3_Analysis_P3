# -*- coding: utf-8 -*-
"""
Created on Sat Feb  5 23:36:53 2022

@author: smhil
"""


import numpy as np
import scipy.io
import pylab as pl
B=-3.36
Bprime=+3.05
alpha=+0.73
#Latpg=0.
#DLong=0.

mat = scipy.io.loadmat('Montoya/Jupiter_calibrated-MatrixU')
fig1,ax1=pl.subplots(1,1,figsize=(6.0,6.0), dpi=150, facecolor="white",
                    sharey=True,sharex=True)

ax1.imshow(mat['mu'])

Latpgtemp=mat['lat']
Latpg=Latpgtemp.flatten()
longtemp=mat['lon']
long=longtemp.flatten()
mu_arr=np.zeros((161,161),dtype=float)
mu0_arr=np.zeros((161,161),dtype=float)
print(Latpg.size,long.size)
for i in range(0,160):
    for j in range(0,160):
        mu_arr[i,j]=np.sin(B*np.pi/180.)*np.sin(Latpg[i]*np.pi/180.) \
            +np.cos(B*np.pi/180.)*np.cos(Latpg[i]*np.pi/180.)*np.cos(long[j]*np.pi/180.)
        mu0_arr[i,j]=np.sin(Bprime*np.pi/180.)*np.sin(Latpg[i]*np.pi/180.) \
            + np.cos(Bprime*np.pi/180.)*np.cos(Latpg[i]*np.pi/180.)*np.cos(long[j]*np.pi/180.)

fig2,ax2=pl.subplots(1,1,figsize=(6.0,6.0), dpi=150, facecolor="white",
                    sharey=True,sharex=True)

ax2.imshow(mu_arr)

fig3,ax3=pl.subplots(1,1,figsize=(6.0,6.0), dpi=150, facecolor="white",
                    sharey=True,sharex=True)

ax3.imshow(mu_arr/mat['mu'])


#print(mu,mu0)
#print(1/mu+1/mu0)
#print(np.arcsin(1/(1/mu+1/mu0))*180/np.pi)
#print(np.arccos(mu)*180/np.pi,np.arccos(mu0)*180/np.pi)


