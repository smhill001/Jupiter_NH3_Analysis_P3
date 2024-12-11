# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 07:19:56 2024

@author: smhil
"""

def plot_map_scatter(patch1,patch2,Real_CM2,LatLims,axscor,PCldlow,PCldhigh,
                 fNH3low,fNH3high,FiveMicron,axis_inv=False):
    import pylab as pl
    import numpy as np
    import copy
    
    #print("000000000: ",patch1.shape,patch2.shape)
    BZ={"SSTB":[-39.6,-36.2],
          "STZ":[-36.2,-32.4],
          "STB":[-32.4,-27.1],
          "STrZ":[-27.1,-19.7],
          "SEB":[-19.7,-7.2],
          "SEZ":[-7.2,0.0],
          "NEZ":[0.0,6.9],
          "NEB":[6.9,17.4],
          "NTrZ":[17.4,24.2],
          "NTB":[24.2,31.4],
          "NTZ":[31.4,35.4],
          "NNTB":[35.4,39.6]}

    BZind=copy.deepcopy(BZ)   
    BZkeys=BZ.keys()
    #patch1=patch1*1000.

    #figcor,axscor=pl.subplots(1,1,figsize=(6.0,4.), dpi=150, facecolor="white",
    #                    sharey=True,sharex=True)          

    for key in BZ.keys():
        #print(key,BZ[key],[90,90]-np.array(BZ[key]),LatLims)
        BZind[key][0]=int(90-BZ[key][0])-LatLims[0]
        BZind[key][1]=int(90-BZ[key][1])-LatLims[0]
        
        if BZind[key][0]<0:
            BZind[key][0]=0
        if BZind[key][1]<0:
            BZind[key][1]=0
        if BZind[key][0]>(LatLims[1]-LatLims[0]):
            BZind[key][0]=LatLims[1]-LatLims[0]
        if BZind[key][1]>(LatLims[1]-LatLims[0]):
            BZind[key][1]=LatLims[1]-LatLims[0]
        
        if BZind[key][0]==BZind[key][1]:
            print("do nothing")
        else:
            print(patch2.shape,patch1.shape)
            axscor.scatter(patch2[BZind[key][1]:BZind[key][0],:],
                           patch1[BZind[key][1]:BZind[key][0],:],
                           marker="o",s=3.0,
                           alpha=0.8,label=key)
     
    axscor.grid(linewidth=0.2)
    axscor.set_ylim(PCldlow,PCldhigh)
    axscor.set_xlim(fNH3low,fNH3high)
    axscor.set_ylabel("Cloud-top Pressure (mb)",fontsize=10)
    if axis_inv:
        axscor.invert_yaxis()
    if FiveMicron:
        axscor.set_xlabel("5um Radiance (Log10(arb. units)",fontsize=10)
    else:
        axscor.set_xlabel("Ammonia Mole Fraction (ppm)",fontsize=10)
                    
    axscor.legend(fontsize=7,ncols=4)
    
    return(BZ)
  
