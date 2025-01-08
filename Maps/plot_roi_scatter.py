# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 07:19:56 2024

@author: smhil
"""

def plot_roi_scatter(patch1,patch2,Real_CM2,LatLims,LonLims,axscor,PCldlow,PCldhigh,
                 fNH3low,fNH3high,FiveMicron,axis_inv=False,ROI=False):
    import pylab as pl
    import numpy as np
    import copy
    
    ###########################################################################
    # LOOP OVER ROIS AND PLOT SCATTER IN APPROPRIATE COLOR
    ###########################################################################
    #print("In ROI",ROI)
    mean1=[]
    stdv1=[]
    mean2=[]
    stdv2=[]
    roilabel=[]
    
    if ROI:
        for R in ROI:
            RLatLims=-LatLims[0]+np.array([ROI[R][0],ROI[R][1]])
            RCM=ROI[R][2]
            RLonRng=ROI[R][3]
            RLonLims=[360-int(RCM+RLonRng),360-int(RCM-RLonRng)]

            RLonLims=np.array(RLonLims)-LonLims[0]

            subpatch1=patch1[RLatLims[0]:RLatLims[1],
                             RLonLims[0]:RLonLims[1]]
            subpatch2=patch2[RLatLims[0]:RLatLims[1],
                             RLonLims[0]:RLonLims[1]]

            axscor.scatter(subpatch2,
                           subpatch1,
                           marker="o",s=3.0,
                           alpha=0.8,label=R)
            mean1.append(np.mean(subpatch1))
            mean2.append(np.mean(subpatch2))
            stdv1.append(np.std(subpatch1))
            stdv2.append(np.std(subpatch2))
            roilabel.append(R)
            #print("ROI",R)
            #print("patch1",mean1,stdv1)
            #print("patch2",mean2,stdv2)
     
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
                    
    axscor.legend(fontsize=7,ncols=3)
    
    
    return(roilabel,mean1,stdv1,mean2,stdv2)
  
