# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 07:19:56 2024

CALLED BY Profile_L2_script.py and Profile_L3_script.py

@author: smhil
"""

def plot_profile_scatter(profile1,profile2,Lats,axscor,PCldlow,PCldhigh,
                 fNH3low,fNH3high,FiveMicron,bands,colors,
                 Level="L2",leg=True,axis_inv=False,date='NA',
                 counter=36,countmax=0):
    
    import pylab as pl
    import numpy as np
    import copy
    from matplotlib.cm import get_cmap
    from matplotlib.colors import Normalize
    
    print("profile1",profile1.shape)
    print("profile2",profile2.shape)
    print("Lats",Lats.shape)
    
    #print("000000000: ",patch1.shape,patch2.shape)
    ###########################################################################
    # SET BELT AND ZONE BOUNDARIES
    ###########################################################################
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
    print("BZkeys",BZkeys)

    ###########################################################################
    # LOOP OVER BELTS AND PLOT SCATTER IN APPROPRIATE COLOR
    ###########################################################################
    idx=0
    print(bands)
    for key in bands:
        print(key)#,BZ[key],[90,90]-np.array(BZ[key]))
        BZind[key][0]=(np.abs(Lats-BZ[key][0])).argmin()
        BZind[key][1]=(np.abs(Lats-BZ[key][1])).argmin()

        """
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
        else:"""
        mean1=np.mean(profile1[BZind[key][0]:BZind[key][1]])
        mean2=np.mean(profile2[BZind[key][0]:BZind[key][1]])
        stdv1=np.std(profile1[BZind[key][0]:BZind[key][1]])
        stdv2=np.std(profile2[BZind[key][0]:BZind[key][1]])
        print("key,mean1,mean2",key,mean1,mean2)
        
        cmap = get_cmap('coolwarm')          # choose colormap
        norm = Normalize(vmin=0, vmax=countmax)
        color = cmap(norm(counter))
        
        if 220601<=int(date)<=221231:
            marker='o'
            clr='C0'
        if 230901<=int(date)<=240331:
            marker='^'
            clr='C1'
        if 240801<=int(date)<=250430:
            marker='s'
            clr='C2'
        if 250901<=int(date)<=260531:
            marker='D'
            clr='C3'
            
        """
        marker='o'
        if key=='SEZ':
            fclr='w'
        elif key=='NEZ':
            fclr=clr
        """
        
        if leg:
            #axscor.errorbar(mean2,mean1,xerr=stdv2,yerr=stdv1,marker=marker,
            #                markersize=7,elinewidth=0.5,ecolor=colors[idx],markeredgewidth=2,
            #                #alpha=counter/countmax,label=key,markerfacecolor=colors[idx],markeredgecolor=color)
            #                alpha=1,label=key,markerfacecolor=colors[idx],markeredgecolor=colors[idx])
            #                #alpha=0.8,label=key,markerfacecolor=colors[idx],markeredgecolor=colors[idx])
            axscor.errorbar(mean2,mean1,xerr=stdv2,yerr=stdv1,marker=marker,
                            markersize=7,elinewidth=0.5,ecolor=colors[idx],markeredgewidth=2,
                            #alpha=counter/countmax,label=key,markerfacecolor=colors[idx],markeredgecolor=color)
                            alpha=0.2+0.8*counter/countmax,label=key,markerfacecolor=colors[idx],markeredgecolor=colors[idx])
                            #alpha=0.8,label=key,markerfacecolor=colors[idx],markeredgecolor=colors[idx])
                            ##alpha=0.8,label=key,markerfacecolor=fclr,markeredgecolor=clr)

        else:
            axscor.errorbar(mean2,mean1,xerr=stdv2,yerr=stdv1,marker=marker,
                            markersize=7,elinewidth=0.5,ecolor=colors[idx],markeredgewidth=2,
                            #alpha=counter/countmax,markerfacecolor=colors[idx],markeredgecolor=color)
                            alpha=0.2+0.8*counter/countmax,markerfacecolor=colors[idx],markeredgecolor=colors[idx])
                            #alpha=0.8,markerfacecolor=colors[idx],markeredgecolor=colors[idx])
                            #alpha=0.8,markerfacecolor=fclr,markeredgecolor=clr)

                       
        #if date!="NA":
        #    axscor.annotate(date,(mean2,mean1),textcoords='offset points',
        #                    xytext=(1,-1),fontsize=6,color=colors[idx])
        
        idx=idx+1
    axscor.grid(linewidth=0.2)
    axscor.set_ylim(PCldlow,PCldhigh)
    axscor.set_xlim(fNH3low,fNH3high)
    if Level=="L3":
        axscor.set_ylabel("Cloud-top Pressure (mb)",fontsize=14)
        if axis_inv:
            axscor.invert_yaxis()
        if FiveMicron:
            axscor.set_xlabel("5um Radiance (Log10(arb. units)",fontsize=12)
        else:
            axscor.set_xlabel("Ammonia Mole Fraction (ppm)",fontsize=12)
    elif Level=="L2":
        axscor.set_ylabel("Methane Equivalent Width (nm)",fontsize=12)
        if axis_inv:
            axscor.invert_yaxis()
        if FiveMicron:
            axscor.set_xlabel("5um Radiance (Log10(arb. units)",fontsize=12)
        else:
            axscor.set_xlabel("Ammonia Equivalent Width (nm)",fontsize=12)
                    
    if leg:
        axscor.legend(fontsize=8,ncols=2)
        
    return(BZ)
  
