"""
Created on Tue Jun  4 12:39:15 2024

@author: smhil
"""
import sys
sys.path.append('c:/Astronomy/Python Play')
import pylab as pl
import numpy as np
from astropy.convolution import convolve
import plot_TEXES_Groups_P3 as PTG
from astropy.convolution import Box1DKernel
sys.path.append('./Maps')
import Profile_L2 as PFL2
from scipy.optimize import curve_fit
from scipy import interpolate

def sorted_reduced(x,y,xmin,xmax):
    """
    Parameters
    ----------
    x : TYPE
        Array of x-values to be ordered by increasing value
    y : TYPE
        Array of y-values to be ordered according to x-values by increasing value
    xmin : TYPE
        Minimum x-value to include.
    xmax : TYPE
        Maximum x-value to include.

    Returns
    -------
    xreduced : TYPE
        Array of x-values order by increasing value and bounded by xmin and xmax.
    yreduced : TYPE
        Array of y-values corresponding to sorted, bounded x-values

    """
    ###VLT22 special treatment
    idx_sort = x.argsort()
    xsort=x[idx_sort]
    ysort=y[idx_sort]

    reducedindices=[i for i, e in enumerate(np.array(xsort)) if (xmin<=e<=xmax)]
    xreduced=xsort[reducedindices]
    yreduced=ysort[reducedindices]

    return(xreduced,yreduced)

def lin(x, a, b):
    return a + (b * x)
def exp(x, a, b):
    return a * np.exp(b * x)
def power(x, a, b):
    return a * x**b
def minnaert(x, a, k):
    return a * (x[0]**k)*(x[1]**(k-1))

def make_fit_plot(NH3Amf_in,NH3Pro_in,CH4Amf_in,CH4Pro_in,Level,fig,ax,path):
    
    from astropy.stats import sigma_clip
    import copy
    import scipy.signal as ss
    
    ###########################################################################
    # Establish cleaned (median filter with 2 sigma threshold) NH3 profile
    ###########################################################################
    NH3Amf,NH3Pro=sorted_reduced(NH3Amf_in,NH3Pro_in,1.0,2.0)
    NH3med=ss.medfilt(np.array(NH3Pro),3)

    tmp=np.array(NH3Pro-NH3med)
    sig=2.0*np.std(tmp)
    clipindices=[i for i, e in enumerate(np.array(tmp)) if (np.abs(e)>sig)]
    print("******clipindices NH3=",clipindices)
    
    NH3pro_clip=copy.deepcopy(NH3Pro)
    if len(clipindices)!=0:
        print(NH3pro_clip.shape,NH3med.shape,)
        NH3pro_clip[clipindices]=NH3med[clipindices]
    
    ###########################################################################
    # Fit NH3 profile, determine R2, and plot
    ###########################################################################
    NH3param, NH3param_cov = curve_fit(power,NH3Amf,NH3pro_clip)
    fityNH3=power(NH3Amf,NH3param[0],NH3param[1])
    tempR=np.corrcoef(NH3pro_clip,fityNH3) 
    NH3R2=tempR[0,1]**2
 
    ax[0].plot(NH3Amf,fityNH3,label="Fit R2="+str(NH3R2)[0:4])
    ax[0].scatter(NH3Amf_in,NH3Pro_in,label="NH3",s=0.5,alpha=0.5)
    ax[0].scatter(NH3Amf,NH3pro_clip,label="NH3 clip",color='k',s=0.5)
    
    
    CH4Amf,CH4Pro=sorted_reduced(CH4Amf_in,CH4Pro_in,1.0,2.0)
    CH4med=ss.medfilt(np.array(CH4Pro),3)
    
    tmp=np.array(CH4Pro-CH4med)
    sig=2.0*np.std(tmp)
    clipindices=[i for i, e in enumerate(np.array(tmp)) if (np.abs(e)>sig)]
    print("******clipindices CH4=",clipindices)

    CH4pro_clip=copy.deepcopy(CH4Pro)
    if len(clipindices)!=0:
        print(NH3pro_clip.shape,NH3med.shape,)
        CH4pro_clip[clipindices]=CH4med[clipindices]
    CH4param, CH4param_cov = curve_fit(power,CH4Amf,CH4pro_clip)
    fity=power(CH4Amf,CH4param[0],CH4param[1])
    tempR=np.corrcoef(CH4pro_clip,fity) 
    CH4R2=tempR[0,1]**2
    
    ax[0].plot(CH4Amf,fity,label="Fit R2="+str(CH4R2)[0:4])
    ax[0].scatter(CH4Amf_in,CH4Pro_in,label="CH4",s=5,alpha=0.5)
    ax[0].scatter(CH4Amf,CH4pro_clip,label="CH4 clip",color='k',s=0.5)

    ###############################################################################
    # SCT23: ADJUSTED SCATTER PLOTS OF NH3 and CH4 VERSUS JOVIAN AIRMASS
    ###############################################################################
    ax[0].scatter(NH3Amf_in,NH3Pro_in*2.6,label="Scaled NH3")
    #ax[0].scatter(NH3Amf_in,NH3Pro_in*NH3Amf_in**0.93,s=5,label="Corr. NH3")
    ax[0].scatter(NH3Amf,NH3Pro*((fityNH3/NH3param[0])**-1.0),s=5,label="Corr. NH3")
    ax[0].scatter(CH4Amf_in,CH4Pro_in*CH4Amf_in**0.43,s=5,label="Corr. CH4")
    
    ax[0].tick_params(axis='both', which='major', labelsize=8)
    ax[0].set_xlabel("One-Way Airmass Factor",fontsize=10)
    ax[0].set_xlim(1,3)
    ax[0].set_xticks(np.linspace(1,3,5), minor=False)
    if Level=="L2":
        ax[0].set_ylim(0.0,2.0)
        ax[0].set_yticks(np.linspace(0.0,2.0,5), minor=False)
        ax[0].set_ylabel("Equivalent Width (nm)",fontsize=10)
    if Level=="L3":
        ax[0].set_ylim(0.0,1200.0)
        ax[0].set_yticks(np.linspace(0.0,1200,7), minor=False)
        ax[0].set_ylabel("PCloud (mb) and fNH3 (ppm)",fontsize=10)
    ax[0].grid(linewidth=0.2)
    ax[0].legend(loc="best",fontsize=7,ncol=3)
    
    ###############################################################################
    # SCT23: SCATTER PLOT OF CH4/NH3 RATIO VERSUS JOVIAN AIRMASS
    ###############################################################################

    Interp=interpolate.interp1d(CH4Amf,CH4pro_clip,kind='linear', 
                                copy=True,bounds_error=False,fill_value="extrapolate")  
    CH4S23ProonGrid=Interp(NH3Amf)
        
    RS23Amf,RS23Pro=NH3Amf,CH4S23ProonGrid/NH3pro_clip
    
    SCT23Rparam, SCT232Rparam_cov = curve_fit(power,RS23Amf,RS23Pro)
    fity=lin(RS23Amf,SCT23Rparam[0],SCT23Rparam[1])
    ax[1].plot(NH3Amf,fity,label="Fit 1")
    ax[1].scatter(RS23Amf,RS23Pro,label="CH4/NH3")
    ratiofity=power(NH3Amf,CH4param[0],CH4param[1])/power(NH3Amf,NH3param[0],NH3param[1])
       
    tempR=np.corrcoef(RS23Pro,ratiofity) 
    RatioR2=tempR[0,1]**2

    ax[1].plot(NH3Amf,ratiofity,label="Fit2 R2="+str(RatioR2)[0:4])
        
    ax[1].tick_params(axis='both', which='major', labelsize=8)
    ax[1].set_xlabel("One-Way Airmass Factor",fontsize=10)
    ax[1].set_xlim(1,3)
    ax[1].set_xticks(np.linspace(1,3,5), minor=False)
    if Level=="L2":
        ax[1].set_ylabel("CH4/NH3 Ratio",fontsize=10)
        ax[1].set_ylim(0.0,6.0)
        ax[1].set_yticks(np.linspace(0.0,6.0,7), minor=False)
    if Level=="L3":
        ax[1].set_ylabel("PCld/fNH3 Ratio",fontsize=10)
        ax[1].set_ylim(0.0,12.0)
        ax[1].set_yticks(np.linspace(0.0,12.0,7), minor=False)
    ax[1].grid(linewidth=0.2)
    ax[1].legend(loc="best",fontsize=7,ncol=1)


    ###############################################################################
    # SCT23: ADJUST SUBPLOTS and SAVE FIGURE
    ###############################################################################
    fig.subplots_adjust(left=0.09, right=0.97, top=0.93, bottom=0.11)
    #!!!! NEED TO MAKE THIS AN ADAPTIVE FILE NAME
    #fig.savefig(path+"Profiles/output/IrwinFig9SCT23.png",dpi=300)

    return(NH3param, NH3param_cov,NH3R2,
           CH4param, CH4param_cov,CH4R2,RatioR2,fig)

