def TestContiguousMap():
    """
    Created on Wed Dec 20 21:02:31 2023

    @author: smhil
    """
    import sys
    drive='c:'
    sys.path.append(drive+'/Astronomy/Python Play')
    sys.path.append(drive+'/Astronomy/Python Play/Util_P3')
    sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy')
    sys.path.append('./Services')

    import os
    #from matplotlib.pyplot import imread
    import pylab as pl
    import numpy as np
    from imageio import imwrite
    #from numpy import inf
    #from re import search
    from astropy.io import fits
    #from astropy.convolution import Gaussian2DKernel
    #from astropy.convolution import convolve
    import RetrievalLibrary as RL
    sys.path.append('./Maps')
    #import get_L2_abs_data as GAOD
    #import make_L3_env_data
    import read_fits_map_L2_L3 as RFM
    
    dataset={'20230815UTa',
             '20230816UTa',
             '20230817UTa',
             '20230818UTa'}
    
    lonlims={'20230815UTa':[130,205],
             '20230816UTa':[285,359],
             '20230817UTa':[40,130],
             '20230818UTa':[205,285]}
    
    outputfNH3=np.zeros([180,360])
    outputPCloud=np.zeros([180,360])
    
    for obskey in dataset:
        print("*******obsdate=",obskey)
        PCloudhdr,PClouddata,fNH3hdr,fNH3data,sza,eza,RGB,RGB_CM2= \
                        RFM.read_fits_map_L2_L3(obskey=obskey,
                                                imagetype="Map",Level="L3")
                        
        amfdata=(1.0/sza+1.0/eza)/2.0
        TestfNH3=fNH3data*amfdata**0.65
        TestPCloud=PClouddata*amfdata**0.25
        print("**********TestfNH3.shape=",TestfNH3.shape)
        
        ll_0=360-lonlims[obskey][0]
        ll_1=360-lonlims[obskey][1]
        outputfNH3[45:135,ll_1:ll_0]= \
            TestfNH3[45:135,ll_1:ll_0]
        outputPCloud[45:135,ll_1:ll_0]= \
            TestPCloud[45:135,ll_1:ll_0]

    fig1,axs1=pl.subplots(2,1,figsize=(8.0,6.0), dpi=150, facecolor="white",
                        sharey=True,sharex=True)

    axs1[0].imshow(outputfNH3*1e6,"jet",vmin=50,vmax=150)
    axs1[0].xlim=[0,360]
    axs1[0].set_xticks(np.linspace(360,0,13), minor=False)
    xticklabels=np.array(np.linspace(360,0,13))
    axs1[0].set_xticklabels(xticklabels.astype(int))
    axs1[1].imshow(outputPCloud,"jet",vmin=400,vmax=900)
      
                        
                        