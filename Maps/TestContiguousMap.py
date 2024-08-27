def TestContiguousMap(axNH3,axCH4,axRGB,collection="20230815-20230818",LonSys='2'):
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
    import plot_patch as PP
    import make_patch_RGB as MPRGB
    import RetrievalLibrary as RL

               
    if collection=="20220810-20220812":
        dataset={'20220810UTa',
                 '20220812UTa'}
        
        lonlims={'20220810UTa':[0,60],
                 '20220812UTa':[240,360]}
        
    elif collection=="20220828-20220901":
        dataset={'20220828UTa',
                 '20220830UTa',
                 '20220901UTa'}
        
        lonlims={'20220828UTa':[0,100],
                 '20220830UTa':[290,360],
                 '20220901UTa':[200,290]}
        
    elif collection=="20220904-20220905":
        dataset={'20220904UTa',
                 '20220905UTa'}
        
        lonlims={'20220904UTa':[0,85],
                 '20220905UTa':[85,210]}
        
    elif collection=="20220912-20220913":
        dataset={'20220912UTa',
                 '20220913UTa'}
        
        lonlims={'20220912UTa':[40,170],
                 '20220913UTa':[170,290]}
        
    elif collection=="20220919-20220919":
        dataset={'20220919UTa',
                 '20220919UTb'}
        
        lonlims={'20220919UTa':[300,360],
                 '20220919UTb':[0,110]}
        
    elif collection=="20221009-20221013":
        dataset={'20221009UTa',
                 '20221009UTb',
                 '20221013UTa'}
        
        lonlims={'20221009UTa':[90,180],
                 '20221009UTb':[180,260],
                 '20221013UTa':[0,90]}
        
    elif collection=="20221019-20221021":
        dataset={'20221019UTa',
                 '20221020UTa',
                 '20221021UTa'}
        
        lonlims={'20221019UTa':[175,260],
                 '20221020UTa':[0,80],
                 '20221021UTa':[80,175]}
        
    elif collection=="20230815-20230818":
        dataset={'20230815UTb',
                 '20230816UTb',
                 '20230817UTa',
                 '20230818UTa'}
        
        lonlims={'20230815UTb':[130,205],
                 '20230816UTb':[285,359],
                 '20230817UTa':[40,130],
                 '20230818UTa':[205,285]}
                
    elif collection=="20230827-20230831":
        dataset={'20230827UTa',
                 '20230827UTb',
                 '20230830UTa',
                 '20230831UTa',
                 '20230831UTb'}
        
        lonlims={'20230827UTa':[115,185],
                 '20230827UTb':[185,205],
                 '20230830UTa':[205,320],
                 '20230831UTa':[0,45],
                 '20230831UTb':[45,115]}
        
    elif collection=="20230905-20230906":
        dataset={'20230905UTa',
                 '20230905UTb',
                 '20230906UTa',
                 '20230906UTb'}
        
        lonlims={'20230905UTa':[20,85],
                 '20230905UTb':[85,150],
                 '20230906UTa':[160,230],
                 '20230906UTb':[230,300]}
        
    elif collection=="20230922-20230929":
        dataset={'20230922UTa',
                 '20230924UTa',
                 '20230929UTa'}
        
        #lonlims={'20230922UTa':[245,335],
        #         '20230924UTa':[10,110],
        #         '20230929UTa':[175,275]}
        lonlims={'20230922UTa':[132,232],
                 '20230924UTa':[283,360],
                 '20230929UTa':[123,223]}
        
    elif collection=="20231005-20231006":
        dataset={'20231005UTa',
                 '20231006UTa',
                 '20231006UTb'}
        
        lonlims={'20231005UTa':[0,70],
                 '20231006UTa':[120,190],
                 '20231006UTb':[190,260]}
        
    elif collection=="20231015-20231019":
        dataset={'20231015UTa',
                 '20231015UTb',
                 '20231017UTa',
                 '20231017UTb',
                 '20231019UTa',
                 '20231019UTb'}
        
        lonlims={'20231015UTa':[60,95],
                 '20231015UTb':[95,150],
                 '20231017UTa':[0,35],
                 '20231017UTb':[35,60],
                 '20231019UTa':[280,335],
                 '20231019UTb':[335,360]}
                
    elif collection=="20231022-20231026":
        dataset={'20231022UTa',
                 '20231023UTa',
                 '20231023UTb',
                 '20231024UTa',
                 '20231024UTb',
                 '20231026UTa'}
        
        lonlims={'20231022UTa':[0,100],
                 '20231023UTa':[130,200],
                 '20231023UTb':[200,270],
                 '20231024UTa':[270,340],
                 '20231024UTb':[340,360],
                 '20231026UTa':[260,290]}
                
    elif collection=="20231103-20231107":
        dataset={'20231103UTa',
                 '20231107UTa',
                 '20231107UTb'}
        
        lonlims={'20231103UTa':[0,100],
                 '20231107UTa':[190,250],
                 '20231107UTb':[250,340]}
        
    elif collection=="20231110-20231110":
        dataset={'20231110UTa',
                 '20231110UTb',
                 '20231110UTc'}
        
        lonlims={'20231110UTa':[290,360],
                 '20231110UTb':[0,30],
                 '20231110UTc':[30,100]}
        
    elif collection=="20231112-20231113":
        dataset={'20231112UTa',
                 '20231112UTb',
                 '20231113UTa',
                 '20231113UTb'}
        
        lonlims={'20231112UTa':[250,300],
                 '20231112UTb':[300,360],
                 '20231113UTa':[0,55],
                 '20231113UTb':[55,110]}
        
    elif collection=="20231115-20231115":
        dataset={'20231115UTa',
                 '20231115UTb'}
        
        lonlims={'20231115UTa':[310,360],
                 '20231115UTb':[0,100]}
        
    elif collection=="20231128-20231129":
        dataset={'20231128UTa',
                 '20231128UTb',
                 '20231129UTa',
                 '20231129UTb',
                 '20231129UTc'}
        
        lonlims={'20231128UTa':[130,200],
                 '20231128UTb':[200,260],
                 '20231129UTa':[260,310],
                 '20231129UTb':[310,360],
                 '20231129UTc':[0,50]}
        
    elif collection=="20231206-20231207":
        dataset={'20231206UTa',
                 '20231206UTb',
                 '20231207UTa',
                 '20231207UTb',
                 '20231207UTc'}
        
        lonlims={'20231206UTa':[230,290],
                 '20231206UTb':[290,360],
                 '20231207UTa':[0,60],
                 '20231207UTb':[60,100],
                 '20231207UTc':[100,160]}
                              
    elif collection=="20231217-20231218":
        dataset={'20231217UTa',
                 '20231217UTb',
                 '20231217UTc',
                 '20231218UTa',
                 '20231218UTb',
                 '20231218UTc',
                 '20231218UTd',
                 '20231218UTd',
                 '20231218UTe',
                 '20231218UTf'}
        
        lonlims={'20231217UTa':[45,80],
                 '20231217UTb':[80,110],
                 '20231217UTc':[110,135],
                 '20231218UTa':[135,195],
                 '20231218UTb':[195,255],
                 '20231218UTc':[255,285],
                 '20231218UTd':[285,315],
                 '20231218UTe':[315,360],
                 '20231218UTf':[0,45]}
                                  
    elif collection=="20231229-20231229":
        dataset={'20231229UTa',
                 '20231229UTb',
                 '20231229UTc'}
        
        lonlims={'20231229UTa':[0,75],
                 '20231229UTb':[75,125],
                 '20231229UTc':[125,200]}
                                
    elif collection=="20240129-20240202":
        dataset={'20240129UTa',
                 '20240129UTb',
                 '20240130UTa',
                 '20240130UTb',
                 '20240130UTc',
                 '20240131UTa',
                 '20240131UTb',
                 '20240202UTa',
                 '20240202UTb'}
        
        lonlims={'20240129UTa':[50,85],
                 '20240129UTb':[85,145],
                 '20240130UTa':[145,200],
                 '20240130UTb':[200,250],
                 '20240130UTc':[250,280],
                 '20240131UTa':[0,30],
                 '20240131UTb':[30,50],
                 '20240202UTa':[280,320],
                 '20240202UTb':[320,360]}
                               
    elif collection=="20240229-20240301":
        dataset={'20240229UTa',
                 '20240229UTb',
                 '20240301UTa',
                 '20240301UTb'}
        
        lonlims={'20240229UTa':[0,53],
                 '20240229UTb':[53,120],
                 '20240301UTa':[130,197],
                 '20240301UTb':[197,265]}

        lonlims={'20240229UTa':[0,90],
                 '20240229UTb':[15,120],
                 '20240301UTa':[130,235],
                 '20240301UTb':[146,265]}
                
    
    outputfNH3=np.zeros([180,360])
    outputPCloud=np.zeros([180,360])
    outputRGB=np.zeros([180,360,3])
    
    stackfNH3=np.zeros([180,360])
    stackPCloud=np.zeros([180,360])
    stackRGB=np.zeros([180,360,3])
    
    #LonSys='1'
    print("dataset=",dataset)
    First=True
    
    for obskey in dataset:
        print("*******obsdate=",obskey)
        PCloudhdr,PClouddata,fNH3hdr,fNH3data,sza,eza,RGB,RGB_CM2,RGBtime= \
                        RFM.read_fits_map_L2_L3(obskey=obskey,LonSys=LonSys,
                                                imagetype="Map",Level="L3")

        amfdata=(1.0/sza+1.0/eza)/2.0
        #TestfNH3=fNH3data*amfdata**0.65
        #TestPCloud=PClouddata*amfdata**0.25
        #TestfNH3=fNH3data*amfdata**0.65
        TestfNH3=fNH3data*(amfdata**0.55)

        TestPCloud=PClouddata*amfdata**0.25
        print("**********TestfNH3.shape=",TestfNH3.shape)
        
        #lats=[30,150]
        lats=[60,120]
        #ll_0=360-lonlims[obskey][0]
        #ll_1=360-lonlims[obskey][1]
        if LonSys=='1':
            ll_0=int(360-(fNH3hdr["CM1"]-50))
            ll_1=int(360-(fNH3hdr["CM1"]+50))
        elif LonSys=='2':
            ll_0=int(360-(fNH3hdr["CM2"]-50))
            ll_1=int(360-(fNH3hdr["CM2"]+50))
        
        #print("***************** ll_0x, ll_1x=",ll_0x,ll_1x)
        if ll_0>360:
            ll_0=360
        if ll_1<0:
            ll_1=0
        print("***************** ll_0, ll_1=",ll_0,ll_1)
        #print("***************** ll_0x, ll_1x=",ll_0x,ll_1x)
        
        outputfNH3[lats[0]:lats[1],ll_1:ll_0]= \
            TestfNH3[lats[0]:lats[1],ll_1:ll_0]
        outputPCloud[lats[0]:lats[1],ll_1:ll_0]= \
            TestPCloud[lats[0]:lats[1],ll_1:ll_0]
        outputRGB[lats[0]:lats[1],ll_1:ll_0,:]= \
            RGB[lats[0]:lats[1],ll_1:ll_0,:]
            
        tempfNH3=np.zeros([180,360])
        tempPCloud=np.zeros([180,360])
        tempfNH3[lats[0]:lats[1],ll_1:ll_0]= \
            TestfNH3[lats[0]:lats[1],ll_1:ll_0]
        tempPCloud[lats[0]:lats[1],ll_1:ll_0]= \
            TestPCloud[lats[0]:lats[1],ll_1:ll_0]

        if First:
            stackfNH3=tempfNH3
            stackfNH3=np.reshape(stackfNH3,(180,360,1))
            stackPCloud=tempPCloud
            stackPCloud=np.reshape(stackPCloud,(180,360,1))
        else:
            tempfNH3=np.reshape(tempfNH3,(180,360,1))
            stackfNH3=np.dstack((stackfNH3,tempfNH3))
            tempPCloud=np.reshape(tempPCloud,(180,360,1))
            stackPCloud=np.dstack((stackPCloud,tempPCloud))
        First=False
        
    print(stackfNH3.shape)
    indzf=np.where(stackfNH3==0)
    stackfNH3[indzf]=np.nan
    blendfNH3=np.nanmean(stackfNH3,axis=2)
    indzP=np.where(stackPCloud==0)
    stackPCloud[indzP]=np.nan
    blendPCloud=np.nanmean(stackPCloud,axis=2)
    print(blendfNH3.shape)


    #6x6 works for 60N-60Z and 360 long
    #8x5 works for 30N-30S and 360 Long
    fig1,axs1=pl.subplots(3,1,figsize=(8.0,5), dpi=150, facecolor="white",
                          sharex=True,sharey=True)
    fig1.suptitle(collection)

    for ix in range(0,2):
        axs1[ix].grid(linewidth=0.2)
        axs1[ix].ylim=[-30.,30.]
        axs1[ix].xlim=[0.,360.]
        axs1[ix].set_xticks(np.linspace(450.,0.,31), minor=False)
        xticklabels=np.array(np.mod(np.linspace(450,0,31),360))
        axs1[ix].set_xticklabels(xticklabels.astype(int))
        axs1[ix].set_yticks(np.linspace(-30,30,7), minor=False)
        axs1[ix].tick_params(axis='both', which='major', labelsize=7)
        axs1[ix].set_ylabel("PG Latitude (deg)")

        axs1[ix].set_adjustable('box') 


    fNH3_patch_mb,vn,vx,tx_fNH3=PP.plot_patch(blendfNH3,lats,[0,360],
                                     180,180,"jet",
                                     axs1[0],'%3.2f',cont=False,n=11,vn=50,vx=150)
    axs1[0].set_title('fNH3 (ppm)',fontsize=10)

    PCld_patch_mb,vn,vx,tx_fNH3=PP.plot_patch(blendPCloud,lats,[0,360],
                                     180,180,"jet",
                                     axs1[1],'%3.2f',cont=False,n=6,vn=600,vx=1100)
    axs1[1].set_title('PCloud (mbar)',fontsize=10)
    
    RGB_patch=MPRGB.make_patch_RGB(outputRGB,lats,[0,360],180,180)
    RGB4Display=np.power(np.array(RGB_patch).astype(float),1.0)
    RGB4Display=RGB4Display/RGB4Display.max()
    show=axs1[2].imshow(RGB4Display,
               extent=[360,0,90-lats[1],
                       90-lats[0]],
                       aspect="equal")
    #temp=RL.make_contours_CH4_patch(axs1[2],outputfNH3,[30.,150.],[0.,360.],
    #                       tx_fNH3,frmt='%3.0f',clr='k')
    axs1[2].tick_params(axis='both', which='major', labelsize=7)
    axs1[2].set_xlabel("Sys. "+LonSys+" Longitude (deg)",fontsize=9)
    axs1[2].set_ylabel("PG Latitude (deg)")

    box = axs1[2].get_position()
    
    #works for 6x6 and 60N-60S and 360 long
    fig1.subplots_adjust(left=0.10, bottom=0.05, right=0.90, top=0.93,
                wspace=0.25, hspace=0.05)     
    axs1[2].set_position([box.x0-0.075, box.y0-0.04, box.width * 1.11, box.height * 1.11])

    #works for 8x5 and 30N-30S and 360 long
    fig1.subplots_adjust(left=0.07, bottom=0.06, right=0.94, top=0.92,
                wspace=0.25, hspace=0.08)     
    axs1[2].set_position([box.x0-0.055, box.y0-0.04, box.width * 1.07, box.height * 1.07])

    pathmapplots="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Studies/maps/"
    fig1.savefig(pathmapplots+collection+"Sys"+LonSys+" map.png",dpi=300)
    
    
    ###########################################################################
    lats=[100,120]
    fNH3_patch_mb,vn,vx,tx_fNH3=PP.plot_patch(blendfNH3,lats,[0,360],
                                     180,180,"jet",
                                     axNH3,'%3.2f',cbarplot=False,cont=False,n=11,vn=50,vx=150)
    axNH3.set_ylabel(collection.replace('-','\n'),rotation='horizontal',fontsize=6)
    axNH3.yaxis.set_label_coords(-0.05,0.15)
    axNH3.tick_params('y', labelleft=False)
    axNH3.tick_params('x', labelsize=8)

    PCld_patch_mb,vn,vx,tx_fNH3=PP.plot_patch(blendPCloud,lats,[0,360],
                                     180,180,"jet",
                                     axCH4,'%3.2f',cbarplot=False,cont=False,n=5,vn=600,vx=1000)
    axCH4.set_ylabel(collection,rotation='horizontal',fontsize=6)
    axCH4.set_ylabel(collection.replace('-','\n'),rotation='horizontal',fontsize=6)
    axCH4.yaxis.set_label_coords(-0.05,0.15)
    axCH4.tick_params('y', labelleft=False)
    axCH4.tick_params('x', labelsize=8)

    RGB_patch=MPRGB.make_patch_RGB(outputRGB,lats,[0,360],180,180)
    RGB4Display=np.power(np.array(RGB_patch).astype(float),1.0)
    RGB4Display=RGB4Display/RGB4Display.max()
    show=axRGB.imshow(RGB4Display,
               extent=[360,0,90-lats[1],
                       90-lats[0]],
                       aspect="equal")
    axRGB.set_ylabel(collection.replace('-','\n'),rotation='horizontal',fontsize=6)
    axRGB.tick_params('y', labelleft=False)
    axRGB.yaxis.set_label_coords(-0.05,0.15)
    axRGB.tick_params('x', labelsize=8)


    return()