def MakeContiguousMap(axNH3,axCH4,axRGB,collection="20241129-20241129",LonSys='2',
                      FiveMicron=False,Five_obskey='',IRTFdataset='',
                      lats=[75,105],LonLims=[0,360],figsz=[6.0,6.0],ROI=False,
                      variance=False,proj='maps',ctbls=['terrain_r','Blues']):
    """
    Created on Wed Dec 20 21:02:31 2023

    axNH3
    axCH4
    axRGB
    collection
    LonSys
    FiveMicron = False, 'png', or 'fits'
    Five_obskey = observation key to which I've attached the 5umfile metadata

    EXAMPLES: 
        1) MakeContiguousMap(False,False,False,collection="20220904-20220905",
                             LonSys='2', FiveMicron='png',
                             Five_obskey='20220904UTa')
        2) MakeContiguousMap(False,False,False,collection="20220904-20220905",
                             LonSys='2', FiveMicron='png',
                             Five_obskey='20220904UTa')     
        3) MakeContiguousMap(False,False,False,collection="20220904-20220905",
                             LonSys='2', FiveMicron='fits',
                             Five_obskey='20220904UTa')   
        4) MakeContiguousMap(False,False,False,collection="20241129-20241129",LonSys='2',
                             FiveMicron=False,Five_obskey='',IRTFdataset='',
                             lats=[60,120],LonLims=[0,360],figsz=[6.0,6.0],ROI=False,
                             variance=False,proj='maps',ctbls=['terrain_r','Blues'])
        !!! Currently FITS version for five micron does not work. I probably
        !!! need to modify read_fits_map_L2_L3.py as a minimum. Then I also
        !!! have to come up with a way to manage get_obs_list.py for 5umfile
        !!! data depending on if it's png or fits                             
    @author: smhil
    """
    import sys
    drive='c:'
    sys.path.append(drive+'/Astronomy/Python Play')
    sys.path.append(drive+'/Astronomy/Python Play/Util_P3')
    sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy')
    sys.path.append('./Services')

    import os
    import pylab as pl
    import numpy as np
    from imageio import imwrite
    from astropy.io import fits
    import scipy.ndimage as ndi
    import RetrievalLibrary as RL
    sys.path.append('./Maps')
    import read_fits_map_L2_L3 as RFM
    import plot_patch as PP
    import make_patch_RGB as MPRGB
    import get_map_collection as gmc
    
    if ctbls[0]=="jet":
        fNH3low=60
        fNH3high=160
        PCldlow=1200
        PCldhigh=2000
    if ctbls[0]=="terrain_r":
        fNH3low=70
        fNH3high=140
        PCldlow=1400
        PCldhigh=2000

    print(lats)
    aspectratio=(LonLims[1]-LonLims[0])/(lats[1]-lats[0])
    if aspectratio==1:
        figsz=[3.0,6.0]
    elif aspectratio==4./3.:
        figsz=[3.5,6.0]
    elif aspectratio==3:
        figsz=[6.0,6.0]
    elif aspectratio==4:
        figsz=[7.05,6.0]
    elif aspectratio==6:
        figsz=[8.5,5.0]
    elif aspectratio==12:
        figsz=[12.,4.0]
    
    
    dataset,lonlims=gmc.get_map_collection(collection)
    ###########################################################################
    # Establish empty arrays for contiguous maps
    ###########################################################################
    outputfNH3=np.zeros([180,360])
    outputPCloud=np.zeros([180,360])
    outputweights=np.zeros([180,360])
    outputRGB=np.zeros([180,360,3])
    
    stackfNH3=np.zeros([180,360])
    stackPCloud=np.zeros([180,360])
    stackweights=np.zeros([180,360])
    stackRGB=np.zeros([180,360,3])
    stackR=np.zeros([180,360])
    stackG=np.zeros([180,360])
    stackB=np.zeros([180,360])
    
    ###########################################################################
    # Loop over observations in each data set and create 3D cubes
    ###########################################################################
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
        
        lonhalfwidth=45
        boxcar=9
        if LonSys=='1':
            ll_0=int(360-(fNH3hdr["CM1"]-lonhalfwidth))
            ll_1=int(360-(fNH3hdr["CM1"]+lonhalfwidth))
            wl_0=int(360-(fNH3hdr["CM1"]-lonhalfwidth+boxcar))
            wl_1=int(360-(fNH3hdr["CM1"]+lonhalfwidth-boxcar))
        elif LonSys=='2':
            ll_0=int(360-(fNH3hdr["CM2"]-lonhalfwidth))
            ll_1=int(360-(fNH3hdr["CM2"]+lonhalfwidth))
            wl_0=int(360-(fNH3hdr["CM2"]-lonhalfwidth+boxcar))
            wl_1=int(360-(fNH3hdr["CM2"]+lonhalfwidth-boxcar))
        
        #print("***************** ll_0x, ll_1x=",ll_0x,ll_1x)
        #!!! Need to fix this so it rolls over 360/0 boundary
        if ll_0>360:
            ll_0=360
        if ll_1<0:
            ll_1=0
        if wl_0>360:
            wl_0=360
        if wl_1<0:
            wl_1=0
        print("***************** ll_0, ll_1=",ll_0,ll_1)
        #print("***************** ll_0x, ll_1x=",ll_0x,ll_1x)
        
        #######################################################################
        # RESET EMPTY FRAMES FOR NEWEST LIMB-CORRECTED REGION
        #######################################################################
        outputfNH3=np.zeros([180,360])
        outputPCloud=np.zeros([180,360])
        outputweights=np.zeros([180,360])
        outputR=np.zeros([180,360])
        outputG=np.zeros([180,360])
        outputB=np.zeros([180,360])
        outputRGB=np.zeros([180,360,3])

        outputfNH3[lats[0]:lats[1],ll_1:ll_0]= \
            TestfNH3[lats[0]:lats[1],ll_1:ll_0]
        outputPCloud[lats[0]:lats[1],ll_1:ll_0]= \
            TestPCloud[lats[0]:lats[1],ll_1:ll_0]
        # Link for uniform_filter1d
        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.uniform_filter1d.html#scipy.ndimage.uniform_filter1d
        outputweights[lats[0]:lats[1],wl_1:wl_0]= \
            outputPCloud[lats[0]:lats[1],wl_1:wl_0]
        outputweights[outputweights>0]=1.0
        outputweights=ndi.uniform_filter1d(outputweights,boxcar,1)
        
        outputR[lats[0]:lats[1],ll_1:ll_0]= \
            RGB[lats[0]:lats[1],ll_1:ll_0,0]
        outputG[lats[0]:lats[1],ll_1:ll_0]= \
            RGB[lats[0]:lats[1],ll_1:ll_0,1]
        outputB[lats[0]:lats[1],ll_1:ll_0]= \
            RGB[lats[0]:lats[1],ll_1:ll_0,2]

        #######################################################################
        #!!! outputRGB is never stacked! It overwrites with latest patch.
        #######################################################################
        outputRGB[lats[0]:lats[1],ll_1:ll_0,:]= \
            RGB[lats[0]:lats[1],ll_1:ll_0,:]
        #######################################################################
        #!!! tempfNH3 and outputfNH3 look to be identical. Same for PCLoud
        #######################################################################
        tempfNH3=np.zeros([180,360])
        tempPCloud=np.zeros([180,360])
        tempfNH3[lats[0]:lats[1],ll_1:ll_0]= \
            TestfNH3[lats[0]:lats[1],ll_1:ll_0]
        tempPCloud[lats[0]:lats[1],ll_1:ll_0]= \
            TestPCloud[lats[0]:lats[1],ll_1:ll_0]

        if First:
            stackfNH3=outputfNH3
            stackfNH3=np.reshape(stackfNH3,(180,360,1))
            stackPCloud=outputPCloud
            stackPCloud=np.reshape(stackPCloud,(180,360,1))
            stackweights=outputweights
            stackweights=np.reshape(stackweights,(180,360,1))
            stackR=outputR
            stackR=np.reshape(stackR,(180,360,1))
            stackG=outputG
            stackG=np.reshape(stackG,(180,360,1))
            stackB=outputB
            stackB=np.reshape(stackB,(180,360,1))
        else:
            tempfNH3=np.reshape(outputfNH3,(180,360,1))
            stackfNH3=np.dstack((stackfNH3,tempfNH3))
            tempPCloud=np.reshape(outputPCloud,(180,360,1))
            stackPCloud=np.dstack((stackPCloud,tempPCloud))
            tempweights=np.reshape(outputweights,(180,360,1))
            stackweights=np.dstack((stackweights,tempweights))
            tempR=np.reshape(outputR,(180,360,1))
            stackR=np.dstack((stackR,tempR))
            tempG=np.reshape(outputG,(180,360,1))
            stackG=np.dstack((stackG,tempG))
            tempB=np.reshape(outputB,(180,360,1))
            stackB=np.dstack((stackB,tempB))
        First=False
        
    ###########################################################################
    # Flatten cubes by taking the mean and ignoring zero values
    #!!!  ADD VARIANCE BELOW!!!
    #!!!  WRITE FITS FILES OF MEAN AND STANDARD DEVIATIONS OF FULL MAPS!!!!
    #
    # Link for taking weighted averages of masked arrays
    # https://numpy.org/doc/2.2/reference/generated/numpy.ma.average.html
    # https://stackoverflow.com/questions/21113384/python-numpy-weighted-average-with-nans
    ###########################################################################
    print(stackfNH3.shape)
    indzf=np.where(stackfNH3==0)
    stackfNH3[indzf]=np.nan
    blendfNH3=np.nanmean(stackfNH3,axis=2)
    
    stackfNH3masked = np.ma.MaskedArray(stackfNH3, mask=np.isnan(stackfNH3))
    blendweightfNH3=np.ma.average(stackfNH3masked, axis=2, weights=stackweights) 
    
    stdvfNH3=np.nanstd(stackfNH3,axis=2)
    fracfNH3=stdvfNH3/blendfNH3
    
    indzP=np.where(stackPCloud==0)
    stackPCloud[indzP]=np.nan
    blendPCloud=np.nanmean(stackPCloud,axis=2)

    stackPCloudmasked = np.ma.MaskedArray(stackPCloud, mask=np.isnan(stackPCloud))
    blendweightPCloud=np.ma.average(stackPCloudmasked, axis=2, weights=stackweights) 
    
    stdvPCloud=np.nanstd(stackPCloud,axis=2)
    fracPCloud=stdvPCloud/blendPCloud

    indzP=np.where(stackR==0)
    stackR[indzP]=np.nan
    blendR=np.nanmean(stackR,axis=2)
    
    indzP=np.where(stackG==0)
    stackG[indzP]=np.nan
    blendG=np.nanmean(stackG,axis=2)
    
    indzP=np.where(stackB==0)
    stackB[indzP]=np.nan
    blendB=np.nanmean(stackB,axis=2)
    
    blendRGB=np.zeros([180,360,3])
    blendRGB[:,:,0]=blendR
    blendRGB[:,:,1]=blendG
    blendRGB[:,:,2]=blendB
    
    print(blendfNH3.shape)

    #6x6 works for 60N-60Z and 360 long
    #8x5 works for 30N-30S and 360 Long
    if FiveMicron=='png' or FiveMicron=='fits':
        rng=[0,1,2,3]
    else:
        rng=[0,1,2]
    fig1,axs1=pl.subplots(rng[-1]+1,1,figsize=(figsz[0],figsz[1]), dpi=150, facecolor="white",
                          sharex=True,sharey=True)
    fig1.suptitle(collection+"\n Average")
    
    for ix in rng:
        axs1[ix].grid(linewidth=0.2)
        axs1[ix].ylim=[-90.,90.]
        axs1[ix].xlim=[0.,360.]
        axs1[ix].set_xticks(np.linspace(450.,0.,31), minor=False)
        xticklabels=np.array(np.mod(np.linspace(450,0,31),360))
        axs1[ix].set_xticklabels(xticklabels.astype(int))
        axs1[ix].set_yticks(np.linspace(-90,90,13), minor=False)
        axs1[ix].tick_params(axis='both', which='major', labelsize=7)
        axs1[ix].set_ylabel("PG Latitude (deg)")
        axs1[ix].set_adjustable('box') 

    print("###############")
    print([360-LonLims[1],360-LonLims[0]])
    fNH3_patch_mb,vn,vx,tx_fNH3=PP.plot_patch(blendweightfNH3,lats,[360-LonLims[1],360-LonLims[0]],
                                     180,180,ctbls[0],
                                     axs1[0],'%3.2f',cont=False,n=11,vn=fNH3low,vx=fNH3high,
                                     cbar_title="")
    #fNH3_patch_mb,vn,vx,tx_fNH3=PP.plot_patch(blendfNH3,lats,[360-LonLims[1],360-LonLims[0]],
    #                                 180,180,ctbls[0],
    #                                 axs1[0],'%3.2f',cont=False,n=11,vn=fNH3low,vx=fNH3high,
    #                                 cbar_title="")
    print("vn,vx,tx_fNH3",vn,vx,tx_fNH3)
    axs1[0].set_title('fNH3 (ppm)',fontsize=10)

    #PCld_patch_mb,vn,vx,tx_fNH3=PP.plot_patch(blendPCloud,lats,[360-LonLims[1],360-LonLims[0]],
    #                                 180,180,ctbls[1],
    #                                 axs1[1],'%3.2f',cont=False,n=5,vn=PCldlow,vx=PCldhigh,
    #                                 cbar_title="",cbar_reverse=True)
    PCld_patch_mb,vn,vx,tx_fNH3=PP.plot_patch(blendweightPCloud,lats,[360-LonLims[1],360-LonLims[0]],
                                     180,180,ctbls[1],
                                     axs1[1],'%3.2f',cont=False,n=5,vn=PCldlow,vx=PCldhigh,
                                     cbar_title="",cbar_reverse=True)
    axs1[1].set_title('PCloud (mbar)',fontsize=10)
    RGBaxs=2

    ###########################################################################
    if variance:
        fig2,axs2=pl.subplots(rng[-1]+1,1,figsize=(figsz[0],figsz[1]), dpi=150, facecolor="white",
                              sharex=True,sharey=True)
        fig2.suptitle(collection+"\n Standard Deviation")
        
        for ix in rng:
            axs2[ix].grid(linewidth=0.2)
            axs2[ix].ylim=[-30.,30.]
            axs2[ix].xlim=[0.,360.]
            axs2[ix].set_xticks(np.linspace(450.,0.,31), minor=False)
            xticklabels=np.array(np.mod(np.linspace(450,0,31),360))
            axs2[ix].set_xticklabels(xticklabels.astype(int))
            axs2[ix].set_yticks(np.linspace(-90,90,13), minor=False)
            axs2[ix].tick_params(axis='both', which='major', labelsize=7)
            axs2[ix].set_ylabel("PG Latitude (deg)")
            axs2[ix].set_adjustable('box') 
    
        print("###############")
        print([360-LonLims[1],360-LonLims[0]])
        #fNH3_patch_mb,vn,vx,tx_fNH3=PP.plot_patch(stdvfNH3,lats,[360-LonLims[1],360-LonLims[0]],
        #                                 180,180,"jet",
        #                                 axs2[0],'%3.2f',cont=False,n=11,vn=0,vx=20,
        #                                 cbar_title="")
        fNH3_patch_mb,vn,vx,tx_fNH3=PP.plot_patch(fracfNH3,lats,[360-LonLims[1],360-LonLims[0]],
                                         180,180,"jet",
                                         axs2[0],'%3.2f',cont=False,n=6,vn=0,vx=0.5,
                                         cbar_title="")
        print("vn,vx,tx_fNH3",vn,vx,tx_fNH3)
        axs2[0].set_title('fNH3 (fractional '+r'$\sigma$'+')',fontsize=10)
    
        #PCld_patch_mb,vn,vx,tx_fNH3=PP.plot_patch(stdvPCloud,lats,[360-LonLims[1],360-LonLims[0]],
        #                                 180,180,"jet",
        #                                 axs2[1],'%3.2f',cont=False,n=5,vn=0,vx=50,
        #                                 cbar_title="")
        PCld_patch_mb,vn,vx,tx_fNH3=PP.plot_patch(fracPCloud,lats,[360-LonLims[1],360-LonLims[0]],
                                         180,180,"jet",
                                         axs2[1],'%3.2f',cont=False,n=6,vn=0,vx=0.25,
                                         cbar_title="")
        axs2[1].set_title('PCloud (fractional '+r'$\sigma$'+')',fontsize=10)

    ###########################################################################
    # READ 5um png file and display patch. Also can be used for other
    # scraped png files from JALPO maps.
    ###########################################################################
    print("FiveMicron",FiveMicron)
    if FiveMicron=='png':
        PCloudhdr,PClouddata,fNH3hdr,fNH3data,sza,eza,Five,Five_CM2,Fivetime= \
                        RFM.read_fits_map_L2_L3(obskey=Five_obskey,LonSys=LonSys,
                                                imagetype="Map",Level="L3",
                                                FiveMicron=FiveMicron)                        

        Five_patch=MPRGB.make_patch_RGB(Five,lats,[360-LonLims[1],360-LonLims[0]],180,180)
        Five4Display=np.power(np.array(Five_patch).astype(float),1.0)
        Five4Display=Five4Display/Five4Display.max()
        show=axs1[2].imshow(Five4Display,
                   extent=[360,0,90-lats[1],
                           90-lats[0]],
                           aspect="equal")
        RGBaxs=3

    ###########################################################################
    # Loop over 5um FITS files to create maps
    ###########################################################################
    if FiveMicron=='fits':
        RGBaxs=3

        pathIRTF_FITS="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/5micron/FITS/output/"
        outputIRTF=np.zeros([180,360])
        stackIRTF=np.zeros([180,360])
   
        print("IRTFdataset=",IRTFdataset)
        FirstIRTF=True
        for IRTFfile in IRTFdataset:
            print("*******obsdate=",IRTFfile)
            
            IRTFhdulist=fits.open(pathIRTF_FITS+IRTFfile)
            IRTFhdulist.info()
            IRTFhdr=IRTFhdulist[0].header
            IRTFdata=IRTFhdulist[0].data
            IRTFhdulist.close()
            
            if LonSys=='1':
                print("LonSys1=",LonSys)
                IRTFroll=IRTFhdr["CM3"]-IRTFhdr["CM1"]
                IRTFdatar=np.roll(IRTFdata,int(IRTFroll),axis=1)
            if LonSys=='2':
                print("LonSys1=",LonSys)
                IRTFroll=IRTFhdr["CM3"]-IRTFhdr["CM2"]
                IRTFdatar=np.roll(IRTFdata,int(IRTFroll),axis=1)

            #lats=[30,150]
            #lats=[60,120]
            #lats=[75,95]
            #ll_0=360-lonlims[obskey][0]
            #ll_1=360-lonlims[obskey][1]
            if LonSys=='1':
                ll_0=int(360-(IRTFhdr["CM1"]-70))
                ll_1=int(360-(IRTFhdr["CM1"]+70))
            elif LonSys=='2':
                ll_0=int(360-(IRTFhdr["CM2"]-70))
                ll_1=int(360-(IRTFhdr["CM2"]+70))
            
            #print("***************** ll_0x, ll_1x=",ll_0x,ll_1x)
            if ll_0>360:
                ll_0=360
            if ll_1<0:
                ll_1=0
            print("***************** ll_0, ll_1=",ll_0,ll_1)
            #print("***************** ll_0x, ll_1x=",ll_0x,ll_1x)
            
            outputIRTF[lats[0]:lats[1],ll_1:ll_0]= \
                IRTFdatar[lats[0]:lats[1],ll_1:ll_0]
                
            tempIRTF=np.zeros([180,360])
            tempIRTF[lats[0]:lats[1],ll_1:ll_0]= \
                IRTFdatar[lats[0]:lats[1],ll_1:ll_0]
    
            if First:
                stackIRTF=tempIRTF
                stackIRTF=np.reshape(stackIRTF,(180,360,1))
            else:
                tempIRTF=np.reshape(tempIRTF,(180,360,1))
                stackIRTF=np.dstack((stackIRTF,tempIRTF))
            FirstIRTF=False
            
        ###########################################################################
        # Flatten cubes by taking the mean and ignoring zero values
        ###########################################################################
        print(stackIRTF.shape)
        indzf=np.where(stackIRTF==0)
        stackIRTF[indzf]=np.nan
        blendIRTF=np.nanmean(stackIRTF,axis=2)
        print(blendIRTF.shape)

        IRTF_patch_mb,vn,vx,tx_fNH3=PP.plot_patch(np.log10(blendIRTF),lats,[360-LonLims[1],360-LonLims[0]],
                                         180,180,"gist_heat",
                                         axs1[2],'%3.2f',cont=False,n=5,vn=1,vx=3.0)
        
    ###########################################################################
    # Done with IRTF branch, now, finally, do the RGB context image
    ###########################################################################

    RGB_patch=MPRGB.make_patch_RGB(blendRGB,lats,[360-LonLims[1],360-LonLims[0]],180,180)
    RGB4Display=np.power(np.array(RGB_patch).astype(float),1.0)
    show=axs1[RGBaxs].imshow(RGB4Display,
               extent=[LonLims[1],LonLims[0],90-lats[1],
                       90-lats[0]],
                       aspect="equal")
    axs1[RGBaxs].set_title('RGB Context',fontsize=10)

    if variance:
        show=axs2[RGBaxs].imshow(RGB4Display,
                   extent=[LonLims[1],LonLims[0],90-lats[1],
                           90-lats[0]],
                           aspect="equal")
        axs2[RGBaxs].set_title('RGB Context',fontsize=10)

    #temp=RL.make_contours_CH4_patch(axs1[2],outputfNH3,[30.,150.],[0.,360.],
    #                       tx_fNH3,frmt='%3.0f',clr='k')
    axs1[RGBaxs].tick_params(axis='both', which='major', labelsize=7)
    axs1[RGBaxs].set_xlabel("Sys. "+LonSys+" Longitude (deg)",fontsize=9)
    axs1[RGBaxs].set_ylabel("PG Latitude (deg)")

    box1 = axs1[RGBaxs].get_position()
    if variance:
        axs2[RGBaxs].tick_params(axis='both', which='major', labelsize=7)
        axs2[RGBaxs].set_xlabel("Sys. "+LonSys+" Longitude (deg)",fontsize=9)
        axs2[RGBaxs].set_ylabel("PG Latitude (deg)")
        box2 = axs2[RGBaxs].get_position()
    
    #Aspect Ratio Customization
    if aspectratio==1:
        fig1.subplots_adjust(left=0.21, bottom=0.07, right=0.79, top=0.88,
                    wspace=0.0, hspace=0.2)     
    if aspectratio==4./3.:
        fig1.subplots_adjust(left=0.20, bottom=0.07, right=0.825, top=0.86,
                    wspace=0.0, hspace=0.2)     
    if aspectratio==3:
        fig1.subplots_adjust(left=0.095, bottom=0.07, right=0.90, top=0.88,
                    wspace=0.25, hspace=0.2)     
    if aspectratio==4:
        fig1.subplots_adjust(left=0.035, bottom=0.07, right=0.94, top=0.88,
                    wspace=0.25, hspace=0.2)     
    if aspectratio==6:
        fig1.subplots_adjust(left=0.04, bottom=0.08, right=0.94, top=0.865,
                    wspace=0.0, hspace=0.2)     
    if aspectratio==12:
        fig1.subplots_adjust(left=0.035, bottom=0.1, right=0.935, top=0.84,
                    wspace=0.0, hspace=0.3)     
    
    if variance:
        if aspectratio==1:
            fig2.subplots_adjust(left=0.21, bottom=0.07, right=0.79, top=0.88,
                        wspace=0.0, hspace=0.2)     
        if aspectratio==4./3.:
            fig2.subplots_adjust(left=0.20, bottom=0.07, right=0.825, top=0.86,
                        wspace=0.0, hspace=0.2)     
        if aspectratio==3:
            fig2.subplots_adjust(left=0.095, bottom=0.07, right=0.90, top=0.88,
                        wspace=0.25, hspace=0.2)     
        if aspectratio==4:
            fig2.subplots_adjust(left=0.035, bottom=0.07, right=0.94, top=0.88,
                        wspace=0.25, hspace=0.2)     
        if aspectratio==6:
            fig2.subplots_adjust(left=0.04, bottom=0.08, right=0.94, top=0.865,
                        wspace=0.0, hspace=0.2)     
        if aspectratio==12:
            fig2.subplots_adjust(left=0.035, bottom=0.1, right=0.935, top=0.84,
                        wspace=0.0, hspace=0.3)     


    if FiveMicron=='png':
        
        if FiveMicron=='png':
            axs1[2].set_title(Fivetime[0:10]+" (5 micron)",fontsize=8)
        axs1[3].set_title("RGB Context Image",fontsize=8)

        box2 = axs1[2].get_position()
        axs1[2].set_position([box2.x0-0.0, box2.y0-0.0, box2.width * 0.95, box2.height * 1.07])
        box3 = axs1[3].get_position()
        axs1[3].set_position([box3.x0-0.01, box3.y0-0.0, box3.width * 1.025, box3.height * 1.07])


    if FiveMicron=='fits':
        box3 = axs1[3].get_position()
        axs1[3].set_position([box3.x0-0.01, box3.y0-0.0, box3.width * 1.025, box3.height * 1.07])

    if ROI:
        for R in ROI:
            for iax in range(0,3):
                axs1[iax].plot(np.array([ROI[R][2]+ROI[R][3],ROI[R][2]-ROI[R][3],
                              ROI[R][2]-ROI[R][3],ROI[R][2]+ROI[R][3],
                              ROI[R][2]+ROI[R][3]]),
                              90.-np.array([ROI[R][0],ROI[R][0],ROI[R][1],
                              ROI[R][1],ROI[R][0]]))
    if variance and ROI:
        for R in ROI:
            for iax in range(0,3):
                axs2[iax].plot(np.array([ROI[R][2]+ROI[R][3],ROI[R][2]-ROI[R][3],
                              ROI[R][2]-ROI[R][3],ROI[R][2]+ROI[R][3],
                              ROI[R][2]+ROI[R][3]]),
                              90.-np.array([ROI[R][0],ROI[R][0],ROI[R][1],
                              ROI[R][1],ROI[R][0]]))


    pathmapplots="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Studies/"+proj+"/"
    fig1.savefig(pathmapplots+collection+" Mean Sys"+LonSys+" map.png",dpi=300)
    if variance:
        fig2.savefig(pathmapplots+collection+" Stdv Sys"+LonSys+" map.png",dpi=300)
        
    ###########################################################################
    # Create Stack Plot subplots on the axes objects passed into the procedure
    ###########################################################################
    if axNH3!=False:
        #lats=[100,120]
        #lats=[80,100]
        fNH3_patch_mb,vn,vx,tx_fNH3=PP.plot_patch(blendweightfNH3,lats,[0,360],
                                         180,180,ctbls[0],
                                         axNH3,'%3.2f',cbarplot=False,cont=False,n=11,vn=fNH3low,vx=fNH3high)
        axNH3.set_ylabel(collection.replace('-','\n'),rotation='horizontal',fontsize=6)
        axNH3.yaxis.set_label_coords(-0.05,0.15)
        axNH3.tick_params('y', labelleft=False)
        axNH3.tick_params('x', labelsize=8)
    
        PCld_patch_mb,vn,vx,tx_fNH3=PP.plot_patch(blendweightPCloud,lats,[0,360],
                                         180,180,ctbls[1],
                                         axCH4,'%3.2f',cbarplot=False,cont=False,n=5,vn=PCldlow,vx=PCldhigh)
        axCH4.set_ylabel(collection,rotation='horizontal',fontsize=6)
        axCH4.set_ylabel(collection.replace('-','\n'),rotation='horizontal',fontsize=6)
        axCH4.yaxis.set_label_coords(-0.05,0.15)
        axCH4.tick_params('y', labelleft=False)
        axCH4.tick_params('x', labelsize=8)
    
        RGB_patch=MPRGB.make_patch_RGB(blendRGB,lats,[0,360],180,180)
        RGB4Display=np.power(np.array(RGB_patch).astype(float),1.0)
        #RGB4Display=RGB4Display/RGB4Display.max()
        show=axRGB.imshow(RGB4Display,
                   extent=[360,0,90-lats[1],
                           90-lats[0]],
                           aspect="equal")
        axRGB.set_ylabel(collection.replace('-','\n'),rotation='horizontal',fontsize=6)
        axRGB.tick_params('y', labelleft=False)
        axRGB.yaxis.set_label_coords(-0.05,0.15)
        axRGB.tick_params('x', labelsize=8)

    print(aspectratio)
    return(lats,blendPCloud,blendfNH3,blendRGB)
    #return(fig1,axs1)