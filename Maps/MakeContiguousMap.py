def MakeContiguousMap(axNH3,axCH4,axRGB,collection="20230815-20230818",LonSys='2',
                      FiveMicron=False,Five_obskey='',IRTFdataset='',LonLims=[0,360]):
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
    import RetrievalLibrary as RL
    sys.path.append('./Maps')
    import read_fits_map_L2_L3 as RFM
    import plot_patch as PP
    import make_patch_RGB as MPRGB
    import get_map_collection as gmc
    
    dataset,lonlims=gmc.get_map_collection(collection)
    ###########################################################################
    # Establish empty arrays for contiguous maps
    ###########################################################################
    outputfNH3=np.zeros([180,360])
    outputPCloud=np.zeros([180,360])
    outputRGB=np.zeros([180,360,3])
    
    stackfNH3=np.zeros([180,360])
    stackPCloud=np.zeros([180,360])
    stackRGB=np.zeros([180,360,3])
    
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
        
        #lats=[30,150]
        lats=[60,120]
        lats=[75,95]
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
        
    ###########################################################################
    # Flatten cubes by taking the mean and ignoring zero values
    ###########################################################################
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
    if FiveMicron=='png' or FiveMicron=='fits':
        rng=[0,1,2,3]
    else:
        rng=[0,1,2]
    fig1,axs1=pl.subplots(rng[-1]+1,1,figsize=(8.0,6), dpi=150, facecolor="white",
                          sharex=True,sharey=True)
    fig1.suptitle(collection)

    for ix in rng:
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

    print("###############")
    print([360-LonLims[1],360-LonLims[0]])
    fNH3_patch_mb,vn,vx,tx_fNH3=PP.plot_patch(blendfNH3,lats,[360-LonLims[1],360-LonLims[0]],
                                     180,180,"jet",
                                     axs1[0],'%3.2f',cont=False,n=11,vn=60,vx=160)
    print("vn,vx,tx_fNH3",vn,vx,tx_fNH3)
    axs1[0].set_title('fNH3 (ppm)',fontsize=10)

    PCld_patch_mb,vn,vx,tx_fNH3=PP.plot_patch(blendPCloud,lats,[360-LonLims[1],360-LonLims[0]],
                                     180,180,"jet",
                                     axs1[1],'%3.2f',cont=False,n=5,vn=600,vx=1000)
    axs1[1].set_title('PCloud (mbar)',fontsize=10)
    RGBaxs=2

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
            lats=[75,95]
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


    RGB_patch=MPRGB.make_patch_RGB(outputRGB,lats,[140,340],180,180)
    RGB4Display=np.power(np.array(RGB_patch).astype(float),1.0)
    RGB4Display=RGB4Display/RGB4Display.max()
    show=axs1[RGBaxs].imshow(RGB4Display,
               extent=[LonLims[1],LonLims[0],90-lats[1],
                       90-lats[0]],
                       aspect="equal")
    #temp=RL.make_contours_CH4_patch(axs1[2],outputfNH3,[30.,150.],[0.,360.],
    #                       tx_fNH3,frmt='%3.0f',clr='k')
    axs1[RGBaxs].tick_params(axis='both', which='major', labelsize=7)
    axs1[RGBaxs].set_xlabel("Sys. "+LonSys+" Longitude (deg)",fontsize=9)
    axs1[RGBaxs].set_ylabel("PG Latitude (deg)")

    box = axs1[RGBaxs].get_position()
    
    #works for 6x6 and 60N-60S and 360 long
    fig1.subplots_adjust(left=0.10, bottom=0.05, right=0.90, top=0.93,
                wspace=0.25, hspace=0.05)     
    axs1[RGBaxs].set_position([box.x0-0.075, box.y0-0.04, box.width * 1.11, box.height * 1.11])

    #works for 8x5 and 30N-30S and 360 long
    fig1.subplots_adjust(left=0.07, bottom=0.06, right=0.94, top=0.92,
                wspace=0.25, hspace=0.08)     
    axs1[RGBaxs].set_position([box.x0-0.055, box.y0-0.04, box.width * 1.07, box.height * 1.07])

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

    pathmapplots="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Studies/maps/"
    fig1.savefig(pathmapplots+collection+"Sys"+LonSys+" map.png",dpi=300)
    
    
    ###########################################################################
    # Create Stack Plot subplots on the axes objects passed into the procedure
    ###########################################################################
    if axNH3!=False:
        #lats=[100,120]
        lats=[80,100]
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

    return(lats)