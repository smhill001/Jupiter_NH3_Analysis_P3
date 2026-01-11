def ReadContiguousMap(LonSys,collection="20220904-20220905",
                      collectionIRTF="20240205-20240205",
                      collection889CH4="20240131-20240201"):
    """
    Retrieve continguous maps from FITS files for plotting and analysis.
    
    !!!Development notes: 
    1) Currently assumes that fNH3, PCld, and RGB 'optical' data is always
    provided. It could be useful to make that optional like the IRTF and 889-nm
    data in order to generalize this function.
    2) A helper function to read the m x n FITS files and perhaps the RGB file
    could dramatically shorten this function and make it more readable.
    

    Parameters
    ----------
    LonSys : String
        Indicates desired longitude system of map plot
    collection : String, optional
        The map collection label for contiguous map FITS files of fNH3, PCld,
        and visual context
    collectionIRTF : String, optional
        The map collection label for contiguous map FITS file of 5-micron data
    collection889CH4 : String, optional
        The map collection label for contiguous map FITS file of 889-nm data

    Returns
    -------
    fNH3data : float array
        Cylindrically mapped mean ammonia mole fraction at each lat-lon in ppm.
    fNH3stdv : float array
        Cylindrically mapped fNH3 standard deviation at each lat-lon in ppm.
    fNH3frac : float array
        Cylindrically mapped fractional uncertainty in fNH3, e.g., 
        fNH3stdv/fNH3data
    fNH3time : float array
        Cylindrically mapped mean fNH3 observation time at each lat-lon in JD(?).
    PClddata : float array
        Cylindrically mapped mean effective cloud pressure at each lat-lon in 
        mbar.
    PCldstdv : float array
        Cylindrically mapped PCld standard deviation at each lat-lon in mbar.
    PCldfrac : float array
        Cylindrically mapped fractional uncertainty in PCld, e.g., 
        PCldstdv/PClddata
    PCldtime : float array
        Cylindrically mapped mean PCld observation time at each lat-lon in JD(?).
    RGBdata : TYPE
        DESCRIPTION.
    RGBstdv : TYPE
        DESCRIPTION.
    RGBfrac : TYPE
        DESCRIPTION.
    RGBtime : TYPE
        DESCRIPTION.
    IRTFdata : float array
        Cylindrically mapped mean effective cloud pressure at each lat-lon in 
        uncalibrated units.
    IRTFstdv : float array
        Cylindrically mapped IRTF standard deviation at each lat-lon in mbar.
    IRTFfrac : float array
        Cylindrically mapped fractional uncertainty in IRTF, e.g., 
        IRTFstdv/IRTFdata
    IRTFtime : float array
        Cylindrically mapped mean IRTF observation time at each lat-lon in JD(?).
    JALPOdata : TYPE
        DESCRIPTION.

    """
    from astropy.io import fits
    import os
    import numpy as np
    
    print("############ IN READ MAP,collectionIRTF=",collectionIRTF)
    
    cm_key = {
        '1': 'CM1',
        '2': 'CM2',
        '3': 'CM3'#,
        #'subobs':180
    }.get(LonSys)
    
    ctype1_key = {
        'Sys. 1 Longitude':'CM1',
        'Sys. 2 Longitude':'CM2',
        'Sys. 3 Longitude':'CM3',
        }
    
    ###########################################################################
    # Identify files to read
    ###########################################################################
    pathinp="C:/Astronomy/Projects/SAS 2021 Ammonia/Data/L4 FITS (cont maps)/"
    contents = os.listdir(pathinp)
    files_in_directory = [item for item in contents if os.path.isfile(os.path.join(pathinp, item))]
    
    fits_in_directory = [item for item in files_in_directory \
                         if '.fits' in item ]  
    collection_in_directory = [item for item in fits_in_directory \
                         if collection in item]  
        
    ###########################################################################
    # Read fNH3 data and 'roll' to desired longitude system (LonSys)
    ###########################################################################
    fNH3file = [item for item in collection_in_directory \
                         if 'L4fNH3' in item]  
    if len(fNH3file)>0:
        fNH3hdulist=fits.open(pathinp+fNH3file[0])
        fNH3hdulist.info()
        fNH3hdr=fNH3hdulist[0].header
        roll=int(fNH3hdr[ctype1_key[fNH3hdr['CTYPE1']]]-fNH3hdr[cm_key])
        fNH3data=np.roll(fNH3hdulist[0].data,roll,axis=1)
        fNH3stdv=np.roll(fNH3hdulist[1].data,roll,axis=1)
        fNH3frac=np.roll(fNH3hdulist[2].data,roll,axis=1)
        fNH3time=np.roll(fNH3hdulist[3].data,roll,axis=1)
        fNH3hdulist.close()       
        
    ###########################################################################
    # Read PCld data and 'roll' to desired longitude system (LonSys)
    ###########################################################################
    PCldfile = [item for item in collection_in_directory \
                         if 'L4PCld' in item]      
    if len(PCldfile)>0:
        PCldhdulist=fits.open(pathinp+PCldfile[0])
        PCldhdulist.info()
        PCldhdr=PCldhdulist[0].header
        roll=int(PCldhdr[ctype1_key[PCldhdr['CTYPE1']]]-PCldhdr[cm_key])
        PClddata=np.roll(PCldhdulist[0].data,roll,axis=1)
        PCldstdv=np.roll(PCldhdulist[1].data,roll,axis=1)
        PCldfrac=np.roll(PCldhdulist[2].data,roll,axis=1)
        PCldtime=np.roll(PCldhdulist[3].data,roll,axis=1)
        PCldhdulist.close()
        
    ###########################################################################
    # Read RGB data and 'roll' to desired longitude system (LonSys)
    ###########################################################################
    RGBfile = [item for item in collection_in_directory \
                         if 'L4RGB' in item]      
    if len(RGBfile)>0:
        RGBhdulist=fits.open(pathinp+RGBfile[0])
        RGBhdulist.info()
        RGBhdr=RGBhdulist[0].header
        roll=int(RGBhdr[ctype1_key[RGBhdr['CTYPE1']]]-RGBhdr[cm_key])
        RGBdata=np.roll(RGBhdulist[0].data,roll,axis=1)
        RGBstdv=np.roll(RGBhdulist[1].data,roll,axis=1)
        RGBfrac=np.roll(RGBhdulist[2].data,roll,axis=1)
        RGBtime=np.roll(RGBhdulist[3].data,roll,axis=1)
        RGBhdulist.close()

    ###########################################################################
    # Read IRTF data and 'roll' to desired longitude system (LonSys)
    ###########################################################################
    if collectionIRTF:        
        collectionIRTF_in_directory = [item for item in fits_in_directory \
                             if collectionIRTF in item]  
        IRTFfile = [item for item in collectionIRTF_in_directory \
                             if 'L4IRTF' in item]      
        if len(IRTFfile)>0:
            print("############### IRTF Map=",IRTFfile)
            IRTFhdulist=fits.open(pathinp+IRTFfile[0])
            IRTFhdulist.info()
            IRTFhdr=IRTFhdulist[0].header
            roll=int(IRTFhdr['CM3']-IRTFhdr[cm_key])
            IRTFdata=np.roll(IRTFhdulist[0].data,roll,axis=1)
            IRTFstdv=np.roll(IRTFhdulist[1].data,roll,axis=1)
            IRTFfrac=np.roll(IRTFhdulist[2].data,roll,axis=1)
            IRTFtime=np.roll(IRTFhdulist[3].data,roll,axis=1)
            IRTFhdulist.close()
        else:
            IRTFdata=np.zeros((180,360))
            IRTFstdv=np.zeros((180,360))
            IRTFfrac=np.zeros((180,360))
            IRTFtime=np.zeros((180,360))       
    else:
        IRTFdata=np.zeros((180,360))
        IRTFstdv=np.zeros((180,360))
        IRTFfrac=np.zeros((180,360))
        IRTFtime=np.zeros((180,360))

    ###########################################################################
    # Read 889-nm JALPO data and 'roll' to desired longitude system (LonSys)
    ###########################################################################
    if collection889CH4:
        
        collection889CH4_in_directory = [item for item in fits_in_directory \
                             if collection889CH4 in item]  
        JALPOfile = [item for item in collection889CH4_in_directory \
                             if 'L4889CH4' in item]      
        if len(JALPOfile)>0:
            JALPOhdulist=fits.open(pathinp+JALPOfile[0])
            JALPOhdulist.info()
            JALPOhdr=JALPOhdulist[0].header
            print("########### ",JALPOhdr['DATE-OBS'])
            print("########### ",JALPOhdr['CM3'])
            print("########### ",cm_key)
    
            roll=int(JALPOhdr['CM3'])-int(JALPOhdr[cm_key])
            JALPOdata=np.roll(JALPOhdulist[0].data,roll,axis=1)
            JALPOhdulist.close()
        else:
            JALPOdata=np.zeros((180,360))
    else:
        JALPOdata=np.zeros((180,360))
    
    return fNH3data,fNH3stdv,fNH3frac,fNH3time,\
        PClddata,PCldstdv,PCldfrac,PCldtime, \
        RGBdata,RGBstdv,RGBfrac,RGBtime, \
        IRTFdata,IRTFstdv,IRTFfrac,IRTFtime, \
        JALPOdata

def plot_optical_maps(LonLims,lats,collection,blendweightfNH3,
                      blendweightPCloud,blendweightTime,blendRGBweight,
                      ctbls,variance,fracfNH3,fracPCloud,RGBaxs):
    """
    

    Parameters
    ----------
    LonLims : TYPE
        DESCRIPTION.
    lats : TYPE
        DESCRIPTION.
    collection : TYPE
        DESCRIPTION.
    blendweightfNH3 : TYPE
        DESCRIPTION.
    blendweightPCloud : TYPE
        DESCRIPTION.
    blendweightTime : TYPE
        DESCRIPTION.
    blendRGBweight : TYPE
        DESCRIPTION.
    ctbls : TYPE
        DESCRIPTION.
    variance : TYPE
        DESCRIPTION.
    fracfNH3 : TYPE
        DESCRIPTION.
    fracPCloud : TYPE
        DESCRIPTION.
    RGBaxs : TYPE
        DESCRIPTION.

    Returns
    -------
    fig1 : TYPE
        DESCRIPTION.
    fig2 : TYPE
        DESCRIPTION.
    axs1 : TYPE
        DESCRIPTION.
    axs2 : TYPE
        DESCRIPTION.
    fNH3_patch_mb : TYPE
        DESCRIPTION.
    PCld_patch_mb : TYPE
        DESCRIPTION.
    RGB_patch : TYPE
        DESCRIPTION.
    tx_fNH3 : TYPE
        DESCRIPTION.
    tx_PCld : TYPE
        DESCRIPTION.
    fNH3low : TYPE
        DESCRIPTION.
    fNH3high : TYPE
        DESCRIPTION.
    PCldlow : TYPE
        DESCRIPTION.
    PCldhigh : TYPE
        DESCRIPTION.
    blendweightTime_patch : TYPE
        DESCRIPTION.
    aspectratio : TYPE
        DESCRIPTION.
    RGB4Display : TYPE
        DESCRIPTION.
    RGBaxs : TYPE
        DESCRIPTION.

    """
    import pylab as pl
    import numpy as np
    import plot_patch as PP
    import make_patch as MP
    
    ###########################################################################
    #BEGIN PLOTTING
    ###########################################################################
    # Determine figure size (inches) based on aspect ratio of data set
    aspectratio=(LonLims[1]-LonLims[0])/(lats[1]-lats[0])
    aspect_ratio_map = {
                        1:     [3.0, 6.0],
                        4/3.:  [3.5, 6.0],
                        2:     [4.5, 6.0],
                        3:     [6.0, 6.0], #e.g. 120x360
                        4:     [7.05, 6.0],
                        6:     [8.5, 5.0],
                        9:     [10, 4.5],
                        12:    [12.0, 4.0] #e.g. 30x360
                        }
    figsz = aspect_ratio_map[aspectratio]
    if RGBaxs==4:
        figsz[1]=figsz[1]*1.36
        
    ###########################################################################
    # Set up default ranges for clouds and ammonia. (If generalized for L2 data
    # will need to add EW ranges for NH3 and CH4). This code could be
    # simplified if I retire the "jet" color table.
    ctbl_settings = {
                    "jet": (70, 140, 1200, 2000),
                    "terrain_r": (60, 160, 1400, 2200),
                    "gray": (60, 160, 1400, 2200)
                }
                
    if ctbls[0] in ctbl_settings:
        fNH3low, fNH3high, PCldlow, PCldhigh = ctbl_settings[ctbls[0]]

    #Set up plots to be 3x1 or 4x1 if there's a 5 micron file
    rng=np.arange(0,RGBaxs+1)

    fig1,axs1=pl.subplots(rng[-1]+1,1,figsize=(figsz[0],figsz[1]), dpi=150, facecolor="white",
                          sharex=True,sharey=True)
    fig1.suptitle(collection+"\n Average")
    
    #Set up plot axes etc.
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

    ###########################################################################    
    #NH3 patch plot    
    print("###############")
    print([360-LonLims[1],360-LonLims[0]])
    fNH3_patch_mb=MP.make_patch(blendweightfNH3,lats,[360-LonLims[1],360-LonLims[0]],
                                     180,180)
    fNH3_patch_mb,vn,vx,tx_fNH3=PP.plot_patch(fNH3_patch_mb,lats,[360-LonLims[1],360-LonLims[0]],
                                     180,180,ctbls[0],
                                     axs1[0],'%3.2f',n=6,
                                     vn=fNH3low,
                                     vx=fNH3high,
                                     cbar_title="")
    
    axs1[0].set_title('fNH3 (ppm)',fontsize=10)
    ###########################################################################
    # Add 6 degree FWHM resolution circle to NH3 plot   
    theta = np.linspace(0, 2*np.pi, 100) 
    r=3 #degrees radius seeing circle (FWHM)
    lonres=LonLims[1]-1.5*r
    print("############################# LonLims[1]-LonLims[0]=",LonLims[1]-LonLims[0])
    x = lonres+r*np.cos(theta)
    y = (90-lats[0])-1.5*r+r*np.sin(theta)
    axs1[0].plot(x,y,'k',clip_on=False)

    ###########################################################################
    #PCloud patch plot
    PCld_patch_mb=MP.make_patch(blendweightPCloud,lats,
                                [360-LonLims[1],360-LonLims[0]],180,180)
    PCld_patch_mb,vn,vx,tx_PCld=PP.plot_patch(PCld_patch_mb,lats,
                                              [360-LonLims[1],360-LonLims[0]],
                                              180,180,ctbls[1],
                                              axs1[1],'%3.2f',n=5,
                                              vn=PCldlow,
                                              vx=PCldhigh,
                                              cbar_title="",cbar_reverse=True)
    axs1[1].set_title('PCloud (mbar)',fontsize=10)
    
    blendweightTime_patch=MP.make_patch(blendweightTime,lats,
                                        [360-LonLims[1],360-LonLims[0]],
                                        180,180)

    ###########################################################################
    #Supports "Find local extrema - Without the NaN it messes up the graphics
    fNH3_patch_mb[np.where(fNH3_patch_mb==0)]=np.nan
    PCld_patch_mb[np.where(PCld_patch_mb==0)]=np.nan


    ###########################################################################
    # Done with IRTF branch, now, finally, do the RGB context image
    ###########################################################################
    RGB_patch=MP.make_patch(blendRGBweight,lats,[360-LonLims[1],360-LonLims[0]],180,180)
    RGB4Display=np.power(np.array(RGB_patch).astype(float),1.0)
    show=axs1[RGBaxs].imshow(RGB4Display,
               extent=[LonLims[1],LonLims[0],90-lats[1],
                       90-lats[0]],
                       aspect="equal")

    ###########################################################################
    #PLOT variance of NH3 and PCld blended maps
    ###########################################################################
    if variance:
        fig2,axs2=pl.subplots(rng[-1]+1,1,figsize=(figsz[0],figsz[1]), dpi=150, facecolor="white",
                              sharex=True,sharey=True)
        fig2.suptitle(collection+"\n Standard Deviation")
        
        for ix in rng:
            axs2[ix].grid(linewidth=0.2)
            axs2[ix].ylim=[-90.,90.]
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
        fracfNH3_patch_mb=MP.make_patch(fracfNH3,lats,[360-LonLims[1],360-LonLims[0]],
                                         180,180)

        fracfNH3_patch_mb,vn,vx,tx_fNH3frac=PP.plot_patch(fracfNH3_patch_mb,lats,[360-LonLims[1],360-LonLims[0]],
                                         180,180,"jet",
                                         axs2[0],'%3.2f',n=6,vn=0,vx=0.5,
                                         cbar_title="")
        print("vn,vx,tx_fNH3",vn,vx,tx_fNH3)
        axs2[0].set_title('fNH3 (fractional '+r'$\sigma$'+')',fontsize=10)
    
        fracPCld_patch_mb=MP.make_patch(fracPCloud,lats,[360-LonLims[1],360-LonLims[0]],
                                        180,180)

        fracPCld_patch_mb,vn,vx,tx_PCldfrac=PP.plot_patch(fracPCld_patch_mb,lats,[360-LonLims[1],360-LonLims[0]],
                                         180,180,"jet",
                                         axs2[1],'%3.2f',n=6,vn=0,vx=0.25,
                                         cbar_title="")
        axs2[1].set_title('PCloud (fractional '+r'$\sigma$'+')',fontsize=10)
        
    else:
        fig2,axs2=False,False
        
    return fig1,fig2,axs1,axs2,fNH3_patch_mb,PCld_patch_mb,RGB_patch,tx_fNH3,tx_PCld,\
        fNH3low,fNH3high,PCldlow,PCldhigh,blendweightTime_patch,aspectratio,\
        RGB4Display,RGBaxs
            
def plot_5u_map(blendweightIRTF,LonLims,lats,axs1,axs2,variance,fracIRTF,IRTFaxs):
    import numpy as np
    import make_patch as MP
    import plot_patch as PP
    ###########################################################################
    #PLOT IRTF data
    ###########################################################################
    print("IN IRTF PLOT ######################")
    IRTF_patch=MP.make_patch(np.log10(blendweightIRTF,where=blendweightIRTF > 0),lats,
                             [360-LonLims[1],360-LonLims[0]],180,180)

    IRTF_patch,vn,vx,tx_IRTF=PP.plot_patch(IRTF_patch,lats,[360-LonLims[1],360-LonLims[0]],
                                     180,180,"gist_heat",
                                     axs1[IRTFaxs],'%3.2f',n=5,vn=1.5,vx=3.5,
                                     cbar_title="test")
    axs1[IRTFaxs].set_title('IRTF',fontsize=10)

    if variance:
        fracIRTF_patch=MP.make_patch(fracIRTF,lats,[360-LonLims[1],360-LonLims[0]],
                                        180,180)

        fracIRTF_patch,vn,vx,tx_PCld=PP.plot_patch(fracIRTF_patch,lats,[360-LonLims[1],360-LonLims[0]],
                                         180,180,"jet",
                                         axs2[IRTFaxs],'%3.2f',n=6,vn=0,vx=0.50,
                                         cbar_title="")
        axs2[IRTFaxs].set_title('IRTF',fontsize=10)
        
    return

def plot_JALPO_map(blendweightCH4889,LonLims,lats,axs1,axs2,variance,CH4889axs):
    import numpy as np
    import make_patch as MP
    import plot_patch as PP
    ###########################################################################
    #PLOT IRTF data
    ###########################################################################
    CH4889_patch=MP.make_patch(blendweightCH4889,lats,
                               [360-LonLims[1],360-LonLims[0]],180,180)

    CH4889_patch,vn,vx,tx_889=PP.plot_patch(CH4889_patch,lats,[360-LonLims[1],360-LonLims[0]],
                                     180,180,"gray",axs1[CH4889axs],'%3.2f',
                                     n=5,vn=1,vx=256.0,cbar_title="test")
    
    axs1[CH4889axs].set_title('889CH4',fontsize=10)

    if variance:
        frac889_patch=np.zeros(CH4889_patch.shape)

        frac889_patch,vn,vx,tx_PCld=PP.plot_patch(frac889_patch,lats,
                                         [360-LonLims[1],360-LonLims[0]],
                                         180,180,"jet",axs2[CH4889axs],'%3.2f',
                                         n=6,vn=0,vx=0.50,cbar_title="")
        axs2[CH4889axs].set_title('889CH4',fontsize=10)
        
    return

def make_bare_map(blendweight,ctbl,low,high,pathmapplots,collection,LonSys,env_data_type):
    """
    Make a PNG bare map from 0-360 longitude and -90 to +90 latitude for 
    projection and animation using WinJUPOS.

    Parameters
    ----------
    blendweight : TYPE
        DESCRIPTION.
    ctbl : TYPE
        DESCRIPTION.
    low : TYPE
        DESCRIPTION.
    high : TYPE
        DESCRIPTION.
    pathmapplots : TYPE
        DESCRIPTION.
    collection : TYPE
        DESCRIPTION.
    LonSys : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    import pylab as pl
    figpngfNH3,axspngfNH3=pl.subplots(1,figsize=(2.4,1.2), dpi=150, facecolor="black")
    axspngfNH3 = figpngfNH3.add_axes([0, 0, 1, 1],facecolor='black')
    axspngfNH3.imshow(blendweight,ctbl,vmin=low,vmax=high)
    axspngfNH3.axis('off')
    figpngfNH3.patch.set_facecolor('black')
    figpngfNH3.savefig(pathmapplots+collection+" "+env_data_type+" Mean Sys"+LonSys+" baremap.png",dpi=150)
    pl.show()
    
def L4_Jup_Map_Plot(collection="20240129-20240202",IRTFcollection='20240205-20240205',
                    CH4889collection="20240131-20240201",LonSys='1',
                      lats=[0,180],LonLims=[0,360],figsz=[6.0,6.0],ROI=False,
                      variance=False,localmax=False,segment=False,
                      proj='AGU2025',ctbls=['terrain_r','Blues'],
                      cont=False,bare_maps=False,cb=False,
                      axNH3=False,axCH4=False,axRGB=False,axIRTF=False,
                      axCH4889=False,axNH3vCloud=False,axNH3vIRTF=False,
                      axIRTFvPCld=False,ax889vPCld=False,counter=0,countmax=0):
    """
    Main program for making plots of blended, contiguous maps, including
    analysis tools and overplotting

    Parameters
    ----------
    collection : TYPE, optional
        DESCRIPTION. The default is "20240129-20240202".
    IRTFcollection : TYPE, optional
        DESCRIPTION. The default is '20240205-20240205'.
    CH4889collection : TYPE, optional
        DESCRIPTION. The default is "20240131-20240201".
    LonSys : TYPE, optional
        DESCRIPTION. The default is '1'.
    lats : TYPE, optional
        DESCRIPTION. The default is [0,180].
    LonLims : TYPE, optional
        DESCRIPTION. The default is [0,360].
    figsz : TYPE, optional
        DESCRIPTION. The default is [6.0,6.0].
    ROI : TYPE, optional
        DESCRIPTION. The default is False.
    variance : TYPE, optional
        DESCRIPTION. The default is False.
    localmax : TYPE, optional
        DESCRIPTION. The default is False.
    segment : TYPE, optional
        DESCRIPTION. The default is False.
    proj : TYPE, optional
        DESCRIPTION. The default is 'AGU2025'.
    ctbls : TYPE, optional
        DESCRIPTION. The default is ['terrain_r','Blues'].
    cont : TYPE, optional
        DESCRIPTION. The default is False.
    bare_maps : TYPE, optional
        DESCRIPTION. The default is False.
    cb : TYPE, optional
        DESCRIPTION. The default is False.
    axNH3 : TYPE, optional
        DESCRIPTION. The default is False.
    axCH4 : TYPE, optional
        DESCRIPTION. The default is False.
    axRGB : TYPE, optional
        DESCRIPTION. The default is False.
    axIRTF : TYPE, optional
        DESCRIPTION. The default is False.
    axCH4889 : TYPE, optional
        DESCRIPTION. The default is False.
    axNH3vCloud : TYPE, optional
        DESCRIPTION. The default is False.
    axNH3vIRTF : TYPE, optional
        DESCRIPTION. The default is False.
    axIRTFvPCld : TYPE, optional
        DESCRIPTION. The default is False.
    ax889vPCld : TYPE, optional
        DESCRIPTION. The default is False.
    counter : TYPE, optional
        DESCRIPTION. The default is 0.
    countmax : TYPE, optional
        DESCRIPTION. The default is 0.

    Returns
    -------
    None.

    """

    import make_patch as MP
    import numpy as np
    import plot_contours_on_patch as PC
    import find_extrema as FX
    import find_blob as FB
    import plot_patch as PP
    import make_lat_lon_str as MLLS
    import sys
    sys.path.append('../Profiles/code')
    import plot_profile_scatter as pps


    if IRTFcollection and CH4889collection:
        RGBaxs=4
        IRTFaxs=2
        CH4889axs=3
    elif IRTFcollection:
        RGBaxs=3
        IRTFaxs=2
    elif CH4889collection:
        RGBaxs=3
        CH4889axs=2
    else:
        RGBaxs=2
    
    blendweightfNH3,fNH3stdv,fNH3frac,blendweightfNH3time,\
        blendweightPCloud,PCldstdv,PCldfrac,blendweightPCldtime,\
        blendRGBweight,RGBstdv,RGBfrac,RGBtime,\
        blendweightIRTF,IRTFstdv,IRTFfrac,IRTFtime,blendweightCH4889=\
        ReadContiguousMap(LonSys,collection=collection,
                              collectionIRTF=IRTFcollection,
                              collection889CH4=CH4889collection)

    fig1,fig2,axs1,axs2,fNH3_patch_mb,PCld_patch_mb,RGB_patch,tx_fNH3,tx_PCld,\
        fNH3low,fNH3high,PCldlow,PCldhigh,blendweightTime_patch,aspectratio,\
            RGB4Display,RGBaxs=\
                plot_optical_maps(LonLims,lats,
                              collection,blendweightfNH3,
                              blendweightPCloud,
                              blendweightfNH3time,
                              blendRGBweight,
                              ctbls,variance,
                              fNH3frac,PCldfrac,RGBaxs)

    if IRTFcollection:
        plot_5u_map(blendweightIRTF,LonLims,lats,axs1,axs2,variance,IRTFfrac,IRTFaxs)

    if CH4889collection:
        plot_JALPO_map(blendweightCH4889,LonLims,lats,axs1,axs2,variance,CH4889axs)

    ###########################################################################
    # Overplot contours
    ###########################################################################
    print("CONTOURS",cont)
    if cont:
        print("CONTOURS",cont)
        print("tx_fNH3",tx_fNH3)
        temp=PC.plot_contours_on_patch(axs1[0],fNH3_patch_mb,lats,[360-LonLims[1],360-LonLims[0]],
                                        tx_fNH3, frmt='%3.0f', clr='k')
        temp=PC.plot_contours_on_patch(axs1[1],PCld_patch_mb,lats,[360-LonLims[1],360-LonLims[0]],
                                        tx_PCld, frmt='%3.0f', clr='r')
        temp=PC.plot_contours_on_patch(axs1[RGBaxs],fNH3_patch_mb,lats,[360-LonLims[1],360-LonLims[0]],
                                        tx_fNH3[-2:], frmt='%3.0f', clr='k')
        temp=PC.plot_contours_on_patch(axs1[RGBaxs],PCld_patch_mb,lats,[360-LonLims[1],360-LonLims[0]],
                                        tx_PCld[:5], frmt='%3.0f', clr='r')
        if IRTFcollection:
            temp=PC.plot_contours_on_patch(axs1[IRTFaxs],fNH3_patch_mb,lats,[360-LonLims[1],360-LonLims[0]],
                                            tx_fNH3[-2:], frmt='%3.0f', clr='k')
            temp=PC.plot_contours_on_patch(axs1[IRTFaxs],PCld_patch_mb,lats,[360-LonLims[1],360-LonLims[0]],
                                            tx_PCld[:5], frmt='%3.0f', clr='r')
        if CH4889collection:
            temp=PC.plot_contours_on_patch(axs1[CH4889axs],fNH3_patch_mb,lats,[360-LonLims[1],360-LonLims[0]],
                                            tx_fNH3[-2:], frmt='%3.0f', clr='k')
            temp=PC.plot_contours_on_patch(axs1[CH4889axs],PCld_patch_mb,lats,[360-LonLims[1],360-LonLims[0]],
                                            tx_PCld[:5], frmt='%3.0f', clr='r')
    
    axs1[RGBaxs].set_title('RGB Context',fontsize=10)
    fig1.text(0.99,0.01,"Dr. Steven Hill, PSI",ha='right',fontsize=10)
    
    ###########################################################################
    # WRITE BARE MAPS FOR ANIMATIONS ETC
    ###########################################################################
    path="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Studies/"+proj+"/"

    if bare_maps:
        temp=make_bare_map(blendweightfNH3,ctbls[0],fNH3low,fNH3high,path,collection,LonSys,"fNH3")
        temp=make_bare_map(blendweightPCloud,ctbls[1],PCldlow,PCldhigh,path,collection,LonSys,"PCld")
        temp=make_bare_map(RGB4Display,ctbls[0],PCldlow,PCldhigh,path,collection,LonSys,"RGB")
        if IRTFcollection:
            temp=make_bare_map(np.log10(blendweightIRTF),'gist_heat',1.0,3.0,path,collection,LonSys,"IRTF")
        

    ###########################################################################

    if variance:
        show=axs2[RGBaxs].imshow(RGB4Display,
                   extent=[LonLims[1],LonLims[0],90-lats[1],
                           90-lats[0]],
                           aspect="equal")
        axs2[RGBaxs].set_title('RGB Context',fontsize=10)
        
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
    adjust_params = {
    1:     {'left': 0.21,  'bottom': 0.07, 'right': 0.79,  'top': 0.88,  'wspace': 0.0,  'hspace': 0.2},
    4/3.:  {'left': 0.20,  'bottom': 0.07, 'right': 0.825, 'top': 0.86,  'wspace': 0.0,  'hspace': 0.2},
    2:     {'left': 0.20,  'bottom': 0.07, 'right': 0.825, 'top': 0.86,  'wspace': 0.0,  'hspace': 0.2},
    3:     {'left': 0.095, 'bottom': 0.07, 'right': 0.90,  'top': 0.88,  'wspace': 0.25, 'hspace': 0.2},
    4:     {'left': 0.035, 'bottom': 0.07, 'right': 0.94,  'top': 0.88,  'wspace': 0.25, 'hspace': 0.2},
    6:     {'left': 0.04,  'bottom': 0.08, 'right': 0.94,  'top': 0.865, 'wspace': 0.0,  'hspace': 0.2},
    9:     {'left': 0.04,  'bottom': 0.10, 'right': 0.94,  'top': 0.865, 'wspace': 0.0,  'hspace': 0.395},
    12:    {'left': 0.035, 'bottom': 0.1,  'right': 0.935, 'top': 0.84,  'wspace': 0.0,  'hspace': 0.3}
            }
    if aspectratio in adjust_params:
        fig1.subplots_adjust(**adjust_params[aspectratio])     

    if variance:
        if aspectratio in adjust_params:
            fig2.subplots_adjust(**adjust_params[aspectratio])     

    if IRTFcollection and CH4889collection:
        
        axs1[3].set_title("CH4 889 nm",fontsize=8)
        axs1[4].set_title("RGB Context Image",fontsize=8)

        #box2 = axs1[2].get_position()
        #axs1[2].set_position([box2.x0-0.0, box2.y0-0.0, box2.width * 0.95, box2.height * 1.07])
        #box3 = axs1[3].get_position()
        #axs1[3].set_position([box3.x0-0.0, box3.y0-0.0, box3.width * 1.0, box3.height * 1.07])
        box4 = axs1[4].get_position()
        axs1[4].set_position([box4.x0+0.085, box4.y0-0.0, box4.width * 1.0, box4.height * 1.07])

    if ROI:
        for R in ROI:
            print("###################### R=",R)
            ##!!!! This color section was custom created for the NEZ ROIs 
            clr='C0'
            if "NH3" in R:
                clr='k'
            if "Plume" in R:
                clr='r'
            if "NEDF" in R:
                clr='y'
            if "NEB ref" in R:
                clr='brown'
            for iax in range(0,3):
                axs1[iax].plot(np.array([ROI[R][2]+ROI[R][3],ROI[R][2]-ROI[R][3],
                              ROI[R][2]-ROI[R][3],ROI[R][2]+ROI[R][3],
                              ROI[R][2]+ROI[R][3]]),
                              90.-np.array([ROI[R][0],ROI[R][0],ROI[R][1],
                              ROI[R][1],ROI[R][0]]),color=clr)
    if variance and ROI:
        for R in ROI:
            for iax in range(0,3):
                axs2[iax].plot(np.array([ROI[R][2]+ROI[R][3],ROI[R][2]-ROI[R][3],
                              ROI[R][2]-ROI[R][3],ROI[R][2]+ROI[R][3],
                              ROI[R][2]+ROI[R][3]]),
                              90.-np.array([ROI[R][0],ROI[R][0],ROI[R][1],
                              ROI[R][1],ROI[R][0]]))

    ###########################################################################
    # WRITE LOCAL MAX AND MINS TO FILE
    ###########################################################################
    latstr,lonstr=MLLS.make_lat_lon_str(lats,LonLims)

    if localmax:
        print("IN LOCAL MAX")
        ###########################################################################
        # Call the function
        # !!!!!!!NEED FIX FOR ADDITIONAL AXES, E.G., 5um and 889CH4

        results_extrema = FX.process_extrema({
            "NH3": fNH3_patch_mb,
            "PCloud": PCld_patch_mb,
            #"RGB": np.mean(RGB_patch,axis=2)
            "RGB":RGB_patch[:,:,0]
        }, blendweightTime_patch, lats, LonLims)

        output_filename=path+collection+" Mean Sys"+LonSys+" "+lonstr+" "+latstr+" extrema.csv"
        FX.export_extrema_to_csv(results_extrema, output_filename)
        FX.extrema_overplot_all(results_extrema,axes = {'axNH3': axs1[0], 'axCH4': axs1[1], 'axRGB': axs1[2]})

    
    ###########################################################################
    # Find 'blob' regions
    # !!!!! NEED FIX FOR ADDITIONAL AXES, E.G., 5um and 889CH4
    ###########################################################################
    if segment:
        #NH3thresh=135 #standard value
        NH3thresh=145 #alt value
        Cloudthresh=1750
        NEDFthresh=1950
        fNH3_mask, labeled_fNH3, props_fNH3= \
            FB.process_blob(fNH3_patch_mb, PCld_patch_mb, lats, LonLims, timearray=blendweightTime_patch,threshold_abs=NH3thresh, mode='max')
        
        Plum_mask, labeled_Plum, props_Plum= \
            FB.process_blob(PCld_patch_mb, fNH3_patch_mb, lats, LonLims, timearray=blendweightTime_patch, threshold_abs=Cloudthresh, mode='min')
        
        NEDF_mask, labeled_NEDF, props_NEDF= \
            FB.process_blob(PCld_patch_mb, fNH3_patch_mb, lats, LonLims, timearray=blendweightTime_patch, threshold_abs=NEDFthresh, mode='max')
        
        FB.export_regions_to_csv(props_fNH3, path+collection+" Mean Sys"+LonSys+" "+lonstr+" "+latstr+" blobs"+" fNH3.csv")
        FB.export_regions_to_csv(props_Plum, path+collection+" Mean Sys"+LonSys+" "+lonstr+" "+latstr+" blobs"+" Plum.csv")
        FB.export_regions_to_csv(props_NEDF, path+collection+" Mean Sys"+LonSys+" "+lonstr+" "+latstr+" blobs"+" NEDF.csv")
        print("$$$$$$$$$$$$$$$$$$$$$$$$$",lats,LonLims)

        FB.plot_regions_on_axis(axs1[2], labeled_fNH3, props_fNH3,lon_lims=LonLims,lats=lats,
                     plot_contours=False, plot_masks=True,plot_labels=False,contour_color='C0')
        FB.plot_regions_on_axis(axs1[2], labeled_Plum, props_Plum,lon_lims=LonLims,lats=lats,
                     plot_contours=False, plot_masks=True,plot_labels=False, contour_color='white')
        FB.plot_regions_on_axis(axs1[2], labeled_NEDF, props_NEDF,lon_lims=LonLims,lats=lats,
                     plot_contours=False, plot_masks=True,plot_labels=False, contour_color='black')

    fig1.savefig(path+collection+" Mean Sys"+LonSys+" "+lonstr+" "+latstr+" map.png",dpi=300)

    if variance:
        fig2.savefig(path+collection+" Stdv Sys"+LonSys+" "+lonstr+" "+latstr+" map.png",dpi=300)

    ###########################################################################
    ###########################################################################
    # Meridional Profiles
    ###########################################################################
    import pylab as pl

    ###########################################################################
    #! Set up belt and zone boundaries (need to make this a common service!)             
    belt={"SSTB":[-39.6,-36.2],
          "STB":[-32.4,-27.1],
          "SEB":[-19.7,-7.2],
          "NEB":[6.9,17.4],
          "NTB":[24.2,31.4],
          "NNTB":[35.4,39.6]}
    
    zone={"STZ":[-36.2,-32.4],
          "STrZ":[-27.1,-19.7],
          "EZ":[-7.2,6.9],
          "NTrZ":[17.4,24.2],
          "NTZ":[31.4,35.4]}

    rngprofileplots=np.arange(0,RGBaxs)
    Coords=np.linspace(-89.5,89.5,180)
        
    arrays=[blendweightfNH3,blendweightPCloud]
    if blendweightIRTF.any():
        arrays=[blendweightfNH3,blendweightPCloud,blendweightIRTF]
        if blendweightCH4889.any():
            arrays=[blendweightfNH3,blendweightPCloud,blendweightIRTF,blendweightCH4889]
    elif blendweightCH4889.any():
        arrays=[blendweightfNH3,blendweightPCloud,blendweightCH4889]

    # 0. Convert zeros to NaNs first
    arrays = [np.where(arr == 0, np.nan, arr) for arr in arrays]
    # 1. Build a combined mask: True where ANY array has a NaN
    mask_any_nan = np.any([np.isnan(arr) for arr in arrays], axis=0)
    # 2. Apply the mask using masked arrays
    masked_arrays = [np.ma.array(arr, mask=mask_any_nan) for arr in arrays]
    # 3. Compute row-wise mean (axis=1), ignoring masked values
    means = [ma.mean(axis=1) for ma in masked_arrays]
    stdvs = [ma.std(axis=1) for ma in masked_arrays]

    figmp,axsmp=pl.subplots(rngprofileplots[-1]+1,1,figsize=(6,6), dpi=150, facecolor="white",
                          sharex=True)
    figmp.suptitle("Multi-parameter Meridional Profiles")
    blendweightfNH3nan=np.where(blendweightfNH3==0.0,np.nan,blendweightfNH3)
    
    axsmp[0].plot(Coords,means[0][::-1],color='C0')
    axsmp[0].fill_between(Coords, means[0][::-1]-stdvs[0][::-1], 
                          means[0][::-1]+stdvs[0][::-1],alpha=0.1,color='C0')
    axsmp[0].set_ylim(60,160)

    blendweightPCldnan=np.where(blendweightPCloud==0.0,np.nan,blendweightPCloud)
    axsmp[1].plot(Coords,means[1][::-1],color='C0')
    axsmp[1].fill_between(Coords, means[1][::-1]-stdvs[1][::-1], 
                          means[1][::-1]+stdvs[1][::-1],alpha=0.1,color='C0')
    axsmp[1].set_ylim(1400,2000)
    
    if blendweightIRTF.any() and IRTFcollection:
        print("&&&&&&&&& IRTFcollection,blendweightIRTF.any()=",IRTFcollection,blendweightIRTF.any())
        axs5um=axsmp[1].twinx()
        axs5um.plot(Coords,np.log10(means[2][::-1]),color='C1')
        axs5um.fill_between(Coords, np.log10(np.maximum(means[2][::-1]-stdvs[2][::-1],0.001)), 
                              np.log10(means[2][::-1]+stdvs[2][::-1]),alpha=0.1,color='C1')
    
        axs5um.set_ylim(1,3)
        axs5um.tick_params(axis='y',labelsize=8)
        axs5um.set_ylabel("5-um Radiance (Log10)",fontsize=8)
    

        blendweightIRTFnan=np.where(blendweightIRTF==0.0,np.nan,blendweightIRTF)
        #axsmp[2].plot(np.nanmean(np.log10(blendweightIRTFnan),axis=1))
        axsmp[2].plot(Coords,np.log10(means[2][::-1]),color='C1')
        axsmp[2].fill_between(Coords, np.log10(np.maximum(means[2][::-1]-stdvs[2][::-1],0.001)), 
                              np.log10(means[2][::-1]+stdvs[2][::-1]),alpha=0.1,color='C1')
        axsmp[2].set_ylim(1,3)

    if blendweightCH4889.any() and CH4889collection: 
        blendweightCH4889nan=np.where(blendweightCH4889==0.0,np.nan,blendweightCH4889)
        #axsmp[3].plot(np.nanmean(blendweightCH4889nan,axis=1))
        
        #!!!!Need to fix indexing of means when IRTF is missing (s/b 2 instead of 3
        axsmp[3].plot(Coords,means[3][::-1],color='C0')
        axsmp[3].fill_between(Coords, means[3][::-1]-stdvs[3][::-1], 
                              means[3][::-1]+stdvs[3][::-1],alpha=0.1,color='C0')
    
        axsmp[3].set_ylim(50,200)
        #axsmp[3].set_ylim(-90,90)
    axsmp[RGBaxs-1].set_xlim(-30,30)
    
    plot_titles=["Ammonia Abundance "+collection,"Cloud Pressure "+collection,
                 "5-um Radiance "+str(IRTFcollection),"Methane CH4 Radiance "+str(CH4889collection)]
    ylabels=["fNH3 (ppm)","Cloud Pressure (mbar)","5-um Radiance (log10)","889nm CH4 Radiance"]
    count=0
    for a in axsmp:
        a.tick_params(axis='x',labelsize=8)
        a.tick_params(axis='y',labelsize=8)
        a.set_ylabel(ylabels[count],fontsize=8)
        a.yaxis.set_label_coords(-0.08,0.5)
        a.set_title(plot_titles[count],fontsize=10,y=0.95)

        for zb in belt:
            print(zb,belt[zb])
            a.fill_between([belt[zb][0],belt[zb][1]],np.array([0.,0.]),
                                    np.array([5000.,5000.]),color="0.5",alpha=0.2)

        count=count+1

    axsmp[RGBaxs-1].set_xlabel("PG Latitude")
    
    yb=51
    yz=51
    for zb in belt:
        axsmp[RGBaxs-1].annotate(zb,xy=[np.mean(belt[zb]),yb],ha="center")
    for zb in zone:
        axsmp[RGBaxs-1].annotate(zb,xy=[np.mean(zone[zb]),yz],ha="center")  
    
    
    figmp.subplots_adjust(bottom=0.07,top=0.92,left=0.10,right=0.92)
    
    figmp.savefig(path+collection+" Mean Sys"+LonSys+" "+lonstr+" "+latstr+" Meridional Profile.png",dpi=300)


    print("###############",means[0].shape)

    ###########################################################################
    ###########################################################################
    # Create Stack Plot subplots on the axes objects passed into the procedure
    ###########################################################################
    if axNH3vCloud!=False: # fNH3 versus PCloud
        fNH3low=100
        fNH3high=150
        PCldlow=1600
        PCldhigh=2000
        #micronlow=0.5
        #micronhigh=3.5

        axNH3vCloud.tick_params(axis='both', which='major', labelsize=10)
        axNH3vCloud.set_title("Cloud-Top Pressure vs Ammonia Abundance: 2022-2025",fontsize=12)
        #bands=["SEB","SEZ","NEZ","NEB"]
        #colors=["C1","C2","C3","C4"]
        #bands=["STrZ","SEB","SEZ","NEZ","NEB","NTrZ"]
        #colors=["C0","C1","C2","C3","C4","C5"]
        #bands=["SEB","NEB"]
        bands=["SEZ","NEZ"]
        colors=["C0","C1"]


        print("##################### counter=",counter)
        if counter==0:
            leg=True
        else:
            leg=False
        pps.plot_profile_scatter(means[1][::-1],means[0][::-1],Coords,axNH3vCloud,PCldlow,PCldhigh,
                         fNH3low,fNH3high,False,bands,colors,Level="L3",
                         leg=leg,axis_inv=True,date=collection[2:8],counter=counter,countmax=countmax)
        
        axNH3vCloud.legend()
        
        xmin,xmax=axNH3vCloud.get_xlim()
        ymin,ymax=axNH3vCloud.get_ylim()
        
        axNH3vCloud.scatter(xmin+0.03*(xmax-xmin),ymax-0.02*(ymax-ymin),marker='o',c='k',s=50)
        axNH3vCloud.scatter(xmin+0.03*(xmax-xmin),ymax-0.08*(ymax-ymin),marker='^',c='k',s=50)
        axNH3vCloud.scatter(xmin+0.03*(xmax-xmin),ymax-0.14*(ymax-ymin),marker='s',c='k',s=50)
        axNH3vCloud.scatter(xmin+0.03*(xmax-xmin),ymax-0.20*(ymax-ymin),marker='D',c='k',s=50)
        
        axNH3vCloud.annotate('2022',xy=(xmin+0.03*(xmax-xmin),ymax-0.01*(ymax-ymin)), xycoords='data',xytext=(xmin+0.05*(xmax-xmin),ymax-0.02*(ymax-ymin)),fontsize=9,verticalalignment='center_baseline')
        axNH3vCloud.annotate('2023',xy=(xmin+0.03*(xmax-xmin),ymax-0.01*(ymax-ymin)), xycoords='data',xytext=(xmin+0.05*(xmax-xmin),ymax-0.08*(ymax-ymin)),fontsize=9,verticalalignment='center_baseline')
        axNH3vCloud.annotate('2024',xy=(xmin+0.03*(xmax-xmin),ymax-0.01*(ymax-ymin)), xycoords='data',xytext=(xmin+0.05*(xmax-xmin),ymax-0.14*(ymax-ymin)),fontsize=9,verticalalignment='center_baseline')
        axNH3vCloud.annotate('2025',xy=(xmin+0.03*(xmax-xmin),ymax-0.01*(ymax-ymin)), xycoords='data',xytext=(xmin+0.05*(xmax-xmin),ymax-0.20*(ymax-ymin)),fontsize=9,verticalalignment='center_baseline')
        
        #######################################################################

    if axNH3vIRTF!=False and blendweightIRTF.any(): # fNH3 versus IRTF
        fNH3low=100
        fNH3high=150
        IRTFlow=1
        IRTFhigh=3
        #micronlow=0.5
        #micronhigh=3.5

        axNH3vIRTF.tick_params(axis='both', which='major', labelsize=10)
        axNH3vIRTF.set_title("IRTF Radiance vs Ammonia Abundance: 2022-2025",fontsize=12)
        #bands=["SEB","SEZ","NEZ","NEB"]
        #colors=["C1","C2","C3","C4"]
        #bands=["STrZ","SEB","SEZ","NEZ","NEB","NTrZ"]
        #colors=["C0","C1","C2","C3","C4","C5"]
        #bands=["SEB","NEB"]
        bands=["SEZ","NEZ"]
        colors=["C0","C1"]


        print("##################### counter=",counter)
        if counter==0:
            leg=True
        else:
            leg=False
        pps.plot_profile_scatter(np.log10(means[2][::-1]),means[0][::-1],Coords,axNH3vIRTF,IRTFlow,IRTFhigh,
                         fNH3low,fNH3high,True,bands,colors,Level="L3",
                         leg=leg,axis_inv=True,date=collection[2:8],counter=counter,countmax=countmax)
        
        axNH3vIRTF.legend()
        
        xmin,xmax=axNH3vIRTF.get_xlim()
        ymin,ymax=axNH3vIRTF.get_ylim()
        
        axNH3vIRTF.scatter(xmin+0.03*(xmax-xmin),ymax-0.02*(ymax-ymin),marker='o',c='k',s=50)
        axNH3vIRTF.scatter(xmin+0.03*(xmax-xmin),ymax-0.08*(ymax-ymin),marker='^',c='k',s=50)
        axNH3vIRTF.scatter(xmin+0.03*(xmax-xmin),ymax-0.14*(ymax-ymin),marker='s',c='k',s=50)
        axNH3vIRTF.scatter(xmin+0.03*(xmax-xmin),ymax-0.20*(ymax-ymin),marker='D',c='k',s=50)
        
        axNH3vIRTF.annotate('2022',xy=(xmin+0.03*(xmax-xmin),ymax-0.01*(ymax-ymin)), xycoords='data',xytext=(xmin+0.05*(xmax-xmin),ymax-0.02*(ymax-ymin)),fontsize=9,verticalalignment='center_baseline')
        axNH3vIRTF.annotate('2023',xy=(xmin+0.03*(xmax-xmin),ymax-0.01*(ymax-ymin)), xycoords='data',xytext=(xmin+0.05*(xmax-xmin),ymax-0.08*(ymax-ymin)),fontsize=9,verticalalignment='center_baseline')
        axNH3vIRTF.annotate('2024',xy=(xmin+0.03*(xmax-xmin),ymax-0.01*(ymax-ymin)), xycoords='data',xytext=(xmin+0.05*(xmax-xmin),ymax-0.14*(ymax-ymin)),fontsize=9,verticalalignment='center_baseline')
        axNH3vIRTF.annotate('2025',xy=(xmin+0.03*(xmax-xmin),ymax-0.01*(ymax-ymin)), xycoords='data',xytext=(xmin+0.05*(xmax-xmin),ymax-0.20*(ymax-ymin)),fontsize=9,verticalalignment='center_baseline')

        #######################################################################
        
        axIRTFvPCld.tick_params(axis='both', which='major', labelsize=10)
        axIRTFvPCld.set_title("Cloud Pressure vs IRTF Radiance: 2022-2025",fontsize=12)
        #bands=["SEB","SEZ","NEZ","NEB"]
        #colors=["C1","C2","C3","C4"]
        #bands=["STrZ","SEB","SEZ","NEZ","NEB","NTrZ"]
        #colors=["C0","C1","C2","C3","C4","C5"]
        #bands=["SEB","NEB"]
        bands=["SEZ","NEZ"]
        colors=["C0","C1"]
        
        print("##################### counter=",counter)
        if counter==0:
            leg=True
        else:
            leg=False
        pps.plot_profile_scatter(means[1][::-1],np.log10(means[2][::-1]),Coords,axIRTFvPCld,
                         PCldlow,PCldhigh,IRTFlow,IRTFhigh,True,bands,colors,Level="L3",
                         leg=leg,axis_inv=True,date=collection[2:8],counter=counter,countmax=countmax)
        
        axNH3vIRTF.legend()
        
        xmin,xmax=axIRTFvPCld.get_xlim()
        ymin,ymax=axIRTFvPCld.get_ylim()
        
        axIRTFvPCld.scatter(xmin+0.03*(xmax-xmin),ymax-0.02*(ymax-ymin),marker='o',c='k',s=50)
        axIRTFvPCld.scatter(xmin+0.03*(xmax-xmin),ymax-0.08*(ymax-ymin),marker='^',c='k',s=50)
        axIRTFvPCld.scatter(xmin+0.03*(xmax-xmin),ymax-0.14*(ymax-ymin),marker='s',c='k',s=50)
        axIRTFvPCld.scatter(xmin+0.03*(xmax-xmin),ymax-0.20*(ymax-ymin),marker='D',c='k',s=50)
        
        axIRTFvPCld.annotate('2022',xy=(xmin+0.03*(xmax-xmin),ymax-0.01*(ymax-ymin)), xycoords='data',xytext=(xmin+0.05*(xmax-xmin),ymax-0.02*(ymax-ymin)),fontsize=9,verticalalignment='center_baseline')
        axIRTFvPCld.annotate('2023',xy=(xmin+0.03*(xmax-xmin),ymax-0.01*(ymax-ymin)), xycoords='data',xytext=(xmin+0.05*(xmax-xmin),ymax-0.08*(ymax-ymin)),fontsize=9,verticalalignment='center_baseline')
        axIRTFvPCld.annotate('2024',xy=(xmin+0.03*(xmax-xmin),ymax-0.01*(ymax-ymin)), xycoords='data',xytext=(xmin+0.05*(xmax-xmin),ymax-0.14*(ymax-ymin)),fontsize=9,verticalalignment='center_baseline')
        axIRTFvPCld.annotate('2025',xy=(xmin+0.03*(xmax-xmin),ymax-0.01*(ymax-ymin)), xycoords='data',xytext=(xmin+0.05*(xmax-xmin),ymax-0.20*(ymax-ymin)),fontsize=9,verticalalignment='center_baseline')

    ###########################################################################

    if ax889vPCld!=False and blendweightCH4889.any():
        CH4889low=100
        CH4889high=200
        IRTFlow=1
        IRTFhigh=3
        #micronlow=0.5
        #micronhigh=3.5

        ax889vPCld.tick_params(axis='both', which='major', labelsize=10)
        ax889vPCld.set_title("Cloud Pressure versus 889-nm Radiance: 2022-2025",fontsize=12)
        #bands=["SEB","SEZ","NEZ","NEB"]
        #colors=["C1","C2","C3","C4"]
        #bands=["STrZ","SEB","SEZ","NEZ","NEB","NTrZ"]
        #colors=["C0","C1","C2","C3","C4","C5"]
        #bands=["SEB","NEB"]
        bands=["SEZ","NEZ"]
        colors=["C0","C1"]


        print("##################### counter=",counter)
        if counter==0:
            leg=True
        else:
            leg=False
        pps.plot_profile_scatter(means[1][::-1],means[3][::-1],Coords,ax889vPCld,PCldlow,PCldhigh,
                         CH4889low,CH4889high,False,bands,colors,Level="L3",
                         leg=leg,axis_inv=True,date=collection[2:8],counter=counter,countmax=countmax)
        
        ax889vPCld.legend()

        xmin,xmax=ax889vPCld.get_xlim()
        ymin,ymax=ax889vPCld.get_ylim()
        
        ax889vPCld.scatter(xmin+0.03*(xmax-xmin),ymax-0.02*(ymax-ymin),marker='o',c='k',s=50)
        ax889vPCld.scatter(xmin+0.03*(xmax-xmin),ymax-0.08*(ymax-ymin),marker='^',c='k',s=50)
        ax889vPCld.scatter(xmin+0.03*(xmax-xmin),ymax-0.14*(ymax-ymin),marker='s',c='k',s=50)
        ax889vPCld.scatter(xmin+0.03*(xmax-xmin),ymax-0.20*(ymax-ymin),marker='D',c='k',s=50)
        
        ax889vPCld.annotate('2022',xy=(xmin+0.03*(xmax-xmin),ymax-0.01*(ymax-ymin)), xycoords='data',xytext=(xmin+0.05*(xmax-xmin),ymax-0.02*(ymax-ymin)),fontsize=9,verticalalignment='center_baseline')
        ax889vPCld.annotate('2023',xy=(xmin+0.03*(xmax-xmin),ymax-0.01*(ymax-ymin)), xycoords='data',xytext=(xmin+0.05*(xmax-xmin),ymax-0.08*(ymax-ymin)),fontsize=9,verticalalignment='center_baseline')
        ax889vPCld.annotate('2024',xy=(xmin+0.03*(xmax-xmin),ymax-0.01*(ymax-ymin)), xycoords='data',xytext=(xmin+0.05*(xmax-xmin),ymax-0.14*(ymax-ymin)),fontsize=9,verticalalignment='center_baseline')
        ax889vPCld.annotate('2025',xy=(xmin+0.03*(xmax-xmin),ymax-0.01*(ymax-ymin)), xycoords='data',xytext=(xmin+0.05*(xmax-xmin),ymax-0.20*(ymax-ymin)),fontsize=9,verticalalignment='center_baseline')

    ###########################################################################
    ###########################################################################
    # STACK MAP PLOTS
    ###########################################################################

    if axNH3!=False:
        fNH3low=60
        fNH3high=160
        PCldlow=1400
        PCldhigh=2200

        #lats=[100,120]
        #lats=[80,100]
        fNH3_patch_mb=MP.make_patch(blendweightfNH3,lats,[360-LonLims[1],360-LonLims[0]],
                                    180,180)
        print("######### cb=",cb)
        #plot_patch(patch,LatLims,LonLims,CM2,LonRng,colorscale,axis,
        #               cbarplot=True,cbar_title="Test",cbar_reverse=False,vn=0.10,vx=0.20,n=6)
        fNH3_patch_mb,vn,vx,tx_fNH3=PP.plot_patch(fNH3_patch_mb,lats,[360-LonLims[1],360-LonLims[0]],
                                         180,180,ctbls[0],
                                         axNH3,cbarplot=cb,n=11,
                                         vn=fNH3low,
                                         vx=fNH3high)
        axNH3.set_ylabel(collection.replace('-','\n'),rotation='horizontal',fontsize=6)
        axNH3.yaxis.set_label_coords(-0.10,0.15)
        axNH3.tick_params('y', labelleft=False)
        axNH3.tick_params('x', labelsize=8)


        PCld_patch_mb=MP.make_patch(blendweightPCloud,lats,[360-LonLims[1],360-LonLims[0]],
                                    180,180)
        PCld_patch_mb,vn,vx,tx_fNH3=PP.plot_patch(PCld_patch_mb,lats,[360-LonLims[1],360-LonLims[0]],
                                         180,180,ctbls[1],
                                         axCH4,cbarplot=cb,
                                         n=5,vn=PCldlow,vx=PCldhigh)
        axCH4.set_ylabel(collection,rotation='horizontal',fontsize=6)
        axCH4.set_ylabel(collection.replace('-','\n'),rotation='horizontal',fontsize=6)
        axCH4.yaxis.set_label_coords(-0.10,0.15)
        axCH4.tick_params('y', labelleft=False)
        axCH4.tick_params('x', labelsize=8)

           
        RGB_patch=MP.make_patch(blendRGBweight,lats,[360-LonLims[1],360-LonLims[0]],180,180)
        RGB4Display=np.power(np.array(RGB_patch).astype(float),1.0)
        #RGB4Display=RGB4Display/RGB4Display.max()
        show=axRGB.imshow(RGB4Display,
                   extent=[LonLims[1],LonLims[0],90-lats[1],
                           90-lats[0]],
                           aspect="equal")
        axRGB.set_ylabel(collection.replace('-','\n'),rotation='horizontal',fontsize=6)
        axRGB.tick_params('y', labelleft=False)
        axRGB.yaxis.set_label_coords(-0.10,0.15)
        axRGB.tick_params('x', labelsize=8)


        IRTF_patch_mb=MP.make_patch(blendweightIRTF,lats,[360-LonLims[1],360-LonLims[0]],
                                    180,180)
        IRTF_patch_mb,vn,vx,tx_fNH3=PP.plot_patch(np.log10(IRTF_patch_mb),lats,[360-LonLims[1],360-LonLims[0]],
                                         180,180,"gist_heat",
                                         axIRTF,cbarplot=cb,
                                         n=5,vn=1.5,vx=3.5)
        if IRTFcollection:
            #axIRTF.set_ylabel(IRTFcollection,rotation='horizontal',fontsize=6)
            axIRTF.set_ylabel(IRTFcollection.replace('-','\n'),rotation='horizontal',fontsize=6)
            axIRTF.yaxis.set_label_coords(-0.10,0.15)
            axIRTF.tick_params('y', labelleft=False)
            axIRTF.tick_params('x', labelsize=8)
        else:
            axIRTF.set_yticklabels([])
            


        CH4889_patch_mb=MP.make_patch(blendweightCH4889,lats,[360-LonLims[1],360-LonLims[0]],
                                    180,180)
        CH4889_patch_mb,vn,vx,tx_fNH3=PP.plot_patch(CH4889_patch_mb,lats,[360-LonLims[1],360-LonLims[0]],
                                         180,180,"gray",
                                         axCH4889,cbarplot=cb,
                                         n=5,vn=80,vx=230)
        
        axCH4889.set_ylabel(None)
        if CH4889collection:
            axCH4889.set_ylabel(CH4889collection,rotation='horizontal',fontsize=6)
            axCH4889.set_ylabel(CH4889collection.replace('-','\n'),rotation='horizontal',fontsize=6)
            axCH4889.yaxis.set_label_coords(-0.10,0.15)
            axCH4889.tick_params('y', labelleft=False)
            axCH4889.tick_params('x', labelsize=8)
        else:
            axCH4889.set_yticklabels([])

        if localmax:
            FX.extrema_overplot_all(results_extrema,axes = {'axNH3': axNH3, 'axCH4': axCH4, 'axRGB': axRGB})

        if segment:
            FB.plot_regions_on_axis(axRGB, labeled_fNH3, props_fNH3,lon_lims=LonLims,lats=lats,
                         plot_contours=False, plot_masks=True,plot_labels=False,contour_color='C0')
            FB.plot_regions_on_axis(axRGB, labeled_Plum, props_Plum,lon_lims=LonLims,lats=lats,
                         plot_contours=False, plot_masks=True,plot_labels=False, contour_color='white')
            FB.plot_regions_on_axis(axRGB, labeled_NEDF, props_NEDF,lon_lims=LonLims,lats=lats,
                         plot_contours=False, plot_masks=True,plot_labels=False, contour_color='black')

    ###########################################################################
    # YET ANOTHER SECTION: LONGITUDINAL CUTS!
    #!!!! Need to make MAPS! then use profiles modules to extract this!
    ###########################################################################
    
    figlc,axslc=pl.subplots(3,1,figsize=(figsz[0],figsz[1]), dpi=150, facecolor="white",
                          sharex=True)
    figlc.suptitle(collection)

    lon_array=np.arange(LonLims[1],LonLims[0],-1)
    fNH3_array0=np.mean(fNH3_patch_mb[4:7,:],axis=0)
    Cld_array0=np.mean(PCld_patch_mb[4:7,:],axis=0)
    fNH3_array1=np.mean(fNH3_patch_mb[5:9,:],axis=0)
    Cld_array1=np.mean(PCld_patch_mb[5:9,:],axis=0)
    fNH3_array2=np.mean(fNH3_patch_mb[9:12,:],axis=0)
    Cld_array2=np.mean(PCld_patch_mb[9:12,:],axis=0)
    #PCld_array=np.mean(PCld_patch_mb[7:9,:],axis=0)
    #NEDF_array=np.mean(PCld_patch_mb[5:7,:],axis=0)

    axslc[0].plot(lon_array,fNH3_array0,color='C2',label='fNH3')
    axslc[0].plot(lon_array,np.roll(fNH3_array0,10),color='C2',linestyle='dashed',linewidth=0.5)
    axslc0=axslc[0].twinx()
    axslc0.plot(lon_array,Cld_array0,color='C0',label='PCld')
    axslc[0].set_ylim(50.,180.)
    axslc0.set_ylim(1500.,2200.)
    axslc0.invert_yaxis()

    axslc[1].plot(lon_array,fNH3_array1,color='C2')
    axslc[1].plot(lon_array,np.roll(fNH3_array1,10),color='C2',linestyle='dashed',linewidth=0.5)
    axslc1=axslc[1].twinx()
    axslc1.plot(lon_array,Cld_array1,color='C0')
    axslc[1].set_ylim(50.,180.)
    axslc1.set_ylim(1500.,2200.)
    axslc1.invert_yaxis()

    axslc[2].plot(lon_array,fNH3_array2,color='C2')
    axslc[2].plot(lon_array,np.roll(fNH3_array2,10),color='C2',linestyle='dashed',linewidth=0.5)
    axslc2=axslc[2].twinx()
    axslc2.plot(lon_array,Cld_array2,color='C0')
    axslc[2].set_ylim(50.,180.)
    axslc2.set_ylim(1500.,2200.)
    axslc2.invert_yaxis()

    #axslc[1].plot(lon_array,PCld_array)
    #axslc[2].plot(lon_array,fNH3_array2)
    
    for i in range(0,3):
        axslc[i].grid(linewidth=0.2)
        #axslc[i].set_xticks(np.linspace(450.,0.,31), minor=False)
        #axslc[i].set_yticks(np.linspace(-90,90,13), minor=False)
        #yticklabels=np.array(np.linspace(-90,90,13))
        #axslc[i].set_yticklabels(yticklabels.astype(int))
        axslc[i].tick_params(axis='both', which='major', labelsize=8)
        axslc[i].set_ylabel("fNH3 (ppm)",color='C2')
        #axslc0.set_ylabel("PCloud (mb)",color='C0')
        #axslc[i].set_adjustable('box')
        axslc[i].set_xlim(360,0.)

    axslc0.set_ylabel("PCloud (mb)",color='C0')
    axslc0.tick_params(axis='both', which='major', labelsize=8)
    axslc1.set_ylabel("PCloud (mb)",color='C0')
    axslc1.tick_params(axis='both', which='major', labelsize=8)
    axslc2.set_ylabel("PCloud (mb)",color='C0')
    axslc2.tick_params(axis='both', which='major', labelsize=8)
    
    axslc[2].set_xlabel("System 3 Longitude (deg)")
    
    axslc[0].set_title("8-11 deg PG Latitude - NEDFs",fontsize=10)
    axslc[1].set_title("6-10 deg PG Latitude - Plumes",fontsize=10)
    axslc[2].set_title("3-6 deg PG Latitude - NH3 Features",fontsize=10)
    
    #xticklabels=np.array(np.linspace(195,15,13))
    #print(xticklabels)
    #axslc[2].set_xticklabels(xticklabels.astype(int))
    #axslc[2].invert_xaxis()
    Rossby=np.sin((lon_array+30)*9*np.pi/180.)*20+100
    #print(lon_array)
    #axslc[0].plot(lon_array,Rossby,color='k')

    figlc.subplots_adjust(bottom=0.07,top=0.92,left=0.10,right=0.90)


    
    figlc.savefig(path+collection+" Wave.png",dpi=150)
    #Autocorrelation (below) doesn't work well due to differing shapes
    """
    from Winds import _xcorr1d_circular
    shift_px, peak, snr = _xcorr1d_circular(1.0/(Cld_array2/np.mean(Cld_array2)),
                                            fNH3_array2/np.mean(fNH3_array2))
    print("%%%%%%%%%%%%%%%%%% shift_px, peak, snr=",shift_px, peak, snr)
    """
    
    
    if segment:
        return(lats,blendweightPCloud,blendweightfNH3,blendRGBweight,
               labeled_fNH3, props_fNH3, 
                labeled_Plum, props_Plum, 
                labeled_NEDF, props_NEDF)
    else:
        return(lats,blendweightPCloud,blendweightfNH3,blendRGBweight)
    
    #return(fig1,axs1)