def grs_zoom_maps(year=2022):
    """
    Created on Fri Jan  5 07:32:06 2024
    
    @author: smhil
    """
    import numpy as np
    import pylab as pl
    import sys
    sys.path.append('./Maps')
    import Map_Jup_Atm_P3 as MJP
    import read_fits_map_L2_L3 as RFM
    import plot_patch as PP
    import RetrievalLibrary as RL
    import make_patch_RGB as MPRGB

    #GRS2022=['20220730UTa','20220810UTa','20220818UTa','20220828UTa','20220904UTa',
    #          '20220919UTa','20220919UTb','20221013UTa','20221020UTa','20230113UTa']

    GRS2022=['20220818UTa','20220828UTa','20220904UTa',
             '20221013UTa','20221020UTa']
    
    #GRS2023=["20230831UTa","20230831UTb","20230905UTa","20230924UTa","20231005UTa",
    #         "20231015UTa","20231017UTa","20231017UTb","20231022UTa","20231103UTa",
    #         "20231110UTb","20231110UTc","20231113UTa","20231113UTb","20231115UTa",
    #         "20231115UTb","20231207UTa"]
    GRS2023=["20230831UTa","20230831UTb","20230905UTa","20230924UTa","20231005UTa",
             "20231017UTa","20231022UTa","20231103UTa",
             "20231113UTa","20231207UTa"]
    #GRS2023=["20230831UTb","20231017UTa","20231022UTa",
    #         "20231110UTb","20231113UTa","20231207UTa"]
    
    if year==2022:
        GRS=GRS2022
        CM=25
    elif year==2023:
        GRS=GRS2023
        CM=45
    LonRng=20
    LonSys='2'
    LonLims=[360-int(CM+LonRng),360-int(CM-LonRng)]
    LatLims=[90,130]

    #coefs=[0.,0.]
    coefs=[0.65,0.25]

    First=True
    for obskey in GRS:
        
        MJP.Map_Jup_Atm_P3(obskey=obskey,imagetype='Map', Smoothing=False,
                            LatLims=[90,130],LonRng=20,CMpref=CM,LonSys=LonSys,
                            showbands=False,coef=coefs,subproj='GRS')

        PCloudhdr,PClouddata,fNH3hdr,fNH3data,sza,eza,RGB,RGB_CM,RGBtime= \
                        RFM.read_fits_map_L2_L3(obskey=obskey,LonSys=LonSys,
                                                imagetype="Map",Level="L3")   
                        
        print("@@@@@@@@@ fNH3data.shape=",fNH3data.shape)
        if First:
            MapAvgArr=fNH3data
            MapAvgArr=np.reshape(MapAvgArr,(180,360,1))
            
            RGB0AvgArr=RGB[:,:,0]
            RGB1AvgArr=RGB[:,:,1]
            RGB2AvgArr=RGB[:,:,2]
            First=False
        else:
            MapAvgArr=np.dstack((MapAvgArr,fNH3data))
            RGB0AvgArr=np.dstack((RGB0AvgArr,RGB[:,:,0]))
            RGB1AvgArr=np.dstack((RGB1AvgArr,RGB[:,:,1]))
            RGB2AvgArr=np.dstack((RGB2AvgArr,RGB[:,:,2]))
    MapAvg=np.mean(MapAvgArr[:,:,:],axis=2)
    MapStd=np.std(MapAvgArr[:,:,:],axis=2)
    
    RGB0Avg=np.mean(RGB0AvgArr[:,:,:],axis=2)
    RGB1Avg=np.mean(RGB1AvgArr[:,:,:],axis=2)
    RGB2Avg=np.mean(RGB2AvgArr[:,:,:],axis=2)
    
    RGBAvg=np.dstack((RGB0Avg,RGB1Avg,RGB2Avg))
    print(MapAvg.shape)
    print(RGB0AvgArr.shape)
    print(RGBAvg.shape)

    
    ###########################################################################
    ## Just RGB and Abundance
    ###########################################################################

    fig1,axs1=pl.subplots(1,2,figsize=(8.0,4.0), dpi=150, facecolor="white",
                        sharey=True,sharex=True)
    #fig1.suptitle(fNH3hdr["DATE-OBS"].replace("_"," ")+", CM"+LonSys+"="
    #              +str(int(fNH3PlotCM))+", Calibration = "+CalModel+
    #              ", Data "+smthtitle,x=0.5,ha='center',color='k')

    for ix in range(0,1):
        axs1[ix].grid(linewidth=0.2)
        axs1[ix].ylim=[-45.,45.]
        axs1[ix].xlim=[360-LonLims[0],360-LonLims[1]]
        axs1[ix].set_xticks(np.linspace(450,0,31), minor=False)
        xticklabels=np.array(np.mod(np.linspace(450,0,31),360))
        axs1[ix].set_xticklabels(xticklabels.astype(int))
        axs1[ix].set_yticks(np.linspace(-45,45,7), minor=False)
        axs1[ix].tick_params(axis='both', which='major', labelsize=9)

        axs1[ix].set_adjustable('box') 

    #TestfNH3=fNH3data*amfdata**coef[0]
    
    fNH3_patch_mb,vn,vx,tx_fNH3=PP.plot_patch(MapAvg,LatLims,LonLims,
                                     CM,LonRng,"jet",
                                     axs1[0],'%3.2f',cont=False,n=6,vn=80,vx=130)

    temp=RL.make_contours_CH4_patch(axs1[0],fNH3_patch_mb,LatLims,LonLims,
                           lvls=tx_fNH3,frmt='%3.0f',clr='k')
    
    axs1[0].set_title(r'$\bar{f_c}(NH3) (ppm)$',fontsize=10)

    gamma=1.3

    #Logic in RGB_patch depends on LonLims and CM being consistent
    RGB_patch=MPRGB.make_patch_RGB(RGBAvg,LatLims,LonLims,CM,LonRng)
    
    RGB4Display=np.power(np.array(RGB_patch).astype(float),gamma)
    RGB4Display=RGB4Display/RGB4Display.max()
    show=axs1[1].imshow(RGB4Display,
               extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],
                       90-LatLims[0]],
                       aspect="equal")
    temp=RL.make_contours_CH4_patch(axs1[1],fNH3_patch_mb,LatLims,LonLims,
                           tx_fNH3,frmt='%3.0f',clr='k')
    box = axs1[1].get_position()
    
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

    showbands=False
    if showbands:
        for zb in belt:
            #print(zb,belt[zb])
            axs1[0].fill_between([360-LonLims[0],360-LonLims[1]],[belt[zb][0],belt[zb][0]],[belt[zb][1],belt[zb][1]],
                                    color="0.5",alpha=0.2)
            axs1[1].fill_between([360-LonLims[0],360-LonLims[1]],[belt[zb][0],belt[zb][0]],[belt[zb][1],belt[zb][1]],
                                    color="0.8",alpha=0.1)
        #axs1[1].annotate(zb,xy=[np.mean(belt[zb]),51],ha="center")
    #for zb in zone:
        #axs1[1].annotate(zb,xy=[np.mean(zone[zb]),51],ha="center")

    axs1[1].tick_params(axis='both', which='major', labelsize=9)
    axs1[1].set_title("RGB Context Image",fontsize=10)

    axs1[0].set_ylabel("Planetographic Latitude (deg)",fontsize=10)
    axs1[0].set_xlabel("Sys. "+LonSys+" Longitude (deg)",fontsize=10)
    axs1[1].set_xlabel("Sys. "+LonSys+" Longitude (deg)",fontsize=10)
    axs1[1].grid(linewidth=0.2)

    fig1.subplots_adjust(left=0.10, bottom=0.03, right=0.98, top=0.95,
                wspace=0.25, hspace=0.05)     
    axs1[1].set_position([box.x0+0.03, box.y0-0.01, box.width * 1.015, box.height * 1.015])
    
    coef=0.0
    if coef==0.0:
        correction='_C0'
    else:
        correction='_C1'
    print("360-LonLims[0],360-LonLims[1]=",
          360-LonLims[0],360-LonLims[1])
    fnskeleton=correction+'_Sys'+LonSys+'_N'+\
                str(90-LatLims[0])+'-S'+str(LatLims[1]-90)+\
                '_Lon'+str(np.mod(360-LonLims[1],360)).zfill(3)+'-'+\
                    str(np.mod(360-LonLims[0],360)).zfill(3)+'.png'
                    
    subproj='GRS'                
    pathmapplots='C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Analysis Data/L3 Plots/'+subproj+'/'
                    
    #fnNH3=fNH3hdr["FILENAME"][:-5]+fnskeleton
    fnNH3=str(int(year))+"_fNH3_Avg5"+fnskeleton
    fig1.savefig(pathmapplots+fnNH3,dpi=300)
    
    return(fig1,MapAvg)