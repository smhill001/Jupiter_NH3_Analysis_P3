# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 09:07:27 2022

@author: smhil
"""

def GRS_NH3_Comparison_2022(LatLims=[70,130],CM2=20,LonRng=30):
    import sys
    drive='c:'
    sys.path.append(drive+'/Astronomy/Python Play')
    sys.path.append(drive+'/Astronomy/Python Play/Util_P3')
    sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy')

    import os
    from matplotlib.pyplot import imread
    import pylab as pl
    import numpy as np
    from imageio import imwrite
    from numpy import inf
    from re import search
    from astropy.io import fits
    from astropy.convolution import Gaussian2DKernel
    from astropy.convolution import convolve
    from astropy.io import fits
    import plot_TEXES_Groups_P3 as PTG


    ###########################################################################
    #  DATA FILES AND METADATA DICTIONARY
    #    !!!!SHOULD MAKE THIS A DATA OBJECT!!!!
    #    !!!!DOUBLE CHECK THAT FITS FILE TIME TAGS ARE ACCURATE BETWEEN CH4 AND NH3
    ###########################################################################
#    sourcedata=obsdate+"_"+imagetype
    sourcefiles={'20220810UT':{'fNH3file':'2022-08-10-1013_0-Jupiter_fNH3Abs647.fits',
                               'RGBfile':'2022-08-10-1030_0-Jupiter_WV-R(AllRED)GB-RGB-WhtBal-Wavelets-Str_CM2_L360_MAP-BARE.png'},
                 '20220818UT':{'fNH3file':'2022-08-18-0733_4-Jupiter_fNH3Abs647.fits',
                               'RGBfile':'2022-08-18-0745_4-Jupiter_AllRED-WV-RGB-WhtBal-Wavelets_CM2_L360_MAP-BARE.png'},
                 '20220828UT':{'fNH3file':'2022-08-28-0608_2-Jupiter_fNH3Abs647.fits',
                               'RGBfile':'2022-08-28-0601_3-Jupiter-RGB-JamesWillinghan-j220828a1_CM2_L360_MAP-BARE.png'},
                 '20220904UT':{'fNH3file':'2022-09-04-0638_2-Jupiter_fNH3Abs647.fits',
                               'RGBfile':'2022-09-04-1644_7-Jupiter-Yamane-j220904j4_CM2_L360_MAP-BARE.png'},
                 '20220919UT':{'fNH3file':'2022-09-19-0453_4-Jupiter_fNH3Abs647.fits',
                               'RGBfile':'2022-09-19-0518_7-Jupiter_WV-R685G550B450-RGB-WhtBal-Wavelets_CM2_L360_MAP-BARE.png'},
                 '20221013UT':{'fNH3file':'2022-10-13-0345_5-Jupiter-fNH3Abs647.fits',
                               'RGBfile':'2022-10-13-0402_0-Jupiter-WV-R685G550B450-RGB-WhtBal-Wavelets_CM2_L360_MAP-BARE.png'},
                 '20221020UT':{'fNH3file':'2022-10-20-0440_4-Jupiter-fNH3Abs647.fits',
                               'RGBfile':'2022-10-20-0422_6-Jupiter-WV-R685G550B450-RGB-WhtBal-Wavelets_CM2_L360_MAP-BARE.png'}}
        
    fig,axs=pl.subplots(2,7,figsize=(12.0,4.0), dpi=150, facecolor="white",
                        sharey=True,sharex=True)
    fig.suptitle('2022 GRS NH3',x=0.5,ha='center',color='k')
    LonLims=[360-int(CM2+LonRng),360-int(CM2-LonRng)]
    target='Jupiter'
    Dates=['20220810UT','20220818UT','20220828UT','20220904UT','20220919UT','20221013UT','20221020UT']
    kernel = Gaussian2DKernel(1)
    #clrmap='gist_heat'
    clrmap='jet'
    
    MeridEWArray=np.zeros((LatLims[1]-LatLims[0],len(Dates)))

    first=True
    firstRGB=True
    for iPlot in range(0,14):
            # Set up
        i=int(iPlot/7)                           #Plot row
        j=np.mod(iPlot,7)                   #Plot column
        axs[i,j].grid(linewidth=0.2)
        axs[i,j].ylim=[-45.,45.]
        axs[i,j].xlim=[360-LonLims[0],360-LonLims[1]]
        print(360-LonLims[0],360-LonLims[1])
        print(LatLims)
        axs[i,j].set_xticks(np.linspace(450,0,31), minor=False)
        xticklabels=np.array(np.mod(np.linspace(450,0,31),360))
        axs[i,j].set_xticklabels(xticklabels.astype(int))
        axs[i,j].set_yticks(np.linspace(-45,45,7), minor=False)
        axs[i,j].tick_params(axis='both', which='major', labelsize=7)

        axs[i,j].set_adjustable('box') 
        path='c:/Astronomy/Projects/Planets/'+target+'/Imaging Data/'+Dates[j][0:10]+'/'
        print(Dates[j])

        fNH3hdulist=fits.open(path+sourcefiles[Dates[j]]['fNH3file'])
        fNH3hdulist.info()
        fNH3hdr=fNH3hdulist[0].header
        fNH3data=fNH3hdulist[0].data
        fNH3hdulist.close()
        
        fNH3_patch=make_patch(fNH3data,LatLims,LonLims,CM2,LonRng)*1.0e6
        fNH3_patch_smooth=make_patch(convolve(fNH3data,kernel,boundary='extend'),LatLims,LonLims,CM2,LonRng)*1.0e6
        if first:
            fNH3_patch_avg=np.zeros(fNH3_patch_smooth.shape)
        if i>>0:
            show=axs[1,j].imshow(fNH3_patch_smooth, clrmap, origin='upper',vmin=100,vmax=160,  
                       extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],
                               90-LatLims[0]],#vmin=0,vmax=1.2,
                               aspect="equal")
            temp=make_contours_CH4_patch(axs[1,j],fNH3_patch_smooth,LatLims,LonLims,
                                   lvls=np.linspace(100,160,num=13,endpoint=True),frmt='%3.0f')
            axs[1,j].set_xlabel("Sys. 2 Longitude (deg)",fontsize=9)
            
            sfile=sourcefiles[Dates[j]]['fNH3file']
            sec=str(int(str(sfile[16:17]))*6)
            sfiletime=(sfile[0:10]+"_"+sfile[11:13]+":"+sfile[13:15]+":"+sec.zfill(2))
            eph=get_WINJupos_ephem(sfiletime)
            ObsCM2=float(eph[1].strip())


            MeridEW=np.flip(np.mean(fNH3_patch_smooth[:,:],axis=1),axis=0)
            MeridEWerror=np.flip(np.std(fNH3_patch_smooth[:,:],axis=1),axis=0)
            Lats=np.linspace(-44.5,44.5,90)
            #print DateCounter; Date
            MeridEWArray[:,j]=MeridEW[:]

            axs[1,j].set_title("CM2 = "+str(ObsCM2),fontsize=9)
            fNH3_patch_avg=fNH3_patch_avg+fNH3_patch_smooth/7.
            first=False
            
        RGB=imread(path+sourcefiles[Dates[j]]['RGBfile'])
        RGB_patch=make_patch(RGB,LatLims,LonLims,CM2,LonRng)
        if firstRGB:
            RGB_patch_avg=np.zeros(RGB_patch.shape)
        if i==0:
            show=axs[0,j].imshow(RGB_patch, vmin=0,vmax=2^16,  
                       extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],
                               90-LatLims[0]],
                               aspect="equal")
            temp=make_contours_CH4_patch(axs[0,j],fNH3_patch_smooth,LatLims,LonLims,
                                   lvls=np.linspace(100,160,num=13,endpoint=True),frmt='%3.0f')
            RGB_patch_avg=RGB_patch_avg+RGB_patch/7.
            firstRGB=False
            axs[0,j].set_title(Dates[j],fontsize=9)

    axs[0,0].set_ylabel("Planetographic Latitude (deg)",fontsize=8)
    axs[1,0].set_ylabel("Planetographic Latitude (deg)",fontsize=8)


    fig.subplots_adjust(left=0.05, bottom=0.09, right=0.98, top=0.90,
                wspace=0.10, hspace=0.10)  
    
    pathout='c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/'

    fig.savefig(pathout+"GRS_NH3_Comparison_2022.png",dpi=300)

    ###########################################################################
    figavg,axsavg=pl.subplots(2,2,figsize=(6.0,6.0), dpi=150, facecolor="white",
                        sharey='all',sharex='col')
    figavg.suptitle('Average 2022 GRS NH3',x=0.5,ha='center',color='k')
    print(fNH3_patch_avg.shape,RGB_patch_avg.shape)
    print(RGB_patch_avg.min(),RGB_patch_avg.max())
    for i in range(0,1):
        for j in range(0,1):
            axsavg[i,j].grid(linewidth=0.2)
            #axsavg[i,j].ylim=[-45.,45.]
            #axsavg[i,j].xlim=[360-LonLims[0],360-LonLims[1]]
            print(360-LonLims[0],360-LonLims[1])
            print(LatLims)
            axsavg[i,j].set_xticks(np.linspace(450,0,31), minor=False)
            xticklabels=np.array(np.mod(np.linspace(450,0,31),360))
            axsavg[i,j].set_xticklabels(xticklabels.astype(int))

            axsavg[i,j].set_yticks(np.linspace(-45,45,7), minor=False)
            axsavg[i,j].tick_params(axis='both', which='major', labelsize=7)

    RGB_patch_avg_sharp=RGB_patch_avg
    RGB_patch_avg_sharp[:,:,0]=sharpen(RGB_patch_avg_sharp[:,:,0])
    RGB_patch_avg_sharp[:,:,1]=sharpen(RGB_patch_avg_sharp[:,:,1])
    RGB_patch_avg_sharp[:,:,2]=sharpen(RGB_patch_avg_sharp[:,:,2])
    show=axsavg[0,0].imshow(RGB_patch_avg_sharp,vmin=0,vmax=512,  
               extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],
                       90-LatLims[0]],#vmin=0,vmax=1.2,
                       aspect="equal")
    temp=make_contours_CH4_patch(axsavg[0,0],fNH3_patch_avg,LatLims,LonLims,
                           lvls=np.linspace(100,160,num=13,endpoint=True),frmt='%3.0f')
    axsavg[0,0].set_title('RGB Context',fontsize=9)

    clrmap='jet'
    show=axsavg[1,0].imshow(fNH3_patch_avg, clrmap, origin='upper',vmin=120,vmax=150,  
               extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],
                       90-LatLims[0]],#vmin=0,vmax=1.2,
                       aspect="equal")
    temp=make_contours_CH4_patch(axsavg[1,0],fNH3_patch_avg,LatLims,LonLims,
                           lvls=np.linspace(120,150,num=7,endpoint=True),frmt='%3.0f')

    cbar = pl.colorbar(show, ticks=[120,130,140,150], 
                       orientation='vertical',cmap='jet',
                       ax=axsavg[1,0],fraction=0.046, pad=0.04)
    cbar.ax.tick_params(labelsize=8,color='k')#if iSession >1:
    axsavg[1,0].set_title('Ammonia Abundance Index (ppm)',fontsize=9)

    #####Fletcher########
    FletcherRGBfile='C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Fletcher-2016-RGB-Screenshot 2022-11-26 232747.jpeg'
    FletcherRGB=imread(FletcherRGBfile)
    FletcherNH3file='C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Fletcher-2016-AmmoniaVMR-Screenshot 2022-11-26 232643.jpeg'
    FletcherNH3=imread(FletcherNH3file)
    show=axsavg[0,1].imshow(FletcherRGB,  
               extent=[148,88,90-LatLims[1],
                       90-LatLims[0]],#vmin=0,vmax=1.2,
                       aspect="equal")
    axsavg[0,1].set_title("RGB Context  \n[Fletcher et al., 2016] (adapted)",fontsize=9,wrap=True)

    show=axsavg[1,1].imshow(fNH3_patch_avg, clrmap, origin='upper',vmin=120,vmax=150,  
               extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],
                       90-LatLims[0]],#vmin=0,vmax=1.2,
                       aspect="equal")
    cbar = pl.colorbar(show, ticks=[122,130.5,139.5,148], 
                       orientation='vertical',cmap='jet',
                       ax=axsavg[1,1],fraction=0.046, pad=0.04)
    cbar.ax.tick_params(labelsize=8,color='k')#if iSession >1:
    cbar.ax.set_yticklabels(['5', '10', '15','20'],color="k")

    #####Fletcher########
        
    show=axsavg[1,1].imshow(FletcherNH3,  
               extent=[148,88,90-LatLims[1],
                       90-LatLims[0]],#vmin=0,vmax=1.2,
                       aspect="equal")
    
    axsavg[1,1].set_title("Ammonia Abundance at 500mb (ppm) \n[Fletcher et al., 2016] (adapted)",fontsize=9,wrap=True)

    #axsavg[1,1].set_xticks(np.linspace(450,0,31), minor=False)
    #xticklabels=np.array(np.mod(np.linspace(450,0,31),360))
    #axsavg[1,1].set_xticklabels(xticklabels.astype(int))
    xticklabels=np.array(np.linspace(135,90,4))
    axsavg[1,1].set_xticks(np.linspace(135,90,4), minor=False)
    axsavg[1,1].set_xticklabels(xticklabels.astype(int))             


    for i in range(0,2):
        for j in range(0,2):
            axsavg[i,j].tick_params(axis='both', which='major', labelsize=9)

    axsavg[0,0].set_ylabel("Planetographic Latitude (deg)",fontsize=9)
    axsavg[1,0].set_ylabel("Planetographic Latitude (deg)",fontsize=9)
    axsavg[1,0].set_xlabel("Sys. 2 Longitude (deg)",fontsize=9)
    axsavg[1,1].set_xlabel("Sys. 3 Longitude (deg)",fontsize=9)
    

    figavg.subplots_adjust(left=0.08, bottom=0.1, right=0.95, top=0.90,
                wspace=0.10, hspace=0.30)  
    
    
    box = axsavg[0,0].get_position()
    axsavg[0,0].set_position([box.x0+0.0, box.y0-0., box.width * 1., box.height * 1.])
    box = axsavg[0,1].get_position()
    axsavg[0,1].set_position([box.x0+0.0, box.y0-0., box.width * 1., box.height * 1.])

    figavg.savefig(pathout+"GRS_NH3_Comparison_AVG_2022.png",dpi=300)
    
    ###########################################################################
    #print DateCounter; Date
    figprof,axsprof=pl.subplots(1,1,figsize=(6.0,4.0), dpi=150, facecolor="white")
    figprof.suptitle("Ammonia Abundance Profile",x=0.5,ha='center',color='k')

    AvgMeridEW=np.mean(MeridEWArray[:,:],axis=1)
    StdMeridEW=np.std(MeridEWArray[:,:],axis=1)

    MeridEW=np.flip(np.mean(fNH3_patch_avg[:,:],axis=1),axis=0)
    MeridEWerror=np.flip(np.std(fNH3_patch_avg[:,:],axis=1),axis=0)
    #Lats=np.linspace(-44.5,44.5,90)
    Lats=np.linspace(89.5-LatLims[1],LatLims[0],LatLims[1]-LatLims[0]) #!!!! Problem! doesn't work for 
    #!!!! anything other than 45-> -45 (45,135)

    print('^^^^^^^^^',Lats)
    print('^^^^^^^^^',MeridEW)
    axsprof.plot(Lats,AvgMeridEW,label='This Work, Ammonia Abundance Index') 
    axsprof.fill_between(Lats, AvgMeridEW-StdMeridEW, AvgMeridEW+StdMeridEW,color='C0',alpha=.2)

    plevel=0.752910
    PTG.plot_TEXES_Groups(axsprof,clr='C2',prs=plevel,mult=1000000.)
    plevel=0.657540
    PTG.plot_TEXES_Groups(axsprof,clr='C3',prs=plevel,mult=1000000.)
    axsprof.set_xlim(-30.,30.)
    axsprof.set_xlabel("Planetographic Latitude (deg)",fontsize=10)
    axsprof.set_ylabel("Ammonia Abundance (ppm)",fontsize=10)
    axsprof.legend(loc=1,ncol=2, borderaxespad=0.,prop={'size':7})
    axsprof.grid(linewidth=0.2)
    axsprof.tick_params(axis='both', which='major', labelsize=8)

    figprof.subplots_adjust(left=0.10, bottom=0.12, right=0.98, top=0.92)  
    figprof.savefig(pathout+"GRS_NH3_AbundanceProfile_AVG_2022.png",dpi=300)

    #axsprof[0].fill_between(Lats, MeridEW-MeridEWerror, MeridEW+MeridEWerror,alpha=.2)
    #axsprof[0].set_title(Date+", CM2="+str(int(CM2deg))+", ASI120MM",fontsize=12)
    #axsprof[0].set_xlabel("Planetographic Latitude (deg)",fontsize=8)


def make_patch(Map,LatLims,LonLims,CM2deg,LonRng):
    """
    Purpose: Make a map patch and handle the case where the data overlap
             the map edges. This is designed for a map with Jovian longitude
             conventions that with the left boundary at 360 ascending from
             the right boundary at 0. In WinJUPOS, the actual map setting
             shows the left boundary at zero, which is of course, also 360.
    """
    import numpy as np
    patch=np.copy(Map[LatLims[0]:LatLims[1],LonLims[0]:LonLims[1]])
    if CM2deg<LonRng:
        patch=np.concatenate((np.copy(Map[LatLims[0]:LatLims[1],LonLims[0]-1:360]),
                              np.copy(Map[LatLims[0]:LatLims[1],0:LonLims[1]-360])),axis=1)
    if CM2deg>360-LonRng:
        patch=np.concatenate((np.copy(Map[LatLims[0]:LatLims[1],360+LonLims[0]:360]),
                              np.copy(Map[LatLims[0]:LatLims[1],0:LonLims[1]])),axis=1)
    patch_pad=np.pad(patch,5,mode='reflect')
    return patch

def make_contours_CH4_patch(ax,CH4Abs_conv,LatLims,LonLims,lvls=[0.71,0.73,0.75,0.77,0.79],frmt='%3.1e'):
    """
    PURPOSE: Overlay countours of NH3 absorption data on Jovian maps.
             Specifically designed for equivalent widths with mean values of
             ~0.55nm
    """
    cs=ax.contour(CH4Abs_conv,origin='upper', 
                  extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],90-LatLims[0]],
                  colors=['w','w','w','w','w'], alpha=0.5,levels=lvls,
                  linewidths=[0.5,0.5,1.0,0.5,0.5],
                  linestyles=['dashed','dashed','solid','dashed','dashed'])
    #ax.clabel(cs,[19.0,19.5,20.0,20.5,21.0],inline=True,fmt='%2.1f',fontsize=8)
    ax.clabel(cs,lvls,inline=True,fmt=frmt,fontsize=8)

def get_WINJupos_ephem(dateobs):
    # Example call: get_WINJupos_ephem('2021-09-05_04:09:00')
    import win32com.shell.shell as shell
    import time
    #shell.ShellExecuteEx(lpVerb='runas', lpFile='cmd.exe', lpParameters='/c '+commands) #run as admin
    ###########################################################################
    # WRITE *.BAT FILE WITH COMMAND SCRIPT FOR WINJUPOS
    ###########################################################################   
    batfile=open("WINJupos_CM.bat",'w')
    Line1='cd "\Program Files\WinJUPOS 12.1.1"\r\n'
    #Line2='WinJUPOS.x64.exe Jupiter /GetCM:2021-09-05_04:09:00 /GeoLong:-104.9 /GeoLat:39.7 /GetAlt >"c:\Astronomy\Projects\SAS 2021 Ammonia\Jupiter_NH3_Analysis_P3\cm.txt"\r\n'
    Line2='WinJUPOS.x64.exe Jupiter /GetCM:'+dateobs+' /GeoLong:-104.9 /GeoLat:39.7 /GetAlt >"c:\Astronomy\Projects\SAS 2021 Ammonia\Jupiter_NH3_Analysis_P3\cm.txt"\r\n'
    #batfile.writelines(["Line1\r\n","Line2\r\n"])
    batfile.writelines([Line1,Line2])
    batfile.close()
    ###########################################################################
    # EXECUTE *.BAT COMMAND FILE FOR WINJUPOS AND WAIT TO READ RESULT FILE
    ###########################################################################
    commands = "WINJupos_CM.bat"  
    shell.ShellExecuteEx(lpFile='cmd.exe', lpParameters='/c '+commands)
    time.sleep(1)
    ephemfile=open("C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/cm.txt",'r')
    Lines=ephemfile.readlines()
    ephemfile.close()
    LineString=str(Lines[0])
    print(LineString)
    ###########################################################################
    # PARSE OUTPUT FILE AND CREATE STRING ARRAY EPH FOR CM1, CM2, CM3, AND ALT
    ###########################################################################
    start=[i for i, letter in enumerate(LineString) if letter == "="]
    end=[i for i, letter in enumerate(LineString) if letter == "°"]
    eph=[]
    for i in range(0,3):
        temp=LineString[int(start[i])+1:int(end[i])]
        #print("CM"+str(i+1)+" = "+LineString[int(start[i])+1:int(end[i])])#,Linestring[start[1]:end[1]],Linestring[start[2]:end[2]])
        print("CM"+str(i+1)+" = "+temp)#,Linestring[start[1]:end[1]],Linestring[start[2]:end[2]])
        eph.extend([temp])        
    temp=LineString[int(start[3])+1:int(end[3])]
    print("Alt =  "+temp)
    eph.extend([temp])        
    return eph

def sharpen(image):
    from astropy.convolution import Gaussian2DKernel
    from astropy.convolution import convolve
    kernel = Gaussian2DKernel(x_stddev=5)
    blurred=convolve(image,kernel)
    tst=image+0.99*(image-blurred) #Need to revalidate that this is correct
    return tst