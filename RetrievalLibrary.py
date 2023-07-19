# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 14:43:04 2023

@author: smhil
"""

def make_patch(Map,LatLims,LonLims,CM2deg,LonRng,pad=True):
    """
    Purpose: Make a map patch and handle the case where the data overlap
             the map edges. This is designed for a map with Jovian longitude
             conventions that with the left boundary at 360 ascending from
             the right boundary at 0. In WinJUPOS, the actual map setting
             shows the left boundary at zero, which is of course, also 360.
    """
    import numpy as np
    patch=np.copy(Map[LatLims[0]:LatLims[1],LonLims[0]:LonLims[1]])
    print("####################### Patch shape",patch.shape)

    if CM2deg<LonRng:
        patch=np.concatenate((np.copy(Map[LatLims[0]:LatLims[1],LonLims[0]-1:360]),
                              np.copy(Map[LatLims[0]:LatLims[1],0:LonLims[1]-360])),axis=1)
    if CM2deg>360-LonRng:
        patch=np.concatenate((np.copy(Map[LatLims[0]:LatLims[1],360+LonLims[0]:360]),
                              np.copy(Map[LatLims[0]:LatLims[1],0:LonLims[1]])),axis=1)
    if pad:
        patch_pad=np.pad(patch,5,mode='reflect')
    print("####################### Patch shape",patch.shape)

    return patch

def make_contours_CH4_patch(ax,CH4Abs_conv,LatLims,LonLims,lvls=[0.71,0.73,0.75,0.77,0.79],frmt='%3.1e',clr='w'):
    """
    PURPOSE: Overlay countours of NH3 absorption data on Jovian maps.
             Specifically designed for equivalent widths with mean values of
             ~0.55nm
    """
    clrs = []
    for i in range(6):
        clrs.append(clr)
    cs=ax.contour(CH4Abs_conv,origin='upper', 
                  extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],90-LatLims[0]],
                  colors=clrs, alpha=0.5,levels=lvls,
                  linewidths=[0.5,0.5,0.5,0.5,0.5,0.5],
                  linestyles=['dashed','dashed','dashed','dashed','dashed'])
    #ax.clabel(cs,[19.0,19.5,20.0,20.5,21.0],inline=True,fmt='%2.1f',fontsize=8)
    ax.clabel(cs,lvls,inline=True,fmt=frmt,fontsize=8)
    

def make_contours_CH4(ax,CH4Abs_conv,lvls=[0.71,0.73,0.75,0.77,0.79],frmt='%3.1e'):
    """
    PURPOSE: Overlay countours of NH3 absorption data on Jovian maps.
             Specifically designed for equivalent widths with mean values of
             ~0.55nm
    """
    cs=ax.contour(CH4Abs_conv,origin='upper', 
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
    ###print(LineString)
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

def patchstats(patch):
    import numpy as np
    from numpy import nanmean,nanstd,nanmin,nanmax

    mean=np.nanmean(patch)
    stdev=np.nanstd(patch)
    minimum=np.nanmin(patch)
    maximum=np.nanmax(patch)
    
    return [mean,stdev,minimum,maximum]