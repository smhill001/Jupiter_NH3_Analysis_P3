# -*- coding: utf-8 -*-
"""
Created on Wed Oct 06 14:36:01 2021

@author: Steven Hill

PURPOSE:    Read a white balanced RGB (656,647,632) image and compute the
            color slope and NH3 absorption images.

EXAMPLE:    filename="2021-09-05-0409_0-Jupiter-R656G647B632-ReAligned-WhtBal.png"
            filename2="2021-09-05-0432_8-Jupiter_NoWV-R940(G)B889-R(G)B-WhtBal.png"
"""

def Jupiter_NH3_CH4_Images_Arrays_P3(filename="2021-09-05-0409_0-Jupiter-R656G647B632-ReAligned-WhtBal.png",
                           filename2="2021-09-05-0432_8-Jupiter_NoWV-R940(G)B889-R(G)B-WhtBal.png",
                           obsdate="20210905UT"):
    import sys
    drive='c:'
    sys.path.append(drive+'/Astronomy/Python Play')
    sys.path.append(drive+'/Astronomy/Python Play/Util_P3')
    sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy')

    import os
    from matplotlib.pyplot import imread
    import pylab as pl
    import numpy as np
    
    ###########################################################################
    # Here's where we filter and select the images to plot. Probably want
    # more than one option, e.g.:
    #   - For posting to ALPO Japan or other public sites
    #   - For basic array of monchromatic derotated images, both with and 
    #     without wavelets processing
    #   - RGB array, e.g., RGB, 889-G-NUV, 656-647-632
    #   - Array corresponding to Maps presented for NH3 analysis
    ###########################################################################

    path='c:/Astronomy/Projects/Planets/Jupiter/Imaging Data/'+obsdate+'/'
    print(path)
    fn=path+filename
    print(fn)
    imageRGB=imread(fn)
    print(path)
    fn2=path+filename2
    print(fn2)
    imageRGB2=imread(fn2)
    
    indices=(imageRGB>0.1)
    mask=np.zeros(imageRGB.shape)
    print(mask.shape)
    print(indices)
    mask[indices]=1.
    
    fig,ax=pl.subplots(3,3,figsize=(6.0,6.0), dpi=150, facecolor="black",
                          sharex=True)
    fig.suptitle("Jovian Methane and Ammonia "+filename[0:10],
                 x=0.05,ha='left',color='w')
    ax[0,0].annotate("C11, FL=2800mm, ST2000XM, Denver, CO; Copyright Steven Hill",[0.10,1.1],
                 ha='left',xycoords='axes fraction',color='w',fontsize=10)
    
    # AMMONIA PANELS (6)
    
    eph=get_WINJupos_ephem('2021-09-05_04:09:00')
    AmmoniaHeader="AMMONIA "+filename[11:15]+"UT; CM1"+eph[0]+"; CM2"+eph[1]+"; CM3"+eph[2]+"; Alt"+eph[3]

    ax[0,0].imshow(imageRGB[:,:,0],'gist_gray')
    ax[0,0].annotate('656nm (Red cont.)',[0.5,0.03],color='white',ha='center',
                    xycoords='axes fraction',fontsize=8)
    ax[0,1].imshow(imageRGB[:,:,1],'gist_gray')
    ax[0,1].annotate('647nm (NH3)',[0.5,0.03],color='white',ha='center',
                    xycoords='axes fraction',fontsize=8)
    ax[0,2].imshow(imageRGB[:,:,2],'gist_gray')
    ax[0,2].annotate('632nm (Blue cont.)',[0.5,0.03],color='white',ha='center',
                    xycoords='axes fraction',fontsize=8)
    clrslp=np.array(imageRGB[:,:,0])/np.array(imageRGB[:,:,2]) #This is a ratio, not a slope
    ax[1,0].imshow(clrslp*mask[:,:,1],'gist_gray',vmin=0.9,vmax=1.1)
    ax[1,0].annotate('656/632 (Color Slope)',[0.5,0.03],color='white',ha='center',
                    xycoords='axes fraction',fontsize=8)
    ax[1,0].annotate('Scale: 0.9 -> 1.1',[0.5,-0.04],color='white',ha='center',
                    xycoords='axes fraction',fontsize=8)
    CNT647=(15./24.)*np.array(imageRGB[:,:,2])+(9./24.)*np.array(imageRGB[:,:,0])
    CNT647a=(12./24.)*np.array(imageRGB[:,:,2])+(12./24.)*np.array(imageRGB[:,:,0])
    ax[1,1].imshow(CNT647,'gist_gray')
    ax[1,1].annotate('Synth. Continuum @ 647nm',[0.5,0.03],color='white',ha='center',
                    xycoords='axes fraction',fontsize=8)
    nh3abs=np.array(imageRGB[:,:,1])/CNT647
    ax[1,2].imshow(nh3abs*mask[:,:,1],'gist_gray',vmin=0.9,vmax=1.1)
    ax[1,2].annotate('647/Synth. Continuum',[0.51,0.03],color='white',ha='center',
                    xycoords='axes fraction',fontsize=8)
    ax[1,2].annotate('Scale: 0.9 -> 1.1',[0.5,-0.04],color='white',ha='center',
                    xycoords='axes fraction',fontsize=8)
    
    ax[0,1].annotate(AmmoniaHeader,[0.5,0.93],color='white',ha='center',
                    xycoords='axes fraction',fontsize=9)
    ax[0,1].set_zorder(1) #required so not blocked by the ax[0,2] image
    #see https://stackoverflow.com/questions/29735743/getting-text-to-display-in-front-of-subplot-images

    # METHANE PANELS (3)

    eph=get_WINJupos_ephem('2021-09-05_04:33:00')
    MethaneHeader="METHANE "+filename2[11:15]+"UT; CM1"+eph[0]+"; CM2"+eph[1]+"; CM3"+eph[2]+"; Alt"+eph[3]

    ax[2,0].imshow(imageRGB2[:,:,0],'gist_gray')
    ax[2,0].annotate('940nm (Continuum)',[0.51,0.03],color='white',ha='center',
                    xycoords='axes fraction',fontsize=8)
    ax[2,1].imshow(imageRGB2[:,:,2],'gist_gray')
    ax[2,1].annotate('889nm (CH4)',[0.5,0.03],color='white',ha='center',
                    xycoords='axes fraction',fontsize=8)
    ch4abs=np.array(imageRGB2[:,:,2])/np.array(imageRGB2[:,:,0])
    ax[2,2].imshow(ch4abs*mask[:,:,1],'gist_gray',vmin=0.5,vmax=2.0)
    ax[2,2].annotate('889/940 Continuum',[0.5,0.03],color='white',ha='center',
                    xycoords='axes fraction',fontsize=8)
    ax[2,2].annotate('Scale: 0.5 -> 2.0',[0.5,-0.04],color='white',ha='center',
                    xycoords='axes fraction',fontsize=8)
    
    ax[2,1].annotate(MethaneHeader,[0.5,0.93],color='white',ha='center',
                    xycoords='axes fraction',fontsize=9)
    ax[2,1].set_zorder(1) #required so not blocked by the ax[2,2] image
    #see https://stackoverflow.com/questions/29735743/getting-text-to-display-in-front-of-subplot-images


    for i in range(0,3):
        for j in range(0,3):
            print(i,j)
            ax[i,j].axis('off')

    fig.subplots_adjust(left=0.00, bottom=0.03, right=1.0, top=0.90,
                wspace=0.001, hspace=0.05)

    fig.savefig(path+"PublicImageArray"+obsdate+".png",dpi=150)
    
def load_png(file_path):
    """
    Purpose: Properly load a 48-bit PNG file
    Read from KITTI .png file
    Args:
        file_path string: file path(absolute)
    Returns:
        data (numpy.array): data of image in (Height, Width, 3) layout
    
    FROM: https://www.programcreek.com/python/example/98900/png.Reader
    """
    import png
    import numpy as np

    flow_object = png.Reader(filename=file_path)
    flow_direct = flow_object.asDirect()
    flow_data = list(flow_direct[2])
    (w, h) = flow_direct[3]['size']

    flow = np.zeros((h, w, 3), dtype=np.float64)
    for i in range(len(flow_data)):
        flow[i, :, 0] = flow_data[i][0::3]
        flow[i, :, 1] = flow_data[i][1::3]
        flow[i, :, 2] = flow_data[i][2::3]

    return flow.astype(np.uint16) 

def get_WINJupos_ephem(dateobs):
    # Example call: get_WINJupos_ephem('2021-09-05_04:09:00')
    import win32com.shell.shell as shell
    import time
    #shell.ShellExecuteEx(lpVerb='runas', lpFile='cmd.exe', lpParameters='/c '+commands) #run as admin
    batfile=open("WINJupos_CM.bat",'w')
    Line1='cd "\Program Files\WinJUPOS 12.0.11"\r\n'
    #Line2='WinJUPOS.x64.exe Jupiter /GetCM:2021-09-05_04:09:00 /GeoLong:-104.9 /GeoLat:39.7 /GetAlt >"c:\Astronomy\Projects\SAS 2021 Ammonia\Jupiter_NH3_Analysis_P3\cm.txt"\r\n'
    Line2='WinJUPOS.x64.exe Jupiter /GetCM:'+dateobs+' /GeoLong:-104.9 /GeoLat:39.7 /GetAlt >"c:\Astronomy\Projects\SAS 2021 Ammonia\Jupiter_NH3_Analysis_P3\cm.txt"\r\n'
    #batfile.writelines(["Line1\r\n","Line2\r\n"])
    batfile.writelines([Line1,Line2])
    batfile.close()
    
    commands = "WINJupos_CM.bat"  
    shell.ShellExecuteEx(lpFile='cmd.exe', lpParameters='/c '+commands)
    time.sleep(1)
    ephemfile=open("C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/cm.txt",'r')
    Lines=ephemfile.readlines()
    ephemfile.close()
    LineString=str(Lines[0])
    print(LineString)
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
    

    
