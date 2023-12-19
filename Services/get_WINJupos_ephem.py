# -*- coding: utf-8 -*-


def get_WINJupos_ephem(dateobs):
    """
    Created on Fri Sep 15 07:57:17 2023

    @author: smhil
    """
    # Example call: get_WINJupos_ephem('2021-09-05_04:09:00')
    import win32com.shell.shell as shell
    import time
    #shell.ShellExecuteEx(lpVerb='runas', lpFile='cmd.exe', lpParameters='/c '+commands) #run as admin
    ###########################################################################
    # WRITE *.BAT FILE WITH COMMAND SCRIPT FOR WINJUPOS
    ###########################################################################   
    batfile=open("WINJupos_CM.bat",'w')
    Line1='cd "\Program Files\WinJUPOS 12.2.6"\r\n'
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
    time.sleep(0.5)
    ephemfile=open("C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/cm.txt",'r')
    Lines=ephemfile.readlines()
    ephemfile.close()
    LineString=str(Lines[0])
    #print(LineString)
    ###########################################################################
    # PARSE OUTPUT FILE AND CREATE STRING ARRAY EPH FOR CM1, CM2, CM3, AND ALT
    ###########################################################################
    start=[i for i, letter in enumerate(LineString) if letter == "="]
    end=[i for i, letter in enumerate(LineString) if letter == "°"]
    eph=[]
    for i in range(0,3):
        temp=LineString[int(start[i])+1:int(end[i])]
        #print("CM"+str(i+1)+" = "+LineString[int(start[i])+1:int(end[i])])#,Linestring[start[1]:end[1]],Linestring[start[2]:end[2]])
        #print("CM"+str(i+1)+" = "+temp)#,Linestring[start[1]:end[1]],Linestring[start[2]:end[2]])
        eph.extend([temp])        
    temp=LineString[int(start[3])+1:int(end[3])]
    #print("Alt =  "+temp)
    eph.extend([temp])        
    return eph