# -*- coding: utf-8 -*-
"""
Created on Thu Jun 07 10:57:55 2018

NAME:       ConfigFiles.py

PURPOSE:    This module provides a base class and various child classes for 
            reading and parsing comma-delimited text configuration files. 
            The base class could also be extended for reading data files.
            
            0CLASS readurllines
            0  init
            1CLASS Observing_Conditions
            1  load_records
            2CLASS readtextfilelines
            2  init
            3CLASS Target_Parameters
            3  INIT
            3  loadtargetparams
            4CLASS measurement_list
            4  load_records
            5CLASS ObsFileNames
            5  GetFileNames
            6CLASS built_path
            6  init
            6  Spectra
            6  Equivalent_Width
            7FUNCTION MakeKeyDate

UPDATE: 2022-01-25 Converted to Python 3 (at least paritally tested)

@author: Steven Hill
"""

drive='f:'
import urllib.request, urllib.error, urllib.parse

class readurllines:
    """
    Base class to read text filtes from a website
    """
    def __init__(self,URLtoRead):
        #print(URLtoRead)
        response = urllib.request.urlopen(URLtoRead)
        temp=response.read()
        self.URLLines=temp.splitlines()
        self.nrecords=len(self.URLLines)
        self.URLtoRead=URLtoRead
        print("Read "+str(self.nrecords)+" URL records")

class Observing_Conditions(readurllines):
    """
    Read observing conditions from the Suomi-net web site
    """
    pass
    def load_records(self,DateUT):
        self.ObsDateUT=[]  #Keyword for star identification
        self.PWV=[]           #Target, e.g., component of a multiple star
        self.Press=[]           #Target, e.g., component of a multiple star
        self.TempC=[]
        self.RelHum=[]
        for recordindex in range(1,self.nrecords):
            RecordDate=self.URLLines[recordindex][5:21]
            if RecordDate[0:10] == DateUT:
                self.ObsDateUT.append(RecordDate+":00UT")
                self.PWV.append(str(self.URLLines[recordindex][22:27]))
                self.Press.append(str(self.URLLines[recordindex][53:61]))
                self.TempC.append(str(self.URLLines[recordindex][62:66]))                
                self.RelHum.append(str(self.URLLines[recordindex][67:72]))                

class readtextfilelines:
    def __init__(self,FiletoRead):
        """
        Base class for reading local text files
        If no path is given, this class reads from the working directory.
        """
        CfgFile=open(FiletoRead,'r')
        self.CfgLines=CfgFile.readlines()
        CfgFile.close()
        self.nrecords=len(self.CfgLines)
        self.FiletoRead=FiletoRead
        print("Read "+str(self.nrecords)+" records")

class Target_Parameters(readtextfilelines):
    """
    This class reads target attribute parameters. It is aimed at deep sky
    objects such as galaxies: position angle, inclination, and classification.  
    """
    pass
    def loadtargetparams(self,TargetID):

        for recordindex in range(1,self.nrecords):
            fields=self.CfgLines[recordindex].split(',')
            print(fields[0], fields[1])
            if fields[0] == TargetID:
                self.TargetID=fields[0]
                self.TargetType=fields[1]
                self.RA2000=float(fields[2])
                self.DE2000=float(fields[3])
                self.Distance=float(fields[4])
                self.DistUnits=fields[5]
                if self.TargetType=="Stars" or self.TargetType=="Planets":
                    self.mV=float(fields[6])
                    self.SpecType=fields[7]
                elif self.TargetType=="Galaxies":
                    self.PA=float(fields[6])
                    self.Incl=float(fields[7])
                    self.AppSize=float(fields[8])

class measurement_list(readtextfilelines):
    """
    This class is used to load spectral observation metadata.
    """
    pass
    
    def load_records(self,MeasTgt="All",DateUTSelect="All",Grating="All"):

        print("Hi in measurement_list>load_records")
        self.MeasTarget=[]  #Keyword for star identification
        self.DataType=[]           #Target, e.g., component of a multiple star
        self.DataTarget=[]           #Target, e.g., component of a multiple star
        self.DateUT=[]           #UT Date of observation: YYYYMMDDUT
        self.Optics=[]       #Instrument code, to be used for aperture
        self.Camera=[]       #Instrument code, to be used for aperture
        self.Grating=[]    #Grating 100lpm or 200lpm or None
        self.FileList=[]         #List of observation image files (FITS)
        for recordindex in range(1,self.nrecords):
            fields=self.CfgLines[recordindex].split(',')
            if MeasTgt=="All" or MeasTgt==str(fields[0]):
                if Grating=="All" or Grating==str(fields[6]):
                    if DateUTSelect=="All" or DateUTSelect==str(fields[3]):
                        self.MeasTarget.extend([str(fields[0])])
                        self.DataType.extend([str(fields[1])])
                        self.DataTarget.extend([str(fields[2])])
                        self.DateUT.extend([str(fields[3])])
                        self.Optics.extend([str(fields[4])])
                        self.Camera.extend([str(fields[5])])
                        self.Grating.extend([str(fields[6])])
                        self.FileList.extend([str(fields[7])])
                        
class meas_extend_test(measurement_list):
    pass

    def load_extra_field(self,MeasTgt="All",DateUTSelect="All",Grating="All"):
        print("Hi in meas_extend_test")
        self.extra_field=[]
        for recordindex in range(1,self.nrecords):
            fields=self.CfgLines[recordindex].split(',')
            if MeasTgt=="All" or MeasTgt==str(fields[0]):
                if Grating=="All" or Grating==str(fields[6]):
                    if DateUTSelect=="All" or DateUTSelect==str(fields[3]):
                        self.extra_field.extend([str(fields[8])])
        

class ObsFileNames(readtextfilelines):
    #Used to be the function GetObsFileNames in SRL
    pass
    def GetFileNames(self):
        self.FNArray=[]
        for recordindex in range(0,self.nrecords):
            fields=str(self.CfgLines[recordindex]).split(',')
            self.FNArray.append(fields[0])

class built_path:
    def __init__(self,TargetParams):
        self.drive="f:"
        self.base_path='/Astronomy/Projects/'
        self.base_path=self.drive+self.base_path+TargetParams.TargetType+'/'+TargetParams.TargetID+'/'
        
    def spectra(self,dateUT):
        self.input_path=self.base_path+'Spectral Data/'+dateUT+'/'
        self.One_D_path=self.base_path+'Spectral Data/1D Spectra/'
        self.EW_path=self.base_path+'Spectral Data/EWs/'
        self.config_path=self.base_path+'Spectral Data/1D Spectra/Configuration Files/'
        self.reference_path=self.drive+'/Astronomy/Python Play/SPLibraries/SpectralReferenceFiles/ReferenceLibrary/'
        
    def Equivalent_Width(self):
        self.input_path=self.base_path+'Spectral Data/1D Spectra/'
        self.output_path=self.base_path+'Spectral Data/EWs/'
        
def MakeKeyDate(FN):
    import datetime

    fields=FN.split('-')
    DT=fields[1]
    Key=FN[0:7]+DT
    print("DT[8]=",DT[8])
    if DT[8]=="T":
        print("In first branch")
        DateTime=datetime.datetime.strptime(DT[0:4]+"-"+DT[4:6]+"-" \
            +DT[6:8]+"T"+DT[9:11]+":"+DT[11:13]+":"+DT[13:15], \
            '%Y-%m-%dT%H:%M:%S')
    else:
        print("In second branch")
        DateTime=datetime.datetime.strptime(DT[0:4]+"-"+DT[4:6]+"-" \
            +DT[6:8]+"T"+DT[8:10]+":"+DT[10:12]+":"+DT[12:14], \
            '%Y-%m-%dT%H:%M:%S')

    return Key,DateTime

class video_metadata_list(readtextfilelines):
    """
    This class is used to load spectral observation metadata.
    """
    pass
    
    def load_records(self):

        #print "Hi in video_metadata_list>load_records"
        self.VideoFile=[]  
        self.MetaFile=[]           
        self.TimeStamp=[]          
        self.Exposure=[]           
        self.Gain=[]      
        self.Gamma=[]      
        self.Binning=[]    
        self.CaptureArea=[]        
        for recordindex in range(1,self.nrecords):
            fields=self.CfgLines[recordindex].split(',')
            self.VideoFile.extend([str(fields[0])])
            self.MetaFile.extend([str(fields[1])])
            self.TimeStamp.extend([str(fields[2])])
            self.Exposure.extend([str(fields[3])])
            self.Gain.extend([str(fields[4])])
            self.Gamma.extend([str(fields[5])])
            self.Binning.extend([str(fields[6])])
            self.CaptureArea.extend([str(fields[7])])
