# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 07:09:10 2024

@author: smhil
"""

import ReadMUSE as RM
import MUSE_Spectrum as MS
MUSEhdr,MUSEdata,MUSEzen,MUSEszen,wavelength,filterdata,path=RM.ReadMUSE("20220730UT")
MS.MUSE_Spectrum("20220730UT",MUSEhdr,MUSEdata,wavelength,filterdata,path)

MUSEhdr,MUSEdata,MUSEzen,MUSEszen,wavelength,filterdata,path=RM.ReadMUSE("20220919UT")
MS.MUSE_Spectrum("20220919UT",MUSEhdr,MUSEdata,wavelength,filterdata,path)

MUSEhdr,MUSEdata,MUSEzen,MUSEszen,wavelength,filterdata,path=RM.ReadMUSE("20230812UT")
MS.MUSE_Spectrum("20230812UT",MUSEhdr,MUSEdata,wavelength,filterdata,path)

MUSEhdr,MUSEdata,MUSEzen,MUSEszen,wavelength,filterdata,path=RM.ReadMUSE("20230815UT")
MS.MUSE_Spectrum("20220815UT",MUSEhdr,MUSEdata,wavelength,filterdata,path)

MUSEhdr,MUSEdata,MUSEzen,MUSEszen,wavelength,filterdata,path=RM.ReadMUSE("20230816UT")
MS.MUSE_Spectrum("20220806UT",MUSEhdr,MUSEdata,wavelength,filterdata,path)

MUSEhdr,MUSEdata,MUSEzen,MUSEszen,wavelength,filterdata,path=RM.ReadMUSE("20230923UT")
MS.MUSE_Spectrum("20220923UT",MUSEhdr,MUSEdata,wavelength,filterdata,path)
