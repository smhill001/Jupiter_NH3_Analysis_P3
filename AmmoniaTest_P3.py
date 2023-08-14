# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 12:38:27 2021

    This code has been significantly overhauled and expanded since its debut
    in March 2021. It is a single script that performs multiple functions, and 
    several of the original functions have mostly been moved to the library
    NH3_Filter_Library_P3.py. The list of functions is as follows:
        
        COMPUTE and PLOT FILTER TRANSMISSIONS CONVOLVED WITH DISK-INTEGRATED 
            ALBEDO AND CONTINUUM
        COMPUTE and PLOT K_eff, l_eff, AND WEIGHTING FUNCTION CALCULATIONS
        COMPUTE and PLOT **RAYLEIGH CANCELING** AND GAS ABSORPTION ONLY
            WEIGHTING FUNCTIONS
        COMPUTE ESTIMATED JUPITER ABSORPTION USING GALILEAN MOONS


@author: Steven Hill
"""
import sys
sys.path.append('c:/Astronomy/Python Play')
sys.path.append('c:/Astronomy/Python Play/Util_P3')
sys.path.append('c:/Astronomy/Python Play/SPLibraries_P3')
sys.path.append('c:/Astronomy/Python Play/SpectroPhotometry/Spectroscopy_P3')
import matplotlib.pyplot as pl
import numpy as np
from copy import deepcopy
import GeneralSpecUtils_P3 as GSU
import NH3_Filter_Library_P3 as NFL
sys.path.append('./Services')
import get_albedo_continua_crossections as gACC

# LOAD JOVIAN DISK-INTEGRATEDALBEDO DATA FROM KARKOSCHKA, 1994 (DATA FROM 1993)
x0,x1,xtks=600.,680.,9
y0,y1,ytks=0.0,0.7,8
Albedo,Continua,CH4,NH3,LONH3=gACC.get_albedo_continua_crossections(x0,x1,xtks,y0,y1,ytks,
                                                              Crossect=True)

# PLOT FILTER TRANSMISSIONS CONVOLVED WITH DISK-INTEGRATED ALBEDO AND CONTINUUM
FilterList=['620','632','647','656']
Model='Piecewise1'
Tele='SCT'
Continuum_Albedo=np.zeros((Continua[Model]['WaveGrid'].size,2))
Continuum_Albedo[:,0]=Continua[Model]['WaveGrid']
Continuum_Albedo[:,1]=Continua[Model]['Albedo']

filterdata,axsFilt=NFL.compute_filter_spectrum(x0,x1,xtks,y0,y1,ytks,FilterList,
                                               Albedo,Continuum_Albedo,
                                               Model,Telescope=Tele)

Jupiterdata = deepcopy(filterdata)
Iodata = deepcopy(filterdata)
Europadata = deepcopy(filterdata)
Ganymededata = deepcopy(filterdata)
Callistodata = deepcopy(filterdata)
MoonsAvgdata = deepcopy(filterdata)

# COMPUTE K_eff, l_eff, AND PLOT VERTICAL TRANSMISSIONS AND WEIGHTING FUNCTION
# PROFILES
fout_sfx='-'+Tele+'-'+Model
P=np.geomspace(10.,1.0e7,num=70,endpoint=True,dtype=float) #in Pascals
Jupiterdata=NFL.compute_vertical_transmission_profiles(Jupiterdata,FilterList,
                                                       CH4,NH3,P,fout_sfx=fout_sfx)

###############################################################################
# BEGIN CALCULATIONS FOR **RAYLEIGH CANCELING** AND GAS ABSORPTION ONLY
#   WEIGHTING FUNCTIONS
###############################################################################

fig_ray,axs_ray=NFL.vert_profile_quad_plot(SupTitle="Two-way Transmission with Rayleigh Subtraction")
fig_raycont,axs_raycont=NFL.vert_profile_quad_plot(SupTitle="Weighting with Rayleigh Subtraction")

filterlistshort=['620','632','647','656']
for filtr in filterlistshort:
    tau_gas=Jupiterdata[filtr]['tau_CH4']+Jupiterdata[filtr]['tau_NH3']
    tmp=NFL.Compute_Transmission(P,Jupiterdata[filtr]['tau_R'],tau_gas,Jupiterdata[filtr]['filtname']+' Ray',axs_ray[0,0],axs_raycont[0,0])
    if filtr=='647':
        tmp=NFL.Compute_Transmission(P,Jupiterdata[filtr]['tau_R']*0.0,Jupiterdata[filtr]['tau_NH3'],Jupiterdata[filtr]['filtname']+' NH3',axs_ray[1,1],axs_raycont[1,1])
    elif filtr=='620':
        tmp=NFL.Compute_Transmission(P,Jupiterdata[filtr]['tau_R']*0.0,Jupiterdata[filtr]['tau_CH4'],Jupiterdata[filtr]['filtname']+' CH4',axs_ray[0,1],axs_raycont[0,1])

###############################################################################
#!!!! This is a simplified approach only looking at the Rayleigh scattering
#!!!! In reality, we're doing the slope that includes both Rayleigh
#!!!! and gas extinction, but we're assuming that gas extinction is minor.
#!!!! So I should be looking at both cases, the ideal and the real.
###############################################################################
tmp=NFL.Compute_Transmission(P,Jupiterdata['620']['tau_R'],tau_gas*0.,Jupiterdata['620']['filtname']+' Ray',axs_ray[1,0],axs_raycont[1,0])
tmp=NFL.Compute_Transmission(P,Jupiterdata['647']['tau_R'],tau_gas*0.,Jupiterdata['647']['filtname']+' Ray',axs_ray[1,0],axs_raycont[1,0])


NH3RaySlp=(Jupiterdata['656']['tau_R']-Jupiterdata['632']['tau_R'])/24.0 
tau_R647_Estimated=15.0*NH3RaySlp+Jupiterdata['632']['tau_R']

tau_All656=Jupiterdata['656']['tau_CH4']+Jupiterdata['656']['tau_NH3']+Jupiterdata['656']['tau_R']
tau_All632=Jupiterdata['632']['tau_CH4']+Jupiterdata['632']['tau_NH3']+Jupiterdata['632']['tau_R']
NH3Slp=(tau_All656-tau_All632)/24.0 
tau_647_Estimated=15.0*NH3Slp+tau_All632

tmp=NFL.Compute_Transmission(P,tau_R647_Estimated,Jupiterdata['647']['tau_NH3']*0.0,"647RayEstimated",axs_ray[1,0],axs_raycont[1,0])
tmp=NFL.Compute_Transmission(P,Jupiterdata['647']['tau_R']-tau_R647_Estimated,Jupiterdata['647']['tau_NH3']*0.0,"647Ray-647RayEst",axs_ray[1,0],axs_raycont[1,0])
tmp=NFL.Compute_Transmission(P,Jupiterdata['647']['tau_R']-tau_R647_Estimated,Jupiterdata['647']['tau_NH3']*0.0,"647Ray-647RayEst",axs_ray[1,1],axs_raycont[1,1])

tmp=NFL.Compute_Transmission(P,tau_647_Estimated,Jupiterdata['647']['tau_NH3']*0.0,"647Estimated",axs_ray[1,0],axs_raycont[1,0])
tmp=NFL.Compute_Transmission(P,Jupiterdata['647']['tau_R']-tau_647_Estimated,Jupiterdata['647']['tau_NH3']*0.0,"647Ray-647Est",axs_ray[1,0],axs_raycont[1,0])
tmp=NFL.Compute_Transmission(P,Jupiterdata['647']['tau_R']-tau_647_Estimated,Jupiterdata['647']['tau_NH3']*0.0,"647Ray-647Est",axs_ray[1,1],axs_raycont[1,1])

tau_R620_Estimated=-12.0*NH3RaySlp+NH3RaySlp+Jupiterdata['632']['tau_R']
tau_620_Estimated=-12.0*NH3Slp+tau_All632

tmp=NFL.Compute_Transmission(P,tau_R620_Estimated,Jupiterdata['620']['tau_CH4']*0.0,"620RayEstimated",axs_ray[1,0],axs_raycont[1,0])
tmp=NFL.Compute_Transmission(P,Jupiterdata['620']['tau_R']-tau_R620_Estimated,Jupiterdata['620']['tau_CH4']*0.0,"620Ray-620RayEst",axs_ray[1,0],axs_raycont[1,0])
tmp=NFL.Compute_Transmission(P,Jupiterdata['620']['tau_R']-tau_R620_Estimated,Jupiterdata['620']['tau_CH4']*0.0,"620Ray-620RayEst",axs_ray[0,1],axs_raycont[0,1])

tmp=NFL.Compute_Transmission(P,tau_620_Estimated,Jupiterdata['620']['tau_CH4']*0.0,"620Estimated",axs_ray[1,0],axs_raycont[1,0])
tmp=NFL.Compute_Transmission(P,Jupiterdata['620']['tau_R']-tau_620_Estimated,Jupiterdata['620']['tau_CH4']*0.0,"620Ray-620Est",axs_ray[1,0],axs_raycont[1,0])
tmp=NFL.Compute_Transmission(P,Jupiterdata['620']['tau_R']-tau_620_Estimated,Jupiterdata['620']['tau_CH4']*0.0,"620Ray-620Est",axs_ray[0,1],axs_raycont[0,1])

axs_ray[0,0].legend(loc=1,ncol=2, borderaxespad=0.,prop={'size':6})
axs_ray[0,1].legend(loc=1,ncol=2, borderaxespad=0.,prop={'size':6})
axs_ray[1,0].legend(loc=1,ncol=2, borderaxespad=0.,prop={'size':6})
axs_ray[1,1].legend(loc=1,ncol=2, borderaxespad=0.,prop={'size':6})
fig_ray.subplots_adjust(left=0.12, right=0.96, top=0.90, bottom=0.09)

axs_raycont[0,0].legend(loc=1,ncol=2, borderaxespad=0.,prop={'size':6})
axs_raycont[0,1].legend(loc=1,ncol=2, borderaxespad=0.,prop={'size':6})
axs_raycont[1,0].legend(loc=1,ncol=2, borderaxespad=0.,prop={'size':6})
axs_raycont[1,1].legend(loc=1,ncol=2, borderaxespad=0.,prop={'size':6})
fig_raycont.subplots_adjust(left=0.12, right=0.96, top=0.90, bottom=0.09)

########## END OF THIRD FUNCTION AND PLOT ##############

###############################################################################
# COMPUTE ESTIMATED JUPITER ABSORPTION USING GALILEAN MOONS
#   !!!! COULD MAKE THIS WHOLE THING LOOP-ABLE AND SELECTABLE WITH A CALL-LIST
#   !!!!   AND/OR DICTIONARY
###############################################################################
#!!!!!!!!Should do leading and trailing edge of Io (and any other varibility I can find?)
Io_Grid=NFL.MoonAlbedos("Io")
EurGrid=NFL.MoonAlbedos('Europa')
GanGrid=NFL.MoonAlbedos('Ganymede')
CalGrid=NFL.MoonAlbedos('Callisto')

fig_moons,axs_moons=pl.subplots(2,1,figsize=(6.0,6.0), dpi=150, facecolor="white",
                      sharex=True)
#Plot Layout Configuration
x0,x1,xtks=600.,680.,9
y0,y1,ytks=0.0,1.2,7
#fig1.suptitle("Ammonia Filters",x=0.5,ha='center',color='k')
axs_moons[0].set_xlim(x0,x1)
# Set x ticks
axs_moons[0].set_xticks(np.linspace(x0,x1,xtks, endpoint=True))
# Set y limits
axs_moons[0].set_ylim(y0,y1)
axs_moons[0].set_yticks(np.linspace(y0,y1,ytks, endpoint=True))
# Set y ticks
axs_moons[0].grid(linewidth=0.2)
axs_moons[0].tick_params(axis='both', which='major', labelsize=8)
axs_moons[0].set_ylabel("Albedo",color="black")

axs_moons[0].set_title("Moons Albedo")

axs_moons[0].plot(CalGrid[:,0],CalGrid[:,1],label='Callisto',linewidth=1,color='b')
axs_moons[0].plot(GanGrid[:,0],GanGrid[:,1],label='Ganymede',linewidth=1,color='g')
axs_moons[0].plot(EurGrid[:,0],EurGrid[:,1],label='Europa',linewidth=1,color='r')
axs_moons[0].plot(Io_Grid[:,0],Io_Grid[:,1],label='Io Leading',linewidth=1,linestyle='dashed',color='k')
#axs_moons[0].plot(Io_trailing1980[:,0]*1000.,Io_trailing1980[:,1],label='Io Trailing',linewidth=1,color='k')
axs_moons[0].legend(fontsize=7)

axs_moons[1].set_xlim(x0,x1)
# Set x ticks
axs_moons[1].set_xticks(np.linspace(x0,x1,xtks, endpoint=True))
# Set y limits
axs_moons[1].set_ylim(y0,y1)
axs_moons[1].set_yticks(np.linspace(y0,y1,ytks, endpoint=True))
# Set y ticks
axs_moons[1].grid(linewidth=0.2)
axs_moons[1].tick_params(axis='both', which='major', labelsize=8)
axs_moons[1].set_ylabel("Albedo x Transmission",color="black")
axs_moons[1].set_xlabel("Wavelength (nm)")
axs_moons[1].set_title("Convolution with Filters")

firstflag=True
for filtr in FilterList:
    Iodata[filtr]['ContProd']=GSU.SpectrumMath(Jupiterdata[filtr]['FiltTrans'],Io_Grid,"Multiply")
    Europadata[filtr]['ContProd']=GSU.SpectrumMath(Jupiterdata[filtr]['FiltTrans'],EurGrid,"Multiply")
    Ganymededata[filtr]['ContProd']=GSU.SpectrumMath(Jupiterdata[filtr]['FiltTrans'],GanGrid,"Multiply")
    Callistodata[filtr]['ContProd']=GSU.SpectrumMath(Jupiterdata[filtr]['FiltTrans'],CalGrid,"Multiply")
    if firstflag:
        axs_moons[1].plot(Iodata[filtr]['ContProd'][:,0],Iodata[filtr]['ContProd'][:,1],linewidth=0.5,color='k',label='Io')
        axs_moons[1].plot(Europadata[filtr]['ContProd'][:,0],Europadata[filtr]['ContProd'][:,1],linewidth=0.5,color='r',label='Europa')
        axs_moons[1].plot(Ganymededata[filtr]['ContProd'][:,0],Ganymededata[filtr]['ContProd'][:,1],linewidth=0.5,color='g',label='Ganymede')
        axs_moons[1].plot(Callistodata[filtr]['ContProd'][:,0],Callistodata[filtr]['ContProd'][:,1],linewidth=0.5,color='b',label='Callisto')
        firstflag=False
    else:
        axs_moons[1].plot(Iodata[filtr]['ContProd'][:,0],Iodata[filtr]['ContProd'][:,1],linewidth=0.5,color='k')
        axs_moons[1].plot(Europadata[filtr]['ContProd'][:,0],Europadata[filtr]['ContProd'][:,1],linewidth=0.5,color='r')
        axs_moons[1].plot(Ganymededata[filtr]['ContProd'][:,0],Ganymededata[filtr]['ContProd'][:,1],linewidth=0.5,color='g')
        axs_moons[1].plot(Callistodata[filtr]['ContProd'][:,0],Callistodata[filtr]['ContProd'][:,1],linewidth=0.5,color='b')

axs_moons[1].legend(fontsize=8)
fig_moons.subplots_adjust(left=0.10, right=0.90, top=0.94, bottom=0.09)

fig_moons.savefig('c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/ColorSlopes.png',dpi=320)

for filt in filterlistshort:
    print('^^^^^^^^^^^^^'+filt)

    Iodata[filt]['Cont_Int'],Iodata[filt]['Absr_Int'],Iodata[filt]['TransInt']= \
        NFL.cont_absorption_calcs(Iodata[filt]['ContProd'],Jupiterdata[filt]['AbsrProd'], \
                                  float(filt)-Jupiterdata[filt]['halfwdth'],\
                                  float(filt)+Jupiterdata[filt]['halfwdth'], \
                                      Jupiterdata[filt]['filtname'],prn=False)
    Europadata[filt]['Cont_Int'],Europadata[filt]['Absr_Int'],Europadata[filt]['TransInt']= \
        NFL.cont_absorption_calcs(Europadata[filt]['ContProd'],Jupiterdata[filt]['AbsrProd'], \
                                  float(filt)-Jupiterdata[filt]['halfwdth'],\
                                  float(filt)+Jupiterdata[filt]['halfwdth'], \
                                      Jupiterdata[filt]['filtname'],prn=False)
    Ganymededata[filt]['Cont_Int'],Ganymededata[filt]['Absr_Int'],Ganymededata[filt]['TransInt']= \
        NFL.cont_absorption_calcs(Ganymededata[filt]['ContProd'],Jupiterdata[filt]['AbsrProd'], \
                                  float(filt)-Jupiterdata[filt]['halfwdth'],\
                                  float(filt)+Jupiterdata[filt]['halfwdth'], \
                                      Jupiterdata[filt]['filtname'],prn=False)
    Callistodata[filt]['Cont_Int'],Callistodata[filt]['Absr_Int'],Callistodata[filt]['TransInt']= \
        NFL.cont_absorption_calcs(Callistodata[filt]['ContProd'],Jupiterdata[filt]['AbsrProd'], \
                                  float(filt)-Jupiterdata[filt]['halfwdth'],\
                                  float(filt)+Jupiterdata[filt]['halfwdth'], \
                                      Jupiterdata[filt]['filtname'],prn=False)


    print(Iodata[filt]['Cont_Int'],Iodata[filt]['Absr_Int'],Iodata[filt]['TransInt'])
    print(Europadata[filt]['Cont_Int'],Europadata[filt]['Absr_Int'],Europadata[filt]['TransInt'])
    print(Ganymededata[filt]['Cont_Int'],Ganymededata[filt]['Absr_Int'],Ganymededata[filt]['TransInt'])
    print(Callistodata[filt]['Cont_Int'],Callistodata[filt]['Absr_Int'],Callistodata[filt]['TransInt'])

########## END OF FOURTH FUNCTION AND PLOT ####################################
###############################################################################
# MODEL 3 FOR SCT DATA:
#  THIS MODEL IS CLOSEST TO THE REAL MEASUREMENT. IT TAKES THE RATIOS OF THE
#  REFERENCE ALBEDOS OF THE MOONS AND JUPITER CONVOLVED THROUGHT THE SCT
#  FILTER PASS BANDS TO PREDICT THE OBSERVED TRANSMISSION.
###############################################################################

IoRelSlope=(Iodata['656']['TransInt']-Iodata['632']['TransInt'])/(656-632)
Io647Slope=Iodata['632']['TransInt']+IoRelSlope*15.
Io620Slope=Iodata['632']['TransInt']+IoRelSlope*(-12.)
IoNH3Abs=Iodata['647']['TransInt']/Io647Slope
IoCH4Abs=Iodata['620']['TransInt']/Io620Slope
print("************ Io Trans",IoNH3Abs)
print("************ Io Trans",IoCH4Abs)
EurRelSlope=(Europadata['656']['TransInt']-Europadata['632']['TransInt'])/(656-632)
Eur647Slope=Europadata['632']['TransInt']+EurRelSlope*15.
Eur620Slope=Europadata['632']['TransInt']+EurRelSlope*(-12.)
EurNH3Abs=Europadata['647']['TransInt']/Eur647Slope
EurCH4Abs=Europadata['620']['TransInt']/Eur620Slope
print("************ Eur Trans",EurNH3Abs)
print("************ Eur Trans",EurCH4Abs)
GanRelSlope=(Ganymededata['656']['TransInt']-Ganymededata['632']['TransInt'])/(656-632)
Gan647Slope=Ganymededata['632']['TransInt']+GanRelSlope*15.
Gan620Slope=Ganymededata['632']['TransInt']+GanRelSlope*(-12.)
GanNH3Abs=Ganymededata['647']['TransInt']/Gan647Slope
GanCH4Abs=Ganymededata['620']['TransInt']/Gan620Slope
print("************ Gan Trans",GanNH3Abs)
print("************ Gan Trans",GanCH4Abs)
CalRelSlope=(Callistodata['656']['TransInt']-Callistodata['632']['TransInt'])/(656-632)
Cal647Slope=Callistodata['632']['TransInt']+CalRelSlope*15.
Cal620Slope=Callistodata['632']['TransInt']+CalRelSlope*(-12.)
CalNH3Abs=Callistodata['647']['TransInt']/Cal647Slope
CalCH4Abs=Callistodata['620']['TransInt']/Cal620Slope
print("************ Cal Trans",CalNH3Abs)
print("************ Cal Trans",CalCH4Abs)

MoonsAvgdata['647']['TransInt'] = np.mean([IoNH3Abs,EurNH3Abs,GanNH3Abs,CalNH3Abs])
MoonsAvgdata['647']['Tau_Albedo']=-np.log(MoonsAvgdata['647']['TransInt'])
MoonsAvgdata['647']['NH3ColDens']=1000.*MoonsAvgdata['647']['Tau_Albedo']/Jupiterdata['647']['keff_NH3']
MoonsAvgdata['647']['CH4ColDens']=1000.*MoonsAvgdata['647']['Tau_Albedo']/Jupiterdata['647']['keff_CH4']

MoonsAvgdata['620']['TransInt'] = np.mean([IoCH4Abs,EurCH4Abs,GanCH4Abs,CalCH4Abs])
MoonsAvgdata['620']['Tau_Albedo']=-np.log(MoonsAvgdata['620']['TransInt'])
MoonsAvgdata['620']['NH3ColDens']=1000.*MoonsAvgdata['620']['Tau_Albedo']/Jupiterdata['620']['keff_NH3']
MoonsAvgdata['620']['CH4ColDens']=1000.*MoonsAvgdata['620']['Tau_Albedo']/Jupiterdata['620']['keff_CH4']

filtereffectivedata = open('c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/filtereffectivedataMod3.csv', 'w')
tmp="Wavelength (nm),Filter Name,k_eff (NH3),l_eff (NH3),k_eff (CH4),l_eff (CH4),Trans,Tau,NH3 (m-atm),CH4 (m-atm)\n"
filtereffectivedata.write(tmp)

for filtr in ['620','647']:
    print(filtr)
    tmp=filtr+","+Jupiterdata[filtr]['filtname']+","+str(Jupiterdata[filtr]['keff_NH3'])+","\
            +str(Jupiterdata[filtr]['leff_NH3'])+","+str(Jupiterdata[filtr]['keff_CH4'])+","\
            +str(Jupiterdata[filtr]['leff_CH4'])+","+str(MoonsAvgdata[filtr]['TransInt'])+","\
            +str(MoonsAvgdata[filtr]['Tau_Albedo'])+","+str(MoonsAvgdata[filtr]['NH3ColDens'])+","\
            +str(MoonsAvgdata[filtr]['CH4ColDens'])+"\n"
    filtereffectivedata.write(tmp)
filtereffectivedata.close()
"""
wv=np.array([610.,620.,632.,647.,656,670.])
line=MUSERelSlope*(wv-632.)+MUSEdata['632']['Absr_Int']
axsFilt.plot(wv,line/20.,label='fit') # divide by 20 for 20 spectral bins in the integral
axsFilt.legend()
"""

###############################################################################
# MODEL FOR THE VLT-MUSE DATA:
#   THE KARKOSCHKA MODEL IS CONVOLVED WITH THE SQUARE SYNTHETIC FILTER BANDS
#   USED WITH THE RADIANCE DATA.
###############################################################################
Tele='VLT'
y0,y1,ytks=0.0,0.7,8
MUSEfilters,axsMUSEFilt= \
    NFL.compute_filter_spectrum(x0,x1,xtks,y0,y1,ytks,FilterList,
                                             Albedo,Continuum_Albedo,Model,
                                             Telescope=Tele)

MUSEdata = deepcopy(MUSEfilters)
for filtr in FilterList[0:4]:
    MUSEdata[filtr]['Cont_Int'],MUSEdata[filtr]['Absr_Int'],MUSEdata[filtr]['TransInt']= \
        NFL.cont_absorption_calcs(MUSEdata[filtr]['ContProd'],MUSEdata[filtr]['AbsrProd'], \
                                  float(filtr)-MUSEdata[filtr]['halfwdth'],\
                                  float(filtr)+MUSEdata[filtr]['halfwdth'], \
                                      MUSEdata[filtr]['filtname'])
    Jupiterdata[filtr]['Tau_Albedo']=-np.log(MUSEdata[filtr]['TransInt'])

    zeros=np.zeros(Jupiterdata[filtr]['FiltTrans'].shape[0])
    
# COMPUTE MUSE TRANSMISSIONS *WITHOUT* ACCOUNTING FOR 6NM WIDTH OF SPLIT HA
#  FILTER       
MUSERelSlope=(MUSEdata['656']['TransInt']-MUSEdata['632']['TransInt'])/(656-632)
MUSE647Slope=MUSEdata['632']['TransInt']+MUSERelSlope*15.
MUSE620Slope=MUSEdata['632']['TransInt']+MUSERelSlope*(-12.)
MUSENH3Abs=MUSEdata['647']['TransInt']/MUSE647Slope
MUSECH4Abs=MUSEdata['620']['TransInt']/MUSE620Slope

print("MUSE *NOT* CONSIDERING 6 NM HA SPLIT FILTER")
print("************ MUSE NH3 Trans",MUSENH3Abs)
print("************ MUSE CH4 Trans",MUSECH4Abs)

# COMPUTE MUSE TRANSMISSIONS *WITH* ACCOUNTING FOR 6NM WIDTH OF SPLIT HA
#  FILTER    
fout_sfx='-'+Tele+'-'+Model
MUSEdata=NFL.compute_vertical_transmission_profiles(MUSEdata,MUSEfilters,CH4,NH3,P,
                                                 fout_sfx=fout_sfx)
slopewt656=np.sum(MUSEdata['656']['FiltTrans'][:,1])/np.sum(MUSEdata['632']['FiltTrans'][:,1])
print("slopewt656=",slopewt656)
MUSERelSlope=(MUSEdata['656']['Absr_Int']/slopewt656-MUSEdata['632']['Absr_Int'])/(656-632)
MUSE647Slope=MUSEdata['632']['Absr_Int']+MUSERelSlope*15.
MUSE620Slope=MUSEdata['632']['Absr_Int']+MUSERelSlope*(-12.)
MUSENH3Abs=MUSEdata['647']['Absr_Int']/MUSE647Slope
MUSECH4Abs=MUSEdata['620']['Absr_Int']/MUSE620Slope

wv=np.array([610.,620.,632.,647.,656,670.])
line=MUSERelSlope*(wv-632.)+MUSEdata['632']['Absr_Int']
axsMUSEFilt.plot(wv,line/20.,label='fit') # divide by 20 for 20 spectral bins in the integral
axsMUSEFilt.legend()

print("")
print("MUSE CONSIDERING 6 NM HA SPLIT FILTER")
print("************ MUSE NH3 Trans",MUSENH3Abs)
print("************ MUSE CH4 Trans",MUSECH4Abs)