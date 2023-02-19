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
#from scipy import interpolate
import GeneralSpecUtils_P3 as GSU
#from numpy import genfromtxt
import NH3_Filter_Library_P3 as NFL

###############################################################################
# LOAD JOVIAN DISK-INTEGRATEDALBEDO DATA FROM KARKOSCHKA, 1994 (DATA FROM 1993)
###############################################################################
#Plot Layout Configuration
ContinuumModel='1'
x0,x1,xtks=600.,1000.,9
y0,y1,ytks=0.0,0.7,8
Albedo,Continuum_Albedo,CH4,NH3=NFL.Get_Albedo_and_Absorption(x0,x1,xtks,y0,y1,ytks,
                                                              ContMod=ContinuumModel,Crossect=True)
########## END OF FIRST FUNCTION AND PLOT ##############

###############################################################################
# RETRIEVE FILTER TRANSMISSIONS FROM MASTER FILTER LIBRARY
#   !!!! NEED TO ADD TELESCOPE THROUGHPUT CALIBRATION HERE
###############################################################################
filterwavelength=['620','632','647','656','658','672','730','889','940']
filterdata={'620':{'transfile':'620CH4/620CH4_Transmission.txt',
                   'filtname':'620CH4','filtwdth':10.},
             '632':{'transfile':'632OI/632OI_Transmission.txt',
                    'filtname':'632OI','filtwdth':10.},
             '647':{'transfile':'647CNT/647CNT_Transmission.txt',
                    'filtname':'647NH3','filtwdth':10.},
             '656':{'transfile':'656HIA/656HIA_Transmission.txt',
                    'filtname':'656HIA','filtwdth':10.},
             '658':{'transfile':'658NII/658NII_Transmission.txt',
                    'filtname':'658NII','filtwdth':5.},
             '672':{'transfile':'672SII/672SII_Transmission.txt',
                    'filtname':'672SII','filtwdth':10.},
             '730':{'transfile':'730OII/730OII_Transmission.txt',
                    'filtname':'730OII','filtwdth':10.},
             '889':{'transfile':'889CH4/889CH4_Transmission.txt',
                    'filtname':'889CH4','filtwdth':10.},
             '940':{'transfile':'940NIR/940NIR_Transmission.txt',
                    'filtname':'940NIR','filtwdth':10.}}

Jupiterdata = deepcopy(filterdata)
Iodata = deepcopy(filterdata)
Europadata = deepcopy(filterdata)
Ganymededata = deepcopy(filterdata)
Callistodata = deepcopy(filterdata)

path='c:/Astronomy/Projects/Techniques/InstrumentPerformance-P3/Filters/'

###############################################################################
# PLOT FILTER TRANSMISSIONS CONVOLVED WITH DISK-INTEGRATED ALBEDO AND CONTINUUM
###############################################################################
fig1,axs1=pl.subplots(1,1,figsize=(6.0,3.5), dpi=150, facecolor="white",
                      sharex=True)
axs1.set_xlim(x0,x1)
axs1.set_xticks(np.linspace(x0,x1,xtks, endpoint=True))
axs1.set_ylim(y0,y1)
axs1.set_yticks(np.linspace(y0,y1,ytks, endpoint=True))
axs1.grid(linewidth=0.2)
axs1.tick_params(axis='both', which='major', labelsize=8)
axs1.set_ylabel("Albedo x Transmission",color="black")
axs1.set_xlabel("Wavelength (nm)")

for filtr in filterwavelength:
    Jupiterdata[filtr]['FiltTrans']=np.loadtxt(path+Jupiterdata[filtr]['transfile'],usecols=range(2))
    Jupiterdata[filtr]['ContProd']=GSU.SpectrumMath(Jupiterdata[filtr]['FiltTrans'],Continuum_Albedo,"Multiply")
    Jupiterdata[filtr]['AbsrProd']=GSU.SpectrumMath(Jupiterdata[filtr]['FiltTrans'],Albedo,"Multiply")

    Jupiterdata[filtr]['Cont_Int'],Jupiterdata[filtr]['Absr_Int'],Jupiterdata[filtr]['TransInt']= \
        NFL.cont_absorption_calcs(Jupiterdata[filtr]['ContProd'],Jupiterdata[filtr]['AbsrProd'], \
                                  float(filtr)-Jupiterdata[filtr]['filtwdth'],\
                                  float(filtr)+Jupiterdata[filtr]['filtwdth'], \
                                      Jupiterdata[filtr]['filtname'])
    Jupiterdata[filtr]['Tau_Albedo']=-np.log(Jupiterdata[filtr]['TransInt'])

    zeros=np.zeros(Jupiterdata[filtr]['FiltTrans'].shape[0])

    if str(filtr)[0:3]=='620':
        axs1.plot(Jupiterdata[filtr]['ContProd'][:,0],Jupiterdata[filtr]['ContProd'][:,1],label='Continuum Albedo',linewidth=1,color='C0')
        axs1.fill_between(Jupiterdata[filtr]['AbsrProd'][:,0],zeros, Jupiterdata[filtr]['AbsrProd'][:,1],label='Jupiter Albedo',color='C0',alpha=0.2)
    else:
        axs1.plot(Jupiterdata[filtr]['ContProd'][:,0],Jupiterdata[filtr]['ContProd'][:,1],linewidth=1,color='C0')
        axs1.fill_between(Jupiterdata[filtr]['AbsrProd'][:,0],zeros, Jupiterdata[filtr]['AbsrProd'][:,1],color='C0',alpha=0.2)

axs1.legend(fontsize=8,loc=2)

axs1.set_title("Convolution with Filters (Continuum Model "+ContinuumModel+")")

fig1.subplots_adjust(left=0.10, right=0.90, top=0.9, bottom=0.14)

fig1.savefig('c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/AmmoniaFilter.png',dpi=320)

########## END OF SECOND FUNCTION AND PLOT ##############

###############################################################################
# BEGIN K_eff, l_eff, AND WEIGHTING FUNCTION CALCULATIONS
###############################################################################
P=np.geomspace(10.,1.0e7,num=70,endpoint=True,dtype=float) #in Pascals

fig_trans,axs_trans=pl.subplots(2,2,figsize=(6.0,6.0), dpi=150, facecolor="white",
                              sharex=True,sharey=True)
for i in range(0,2):
    for j in range(0,2):
        axs_trans[i,j].set_ylim(100.,5000000.)
        axs_trans[i,j].set_yscale('log')
        axs_trans[i,j].set_ylim(axs_trans[i,j].get_ylim()[::-1]) #reverse y-axis
        axs_trans[i,j].set_xlim(0.,1.) #reverse y-axis
        axs_trans[i,j].set_xscale('linear')
        axs_trans[i,j].grid(which='both')
        if i==1:
            axs_trans[i,j].set_xlabel("Transmission")
        if j==0:
            axs_trans[i,j].set_ylabel("Pressue (Pa)")
axs_trans[0,0].set_title("Gas+Rayleigh")
axs_trans[0,1].set_title("CH4")
axs_trans[1,0].set_title("Rayleigh")
axs_trans[1,1].set_title("NH3")

fig_Keff,axs_Keff=pl.subplots(2,2,figsize=(6.0,6.0), dpi=150, facecolor="white",
                              sharex=True,sharey=True)
for i in range(0,2):
    for j in range(0,2):

        axs_Keff[i,j].set_ylim(1000.,1000000.)
        axs_Keff[i,j].set_yscale('log')
        axs_Keff[i,j].set_ylim(axs_Keff[i,j].get_ylim()[::-1]) #reverse y-axis
        axs_Keff[i,j].set_xlim(0.,1.) #reverse y-axis
        axs_Keff[i,j].set_xscale('linear')
        axs_Keff[i,j].grid(which='both')
        if i==1:
            axs_Keff[i,j].set_xlabel("Normalized Weight")
        if j==0:
            axs_Keff[i,j].set_ylabel("Pressue (Pa)")

axs_Keff[0,0].set_title("Gas+Rayleigh")
axs_Keff[0,1].set_title("CH4")
axs_Keff[1,0].set_title("Rayleigh")
axs_Keff[1,1].set_title("NH3")

###############################################################################
# Write data for each filter to a csv file
###############################################################################
filtereffectivedata = open('c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/filtereffectivedata.csv', 'w')
tmp="Wavelength (nm),Filter Name,k_eff (NH3),l_eff (NH3),k_eff (CH4),l_eff (CH4),Trans,Tau,NH3 (m-atm),CH4 (m-atm)\n"
filtereffectivedata.write(tmp)
for filtr in filterwavelength:
    Jupiterdata[filtr]['keff_CH4'],Jupiterdata[filtr]['leff_CH4']=NFL.K_eff(P,Jupiterdata[filtr]['FiltTrans'],CH4,\
                  float(filtr)-Jupiterdata[filtr]['filtwdth'],float(filtr)+Jupiterdata[filtr]['filtwdth'],Jupiterdata[filtr]['filtname'],axs_Keff)
    Jupiterdata[filtr]['keff_NH3'],Jupiterdata[filtr]['leff_NH3']=NFL.K_eff(P,Jupiterdata[filtr]['FiltTrans'],NH3,\
                  float(filtr)-Jupiterdata[filtr]['filtwdth'],float(filtr)+Jupiterdata[filtr]['filtwdth'],Jupiterdata[filtr]['filtname'],axs_Keff)
    Jupiterdata[filtr]['tau_CH4']=NFL.tau_gas_versus_P(P,Jupiterdata[filtr]['keff_CH4'],Jupiterdata[filtr]['filtname'],axs_Keff,gas='CH4')
    Jupiterdata[filtr]['tau_NH3']=NFL.tau_gas_versus_P(P,Jupiterdata[filtr]['keff_NH3'],Jupiterdata[filtr]['filtname'],axs_Keff,gas='NH3')
    Jupiterdata[filtr]['tau_R']=NFL.tau_rayleigh_versus_P(P,Jupiterdata[filtr]['leff_CH4'],Jupiterdata[filtr]['filtname'],axs_Keff)
    tau_gas=Jupiterdata[filtr]['tau_CH4']+Jupiterdata[filtr]['tau_NH3']
    tmp=NFL.Compute_Transmission(P,Jupiterdata[filtr]['tau_R'],tau_gas,Jupiterdata[filtr]['filtname'],axs_trans[0,0],axs_Keff[0,0])
    tmp=NFL.Compute_Transmission(P,Jupiterdata[filtr]['tau_R']*0.0,Jupiterdata[filtr]['tau_CH4'],Jupiterdata[filtr]['filtname']+' CH4',axs_trans[0,1],axs_Keff[0,1])
    tmp=NFL.Compute_Transmission(P,Jupiterdata[filtr]['tau_R'],tau_gas*0.0,Jupiterdata[filtr]['filtname']+' Ray',axs_trans[1,0],axs_Keff[1,0])
    tmp=NFL.Compute_Transmission(P,Jupiterdata[filtr]['tau_R']*0.0,Jupiterdata[filtr]['tau_NH3'],Jupiterdata[filtr]['filtname']+' NH3',axs_trans[1,1],axs_Keff[1,1])
    #print(filtr,Jupiterdata[filtr]['filtname'],Jupiterdata[filtr]['keff_NH3'],Jupiterdata[filtr]['leff_NH3'],\
    #      Jupiterdata[filtr]['keff_CH4'],Jupiterdata[filtr]['leff_CH4'])
    Jupiterdata[filtr]['NH3ColDens']=1000.*Jupiterdata[filtr]['Tau_Albedo']/Jupiterdata[filtr]['keff_NH3']
    Jupiterdata[filtr]['CH4ColDens']=1000.*Jupiterdata[filtr]['Tau_Albedo']/Jupiterdata[filtr]['keff_CH4']

    tmp=filtr+","+Jupiterdata[filtr]['filtname']+","+str(Jupiterdata[filtr]['keff_NH3'])+","\
            +str(Jupiterdata[filtr]['leff_NH3'])+","+str(Jupiterdata[filtr]['keff_CH4'])+","\
            +str(Jupiterdata[filtr]['leff_CH4'])+","+str(Jupiterdata[filtr]['TransInt'])+","\
            +str(Jupiterdata[filtr]['Tau_Albedo'])+","+str(Jupiterdata[filtr]['NH3ColDens'])+","\
            +str(Jupiterdata[filtr]['CH4ColDens'])+"\n"
    filtereffectivedata.write(tmp)

filtereffectivedata.close()

axs_trans[0,0].legend(loc=1,ncol=3, borderaxespad=0.,prop={'size':6})
fig_trans.subplots_adjust(left=0.12, right=0.96, top=0.94, bottom=0.09)
        
axs_Keff[0,0].legend(loc=1,ncol=3, borderaxespad=0.,prop={'size':6})
fig_Keff.subplots_adjust(left=0.12, right=0.96, top=0.94, bottom=0.09)


fig_trans.savefig('c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/TransmissionFunctions.png',dpi=320,bbox_inches = 'tight')
fig_Keff.savefig('c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/ContributionFunctions.png',dpi=320,bbox_inches = 'tight')

###############################################################################
# BEGIN CALCULATIONS FOR **RAYLEIGH CANCELING** AND GAS ABSORPTION ONLY
#   WEIGHTING FUNCTIONS
###############################################################################
fig_ray,axs_ray=pl.subplots(2,2,figsize=(6.0,6.0), dpi=150, facecolor="white",
                      sharex=True,sharey=True)
for i in range(0,2):
    for j in range(0,2):
        axs_ray[i,j].set_title("Rayleigh Compare (2-way)")
        axs_ray[i,j].set_ylim(100.,5000000.)
        axs_ray[i,j].set_yscale('log')
        axs_ray[i,j].set_ylim(axs_ray[i,j].get_ylim()[::-1]) #reverse y-axis
        axs_ray[i,j].set_xscale('linear')
        axs_ray[i,j].grid(which='both')
        if i==1:
            axs_ray[i,j].set_xlabel("Transmission")
        if j==0:
            axs_ray[i,j].set_ylabel("Pressue (Pa)")
            axs_ray[i,j].set_xlim(0.,1.) #reverse y-axis
        elif j==1:
            axs_ray[i,j].set_xlim(0.,1.) #reverse y-axis


axs_ray[0,0].set_title("Gas+Rayleigh")
axs_ray[0,1].set_title("CH4")
axs_ray[1,0].set_title("Rayleigh")
axs_ray[1,1].set_title("NH3")

fig_raycont,axs_raycont=pl.subplots(2,2,figsize=(6.0,6.0), dpi=150, facecolor="white",
                      sharex=True)
for i in range(0,2):
    for j in range(0,2):
        axs_raycont[i,j].set_title("Rayleigh Compare (2-way)")
        axs_raycont[i,j].set_ylim(1000.,50000000.)
        axs_raycont[i,j].set_yscale('log')
        axs_raycont[i,j].set_ylim(axs_raycont[i,j].get_ylim()[::-1]) #reverse y-axis
        axs_raycont[i,j].set_xlim(0.,1.) #reverse y-axis
        axs_raycont[i,j].set_xscale('linear')
        axs_raycont[i,j].grid(which='both')
        if i==1:
            axs_raycont[i,j].set_xlabel("Weighting Functions (2-way)")
        if j==0:
            axs_raycont[i,j].set_ylabel("Pressue (Pa)")
            axs_raycont[i,j].set_xlim(0.,1.) #reverse y-axis
        elif j==1:
            axs_raycont[i,j].set_xlim(0.,1.) #reverse y-axis

axs_raycont[0,0].set_title("Gas+Rayleigh")
axs_raycont[0,1].set_title("CH4")
axs_raycont[1,0].set_title("Rayleigh")
axs_raycont[1,1].set_title("NH3")

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
fig_ray.subplots_adjust(left=0.12, right=0.96, top=0.94, bottom=0.09)

axs_raycont[0,0].legend(loc=1,ncol=2, borderaxespad=0.,prop={'size':6})
axs_raycont[0,1].legend(loc=1,ncol=2, borderaxespad=0.,prop={'size':6})
axs_raycont[1,0].legend(loc=1,ncol=2, borderaxespad=0.,prop={'size':6})
axs_raycont[1,1].legend(loc=1,ncol=2, borderaxespad=0.,prop={'size':6})
fig_raycont.subplots_adjust(left=0.12, right=0.96, top=0.94, bottom=0.09)

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
x0,x1,xtks=600.,700.,11
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
for filtr in filterwavelength:
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

filterloop=['620','632','647','656']
for filt in filterloop:
    print('^^^^^^^^^^^^^'+filt)

    Iodata[filt]['Cont_Int'],Iodata[filt]['Absr_Int'],Iodata[filt]['TransInt']= \
        NFL.cont_absorption_calcs(Iodata[filt]['ContProd'],Jupiterdata[filt]['AbsrProd'], \
                                  float(filt)-Jupiterdata[filt]['filtwdth'],\
                                  float(filt)+Jupiterdata[filt]['filtwdth'], \
                                      Jupiterdata[filt]['filtname'],prn=False)
    Europadata[filt]['Cont_Int'],Europadata[filt]['Absr_Int'],Europadata[filt]['TransInt']= \
        NFL.cont_absorption_calcs(Europadata[filt]['ContProd'],Jupiterdata[filt]['AbsrProd'], \
                                  float(filt)-Jupiterdata[filt]['filtwdth'],\
                                  float(filt)+Jupiterdata[filt]['filtwdth'], \
                                      Jupiterdata[filt]['filtname'],prn=False)
    Ganymededata[filt]['Cont_Int'],Ganymededata[filt]['Absr_Int'],Ganymededata[filt]['TransInt']= \
        NFL.cont_absorption_calcs(Ganymededata[filt]['ContProd'],Jupiterdata[filt]['AbsrProd'], \
                                  float(filt)-Jupiterdata[filt]['filtwdth'],\
                                  float(filt)+Jupiterdata[filt]['filtwdth'], \
                                      Jupiterdata[filt]['filtname'],prn=False)
    Callistodata[filt]['Cont_Int'],Callistodata[filt]['Absr_Int'],Callistodata[filt]['TransInt']= \
        NFL.cont_absorption_calcs(Callistodata[filt]['ContProd'],Jupiterdata[filt]['AbsrProd'], \
                                  float(filt)-Jupiterdata[filt]['filtwdth'],\
                                  float(filt)+Jupiterdata[filt]['filtwdth'], \
                                      Jupiterdata[filt]['filtname'],prn=False)

    print(Iodata[filt]['Cont_Int'],Iodata[filt]['Absr_Int'],Iodata[filt]['TransInt'])
    print(Europadata[filt]['Cont_Int'],Europadata[filt]['Absr_Int'],Europadata[filt]['TransInt'])
    print(Ganymededata[filt]['Cont_Int'],Ganymededata[filt]['Absr_Int'],Ganymededata[filt]['TransInt'])
    print(Callistodata[filt]['Cont_Int'],Callistodata[filt]['Absr_Int'],Callistodata[filt]['TransInt'])

########## END OF FOURTH FUNCTION AND PLOT ##############

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

