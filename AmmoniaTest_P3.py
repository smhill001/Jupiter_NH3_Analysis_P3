# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 12:38:27 2021

This code creates two plots, each with two subplots
    Plot 1:     Presents information on Jupiter reflectance and filters
        Subplot 1:  Karkoschka, 1994 disk-integrated albedo plus a spline fit that removes
                    the NH3 absorption
        Subplot 2:  Filter transmissions computed (647, 656, 658, 672) and convovled
                    with reference albedo (and spline fit for 647)

    Plot 2:     Presents information on potential color biases due to differences in 
                regional/feature albedos on Jupiter and Galilean moon colors
        Subplot 1:  Plots regional reflectivities (I/F) for NEB, EZ, and SEB from Dahl, 2021
        Subplot 2:  Plots Galilean moon reflectivites from Clark & McCord, 1980

UPDATE 2022-01-25:  Converted the code to Python 3 on the new Astronomy laptop
UPDATE 2022-06-19:  

@author: Steven Hill
"""
import sys
sys.path.append('c:/Astronomy/Python Play')
sys.path.append('c:/Astronomy/Python Play/Util_P3')
sys.path.append('c:/Astronomy/Python Play/SPLibraries_P3')
sys.path.append('c:/Astronomy/Python Play/SpectroPhotometry/Spectroscopy_P3')
import matplotlib.pyplot as pl
import numpy as np
from scipy import interpolate
import GeneralSpecUtils_P3 as GSU
from numpy import genfromtxt
import NH3_Filter_Library_P3 as NFL

###### Get reference disk-integrated albedo from Karkoschka, 1994

Jupiter_Karkoschka1993 = np.fromfile(file="c:/Astronomy/Projects/Planets/Saturn/Spectral Data/Karkoschka/1993.tab.txt", dtype=float, count=-1, sep=" ")
kark1993nrows=int(Jupiter_Karkoschka1993.size/8)
Jupiter_Karkoschka1993=np.reshape(Jupiter_Karkoschka1993,[kark1993nrows,8])

Jupiter_KarkRef1993=np.zeros((kark1993nrows,2))
Jupiter_KarkRef1993[:,0]=Jupiter_Karkoschka1993[:,0]
Jupiter_KarkRef1993[:,1]=Jupiter_Karkoschka1993[:,3] #Albedo

WaveGrid,SignalonGrid=GSU.uniform_wave_grid(Jupiter_KarkRef1993[:,0],Jupiter_KarkRef1993[:,1],
                                        Extend=False,Fine=False)

JK=np.zeros((WaveGrid.size,2))
JK[:,0]=WaveGrid
JK[:,1]=SignalonGrid
#### CODE THIS INLINE WITH THE NEW SCIPY PYTHON 3 IMPLEMENTATION OF
#### SPLINE FITTING
SplineWV1= np.array([560.0, 580.0, 600.0, 635.0, 660.0, 675.0, 690., 714.0,
                     745.0, 830.0, 945.0])
SplineWV2 = np.array([560.0, 580.0, 600.0, 677.0, 690., 
                     745.0, 830.0, 945.0,1050.0])
SplineWV=SplineWV2
SplineMag=np.ones(SplineWV.size)
for i in range(0,SplineWV.size):
    Start=SplineWV[i]-.0000001
    End=SplineWV[i]+.0000001
    SplineWVIndices=np.where((JK[:,0] >Start) & \
         (JK[:,0] < End))
    print("i= ",i,SplineWVIndices)
    SplineMag[i]=np.log10(JK[SplineWVIndices[0],1])

print(SplineMag)
SplineMag[SplineMag.size-1]=np.log10(0.42)
x = np.arange(0, 2*np.pi+np.pi/4, 2*np.pi/8)
y = np.sin(x)
tck = interpolate.splrep(SplineWV, SplineMag, s=0)
#xnew = np.arange(0, 2*np.pi, np.pi/50)
Temp = 10**interpolate.splev(WaveGrid, tck, der=0)



#print(SplineWV)
#Temp=SPL.log_spline_smoother(SplineWV,JK)

Continuum_Albedo=np.zeros((WaveGrid.size,2))
Continuum_Albedo[:,0]=WaveGrid
Continuum_Albedo[:,1]=Temp

CH4_KarkRef1993=np.zeros((kark1993nrows,2))
CH4_KarkRef1993[:,0]=Jupiter_Karkoschka1993[:,0]
CH4_KarkRef1993[:,1]=Jupiter_Karkoschka1993[:,2] #CH4 Coef
print("***",CH4_KarkRef1993.shape)
WaveGrid,SignalonGrid=GSU.uniform_wave_grid(CH4_KarkRef1993[:,0],CH4_KarkRef1993[:,1],
                                        Extend=False,Fine=False)
CH4=np.zeros((WaveGrid.size,2))
CH4[:,0]=WaveGrid
CH4[:,1]=SignalonGrid

fn='Lutz&Owen1980_Figure5_AmmoniaCrossSection.csv'
pth="c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/"
NH3_Lutz_Owen_1980 = np.array(genfromtxt(pth+fn, delimiter=','))


fn='Exomol_NH3.csv'
pth="c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/"
NH3_Exomol = np.array(genfromtxt(pth+fn, delimiter=','))
print("***",NH3_Exomol.shape)
wv,csec=GSU.uniform_wave_grid(NH3_Exomol[:,0]*1000.,NH3_Exomol[:,1],Extend=False,Fine=False)
NH3_Exomol_regrid=np.transpose(np.array([wv,csec]))

NH3=NH3_Exomol_regrid
print("***",NH3_Exomol_regrid.shape,NH3_Exomol_regrid.max(1))

###### Get reference regional I/F reflectivities from Dahl, 2021?

Dahl_NEB="c:/Astronomy/Projects/SAS 2021 Ammonia/Dahl Spectra-NEB.txt"
NEB = genfromtxt(Dahl_NEB, delimiter=',')
Dahl_SEB="c:/Astronomy/Projects/SAS 2021 Ammonia/Dahl Spectra-SEB.txt"
SEB = genfromtxt(Dahl_SEB, delimiter=',')
Dahl_EZ="c:/Astronomy/Projects/SAS 2021 Ammonia/Dahl Spectra-EZ.txt"
EZ = genfromtxt(Dahl_EZ, delimiter=',')

###### Plot disk-integrated reference albedo and simulated ammonia-free albedo
fig1,axs1=pl.subplots(2,1,figsize=(6.0,6.0), dpi=150, facecolor="white",
                      sharex=True)
#fig1.suptitle("Ammonia Filters",x=0.5,ha='center',color='k')

#Plot Layout Configuration
x0,x1,xtks=600.,1000.,9
y0,y1,ytks=0.0,0.7,8

# Set x limits
axs1[0].set_xlim(x0,x1)
# Set x ticks
axs1[0].set_xticks(np.linspace(x0,x1,xtks, endpoint=True))
# Set y limits
axs1[0].set_ylim(y0,y1)
axs1[0].set_yticks(np.linspace(y0,y1,ytks, endpoint=True))
# Set y ticks
axs1[0].grid(linewidth=0.2)
axs1[0].tick_params(axis='both', which='major', labelsize=8)
axs1[0].set_ylabel("Albedo",color="black")

axs1[0].plot(Continuum_Albedo[:,0],Continuum_Albedo[:,1],label='Continuum Albedo',
             linewidth=1,linestyle='--',color='C0')
axs1[0].plot(JK[:,0],JK[:,1],label='Jupiter Albedo',linewidth=1.0,color='C0')
axs1[0].set_title("Jovian Albedo and Methane and Ammonia Absorption")

axs1b = axs1[0].twinx()  # instantiate a second axes that shares the same x-axis
axs1b.ticklabel_format(axis='y')
axs1b.tick_params(axis='y', which='major', labelsize=8)
axs1b.set_yscale('log')
axs1b.set_ylim(1e-4,1e3)
axs1b.set_ylabel("Absorption Coefficient (km-atm)")#,color="green")


axs1b.plot(CH4_KarkRef1993[:,0],CH4_KarkRef1993[:,1],label='CH4 Abs. Coef. ',linewidth=0.5,color='C2')
axs1b.plot(NH3[:,0],NH3[:,1],label='NH3 Abs. Coef. ',linewidth=1.0,color='C3')
axs1b.plot(NH3_Lutz_Owen_1980[:,0],NH3_Lutz_Owen_1980[:,1],label='NH3 Abs. Coef. L&O ',linewidth=0.5,color='C3')

#axs1b.plot(CH4[:,0],CH4[:,1],label='CH4 Abs. Coef. ',linewidth=1,color='k')
axs1[0].legend(fontsize=8, loc=2)
axs1b.legend(fontsize=8, loc=1)

path='c:/Astronomy/Projects/Techniques/InstrumentPerformance-P3/Filters/'


###### Retrieve filter transmissions and convovle with disk integrated albedoes


Transmission620 = np.loadtxt(path+'620CH4/620CH4_Transmission.txt',usecols=range(2))
Transmission632 = np.loadtxt(path+'632OI/632OI_Transmission.txt',usecols=range(2))
Transmission647 = np.loadtxt(path+'647CNT/647CNT_Transmission.txt',usecols=range(2))
Transmission656 = np.loadtxt(path+'656HIA/656HIA_Transmission.txt',usecols=range(2))
Transmission658 = np.loadtxt(path+'658NII/658NII_Transmission.txt',usecols=range(2))
Transmission672 = np.loadtxt(path+'672SII/672SII_Transmission.txt',usecols=range(2))
Transmission730 = np.loadtxt(path+'730OII/730OII_Transmission.txt',usecols=range(2))
Transmission889 = np.loadtxt(path+'889CH4/889CH4_Transmission.txt',usecols=range(2))
Transmission940 = np.loadtxt(path+'940NIR/940NIR_Transmission.txt',usecols=range(2))


##########

ContinuumProduct940=GSU.SpectrumMath(Transmission940,Continuum_Albedo,"Multiply")
ContinuumProduct889=GSU.SpectrumMath(Transmission889,Continuum_Albedo,"Multiply")
ContinuumProduct730=GSU.SpectrumMath(Transmission730,Continuum_Albedo,"Multiply")
ContinuumProduct672=GSU.SpectrumMath(Transmission672,Continuum_Albedo,"Multiply")
ContinuumProduct658=GSU.SpectrumMath(Transmission658,Continuum_Albedo,"Multiply")
ContinuumProduct656=GSU.SpectrumMath(Transmission656,Continuum_Albedo,"Multiply")
ContinuumProduct647=GSU.SpectrumMath(Transmission647,Continuum_Albedo,"Multiply")
ContinuumProduct632=GSU.SpectrumMath(Transmission632,Continuum_Albedo,"Multiply")
ContinuumProduct620=GSU.SpectrumMath(Transmission620,Continuum_Albedo,"Multiply")

AbsorptionProduct940=GSU.SpectrumMath(Transmission940,JK,"Multiply")
AbsorptionProduct889=GSU.SpectrumMath(Transmission889,JK,"Multiply")
AbsorptionProduct730=GSU.SpectrumMath(Transmission730,JK,"Multiply")
AbsorptionProduct672=GSU.SpectrumMath(Transmission672,JK,"Multiply")
AbsorptionProduct658=GSU.SpectrumMath(Transmission658,JK,"Multiply")
AbsorptionProduct656=GSU.SpectrumMath(Transmission656,JK,"Multiply")
AbsorptionProduct647=GSU.SpectrumMath(Transmission647,JK,"Multiply")
AbsorptionProduct632=GSU.SpectrumMath(Transmission632,JK,"Multiply")
AbsorptionProduct620=GSU.SpectrumMath(Transmission620,JK,"Multiply")

###### Plot filter transmissions convolved with disk-integrated albedos

# Set x limits
axs1[1].set_xlim(x0,x1)
# Set x ticks
axs1[1].set_xticks(np.linspace(x0,x1,xtks, endpoint=True))
# Set y limits
axs1[1].set_ylim(y0,y1)
axs1[1].set_yticks(np.linspace(y0,y1,ytks, endpoint=True))
# Set y ticks
axs1[1].grid(linewidth=0.2)
axs1[1].tick_params(axis='both', which='major', labelsize=8)
axs1[1].set_ylabel("Albedo x Transmission",color="black")
axs1[1].set_xlabel("Wavelength (nm)")
axs1[1].plot(ContinuumProduct940[:,0],ContinuumProduct940[:,1],label='Continuum Albedo',linewidth=1,color='b')
axs1[1].plot(ContinuumProduct889[:,0],ContinuumProduct889[:,1],linewidth=1,color='b')
axs1[1].plot(ContinuumProduct730[:,0],ContinuumProduct730[:,1],linewidth=1,color='b')
axs1[1].plot(ContinuumProduct672[:,0],ContinuumProduct672[:,1],linewidth=1,color='b')
axs1[1].plot(ContinuumProduct658[:,0],ContinuumProduct658[:,1],linewidth=1,color='b')
axs1[1].plot(ContinuumProduct656[:,0],ContinuumProduct656[:,1],linewidth=1,color='b')
axs1[1].plot(ContinuumProduct647[:,0],ContinuumProduct647[:,1],linewidth=1,color='b')
axs1[1].plot(ContinuumProduct632[:,0],ContinuumProduct632[:,1],linewidth=1,color='b')
axs1[1].plot(ContinuumProduct620[:,0],ContinuumProduct620[:,1],linewidth=1,color='b')

axs1[1].plot(AbsorptionProduct940[:,0],AbsorptionProduct940[:,1],label='Jupiter Albedo',linewidth=0.5,color='r')
axs1[1].plot(AbsorptionProduct889[:,0],AbsorptionProduct889[:,1],linewidth=0.5,color='r')
axs1[1].plot(AbsorptionProduct730[:,0],AbsorptionProduct730[:,1],linewidth=0.5,color='r')
axs1[1].plot(AbsorptionProduct672[:,0],AbsorptionProduct672[:,1],linewidth=0.5,color='r')
axs1[1].plot(AbsorptionProduct658[:,0],AbsorptionProduct658[:,1],linewidth=0.5,color='r')
axs1[1].plot(AbsorptionProduct656[:,0],AbsorptionProduct656[:,1],linewidth=0.5,color='r')
axs1[1].plot(AbsorptionProduct647[:,0],AbsorptionProduct647[:,1],linewidth=0.5,color='r')
axs1[1].plot(AbsorptionProduct632[:,0],AbsorptionProduct632[:,1],linewidth=0.5,color='r')
axs1[1].plot(AbsorptionProduct620[:,0],AbsorptionProduct620[:,1],linewidth=0.5,color='r')

axs1[1].legend(fontsize=8,loc=2)
axs1[1].set_title("Convolution with Filters")

fig1.subplots_adjust(left=0.10, right=0.90, top=0.94, bottom=0.09)

fig1.savefig('c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/AmmoniaFilter.png',dpi=320)

###############################################################################
# Begin K_eff and weighting function calculations
fig_trans,axs_trans=pl.subplots(1,1,figsize=(6.0,6.0), dpi=150, facecolor="white",
                              sharex=True)
axs_trans.set_title("Transmission (2-way)")
axs_trans.set_ylim(100.,5000000.)
axs_trans.set_yscale('log')
axs_trans.set_ylim(axs_trans.get_ylim()[::-1]) #reverse y-axis
axs_trans.set_xlim(0.,1.) #reverse y-axis
axs_trans.set_xscale('linear')
axs_trans.grid(which='both')
axs_trans.set_xlabel("Transmission Functions (2-way)")
axs_trans.set_ylabel("Pressue (Pa)")

fig_Keff,axs_Keff=pl.subplots(1,1,figsize=(6.0,4.0), dpi=150, facecolor="white",
                              sharex=True)
axs_Keff.set_title("Weighting Functions (2-way)", fontsize=18)
P=np.geomspace(10.,1.0e7,num=70,endpoint=True,dtype=float) #in Pascals

#axs_Keff = pl.gca() 
axs_Keff.set_ylim(1000.,1000000.)
axs_Keff.set_yscale('log')
axs_Keff.set_ylim(axs_Keff.get_ylim()[::-1]) #reverse y-axis
axs_Keff.set_xlim(0.,1.) #reverse y-axis
axs_Keff.set_xscale('linear')
axs_Keff.grid(which='both')
axs_Keff.set_xlabel("Normalized Weight",fontsize=14)
axs_Keff.set_ylabel("Pressue (Pa)",fontsize=14)

"""
keff_CH4,leff_CH4=NFL.K_eff(P,Transmission940,CH4,930.,950.,"940NIR",axs_Keff)
tau_CH4=NFL.Weighting_Function(P,keff_CH4,"940NIR",axs_Keff)
tau_R=NFL.Rayleigh_Function(P,leff_CH4,"940NIR",axs_Keff)
tmp=NFL.Compute_Transmission(P,tau_R,tau_CH4,"940NIR",axs_trans,axs_Keff)
"""
keff_CH4,leff_CH4=NFL.K_eff(P,Transmission889,CH4,879.,899.,"889CH4",axs_Keff)
tau_CH4=NFL.Weighting_Function(P,keff_CH4,"889CH4",axs_Keff)
tau_R=NFL.Rayleigh_Function(P,leff_CH4,"889CH4",axs_Keff)
tmp=NFL.Compute_Transmission(P,tau_R,tau_CH4,"889CH4",axs_trans,axs_Keff)

keff_CH4,leff_CH4=NFL.K_eff(P,Transmission730,CH4,720.,740.,"730OII",axs_Keff)
tau_CH4=NFL.Weighting_Function(P,keff_CH4,"730OII",axs_Keff)
tau_R=NFL.Rayleigh_Function(P,leff_CH4,"730OII",axs_Keff)
tmp=NFL.Compute_Transmission(P,tau_R,tau_CH4,"730OII",axs_trans,axs_Keff)
"""
keff_CH4,leff_CH4=NFL.K_eff(P,Transmission672,CH4,662.,682.,"672SII",axs_Keff)
tau_CH4=NFL.Weighting_Function(P,keff_CH4,"672SII",axs_Keff)
tau_R=NFL.Rayleigh_Function(P,leff_CH4,"672SII",axs_Keff)
tmp=NFL.Compute_Transmission(P,tau_R,tau_CH4,"672SII",axs_trans,axs_Keff)

keff_CH4,leff_CH4=NFL.K_eff(P,Transmission658,CH4,653.,663.,"658NII",axs_Keff)
tau_CH4=NFL.Weighting_Function(P,keff_CH4,"658NII",axs_Keff)
tau_R=NFL.Rayleigh_Function(P,leff_CH4,"658NII",axs_Keff)
tmp=NFL.Compute_Transmission(P,tau_R,tau_CH4,"658NII",axs_trans,axs_Keff)
"""
keff_CH4,leff_CH4=NFL.K_eff(P,Transmission656,CH4,646.,666.,"656HIA",axs_Keff)
tau_CH4=NFL.Weighting_Function(P,keff_CH4,"656HIA",axs_Keff)
tau_R=NFL.Rayleigh_Function(P,leff_CH4,"656HIA",axs_Keff)
tmp=NFL.Compute_Transmission(P,tau_R,tau_CH4,"656HIA",axs_trans,axs_Keff)

"""
keff_CH4,leff_CH4=NFL.K_eff(P,Transmission647,CH4,637.,657.,"647CNT",axs_Keff)
tau_CH4=NFL.Weighting_Function(P,keff_CH4,"647CNT",axs_Keff)
tau_R=NFL.Rayleigh_Function(P,leff_CH4,"647CNT",axs_Keff)
tmp=NFL.Compute_Transmission(P,tau_R,tau_CH4,"647CNT",axs_trans,axs_Keff)
"""
keff_NH3,leff_NH3=NFL.K_eff(P,Transmission647,NH3,637.,657.,"647CNT",axs_Keff)
tau_NH3=NFL.Weighting_Function(P,keff_NH3,"647CNT",axs_Keff)
tau_R=NFL.Rayleigh_Function(P,leff_NH3,"647CNT",axs_Keff)
tmp=NFL.Compute_Transmission(P,tau_R,tau_NH3,"647CNT",axs_trans,axs_Keff)


keff_CH4,leff_CH4=NFL.K_eff(P,Transmission632,CH4,622.,642.,"632OI",axs_Keff)
tau_CH4=NFL.Weighting_Function(P,keff_CH4,"632OI",axs_Keff)
tau_R=NFL.Rayleigh_Function(P,leff_CH4,"632OI",axs_Keff)
tmp=NFL.Compute_Transmission(P,tau_R,tau_CH4,"632OI",axs_trans,axs_Keff)

keff_CH4,leff_CH4=NFL.K_eff(P,Transmission620,CH4,610.,630.,"620CH4",axs_Keff)
tau_CH4=NFL.Weighting_Function(P,keff_CH4,"620CH4",axs_Keff)
tau_R=NFL.Rayleigh_Function(P,leff_CH4,"620CH4",axs_Keff)
tmp=NFL.Compute_Transmission(P,tau_R,tau_CH4,"620CH4",axs_trans,axs_Keff)

axs_Keff.legend(loc=1,ncol=3, borderaxespad=0.,prop={'size':10})
fig_Keff.subplots_adjust(left=0.12, right=0.96, top=0.94, bottom=0.09)
axs_trans.legend(loc=1,ncol=4, borderaxespad=0.,prop={'size':7})
fig_trans.subplots_adjust(left=0.12, right=0.96, top=0.94, bottom=0.09)

fig_trans.savefig('c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/TransmissionFunctions.png',dpi=320,bbox_inches = 'tight')
fig_Keff.savefig('c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/ContributionFunctions.png',dpi=320,bbox_inches = 'tight')

###### Plot regional and moon reflectivities
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

Callisto1980 = np.fromfile(file="c:/Astronomy/Projects/Planets/JovianMoons/References/callisto_no_header.txt", dtype=float, count=-1, sep=" ")    
Callisto1980=np.reshape(Callisto1980,[int(Callisto1980.size/3),3])

WaveGrid,SignalonGrid=GSU.uniform_wave_grid(Callisto1980[:,0]*1000.,Callisto1980[:,1],
                                        Extend=False,Fine=False)
CalGrid=np.zeros((WaveGrid.size,2))
CalGrid[:,0]=WaveGrid
CalGrid[:,1]=SignalonGrid

CallistoProduct672=GSU.SpectrumMath(Transmission672,CalGrid,"Multiply")
CallistoProduct658=GSU.SpectrumMath(Transmission658,CalGrid,"Multiply")
CallistoProduct656=GSU.SpectrumMath(Transmission656,CalGrid,"Multiply")
CallistoProduct647=GSU.SpectrumMath(Transmission647,CalGrid,"Multiply")
CallistoProduct632=GSU.SpectrumMath(Transmission632,CalGrid,"Multiply")
CallistoProduct620=GSU.SpectrumMath(Transmission620,CalGrid,"Multiply")


Ganymede1980 = np.fromfile(file="c:/Astronomy/Projects/Planets/JovianMoons/References/ganymede_no_header.txt", dtype=float, count=-1, sep=" ")    
Ganymede1980=np.reshape(Ganymede1980,[int(Ganymede1980.size/3),3])

WaveGrid,SignalonGrid=GSU.uniform_wave_grid(Ganymede1980[:,0]*1000.,Ganymede1980[:,1],
                                        Extend=False,Fine=False)
GanGrid=np.zeros((WaveGrid.size,2))
GanGrid[:,0]=WaveGrid
GanGrid[:,1]=SignalonGrid

GanymedeProduct672=GSU.SpectrumMath(Transmission672,GanGrid,"Multiply")
GanymedeProduct658=GSU.SpectrumMath(Transmission658,GanGrid,"Multiply")
GanymedeProduct656=GSU.SpectrumMath(Transmission656,GanGrid,"Multiply")
GanymedeProduct647=GSU.SpectrumMath(Transmission647,GanGrid,"Multiply")
GanymedeProduct632=GSU.SpectrumMath(Transmission632,GanGrid,"Multiply")
GanymedeProduct620=GSU.SpectrumMath(Transmission620,GanGrid,"Multiply")

Europa1980 = np.fromfile(file="c:/Astronomy/Projects/Planets/JovianMoons/References/europa_no_header.txt", dtype=float, count=-1, sep=" ")    
Europa1980=np.reshape(Europa1980,[int(Europa1980.size/3),3])

WaveGrid,SignalonGrid=GSU.uniform_wave_grid(Europa1980[:,0]*1000.,Europa1980[:,1],
                                        Extend=False,Fine=False)
EurGrid=np.zeros((WaveGrid.size,2))
EurGrid[:,0]=WaveGrid
EurGrid[:,1]=SignalonGrid

EuropaProduct672=GSU.SpectrumMath(Transmission672,EurGrid,"Multiply")
EuropaProduct658=GSU.SpectrumMath(Transmission658,EurGrid,"Multiply")
EuropaProduct656=GSU.SpectrumMath(Transmission656,EurGrid,"Multiply")
EuropaProduct647=GSU.SpectrumMath(Transmission647,EurGrid,"Multiply")
EuropaProduct632=GSU.SpectrumMath(Transmission632,EurGrid,"Multiply")
EuropaProduct620=GSU.SpectrumMath(Transmission620,EurGrid,"Multiply")

Io_leading1980 = np.fromfile(file="c:/Astronomy/Projects/Planets/JovianMoons/References/io.leading_no_header.txt", dtype=float, count=-1, sep=" ")    
Io_leading1980=np.reshape(Io_leading1980,[int(Io_leading1980.size/3),3])

Io_trailing1980 = np.fromfile(file="c:/Astronomy/Projects/Planets/JovianMoons/References/io.trailing_no_header.txt", dtype=float, count=-1, sep=" ")    
Io_trailing1980=np.reshape(Io_trailing1980,[int(Io_trailing1980.size/3),3])

WaveGrid,SignalonGrid=GSU.uniform_wave_grid(Io_trailing1980[:,0]*1000.,Io_trailing1980[:,1],
                                        Extend=False,Fine=False)
Io_Grid=np.zeros((WaveGrid.size,2))
Io_Grid[:,0]=WaveGrid
Io_Grid[:,1]=SignalonGrid

IoProduct672=GSU.SpectrumMath(Transmission672,Io_Grid,"Multiply")
IoProduct658=GSU.SpectrumMath(Transmission658,Io_Grid,"Multiply")
IoProduct656=GSU.SpectrumMath(Transmission656,Io_Grid,"Multiply")
IoProduct647=GSU.SpectrumMath(Transmission647,Io_Grid,"Multiply")
IoProduct632=GSU.SpectrumMath(Transmission632,Io_Grid,"Multiply")
IoProduct620=GSU.SpectrumMath(Transmission620,Io_Grid,"Multiply")

axs_moons[0].plot(Callisto1980[:,0]*1000.,Callisto1980[:,1],label='Callisto',linewidth=1,color='b')
axs_moons[0].plot(Ganymede1980[:,0]*1000.,Ganymede1980[:,1],label='Ganymede',linewidth=1,color='g')
axs_moons[0].plot(Europa1980[:,0]*1000.,Europa1980[:,1],label='Europa',linewidth=1,color='r')
axs_moons[0].plot(Io_leading1980[:,0]*1000.,Io_leading1980[:,1],label='Io Leading',linewidth=1,linestyle='dashed',color='k')
axs_moons[0].plot(Io_trailing1980[:,0]*1000.,Io_trailing1980[:,1],label='Io Trailing',linewidth=1,color='k')
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

axs_moons[1].plot(CallistoProduct672[:,0],CallistoProduct672[:,1],label='Callisto',linewidth=0.5,color='b')
axs_moons[1].plot(CallistoProduct658[:,0],CallistoProduct658[:,1],linewidth=0.5,color='b')
axs_moons[1].plot(CallistoProduct656[:,0],CallistoProduct656[:,1],linewidth=0.5,color='b')
axs_moons[1].plot(CallistoProduct647[:,0],CallistoProduct647[:,1],linewidth=0.5,color='b')
axs_moons[1].plot(CallistoProduct632[:,0],CallistoProduct632[:,1],linewidth=0.5,color='b')
axs_moons[1].plot(CallistoProduct620[:,0],CallistoProduct620[:,1],linewidth=0.5,color='b')

axs_moons[1].plot(GanymedeProduct672[:,0],GanymedeProduct672[:,1],label='Ganymede',linewidth=0.5,color='g')
axs_moons[1].plot(GanymedeProduct658[:,0],GanymedeProduct658[:,1],linewidth=0.5,color='g')
axs_moons[1].plot(GanymedeProduct656[:,0],GanymedeProduct656[:,1],linewidth=0.5,color='g')
axs_moons[1].plot(GanymedeProduct647[:,0],GanymedeProduct647[:,1],linewidth=0.5,color='g')
axs_moons[1].plot(GanymedeProduct632[:,0],GanymedeProduct632[:,1],linewidth=0.5,color='g')
axs_moons[1].plot(GanymedeProduct620[:,0],GanymedeProduct620[:,1],linewidth=0.5,color='g')

axs_moons[1].plot(EuropaProduct672[:,0],EuropaProduct672[:,1],label='Europa',linewidth=0.5,color='r')
axs_moons[1].plot(EuropaProduct658[:,0],EuropaProduct658[:,1],linewidth=0.5,color='r')
axs_moons[1].plot(EuropaProduct656[:,0],EuropaProduct656[:,1],linewidth=0.5,color='r')
axs_moons[1].plot(EuropaProduct647[:,0],EuropaProduct647[:,1],linewidth=0.5,color='r')
axs_moons[1].plot(EuropaProduct632[:,0],EuropaProduct632[:,1],linewidth=0.5,color='r')
axs_moons[1].plot(EuropaProduct620[:,0],EuropaProduct620[:,1],linewidth=0.5,color='r')

axs_moons[1].plot(IoProduct672[:,0],IoProduct672[:,1],label='Io Trailing',linewidth=0.5,color='k')
axs_moons[1].plot(IoProduct658[:,0],IoProduct658[:,1],linewidth=0.5,color='k')
axs_moons[1].plot(IoProduct656[:,0],IoProduct656[:,1],linewidth=0.5,color='k')
axs_moons[1].plot(IoProduct647[:,0],IoProduct647[:,1],linewidth=0.5,color='k')
axs_moons[1].plot(IoProduct632[:,0],IoProduct632[:,1],linewidth=0.5,color='k')
axs_moons[1].plot(IoProduct620[:,0],IoProduct620[:,1],linewidth=0.5,color='k')

axs_moons[1].legend(fontsize=8)
fig_moons.subplots_adjust(left=0.10, right=0.90, top=0.94, bottom=0.09)

fig_moons.savefig('c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/ColorSlopes.png',dpi=320)

###############################################################################


ContimIntegral620,AbsorpIntegral620=NFL.cont_absorption_calcs(ContinuumProduct620,AbsorptionProduct620,615.,635.,"620CH4")
ContimIntegral632,AbsorpIntegral632=NFL.cont_absorption_calcs(ContinuumProduct632,AbsorptionProduct632,622.,642.,"632OI")
ContimIntegral647,AbsorpIntegral647=NFL.cont_absorption_calcs(ContinuumProduct647,AbsorptionProduct647,635.,659.,"647CNT")
ContimIntegral656,AbsorpIntegral656=NFL.cont_absorption_calcs(ContinuumProduct656,AbsorptionProduct656,645.,672.,"656HIA")
ContimIntegral658,AbsorpIntegral658=NFL.cont_absorption_calcs(ContinuumProduct658,AbsorptionProduct658,653.,663.,"658NII")
ContimIntegral672,AbsorpIntegral672=NFL.cont_absorption_calcs(ContinuumProduct672,AbsorptionProduct672,662.,682.,"672SII")
ContimIntegral730,AbsorpIntegral730=NFL.cont_absorption_calcs(ContinuumProduct730,AbsorptionProduct730,720.,740.,"730OII")
ContimIntegral889,AbsorpIntegral889=NFL.cont_absorption_calcs(ContinuumProduct889,AbsorptionProduct889,877.,902.,"889CH4")
ContimIntegral940,AbsorpIntegral940=NFL.cont_absorption_calcs(ContinuumProduct940,AbsorptionProduct940,926.,954.,"940NIR")

print()
print("########### Jupiter 647/656 ratios")
print("647/656 Continuum=",ContimIntegral647/ContimIntegral656)
print("647/656 with NH3=",AbsorpIntegral647/AbsorpIntegral656)
print("(647/656 with NH3)/(647/656 Continuum)=",(AbsorpIntegral647/AbsorpIntegral656)/(ContimIntegral647/ContimIntegral656))
print()

print()
print("########### Jupiter 889/940 ratios")
print("889/940 Continuum=",ContimIntegral889/ContimIntegral940)
print("889/940 with NH3=",AbsorpIntegral889/AbsorpIntegral940)
print("(889/940 with NH3)/(889/940 Continuum)=",(AbsorpIntegral889/AbsorpIntegral940)/(ContimIntegral889/ContimIntegral940))
print()

###############################################################################
StartIndex647=np.where(CallistoProduct647[:,0]==635.)
EndIndex647=np.where(CallistoProduct647[:,0]==659.)

CalIntegral647=sum(CallistoProduct647[StartIndex647[0][0]:EndIndex647[0][0],1])
GanIntegral647=sum(GanymedeProduct647[StartIndex647[0][0]:EndIndex647[0][0],1])
EurIntegral647=sum(EuropaProduct647[StartIndex647[0][0]:EndIndex647[0][0],1])
IoIntegral647=sum(IoProduct647[StartIndex647[0][0]:EndIndex647[0][0],1])

StartIndex656=np.where(CallistoProduct647[:,0]==645.)
EndIndex656=np.where(CallistoProduct647[:,0]==672.)

CalIntegral656=sum(CallistoProduct656[StartIndex656[0][0]:EndIndex656[0][0],1])
GanIntegral656=sum(GanymedeProduct656[StartIndex656[0][0]:EndIndex656[0][0],1])
EurIntegral656=sum(EuropaProduct656[StartIndex656[0][0]:EndIndex656[0][0],1])
IoIntegral656=sum(IoProduct656[StartIndex656[0][0]:EndIndex656[0][0],1])
print()
print("########### Moons 647/656")
print("Callisto=",CalIntegral647/CalIntegral656)
print("Ganymede=",GanIntegral647/GanIntegral656)
print("Europa=",EurIntegral647/EurIntegral656)
print("Io=",IoIntegral647/IoIntegral656)
print()
print("########### Moons Correction Factor to Jupiter Continuum Estimate")
print("Callisto=",(ContimIntegral647/ContimIntegral656)/(CalIntegral647/CalIntegral656))
print("Ganymede=",(ContimIntegral647/ContimIntegral656)/(GanIntegral647/GanIntegral656))
print("Europa=",(ContimIntegral647/ContimIntegral656)/(EurIntegral647/EurIntegral656))
print("Io=",(ContimIntegral647/ContimIntegral656)/(IoIntegral647/IoIntegral656))

