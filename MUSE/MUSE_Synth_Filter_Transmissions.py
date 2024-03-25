# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 12:38:27 2021

This code creates the transmission files for synthetic filters used for
analysis of the MUSE-VLT data from 9/19/2023. 

@author: Steven Hill
"""
import sys
sys.path.append('c:/Astronomy/Python Play')
sys.path.append('c:/Astronomy/Python Play/Util_P3')
sys.path.append('c:/Astronomy/Python Play/SPLibraries_P3')
sys.path.append('c:/Astronomy/Python Play/SpectroPhotometry/Spectroscopy_P3')
import matplotlib.pyplot as pl
import numpy as np
import GeneralSpecUtils_P3 as GSU

path='c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/MUSE/'

WaveGrid=np.arange(115,1062.5,0.5,dtype=float)
F620_wvs=[615.,625.]
F620_idx=[np.argmin(abs(F620_wvs[0]-WaveGrid)),np.argmin(abs(F620_wvs[1]-WaveGrid))]

Filter620=np.zeros(len(WaveGrid))
Filter620[F620_idx[0]:F620_idx[1]]=1.0
#Filter620[Filter620==0]=np.NaN
print(np.column_stack((WaveGrid,Filter620)).shape)
Transmission620=np.column_stack((WaveGrid,Filter620))
np.savetxt(path+'620CH4_MUSE_Transmission.txt',Transmission620,delimiter=" ",
           fmt="%10.3F %10.7F")

F632_wvs=[627.,637.]
F632_idx=[np.argmin(abs(F632_wvs[0]-WaveGrid)),np.argmin(abs(F632_wvs[1]-WaveGrid))]

Filter632=np.zeros(len(WaveGrid))
Filter632[F632_idx[0]:F632_idx[1]]=1.0
#Filter632[Filter632==0]=np.NaN
print(np.column_stack((WaveGrid,Filter632)).shape)
Transmission632=np.column_stack((WaveGrid,Filter632))
np.savetxt(path+'632OI_MUSE_Transmission.txt',Transmission632,delimiter=" ",
           fmt="%10.3F %10.7F")

F647_wvs=[642.,652.]
F647_idx=[np.argmin(abs(F647_wvs[0]-WaveGrid)),np.argmin(abs(F647_wvs[1]-WaveGrid))]

Filter647=np.zeros(len(WaveGrid))
Filter647[F647_idx[0]:F647_idx[1]]=1.0
#Filter647[Filter647==0]=np.NaN
print(np.column_stack((WaveGrid,Filter647)).shape)
Transmission647=np.column_stack((WaveGrid,Filter647))
np.savetxt(path+'647NH3_MUSE_Transmission.txt',Transmission647,delimiter=" ",
           fmt="%10.3F %10.7F")

F656_wvsBLU=[652.0,655.0]
F656_wvsRED=[658.0,661.0]
F656_idxRED=[np.argmin(abs(F656_wvsRED[0]-WaveGrid)),np.argmin(abs(F656_wvsRED[1]-WaveGrid))]
F656_idxBLU=[np.argmin(abs(F656_wvsBLU[0]-WaveGrid)),np.argmin(abs(F656_wvsBLU[1]-WaveGrid))]

Filter656=np.zeros(len(WaveGrid))
Filter656[F656_idxRED[0]:F656_idxRED[1]]=1.0
Filter656[F656_idxBLU[0]:F656_idxBLU[1]]=1.0
#Filter656[Filter656==0]=np.NaN
print(np.column_stack((WaveGrid,Filter656)).shape)
Transmission656=np.column_stack((WaveGrid,Filter656))
np.savetxt(path+'656HIA_MUSE_Transmission.txt',Transmission656,delimiter=" ",
           fmt="%10.3F %10.7F")

###### Plot filter transmissions convolved with disk-integrated albedos
pl.figure(figsize=(6.5, 4.0), dpi=150, facecolor="white")

x0=600.
x1=680.
xtks=9
y0=0.0
y1=1.1
ytks=12

# Set x limits
pl.xlim(x0,x1)
# Set x ticks
pl.xticks(np.linspace(x0,x1,xtks, endpoint=True))
# Set y limits
pl.ylim(y0,y1)
pl.yticks(np.linspace(y0,y1,ytks, endpoint=True))
# Set y ticks
pl.grid(linewidth=0.2)
pl.tick_params(axis='both', which='major', labelsize=8)
pl.title("VLT-MUSE Synthetic Filters for Ammonia Investigation",fontsize=12)
pl.ylabel("Filter Transmission",fontsize=10,color="black")
pl.xlabel("Wavelength (nm)",fontsize=10)

zeros=np.zeros(Filter620.shape[0])

pl.plot(Transmission656[:,0],Transmission656[:,1],linewidth=1,color='C2',label='656HIA')
pl.fill_between(Transmission656[:,0], zeros, Transmission656[:,1],color='C0',alpha=.2,label='Continuum')
pl.plot(Transmission647[:,0],Transmission647[:,1],linewidth=1,color='C3',label='647NH3')
pl.fill_between(Transmission647[:,0], zeros, Transmission647[:,1],color='C3',alpha=.2,label='Ammonia')
pl.plot(Transmission632[:,0],Transmission632[:,1],linewidth=1,color='C4',label='632OI')
pl.fill_between(Transmission632[:,0], zeros, Transmission632[:,1],color='C0',alpha=.2)
pl.plot(Transmission620[:,0],Transmission620[:,1],linewidth=1,color='C8',label='620CH4')
pl.fill_between(Transmission620[:,0], zeros, Transmission620[:,1],color='C8',alpha=.2,label='Methane')

pl.legend(fontsize=7)
pl.subplots_adjust(left=0.09, bottom=0.12, right=0.97, top=0.92)  


pl.savefig('c:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/MUSE/MUSE_Synth_Filter_Transmissions.png',dpi=320)

