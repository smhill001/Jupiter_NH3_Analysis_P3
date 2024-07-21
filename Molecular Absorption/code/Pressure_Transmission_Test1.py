# -*- coding: utf-8 -*-
"""
Created on Sun Jul 14 15:05:48 2024

@author: smhil
"""

import numpy as np

Keff=0.427
Pbar=10 #bar
print("Pbar=",Pbar)
#From NH3_Filter_Library_P3.py
amagat=2.69e24 #Lodschmits number?
#gravity=2228.0
gravity=2479.0 #cm/s^2
mean_mol_wt=3.85e-24 #gm/molecule, which is 2.22 gm/mole
fCH4=2.04e-3#1.81e-3
STP=1.01325e6

P=Pbar*1e5 #Pascals
kmatm=(P/1.0e5)*STP*fCH4/(amagat*mean_mol_wt*gravity)

print("kmatm=",kmatm)

CH4_tau=kmatm*Keff
print("CH4_tau=",CH4_tau)
#tau_mend=22.4e4*(P/1.0e5)*fgas*Keff/(2.22*gravity)
#print("tau_gas/tau_mend = ",tau_gas/tau_mend)
trans=np.exp(-2.0*CH4_tau)
print("trans=",trans)




###########################################################################
#  SET NECESSARY CONSTANTS from make_L3_env_data.py
###########################################################################
amagat=2.69e24 #Lodschmits number. (cm-2) This is really km-amagat
#gravity=2228.0 #cm/s^2
gravity=2479.0 #cm/s^2
mean_mol_wt=3.85e-24 #gm/molecule, which is 2.22 gm/mole
fCH4=2.04e-3#1.81e-3
STP=1.01e6  #dyne/cm^2 [(g-cm/s^2)/cm^2]`

CH4_Ncol=CH4_tau/Keff
print("CH4_Ncol=",CH4_Ncol)
CH4_Cloud_Press=(CH4_Ncol)*amagat*gravity*mean_mol_wt/(fCH4*STP)
print("CH4_Cloud_Press=",CH4_Cloud_Press)
print()
###########################################################################
#From Hill et al paper
###########################################################################

R=8.31446261815324 #J K^-1 mol-1
kb=1.380649e-23
print("R/kb=",R/kb)
gravity=24.79 #m/s^2

mean_mol_wt=0.00222 #kg/mole
NA=6.02214097e23 #Avogadro's number (mol^-1), recipricol moles
amagat=2.68678e25 #molecules per cubic meter at STP
ratio=amagat/NA
print("ratio km-amagat to mol/m**2=",1000*ratio)

P=500000. #Pascals
eta=2.0 #nadir viewing

#eq 5 - The one way column density at a given pressure
X=fCH4*P*NA/(mean_mol_wt*gravity)
print("X (molecules/m^2)=",X)
NCH4=X/(NA*44615)
print("NCH4=",NCH4) #km-amagat


nadirtau=eta*NCH4*Keff
print("nadirtau=",nadirtau)
nadirtrans=np.exp(-nadirtau)
print("nadirtrans=",nadirtrans)


#eq 6a
Pcl1=X*mean_mol_wt*gravity/(fCH4*NA)
print("Pcl1=",Pcl1)
