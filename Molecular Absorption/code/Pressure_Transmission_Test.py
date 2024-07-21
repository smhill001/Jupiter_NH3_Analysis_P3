###########################################################################
#From Hill et al paper
###########################################################################
import numpy as np

###########################################################################
#Inputs
###########################################################################
P=5e5 #Pascals
print("Pressure (Pa)=",P)
eta=2.0 #nadir viewing
print("Jovian Airmass Factor, eta=",eta)

###########################################################################
#Constants & Checks
###########################################################################

R=8.31446261815324 #J K^-1 mol-1
kb=1.380649e-23
print("R/kb=",R/kb)
gravity=24.79 #m/s^2
fCH4=2.04e-3
Keff=0.427 #(km-amagat)^-1

mean_mol_wt=0.00222 #kg/mole
NA=6.02214097e23 #Avogadro's number (mol^-1), recipricol moles
amagat=2.68678e25 #molecules per cubic meter at STP
ratio=amagat/NA
print("ratio km-amagat to mol/m**2=",1000*ratio)


###########################################################################
#eq 5 - The one way column density at a given pressure
###########################################################################
X=fCH4*P*NA/(mean_mol_wt*gravity)
print("X (molecules/m^2)=",X)
###########################################################################
#eq 6 - The one way column density at a given pressure
###########################################################################
#Converting to km-amagat
NCH4=X/(NA*44615)
print("NCH4 (km-amagat)=",NCH4) #km-amagat

###########################################################################
#eq 3 - Compute tau and transmission for the given airmass factor, eta
###########################################################################

nadirtau=eta*NCH4*Keff
print("nadirtau=",nadirtau)
nadirtrans=np.exp(-nadirtau)
print("nadirtrans=",nadirtrans)
