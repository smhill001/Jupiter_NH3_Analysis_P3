# -*- coding: utf-8 -*-
"""
Created on Sun Mar 21 23:45:33 2021

PURPOSE:    Compute count rates for multiple targets in an image. This version
            is particularly focused on Jupiter and it's moons and so only
            has two fixed/configurable apertures, one for 'moons' and 
            one for 'Jupiter'

UPDATE 2022-01-25:  Converted to Python 3
@author: Steven Hill
"""

def ComputeNetRateJupiter(scidata,header,TargetIDs,SessionID,positions,radii,figure,axis):
    from photutils import CircularAperture
    from photutils import aperture_photometry
    from photutils import CircularAnnulus
    import Meta_and_Control_Data_Operations_P3 as Meta
    from astropy.table import Table, hstack
    import pylab as pl

    #Create aperture objects    
    apertures = CircularAperture(positions, r=radii[0])
    annulus_apertures = CircularAnnulus(positions, r_in=radii[1], r_out=radii[2])

    #Compute raw fluxes
    rawflux_table = aperture_photometry(scidata, apertures)
    bkgflux_table = aperture_photometry(scidata, annulus_apertures)
    
    phot_table = hstack([rawflux_table, bkgflux_table], table_names=['raw', 'bkg'])
    print("phot_table['aperture_sum_bkg']=",phot_table['aperture_sum_bkg'])
    print("annulus_apertures.area()=",annulus_apertures.area)
    bkg_mean = phot_table['aperture_sum_bkg'] / annulus_apertures.area
    #print "bkg_mean=",bkg_mean
    bkg_sum = bkg_mean * apertures.area
    #print "bkg_sum=",bkg_sum
    final_sum = phot_table['aperture_sum_raw'] - bkg_sum
    #print "final_sum=",final_sum
    rate=final_sum/header['EXPTIME']
    phot_table['net_count_rate'] = rate
    #print "Raw=",rawflux_table
    #print 'Bkg=',bkgflux_table
    phot_table['Target']=TargetIDs
    phot_table['Filter']=header['Filter']
    phot_table['Date-Obs']=header['MIDPOINT']
    phot_table['SessionID']=SessionID
    phot_table.remove_column('id_bkg')
    phot_table.remove_column('xcenter_bkg')
    phot_table.remove_column('ycenter_bkg')
    #print phot_table
    Filter=Meta.FilterParameters(header['FILTER'])
    WVCenter=Filter.CenterWV###Testing Area
    
    
    #Code to display diagnostic image plots with photometric annuli overlaid
    
    if axis!='None':
            
        if TargetIDs[0]=='0_Jupiter':
            vmx=5.0e5
        elif TargetIDs[0]!='0_Jupiter':
            vmx=5.0e4
            
        axis.imshow(scidata,vmin=0,vmax=vmx)
    
        pl.axes(axis)
        ap_patches = apertures.plot(color='white', lw=0.5,
                               label='Photometry aperture')
        ann_patches = annulus_apertures.plot(color='red', lw=0.5,
                                        label='Background annulus')
        
        title=header['FILTER']
        axis.set_title(title,fontsize=7)
        xlbl=TargetIDs[0][2:]
        for ixlbl in range(1,len(TargetIDs)):
            xlbl=xlbl+" "+TargetIDs[ixlbl][2:]
        axis.set_xlabel(xlbl,fontsize=7)
        axis.tick_params(axis='both', which='major', labelsize=7)
    
        axis.set_adjustable('box') 

    #labels = (ap_patches[0], ann_patches[0])
    #pl.legend(font=10)
    
    return rate,WVCenter,phot_table

def uniform_lat_grid(Latitude,Signal,Fine=False):
    """
    Takes an existing latitude profile on a non-standard or even irregular
    grid and performs linear interpolation to place the data
    on one of two uniform grids:
        1) -90 to 90 with 1 deg bins
        2) -90 to 90 with 2 deg bins
    """
    import numpy as np
    from scipy import interpolate

    if Fine:    #Set grid interval
        dlat=1.0
    else:
        dlat=2.0

    LatGrid=np.arange(-90.,90.1,dlat,dtype=float)
    #print Wavelength.size,Signal.size    
    Interp=interpolate.interp1d(Latitude,Signal,kind='linear', 
                                copy=True,bounds_error=False, 
                                fill_value=np.NaN,axis=0)  
    SignalonGrid=Interp(LatGrid)

    return LatGrid,SignalonGrid