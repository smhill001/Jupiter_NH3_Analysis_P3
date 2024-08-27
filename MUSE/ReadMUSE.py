def ReadMUSE(date):
    """
    Created on Thu Jan 12 14:21:05 2023
    
    PROGRAM: ReadMUSE
    PURPOSE: Read MUSE fits file for Jupiter observation on 9/19/2022 and
                1) Read MUSE FITS cubes
                2) Set up filter metadata
    @author: smhil
    """
    
    import sys
    drive='c:'
    sys.path.append(drive+'/Astronomy/Python Play')
    sys.path.append(drive+'/Astronomy/Python Play/Util_P3')
    sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy_P3')
    
    import os
    import pylab as pl
    import numpy as np
    from imageio import imwrite
    from astropy.io import fits
    import GeneralSpecUtils_P3 as GSU
    
    ###########################################################################
    # Read MUSE FITS cube and set wavelength
    ###########################################################################
    path='c:/Astronomy/Projects/SAS 2021 Ammonia/VLT MUSE/'
    if date=='20220730UT':
        MUSEfile='2022-07-30_P109_006_obs9_smooth.fits'
        delta=0
        wavelength=np.linspace(471.+delta,935.+delta,465)
    elif date=='20220919UT':
        MUSEfile='2022-09-19_obs18_proj.fits'
        delta=-0.0
        wavelength=np.linspace(470.+delta,935.13+delta,3722)
    if date=='20230812UT':
        MUSEfile='23_08_12_obs45smooth_proj_DATACUBE_FINAL.fits'
        delta=0
        wavelength=np.linspace(471.+delta,935.+delta,465)
    if date=='20230815UT':
        MUSEfile='23_08_15_obs49smooth_proj_DATACUBE.fits'
        delta=0
        wavelength=np.linspace(471.+delta,935.+delta,465)
    if date=='20230816UT':
        MUSEfile='23_08_16_obs55_smooth_proj_DATACUBE.fits'
        delta=0
        wavelength=np.linspace(471.+delta,935.+delta,465)
    if date=='20230923UT':
        MUSEfile='23_09_23_obs90_smooth_proj_DATACUBE.fits'
        delta=0
        wavelength=np.linspace(471.+delta,935.+delta,465)

        
    MUSEhdulist=fits.open(path+MUSEfile)
    MUSEhdulist.info()
    MUSEhdr=MUSEhdulist[0].header
    MUSEdata=MUSEhdulist[1].data
    MUSEhdr1=MUSEhdulist[1].header
    MUSEzen=MUSEhdulist[5].data
    MUSEszen=MUSEhdulist[6].data
    MUSEhdulist.close()
    
    ###########################################################################
    # Load filter metadata for MUSE
    ########################################################################### 
    filterwavelength=['620','632','647','656','658','672','730','889','940']
    """
    #Filter wv lims geared toward the core of the filter band for direct integration
    #without filter convolution
    filterdata={'620':{'transfile':'620CH4/620CH4_Transmission.txt',
                       'filtname':'620CH4','filtwdth':10.,'wvs':[615.,625.]},
                 '632':{'transfile':'632OI/632OI_Transmission.txt',
                        'filtname':'632OI','filtwdth':10.,'wvs':[627.,637.]},
                 '647':{'transfile':'647CNT/647CNT_Transmission.txt',
                        'filtname':'647NH3','filtwdth':10.,'wvs':[642.,652.]},
                 '656':{'transfile':'656HIA/656HIA_Transmission.txt',
                        'filtname':'656HIA','filtwdth':10.,
                        'wvsRED':[658.0,660.5],'wvsBLU':[651.0,654.0]},
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
    """
    filterdata={'620':{'transfile':'620CH4/620CH4_Transmission.txt',
                       'filtname':'620CH4','filtwdth':10.,'wvs':[610.,630.]},
                 '632':{'transfile':'632OI/632OI_Transmission.txt',
                        'filtname':'632OI','filtwdth':10.,'wvs':[622.,642.]},
                 '647':{'transfile':'647CNT/647CNT_Transmission.txt',
                        'filtname':'647NH3','filtwdth':10.,'wvs':[637.,657.]},
                 '656':{'transfile':'656HIA/656HIA_Transmission.txt',
                        'filtname':'656HIA','filtwdth':10.,
                        'wvsRED':[658.0,666.0],'wvsBLU':[646.0,654.0]},
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
    return(MUSEhdr,MUSEhdr1,MUSEdata,MUSEzen,MUSEszen,wavelength,filterdata,path)