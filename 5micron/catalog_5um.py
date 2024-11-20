def catalog_5um():
    
    import sys
    drive='c:'
    sys.path.append(drive+'/Astronomy/Python Play')
    sys.path.append(drive+'/Astronomy/Python Play/Util_P3')
    sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy')
    sys.path.append('./Services')

    from astropy.io import fits
    sys.path.append('./Services')
    import os
    from astropy.table import Table
    from astropy.io import ascii


    import numpy as np

    path5um="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/5micron/FITS/output/"
    fitsfilelist=os.listdir(path5um)
    print(fitsfilelist)
    
    filename="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/5micron/5umCatalog.csv"

    t = Table(names=('ObsID', 'Time','CM1','CM2','CM3',
                     'Telescope','FL','Camera','Seeing','Transparency','File Name'),
              dtype=('S10','S10','S10','S10','S10',
                     'S10','S10','S10','S10','S10','S20'))
    
    for fn in fitsfilelist:
        hdulist=fits.open(path5um+fn)
        hdulist.info()
        hdr=hdulist[0].header
        hdulist.close()

        #for item in sourcefiles[ObsID]:
        #    print(item)
        print("%%",hdr["DATE-OBS"])
        print("%%",hdr["DATE-OBS"][17:19])
        sec=str(int(str(hdr["DATE-OBS"][17:19])))
        time=hdr["DATE-OBS"]
        CM1=float(hdr["CM1"])
        CM2=float(hdr["CM2"])
        CM3=float(hdr["CM3"])

        print("ObsID",hdr["TELESCOP"],"FL",hdr["INSTRUME"],"Seeing","Transparency",fn)
        t.add_row(("ObsID",(time.replace("T"," ")).replace("Z",""),str(CM1),str(CM2),str(CM3),hdr["TELESCOP"],
                   "FL",hdr["INSTRUME"],"Seeing","Transparency",fn))
        
    ascii.write(t,filename,format='basic',overwrite=True,delimiter=',')


