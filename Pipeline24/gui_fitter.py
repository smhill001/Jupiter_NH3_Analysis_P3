def gui_fitter(dataset,planet='Saturn'):

    import glob
    import planetmapper
    
    print(dataset[0])
    print()
    print(dataset[1])
    
    file_list = dataset[1]
    path = dataset[0]
    path_file_list = [path + x for x in file_list]
    #path_file_list = [x + string for x in my_list]
    
    fnclip=[file_list[0]]
    print("in gui_fitter")
    print(fnclip)
    #for fn in file_list:
    RGB=["685NIR","550GRN","450BLU"]
    science=["656HIA","632OI","620CH4","647NH3"]
    
    #for k in range(0,fn_bases.nrecords-1):
    #    fn_base=fn_bases.VideoFile[k].split('.')[0]
    #    print(fn_base)
    First=True
    for filtr in RGB:
        print(filtr)
        for fn in file_list:
            if filtr in fn:
                time=fn[0:10]+"T"+fn[11:13]+":"+fn[13:15]
                print(time)
                # Running from Python allows you to customise SPICE settings like the aberration correction
                observation = planetmapper.Observation(path+fn,target=planet,utc=time)
                
                #observation.set_plate_scale_arcsec(0.11141)

                #observation.get_disc_params(246.6,377.5)
                #observation.set_x0(246.6)
                #observation.set_y0(377.5)
                
                #params=observation.get_disc_params()
                #print("######### params2=",params)
                
                print("1##########observation.backplanes=",list(observation.backplanes.keys()))
                
                del observation.backplanes['DOPPLER']
                del observation.backplanes['LON-CENTRIC']
                del observation.backplanes['LAT-CENTRIC']
                del observation.backplanes['RA']
                del observation.backplanes['DEC']
                del observation.backplanes['KM-X']
                del observation.backplanes['KM-Y']
                del observation.backplanes['RING-RADIUS']
                del observation.backplanes['RING-LON-GRAPHIC']
                del observation.backplanes['RING-DISTANCE']
                del observation.backplanes['LIMB-LON-GRAPHIC']
                del observation.backplanes['LIMB-DISTANCE']
                del observation.backplanes['RADIAL-VELOCITY']
                del observation.backplanes['LOCAL-SOLAR-TIME']
                del observation.backplanes['AZIMUTH']
                del observation.backplanes['PHASE']
                del observation.backplanes['LIMB-LAT-GRAPHIC']
                del observation.backplanes['ANGULAR-X']
                del observation.backplanes['ANGULAR-Y']
                del observation.backplanes['PIXEL-X']
                del observation.backplanes['PIXEL-Y']
                del observation.backplanes['DISTANCE']
                print("******************")
                print("2##########observation.backplanes=",observation.backplanes.keys())
                
                #observation.save_observation(path+fn.replace("png","fits"))
            
                # Run some custom setup
                #observation.add_other_bodies_of_interest('Io', 'Europa', 'Ganymede', 'Callisto')
                #observation.set_plate_scale_arcsec(42) # set a custom plate scale
                #observation.rotation_from_wcs() # get the disc rotation from the header's WCS info
            
                # Run the GUI to fit the observation interactively
                # This will open a GUI window every loop
                if First:
                    coords = observation.run_gui()
                    #print("coords",coords)
                    params=observation.get_disc_params()
                    print("######### params1=",params)
                else:
                    observation.set_disc_params(params[0],params[1],params[2],params[3])

                #observation.set_coords(coords)
                observation.save_observation(path+fn.replace(".png",".fits"))
                observation.save_mapped_observation(path+fn.replace(".png","map.fits"))
                First=False

                # More custom code can go here to use the fitted observation...
                # for example, we can print some values for the last click location
                #if coords:
                #    x, y = coords[-1]
                #    print(observation.xy2lonlat(x, y))
    for filtr in science:
        print(filtr)
        for fn in file_list:
            if filtr in fn:
                time=fn[0:10]+"T"+fn[11:13]+":"+fn[13:15]
                print(time)
                # Running from Python allows you to customise SPICE settings like the aberration correction
                observation = planetmapper.Observation(path+fn,target=planet,utc=time)
                del observation.backplanes['DOPPLER']
                del observation.backplanes['LON-CENTRIC']
                del observation.backplanes['LAT-CENTRIC']
                del observation.backplanes['RA']
                del observation.backplanes['DEC']
                del observation.backplanes['KM-X']
                del observation.backplanes['KM-Y']
                del observation.backplanes['RING-RADIUS']
                del observation.backplanes['RING-LON-GRAPHIC']
                del observation.backplanes['RING-DISTANCE']
                del observation.backplanes['LIMB-LON-GRAPHIC']
                del observation.backplanes['LIMB-DISTANCE']
                del observation.backplanes['RADIAL-VELOCITY']
                del observation.backplanes['LOCAL-SOLAR-TIME']
                del observation.backplanes['AZIMUTH']
                del observation.backplanes['PHASE']
                del observation.backplanes['LIMB-LAT-GRAPHIC']
                del observation.backplanes['ANGULAR-X']
                del observation.backplanes['ANGULAR-Y']
                del observation.backplanes['PIXEL-X']
                del observation.backplanes['PIXEL-Y']
                del observation.backplanes['DISTANCE']

                observation.set_disc_params(params[0],params[1],params[2],params[3])
                observation.save_observation(path+fn.replace(".png",".fits"))
                observation.save_mapped_observation(path+fn.replace(".png","map.fits"))

                
                #observation.set_plate_scale_arcsec(0.22282)
    
            
                # Run some custom setup
                #observation.add_other_bodies_of_interest('Io', 'Europa', 'Ganymede', 'Callisto')
                #observation.set_plate_scale_arcsec(42) # set a custom plate scale
                #observation.rotation_from_wcs() # get the disc rotation from the header's WCS info
            
                # Run the GUI to fit the observation interactively
                # This will open a GUI window every loop
                #coords = observation.run_gui()
            
                # More custom code can go here to use the fitted observation...
                # for example, we can print some values for the last click location
                #if coords:
                #    x, y = coords[-1]
                #    print(observation.xy2lonlat(x, y))
        
