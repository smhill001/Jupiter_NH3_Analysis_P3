def get_L0A_png_list(obskey,planet="Saturn"):
    import sys
    drive='c:'
    sys.path.append(drive+'/Astronomy/Python Play')
    sys.path.append(drive+'/Astronomy/Python Play/Util_P3')
    sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy')
    sys.path.append(drive+'/Astronomy/Python Play/SPLibraries_P3')
    import os
    import ConfigFiles as CF

    path="C:/Astronomy/Projects/Planets/"+planet+"/Imaging Data/"+obskey[:-1]+"/"
    
    dirlist=os.listdir(path)
    pngs=[i for i in dirlist if '.png' in i]
    stacked=[i for i in pngs if 'FlatStack' in i]
    aligned=[i for i in stacked if 'Aligned' in i]
    stretched=[i for i in aligned if 'Str0' in i]
    #for string in stretched:
    #    stacked.remove(string)
    WV=[i for i in aligned if 'WV' in i]
    for string in WV:
        stacked.remove(string)
    
    
    print("###################")
    #print(pngs)
    print("################### NO WAVELETS")
    #for i in stacked:
    #    print(i)
    print("################### WAVELETS")
    #for i in WV:
        #print(i)
    
    
    RGB=["685NIR","550GRN","450BLU"]
    science=["656HIA","632OI","620CH4","647NH3"]
    fn_bases=CF.video_metadata_list(path+obskey+".csv") #read the metadata
    fn_bases.load_records()
    #print(fn_bases.VideoFile)
    
    fn_L0A_temp=[]
    #for k in range(0,fn_bases.nrecords-1):
    #    fn_base=fn_bases.VideoFile[k].split('.')[0]
    #    print(fn_base)
    for filtr in RGB:
        #print(filtr)
        #print(WV)
        ftemp=[i for i in WV if filtr in i]
        fn_L0A_temp.extend(ftemp)
    for filtr in science:
        #print(filtr)
        #print(aligned)
        ftemp=[i for i in aligned if filtr in i]
        fn_L0A_temp.extend(ftemp)
    #print("################### L0A list")
    #for i in fn_L0A_temp:
        #print(i)
    #print("###################")
    fn_L0A_list=[]
    #print(fn_bases.VideoFile)
    for fnbase in fn_bases.VideoFile:
        ftemp=[i for i in fn_L0A_temp if fnbase.split('.')[0] in i]
        fn_L0A_list.extend(ftemp)

    #!Need to add one further check that the selected files are part of the
    #! specific observation (obskey) by comparing with that list. This is
    #! crucial when there is more than one observation contained in a folder.
    
    #print(path,fn_L0A_list)
    return(path,fn_L0A_list)