def ZonalHSTWinds(HSTkey1,HSTkey2):
    
    import os
    import sys
    sys.path.append('./Hubble')
    sys.path.append('././Hubble')

    from Winds import estimate_zonal_winds
    from read_HST import getFilename,read_HST
    from matplotlib import pyplot as pl
    from datetime import datetime
    import numpy as np
    
    pathHST='C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Hubble/Data/'
    filenames=os.listdir(pathHST)
    
    filename1=getFilename(filenames,HSTkey1)
    mapHST1,dateobs1=read_HST(obskey=filename1,LonSys='3')
    filename2=getFilename(filenames,HSTkey2)
    mapHST2,dateobs2=read_HST(obskey=filename2,LonSys='3')

    t1 = datetime.strptime(dateobs1, "%Y-%m-%dT%H:%M:%S")
    t2 = datetime.strptime(dateobs2, "%Y-%m-%dT%H:%M:%S")
    dt=(t2-t1).total_seconds()

    lats,u=estimate_zonal_winds(mapHST1,mapHST2,dt,deg_per_px=0.1,ref_latitude_deg=0.0)
    
    kernel_size = 10
    kernel = np.ones(kernel_size) / kernel_size
    u_convolved = np.convolve(u, kernel, mode='same')
    fig,axs=pl.subplots(1,figsize=(6,6), dpi=150, facecolor="white")
    #axs.plot(lats,u)
    axs.plot(u_convolved,lats,color='k',label="Global (HST)")
    axs.set_title("Zonal Wind Profile")
    axs.set_ylabel("PG Latitude (deg)",fontsize=10)
    axs.set_xlabel("Zonal Velocity (m/s)",fontsize=10)
    axs.tick_params('x', labelsize=8)
    axs.set_xlim(0,150)
    axs.set_ylim(0,15)
    axs.grid(linewidth=0.2)
    
    return fig,axs

def overplot_NEZ_feature_speeds(axs):
    import numpy as np
    pth="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Studies/NEZ/"
    config={"NH3":{"file":"NH3speed_by_label.csv","color":'C2'},
             "Plume":{"file":"Plumespeed_by_label.csv","color":'C0'},   
             "NEDF":{"file":"NEDFspeed_by_label.csv","color":'r'}}


    for item in config:
        print(item,config[item]["file"])
        infile=config[item]["file"]
        data_array = np.genfromtxt(pth+infile, delimiter=',',dtype=None)
        y=np.array(data_array[1:,:],dtype=float)
        print(data_array[1:,3])
        axs.scatter(y[:,3],y[:,1],color=config[item]["color"],label=item)
        
        
def NEZ_feature_quiver_plot():
    import numpy as np
    import matplotlib.pyplot as pl
    pth="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Studies/NEZ/"
    config={"NH3":{"file":"NH3speed_by_label.csv","color":'C2'},
             "Plume":{"file":"Plumespeed_by_label.csv","color":'C0'},   
             "NEDF":{"file":"NEDFspeed_by_label.csv","color":'r'}}

    NH3_array = np.genfromtxt(pth+config["NH3"]["file"], delimiter=',',dtype=str)
    NH3=np.array(NH3_array[1:,:],dtype=float)

    Plume_array = np.genfromtxt(pth+config["Plume"]["file"], delimiter=',',dtype=str)
    Plume=np.array(Plume_array[1:,:],dtype=float)

    NEDF_array = np.genfromtxt(pth+config["NEDF"]["file"], delimiter=',',dtype=str)
    NEDF=np.array(NEDF_array[1:,:],dtype=float)

    avglon=[]
    labels=[]
    NH3_indices=[]
    Plume_indices=[]
    NEDF_indices=[]
    fig,axs=pl.subplots(1,figsize=(8,6), dpi=150, facecolor="white")
    fig1,axs1=pl.subplots(1,figsize=(10,5), dpi=150, facecolor="white")
    First=True
    for k in range(1,10):
        flabel=float(k)
        if (flabel in NH3[:,0]) and (flabel in Plume[:,0]) and (flabel in NEDF[:,0]):
            NH3idx=np.where(NH3[:,0]==flabel)
            Plumeidx=np.where(Plume[:,0]==flabel)
            NEDFidx=np.where(NEDF[:,0]==flabel)
            
            avglon=np.mean([NH3[NH3idx,2],Plume[Plumeidx,2],NEDF[NEDFidx,2]])
            print("Yay",k,NEDF[NEDFidx,2],avglon)
            #avglon.append(np.mean([NH3[i,2],Plume[i,2],NEDF[i,2]]))
            labels.append(k)  
            NH3_indices.append(NH3idx)  
            Plume_indices.append(Plumeidx)  
            NEDF_indices.append(NEDFidx)  
            
            if First:
                axs.quiver([NH3[NH3idx[0][0],2]]-avglon,[NH3[NH3idx[0][0],1]],
                         [NH3[NH3idx[0][0],3]],[0.],pivot='mid',color='C2',label='NH3')
                axs.quiver([Plume[Plumeidx[0][0],2]]-avglon,[Plume[Plumeidx[0][0],1]],
                         [Plume[Plumeidx[0][0],3]],[0.],pivot='mid',color='C0',label="Plume")
                axs.quiver([NEDF[NEDFidx[0][0],2]]-avglon,[NEDF[NEDFidx[0][0],1]],
                         [NEDF[NEDFidx[0][0],3]],[0.],pivot='mid',color='r',label="NEDF")
                First=False
            else:
                axs.quiver([NH3[NH3idx[0][0],2]]-avglon,[NH3[NH3idx[0][0],1]],
                         [NH3[NH3idx[0][0],3]],[0.],pivot='mid',color='C2')
                axs.quiver([Plume[Plumeidx[0][0],2]]-avglon,[Plume[Plumeidx[0][0],1]],
                         [Plume[Plumeidx[0][0],3]],[0.],pivot='mid',color='C0')
                axs.quiver([NEDF[NEDFidx[0][0],2]]-avglon,[NEDF[NEDFidx[0][0],1]],
                         [NEDF[NEDFidx[0][0],3]],[0.],pivot='mid',color='r')

            
            axs.annotate(str(k),(NH3[NH3idx[0][0],2]-avglon,NH3[NH3idx[0][0],1]),
                          xytext=(NH3[NH3idx[0][0],2]-avglon,NH3[NH3idx[0][0],1]+0.5),
                     color='C2')
            axs.annotate(str(k),(Plume[Plumeidx[0][0],2]-avglon,Plume[Plumeidx[0][0],1]),
                          xytext=(Plume[Plumeidx[0][0],2]-avglon,Plume[Plumeidx[0][0],1]+0.5),
                     color='C0')
            axs.annotate(str(k),(NEDF[NEDFidx[0][0],2]-avglon,NEDF[NEDFidx[0][0],1]),
                          xytext=(NEDF[NEDFidx[0][0],2]-avglon,NEDF[NEDFidx[0][0],1]+0.5),
                     color='r')
            axs.legend()

            axs1.quiver([NH3[NH3idx[0][0],2]],[NH3[NH3idx[0][0],1]],
                     [NH3[NH3idx[0][0],3]],[0.],pivot='mid',color='C2')
            axs1.quiver([Plume[Plumeidx[0][0],2]],[Plume[Plumeidx[0][0],1]],
                     [Plume[Plumeidx[0][0],3]],[0.],pivot='mid',color='C0')
            axs1.quiver([NEDF[NEDFidx[0][0],2]],[NEDF[NEDFidx[0][0],1]],
                     [NEDF[NEDFidx[0][0],3]],[0.],pivot='mid',color='r')
            
            axs1.annotate(str(k),(NH3[NH3idx[0][0],2],NH3[NH3idx[0][0],1]),
                          xytext=(NH3[NH3idx[0][0],2],NH3[NH3idx[0][0],1]+0.5),
                     color='C2')
            axs1.annotate(str(k),(Plume[Plumeidx[0][0],2],Plume[Plumeidx[0][0],1]),
                          xytext=(Plume[Plumeidx[0][0],2],Plume[Plumeidx[0][0],1]+0.5),
                     color='C0')
            axs1.annotate(str(k),(NEDF[NEDFidx[0][0],2],NEDF[NEDFidx[0][0],1]),
                          xytext=(NEDF[NEDFidx[0][0],2],NEDF[NEDFidx[0][0],1]+0.5),
                     color='r')
            print("Done with IF")
            
        axs.set_title("Speed and Location")
        axs.set_xlabel("Longitude (deg)",fontsize=10)
        axs.set_ylabel("PG Latitude (deg)",fontsize=10)
        axs.set_xlim(-20,20)
        axs.set_ylim(0,15)
        axs.grid(linewidth=0.2)
        axs.invert_xaxis()
        
        axs1.set_title("Speed and Location")
        axs1.set_xlabel("Longitude (deg)",fontsize=10)
        axs1.set_ylabel("PG Latitude (deg)",fontsize=10)
        axs1.set_xlim(15,195)
        axs1.set_ylim(0,15)
        axs1.grid(linewidth=0.2)
        axs1.invert_xaxis()
        
    pathmapplots="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Studies/NEZ/"  

    fig.savefig(pathmapplots+"Speed and Location.png",dpi=300)

    return labels,NH3_indices,Plume_indices,NEDF_indices
    
fig,axs=ZonalHSTWinds("2024c_f631n","2024d_f631n")
overplot_NEZ_feature_speeds(axs)
axs.legend()
pathmapplots="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Studies/NEZ/"  

fig.savefig(pathmapplots+"WindProfile.png",dpi=300)
