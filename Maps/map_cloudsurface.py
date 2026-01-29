def surface_helper(fig3d,ax3d,LonSys,Lons,Lats,zmin,zmax,title):
    from matplotlib import cm, colors
    import pylab as pl
    import numpy as np

    ax3d.set_zlim(zmin,zmax)
    ax3d.invert_xaxis()
    ax3d.set_box_aspect((np.ptp(Lons), np.ptp(Lats), 30))
    ax3d.view_init(45, -60, 0) 
    ax3d.set_xlabel('Sys '+LonSys+' Longitude deg')
    ax3d.set_ylabel('PG Latitude (deg)')
    ax3d.set_title(title)

    if "on Pressure Surface" in title:
        ax3d.set_zlabel('Pressure (mb)')
        ax3d.invert_zaxis()
    elif "on Ammonia Surface" in title:
        ax3d.set_zlabel('fNH3 (ppm)')
        
    fig3d.subplots_adjust(left=0.05, bottom=0.03, right=0.98, top=0.95)     


def map_cloudsurface(PCld_patch,fNH3_patch_mb,RGB4Display,
                     HDRCld,HDRfNH3,timeRGB,
                     LonSys,LatLims,LonLims,LonRng,PlotCM,pathout):

    import sys
    drive='c:'
    sys.path.append(drive+'/Astronomy/Python Play')
    sys.path.append(drive+'/Astronomy/Python Play/Util_P3')
    sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy')
    sys.path.append('./Services')

    from matplotlib import cm, colors
    import pylab as pl
    import numpy as np

    fig3dCld, ax3dCld = pl.subplots(1,figsize=(12,8),subplot_kw={"projection": "3d"})
    fig3dboth, ax3dboth = pl.subplots(1,figsize=(12,8),subplot_kw={"projection": "3d"})
    fig3dRGBPCld, ax3dRGBPCld = pl.subplots(1,figsize=(12,8),subplot_kw={"projection": "3d"})
    fig3dRGBfNH3, ax3dRGBfNH3 = pl.subplots(1,figsize=(12,8),subplot_kw={"projection": "3d"})
    fig3dNH3, ax3dNH3 = pl.subplots(1,figsize=(12,8),subplot_kw={"projection": "3d"})
    
    # Make data.
    Lons = np.arange(360-LonLims[0],360-LonLims[1],-1.0)
    Lats = np.arange(90-LatLims[0],90-LatLims[1], -1.0)
    zmin=1000.
    zmax=2500.
    print(Lats)
    X, Y = np.meshgrid(Lons, Lats)
    Z = PCld_patch
    W = fNH3_patch_mb
    Img=RGB4Display
    
    fnskeleton='_Sys'+LonSys+'_N'+\
                str(90-LatLims[0])+'-S'+str(LatLims[1]-90)+\
                '_Lon'+str(np.mod(360-LonLims[1],360)).zfill(3)+'-'+\
                    str(np.mod(360-LonLims[0],360)).zfill(3)

    print("############# fnskel=",fnskeleton)
    ###########################################################################
    # Plot the Cloud surface.
    surfPCld = ax3dCld.plot_surface(X, Y, Z, cmap="Blues",linewidth=0, 
                                    antialiased=False,
                                    vmin=1400,vmax=2200)
    title='Cloud Pressure Plotted on Pressure Surface'
    if HDRCld:
        title=title+'\n'+HDRCld['DATE-OBS']
    surface_helper(fig3dCld,ax3dCld,LonSys,Lons,Lats,zmin,zmax,title)

    cbarCld=fig3dCld.colorbar(surfPCld, shrink=0.35, aspect=15, label='Pressure (mb)',
                       pad=0.2)
    cbarCld.set_ticks(ticks=[1400,1600,1800,2000,2200], labels=[1400,1600,1800,2000,2200])
    cbarCld.ax.invert_yaxis()
    if HDRCld:
        fig3dCld.savefig(pathout+HDRCld['FILENAME'][0:26]+'L3Surf_PonP'+fnskeleton+'.png',dpi=300)

    ###########################################################################
    # Create fNH3 facecolors as 4th data set
    norm = colors.Normalize(vmin=60, vmax=160)
    cmap = cm.get_cmap('terrain_r') 
    facecolors = cmap(norm(W))

    ###########################################################################
    # Plot fNH3 on Pressure surface
    surfboth = ax3dboth.plot_surface(X, Y, Z, facecolors=facecolors, 
                                     rstride=1, cstride=1, antialiased=False)
    title='fNH3 Plotted on Pressure Surface'
    if HDRCld and HDRfNH3:
        title='fNH3 Plotted on Pressure Surface'
        title=title+'\nPressure:'+HDRCld['DATE-OBS']+' | fNH3: '+HDRfNH3['DATE-OBS']
    surface_helper(fig3dboth,ax3dboth,LonSys,Lons,Lats,zmin,zmax,title)
    cbarboth=fig3dboth.colorbar(surfboth, shrink=0.35, aspect=15, label='fNH3 (ppm)',
                       pad=0.2)
    cbarboth.set_ticks(ticks=[0,0.25,0.5,0.75,1.0], labels=[160,135,110,85,60])
    cbarboth.ax.invert_yaxis()
    #!! Need to fix this colorbar
    if HDRCld:
        fig3dboth.savefig(pathout+HDRCld['FILENAME'][0:26]+'L3Surf_NonP'+fnskeleton+'.png',dpi=300)

    ###########################################################################
    # Plot Context Image (RGB) on Pressure surface
    surf = ax3dRGBPCld.plot_surface(X, Y, Z, facecolors=Img,
                                    rstride=1, cstride=1, antialiased=False)
    title='Context Image Plotted on Pressure Surface'
    if HDRCld and timeRGB:
        title=title+'\nPressure:'+HDRCld['DATE-OBS']+' | RGB: '+timeRGB+'Z'
    surface_helper(fig3dRGBPCld,ax3dRGBPCld,LonSys,Lons,Lats,zmin,zmax,title)
    cbarRGBPCld=fig3dCld.colorbar(surf, shrink=0.35, aspect=15, label='Pressure (mb)',
                       pad=0.2)
    cbarRGBPCld.ax.set_visible(False)
    if HDRCld:
        fig3dRGBPCld.savefig(pathout+HDRCld['FILENAME'][0:26]+'L3Surf_ConP'+fnskeleton+'.png',dpi=300)

    ###########################################################################
    # Plot Context Image on fNH3 surface
    surf = ax3dRGBfNH3.plot_surface(X, Y, W, facecolors=Img, 
                                    rstride=1, cstride=1, antialiased=False)
    title='Context Image Plotted on Ammonia Surface'
    if HDRfNH3 and timeRGB:
        title=title+'\nfNH3:'+HDRfNH3['DATE-OBS']+' | RGB: '+timeRGB+'Z'
    surface_helper(fig3dRGBfNH3,ax3dRGBfNH3,LonSys,Lons,Lats,0,200,title)
    cbarRGBfNH3=fig3dCld.colorbar(surf, shrink=0.35, aspect=15, label='Pressure (mb)',
                       pad=0.2)
    cbarRGBfNH3.ax.set_visible(False)
    if HDRCld:
        fig3dRGBfNH3.savefig(pathout+HDRCld['FILENAME'][0:26]+'L3Surf_ConN'+fnskeleton+'.png',dpi=300)

    ###########################################################################
    # Plot the NH3 surface.
    surffNH3 = ax3dNH3.plot_surface(X, Y, W, cmap="terrain_r",
                           linewidth=0, antialiased=False,vmin=60,vmax=160)        
    title='fNH3 Plotted on Ammonia Surface'
    if HDRfNH3:
        title=title+'\n'+HDRfNH3['DATE-OBS']
    surface_helper(fig3dNH3,ax3dNH3,LonSys,Lons,Lats,0,200,title)
    cbarNH3=fig3dboth.colorbar(surffNH3, shrink=0.35, aspect=15, label='fNH3 (ppm)',
                       pad=0.2)
    if HDRCld:
        fig3dNH3.savefig(pathout+HDRCld['FILENAME'][0:26]+'L3Surf_NonN'+fnskeleton+'.png',dpi=300)
