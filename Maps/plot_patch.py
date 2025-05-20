def plot_patch(patch,LatLims,LonLims,CM2,LonRng,colorscale,axis,frmt,
               cont=True,cbarplot=True,cbar_title="Test",cbar_reverse=False,vn=0.10,vx=0.20,n=6):
    """
    Created on Thu Apr 11 18:54:35 2024
    
    @author: smhil
    """    
    import numpy as np
    import pylab as pl
    #import make_patch as MP
    import plot_contours_on_patch as PC

    print("@@@@@@@@@@@@ LatLims, LonLims, CM2, LonRng",LatLims, LonLims, CM2, LonRng)
    #patch=MP.make_patch(fullmap,LatLims,LonLims,CM2,LonRng)
    np.nan_to_num(patch, copy=False, nan=0.0, posinf=0.0, neginf=0.0)
    tx=np.linspace(vn,vx,n,endpoint=True)

    show=axis.imshow(patch, colorscale, origin='upper',vmin=vn,vmax=vx,  
               extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],
                       90-LatLims[0]],
                       aspect="equal")
    #if cont:
    #    temp=PC.plot_contours_on_patch(axis,patch,LatLims,LonLims,
    #                       lvls=tx,frmt=frmt,clr='k')

    im_ratio = patch.shape[0]/patch.shape[1]
    if cbarplot:
        cbar = pl.colorbar(show, ticks=tx, 
                   orientation='vertical',cmap='gist_heat',
                   ax=axis,fraction=0.046*im_ratio, pad=0.04)
        cbar.ax.set_yticklabels(np.around(tx,3))
        cbar.ax.tick_params(labelsize=6,color="k")#if iSession >1:
        cbar.ax.set_ylabel(cbar_title,size=8)#,loc="top")
        cbar.ax.yaxis.set_label_coords(-0.7, 0.5)
        if cbar_reverse:
            cbar.ax.invert_yaxis()
    #if colorscale=="Greys":
    #    cbar.set_label('Cloud Top Pressure (mb)',fontsize=7)

    return patch,vn,vx,tx



