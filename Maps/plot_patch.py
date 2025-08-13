def plot_patch(patch,LatLims,LonLims,CM2,LonRng,colorscale,axis,
               cbarplot=True,cbar_title="Test",cbar_reverse=False,vn=0.10,vx=0.20,n=6):
    """
    Purpose:
        To plot a map patch with appropriate latitude and longitude scales,
        and, optionally, a color bar

    Parameters
    ----------
    patch : TYPE
        DESCRIPTION.
    LatLims : TYPE
        DESCRIPTION.
    LonLims : TYPE
        DESCRIPTION.
    CM2 : TYPE
        DESCRIPTION.
    LonRng : TYPE
        DESCRIPTION.
    colorscale : TYPE
        DESCRIPTION.
    axis : TYPE
        DESCRIPTION.
    cbarplot : TYPE, optional
        DESCRIPTION. The default is True.
    cbar_title : TYPE, optional
        DESCRIPTION. The default is "Test".
    cbar_reverse : TYPE, optional
        DESCRIPTION. The default is False.
    vn : TYPE, optional
        DESCRIPTION. The default is 0.10.
    vx : TYPE, optional
        DESCRIPTION. The default is 0.20.
    n : TYPE, optional
        DESCRIPTION. The default is 6.

    Returns
    -------
    None.

    """
    import numpy as np
    import pylab as pl

    np.nan_to_num(patch, copy=False, nan=0.0, posinf=0.0, neginf=0.0)
    tx=np.linspace(vn,vx,n,endpoint=True)

    show=axis.imshow(patch, colorscale, origin='upper',vmin=vn,vmax=vx,  
               extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],
                       90-LatLims[0]],
                       aspect="equal")

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

    return patch,vn,vx,tx



