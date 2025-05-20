# -*- coding: utf-8 -*-
"""
Created on Tue May 13 09:10:36 2025

@author: smhil
"""

def plot_contours_on_patch(ax,CH4Abs_conv,LatLims,LonLims,lvls=[0.71,0.73,0.75,0.77,0.79],frmt='%3.1e',clr='w'):
    """
    PURPOSE: Overlay countours of NH3 absorption data on Jovian maps.
             Specifically designed for equivalent widths with mean values of
             ~0.55nm
    """
    clrs = []
    for i in range(6):
        clrs.append(clr)
    cs=ax.contour(CH4Abs_conv,origin='upper', 
                  extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],90-LatLims[0]],
                  colors=clrs, alpha=0.8,levels=lvls,
                  #linewidths=[0.5,0.5,0.5,0.5,0.5,0.5],
                  linewidths=[1.0,1.0,1.0,1.0,1.0,1.0],
                  linestyles=['dashed','dashed','dashed','dashed','dashed'])
    #ax.clabel(cs,[19.0,19.5,20.0,20.5,21.0],inline=True,fmt='%2.1f',fontsize=8)
    #print(lvls)
    ax.clabel(cs,lvls,inline=True,fmt=frmt,fontsize=9)
    