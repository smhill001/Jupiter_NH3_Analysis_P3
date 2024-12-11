# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 11:48:42 2024

@author: smhil
"""

import L3_Jup_Map_Plot as L3M
import L2_Jup_Map_Plot as L2M

L2M.L2_Jup_Map_Plot(obskey="20220904UTa",imagetype='Map',target="Jupiter",
                    Smoothing=False,LatLims=[45,135],LonRng=45,CMpref='subobs',
                    LonSys='2',showbands=False,coef=[0.9,0.45],subproj='',
                    figxy=[8.0,4.0],FiveMicron='fits')

#Figure 1

L3M.L3_Jup_Map_Plot(obskey="20220904UTa",imagetype='Map',target="Jupiter",
                Smoothing=False,LatLims=[45,135],LonRng=45,CMpref='subobs',
                LonSys='2',showbands=False,coef=[0.75,0.45],subproj='',
                figxy=[8.0,4.0],FiveMicron='fits')


###############################################################################
#GRS Figures
###############################################################################

#7/30/2022

#L3M.L3_Jup_Map_Plot(obskey="20220730UTa",imagetype='Map',target="Jupiter",
#                    Smoothing=False,LatLims=[90,130],LonRng=20,CMpref=20,
#                    LonSys='2',showbands=False,coef=[0.75,0.45],subproj='',
#                    figxy=[8.0,4.0],FiveMicron='fits')

#8/18/2022
L3M.L3_Jup_Map_Plot(obskey="20220818UTa",imagetype='Map',target="Jupiter",
                    Smoothing=False,LatLims=[90,130],LonRng=20,CMpref=25,
                    LonSys='2',showbands=False,coef=[0.75,0.45],subproj='',
                    figxy=[8.0,4.0],FiveMicron='fits')

#9/04/2022
L3M.L3_Jup_Map_Plot(obskey="20220904UTa",imagetype='Map',target="Jupiter",
                    Smoothing=False,LatLims=[90,130],LonRng=20,CMpref=25,
                    LonSys='2',showbands=False,coef=[0.75,0.45],subproj='',
                    figxy=[8.0,4.0],FiveMicron='fits')

#10/17/2023
L3M.L3_Jup_Map_Plot(obskey="20231017UTa",imagetype='Map',target="Jupiter",
                    Smoothing=False,LatLims=[90,130],LonRng=20,CMpref=45,
                    LonSys='2',showbands=False,coef=[0.75,0.45],subproj='',
                    figxy=[8.0,4.0],FiveMicron='fits')

###############################################################################
#NEDFs
###############################################################################

#8/18/2022
L3M.L3_Jup_Map_Plot(obskey="20220818UTa",imagetype='Map',target="Jupiter",
                    Smoothing=False,LatLims=[75,95],LonRng=25,CMpref=90,
                    LonSys='1',showbands=False,coef=[0.75,0.45],subproj='',
                    figxy=[8.0,2.5],FiveMicron='fits')

#9/04/2022
L3M.L3_Jup_Map_Plot(obskey="20220904UTa",imagetype='Map',target="Jupiter",
                    Smoothing=False,LatLims=[75,95],LonRng=25,CMpref=230,
                    LonSys='1',showbands=False,coef=[0.75,0.45],subproj='',
                    figxy=[8.0,2.5],FiveMicron='fits')

#10/17/2022
L3M.L3_Jup_Map_Plot(obskey="20231017UTa",imagetype='Map',target="Jupiter",
                    Smoothing=False,LatLims=[75,95],LonRng=25,CMpref=120,
                    LonSys='1',showbands=False,coef=[0.75,0.45],subproj='',
                    figxy=[8.0,2.5],FiveMicron='fits')

###############################################################################
#Other Features
###############################################################################

#7/30/2022
L3M.L3_Jup_Map_Plot(obskey="20220730UTa",imagetype='Map',target="Jupiter",
                    Smoothing=False,LatLims=[45,65],LonRng=25,CMpref=50,
                    LonSys='2',showbands=False,coef=[0.75,0.45],subproj='',
                    figxy=[8.0,2.5],FiveMicron='fits')

#7/30/2022
L3M.L3_Jup_Map_Plot(obskey="20220730UTa",imagetype='Map',target="Jupiter",
                    Smoothing=False,LatLims=[115,135],LonRng=25,CMpref=60,
                    LonSys='2',showbands=False,coef=[0.75,0.45],subproj='',
                    figxy=[8.0,2.5],FiveMicron='fits')

#9/04/2022
L3M.L3_Jup_Map_Plot(obskey="20220904UTa",imagetype='Map',target="Jupiter",
                    Smoothing=False,LatLims=[115,135],LonRng=25,CMpref=45,
                    LonSys='2',showbands=False,coef=[0.75,0.45],subproj='',
                    figxy=[8.0,2.5],FiveMicron='fits')
