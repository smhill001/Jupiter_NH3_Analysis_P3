# -*- coding: utf-8 -*-
"""
Created on Thu Nov 20 17:43:30 2025

@author: smhil
"""

import L4_Jup_Map_Plot as L4MP
"""
#20220730
L4MP.L4_Jup_Map_Plot(lats=[0,180],LonLims=[0,360],LonSys='3',
                     collection="20220730-20220730",
                     IRTFcollection="20220725-20220726",
                     CH4889collection="20220802-20220803",
                     cont=False)

#20220818
L4MP.L4_Jup_Map_Plot(lats=[0,180],LonLims=[0,360],LonSys='3',
                     collection="20220818-20220818",
                     IRTFcollection="20220817-20220818",
                     CH4889collection="20220816-20220816",
                     cont=False)

#20220925
L4MP.L4_Jup_Map_Plot(lats=[0,180],LonLims=[0,360],
                     LonSys='3',collection="20220925-20220925",
                     IRTFcollection="20220929-20220929",
                     CH4889collection="20220929-20220929",
                     cont=False)

#20231015-19
L4MP.L4_Jup_Map_Plot(lats=[0,180],LonLims=[0,360],LonSys='3',
                     collection="20231015-20231019",
                     IRTFcollection="20231014-20231015",
                     CH4889collection="20231017-20231017",
                     cont=False)

#20240129-0202
L4MP.L4_Jup_Map_Plot(lats=[0,180],LonLims=[0,360],LonSys='3',
                     collection="20240129-20240202",
                     IRTFcollection="20240205-20240205",
                     CH4889collection="20240131-20240201",
                     cont=False)

#20240919
L4MP.L4_Jup_Map_Plot(lats=[0,180],LonLims=[0,360],LonSys='3',
                     collection="20240919-20240919",  
                     IRTFcollection="20240920-20240920",
                     CH4889collection="NA",cont=False)

#20240925-29
L4MP.L4_Jup_Map_Plot(lats=[0,180],LonLims=[0,360],LonSys='3',
                     collection="20240925-20240929",  
                     IRTFcollection="20240920-20240920",
                     CH4889collection="NA",cont=False)

#20241022-23
L4MP.L4_Jup_Map_Plot(lats=[0,180],LonLims=[0,360],LonSys='3',
                     collection="20241022-20241023",
                     IRTFcollection="20241022-20241022",
                     CH4889collection="NA",cont=False)

#20250308
L4MP.L4_Jup_Map_Plot(lats=[0,180],LonLims=[0,360],LonSys='3',
                     collection="20250309-20250309",
                     IRTFcollection="20250308-20250308",
                     CH4889collection="20250310-20250310",cont=False)

#20250919
L4MP.L4_Jup_Map_Plot(lats=[0,180],LonLims=[0,360],LonSys='3',
                     collection="20250919-20250919",
                     IRTFcollection="20250913-20250913",
                     CH4889collection="20250919-20250920",cont=False)
"""

import matplotlib.pyplot as pl
import numpy as np

path='C:/Astronomy/Projects/SAS 2021 Ammonia/HST GO 18055/'
fNH3file="Picture1a.png"
PCldfile="Picture2a.png"
RGBfile="Picture3a.png"
fNH3=pl.imread(path+fNH3file)
PCld=pl.imread(path+PCldfile)
RGB=pl.imread(path+RGBfile)

fNH3a=np.mean(fNH3[:,:,0:2],axis=2)
PClda=np.mean(PCld[:,:,0:2],axis=2)

figfNH3,axsfNH3=pl.subplots(1,figsize=(6,4), dpi=150, facecolor="white")
figPCld,axsPCld=pl.subplots(1,figsize=(6,4), dpi=150, facecolor="white")
figRGB,axsRGB=pl.subplots(1,figsize=(6,4), dpi=150, facecolor="white")

axsfNH3.imshow(fNH3a,cmap='terrain_r')
axsfNH3.axis('off')
axsPCld.imshow(PClda,cmap='Blues_r')
axsPCld.axis('off')
axsRGB.imshow(RGB)
axsRGB.axis('off')


"""
#20251016
L4MP.L4_Jup_Map_Plot(lats=[75,105],LonLims=[0,60],LonSys='3',
                     collection="20251016-20251017",
                     IRTFcollection=False,
                     CH4889collection=False,cont=False)
"""