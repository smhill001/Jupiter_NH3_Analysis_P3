"""
Created on Sun Aug 25 15:04:22 2024

@author: smhil
"""
import pylab as pl
import MakeContiguousMap as MCM

LonSys='2'
lats=[95,125]
maps2022=["20220810-20220812","20220828-20220901","20220904-20220905",
          "20220912-20220913","20220919-20220919","20221009-20221013",
          "20221019-20221021"]

maps2023=["20230815-20230818","20230827-20230831",
          "20230905-20230906","20230922-20230929","20231005-20231006",
          "20231015-20231019","20231022-20231026","20231103-20231107",
          "20231110-20231110","20231112-20231113","20231115-20231115",
          "20231128-20231129","20231206-20231207","20231217-20231218",
          "20231229-20231229","20240129-20240202","20240229-20240301"]

maps2024=["20241006-20241010","20241027-20241027","20241105-20241105","20241115-20241115",
          "20241118-20241118","20241128-20241129","20241202-20241203",
          "20241205-20241205","20250129-20250129"]

###############################################################################
# 2022
###############################################################################
fig22NH3,axs22NH3=pl.subplots(7,1,figsize=(6.0,6.0), dpi=150, facecolor="white",
                      sharex=True,sharey=True)   
fig22CH4,axs22CH4=pl.subplots(7,1,figsize=(6.0,6.0), dpi=150, facecolor="white",
                      sharex=True,sharey=True)   
fig22RGB,axs22RGB=pl.subplots(7,1,figsize=(6.0,6.0), dpi=150, facecolor="white",
                      sharex=True,sharey=True)   
counter=0
for mp in maps2022:
    MCM.MakeContiguousMap(axs22NH3[counter],axs22CH4[counter],axs22RGB[counter],
                          collection=mp,LonSys=LonSys,lats=lats)
    #pl.show()
    counter=counter+1

###############################################################################
# 2023
###############################################################################   
fig23NH3,axs23NH3=pl.subplots(17,1,figsize=(6.0,6.0), dpi=150, facecolor="white",
                      sharex=True,sharey=True)   
fig23NH3.subplots_adjust(left=0.10, bottom=0.08, right=0.98, top=0.94,
            wspace=0.25, hspace=0.08)     
fig23NH3.suptitle("Ammonia Abundance (ppm)")
axs23NH3[16].set_xlabel("System "+LonSys+" Longitude (deg)",fontsize=8)



fig23CH4,axs23CH4=pl.subplots(17,1,figsize=(6.0,6.0), dpi=150, facecolor="white",
                      sharex=True,sharey=True)   
fig23CH4.subplots_adjust(left=0.10, bottom=0.08, right=0.98, top=0.94,
            wspace=0.25, hspace=0.08)     
fig23CH4.suptitle("Effective Cloud-Top Pressure (mb)")
axs23CH4[16].set_xlabel("System "+LonSys+" Longitude (deg)",fontsize=8)


fig23RGB,axs23RGB=pl.subplots(17,1,figsize=(6.0,6.0), dpi=150, facecolor="white",
                      sharex=True,sharey=True)   
fig23RGB.subplots_adjust(left=0.10, bottom=0.08, right=0.98, top=0.94,
            wspace=0.25, hspace=0.08)     
fig23RGB.suptitle("Visual Context")
axs23RGB[16].set_xlabel("System "+LonSys+" Longitude (deg)",fontsize=8)


counter=0
for mp in maps2023:
    MCM.MakeContiguousMap(axs23NH3[counter],axs23CH4[counter],axs23RGB[counter],
                                                               collection=mp,LonSys=LonSys,lats=lats)
    #pl.show()
    counter=counter+1

###############################################################################
# 2024
###############################################################################
fig24NH3,axs24NH3=pl.subplots(9,1,figsize=(6.0,6.0), dpi=150, facecolor="white",
                      sharex=True,sharey=True)   
fig24CH4,axs24CH4=pl.subplots(9,1,figsize=(6.0,6.0), dpi=150, facecolor="white",
                      sharex=True,sharey=True)   
fig24RGB,axs24RGB=pl.subplots(9,1,figsize=(6.0,6.0), dpi=150, facecolor="white",
                      sharex=True,sharey=True)   
counter=0
for mp in maps2024:
    MCM.MakeContiguousMap(axs24NH3[counter],axs24CH4[counter],axs24RGB[counter],
                          collection=mp,LonSys=LonSys,lats=lats)
    #pl.show()
    counter=counter+1


if int(lats[0])<90:
    latstr=str(90-lats[0])+"N"
if int(lats[0])==90:
    latstr=str(90-lats[0])
if int(lats[0])>90:
    latstr=str(lats[0]-90)+"S"
    
if int(lats[1])<90:
    latstr=latstr+"-"+str(90-lats[1])+"N"
if int(lats[1])==90:
    latstr=latstr+"-"+str(90-lats[1])
if int(lats[1])>90:
    latstr=latstr+"-"+str(lats[1]-90)+"S"
    
pathmapplots="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Studies/maps/"
fig22NH3.savefig(pathmapplots+"2022 NH3 Stack Sys"+LonSys+" "+latstr+"_map.png",dpi=300)
fig22CH4.savefig(pathmapplots+"2022 CH4 Stack Sys"+LonSys+" "+latstr+"_map.png",dpi=300)
fig22RGB.savefig(pathmapplots+"2022 RGB Stack Sys"+LonSys+" "+latstr+"_map.png",dpi=300)

fig23NH3.savefig(pathmapplots+"2023 NH3 Stack Sys"+LonSys+" "+latstr+"_map.png",dpi=300)
fig23CH4.savefig(pathmapplots+"2023 CH4 Stack Sys"+LonSys+" "+latstr+"_map.png",dpi=300)
fig23RGB.savefig(pathmapplots+"2023 RGB Stack Sys"+LonSys+" "+latstr+"_map.png",dpi=300)

fig24NH3.savefig(pathmapplots+"2024 NH3 Stack Sys"+LonSys+" "+latstr+"_map.png",dpi=300)
fig24CH4.savefig(pathmapplots+"2024 CH4 Stack Sys"+LonSys+" "+latstr+"_map.png",dpi=300)
fig24RGB.savefig(pathmapplots+"2024 RGB Stack Sys"+LonSys+" "+latstr+"_map.png",dpi=300)
