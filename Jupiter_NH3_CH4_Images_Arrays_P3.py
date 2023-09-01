# -*- coding: utf-8 -*-
"""
Created on Wed Oct 06 14:36:01 2021

@author: Steven Hill

PURPOSE:    Read a white balanced RGB (656,647,632) image and compute the
            color slope and NH3 absorption images.

EXAMPLE:    Jupiter_NH3_CH4_Images_Arrays_P3(obsdate="20210905UT")
            
"""

def Jupiter_NH3_CH4_Images_Arrays_P3(obsdate="20220919UTb",target="Jupiter",
                                     imagetype='Map'):
    import sys
    drive='c:'
    sys.path.append(drive+'/Astronomy/Python Play')
    sys.path.append(drive+'/Astronomy/Python Play/Util_P3')
    sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Spectroscopy')

    import os
    from matplotlib.pyplot import imread
    import pylab as pl
    import numpy as np
    from imageio import imwrite
    from numpy import inf
    from re import search
    from astropy.io import fits
 
    ###########################################################################
    #  DATA FILES AND METADATA DICTIONARY
    #    !!!!SHOULD MAKE THIS A DATA OBJECT!!!!
    ###########################################################################
    sourcedata=obsdate+"_"+imagetype
    sourcefiles={'20210622UT':{'NH3file':'2021-06-22-1050_1-Jupiter_NoWV-R656G647B656-RGB-WhtBal.png',
                               'CH4file':['2021-06-22-1054_6-Jupiter_NoWV-R656(G)B889-R(G)B-WhtBal.png'],
                               'Context':{'NUVfile':'2021-06-22-1035_3-Jupiter_380NUV-Derot-Stack.png',
                                          'RGBfile':'2021-06-22-1043_2-Jupiter-(685)GB-RGB-WhtBal-Wavelets.png'}},
                 
                 '20210708UT':{'NH3file':'2021-07-08-1044_4-Jupiter-NoWV-DR-ST-R656G647B632-WhtBal.png',
                               'CH4file':['2021-07-08-1049_5-Jupiter_NoWV-R656(G)B889-R(G)B-WhtBal.png'],
                               'Context':{'NUVfile':'2021-07-08-1041_5-Jupiter_380NUV_-Derot-Stack-Wavelets4x10+5x20.png',
                                          'RGBfile':'2021-07-08-1037_5-Jupiter-AllREDGB-RGB-WhtBal-ClrSmth-Smth-Wavelets.png'}},
                 
                 '20210719UT':{'NH3file':'2021-07-19-1053_0-Jupiter_NoWV-R656G647B630-RGB-WhtBal.png',
                               'CH4file':['2021-07-19-1055_4-Jupiter_NoWV-R656(G)B889-R(G)B-WhtBal.png'],
                               'Context':{'NUVfile':'2021-07-19-1043_7-Jupiter_380NUV-Derot-Stack.png',
                                          'RGBfile':'2021-07-19-1056_7-Jupiter-Wavelets-RGB-WhtBal-ClrSmth-Smth-Wavelets.png'}},
                 
                 '20210720UT':{'NH3file':'2021-07-20-1108_1-Jupiter-DR-ST-R656G647B632-WhtBal.png',
                               'CH4file':['2021-07-20-1110_2-Jupiter_NoWV-R656(G)B889-R(G)B-WhtBal.png'],
                               'Context':{'NUVfile':'2021-07-20-1056_7-Jupiter_380NUV-Derot-Stack.png',
                                          'RGBfile':'2021-07-20-1053_5-Jupiter-Wavelets-RGB-WhtBal-ClrSmth-Smth-Wavelets.png'}},
                 
                 '20210905UT':{'NH3file':'2021-09-05-0409_0-Jupiter-R656G647B632-ReAligned-WhtBal.png',
                               'CH4file':['2021-09-05-0432_8-Jupiter_NoWV-R940(G)B889-R(G)B-WhtBal.png'],
                               'Context':{'NUVfile':'2021-09-05-0431_9-Jupiter-380NUV-Boudreau.png',
                                          'RGBfile':'2021-09-05-0423_5-Jupiter-RGB-Boudreau.png'}},
                 
                 '20210910UT':{'NH3file':'2021-09-10-0414_7-Jupiter-NoWV-DR-ST-R656G647B632-WhtBal.png',
                               'CH4file':['2021-09-10-0429_4-Jupiter_NoWV-R656(G)B889-R(G)B-WhtBal.png'],
                               'Context':{'NUVfile':'2021-09-10-0430_6-Jupiter_380NUV-Wavelets-Derot-Stack.png',
                                          'RGBfile':'2021-09-10-0421_0-Jupiter-Wavelets-R(AllRED)GB-RGB-WhtBal-ClrSmth-Smth-Wavelets.png'}},
                 
                 '20210913UT':{'NH3file':'2021-09-13-0443_3-Jupiter_NoWV-R656G647B630-RGB-WhtBal.png',
                               'CH4file':['2021-09-13-0458_3-Jupiter_NoWV-R656(G)B889-R(G)B-WhtBal.png'],
                               'Context':{'NUVfile':'2021-09-13-0504_6-Jupiter_380NUV-Wavelets-Derot-Stack.png',
                                          'RGBfile':'2021-09-13-0450_4-Jupiter-Wavelets-R(AllRED)GB-RGB-WhtBal-ClrSmth-Smth-Wavelets.png'}},
                 
                 '20210915UT':{'NH3file':'2021-09-15-0349_0-Jupiter-NoWV-DR-ST-R656G647B632-WhtBal.png',
                               'CH4file':['2021-09-15-0401_4-Jupiter_NoWV-R656(G)B889-R(G)B-WhtBal.png'],
                               'Context':{'NUVfile':'2021-09-15-0404_9-Jupiter_380NUV-Wavelets-Derot-Stack.png',
                                          'RGBfile':'2021-09-15-0355_8-Jupiter-Wavelets-R(AllRED)GB-RGB-WhtBal-ClrSmth-Smth-Wavelets.png'}},
 
                 '20210919UT':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'X/10','Transparency':'X/10'},
                               'NH3file':'2021-09-19-0508_9-Jupiter-Wavelets-R(656)G(647)B(632)-RGB-WhtBal.png',
                               'CH4file':['NA'],
                               'CH4labels':[],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2021-09-19-0505_4-Jupiter-Wavelets-R(AllRed)GB-RGB-WhtBal-ClrSmth-Smth-Wavelets.png'}},

                 '20210923UT':{'NH3file':'2021-09-23-0434_5-Jupiter_NoWV-R656G647B630-RGB-WhtBal.png',
                               'CH4file':['2021-09-23-0452_1-Jupiter_NoWV-R656(G)B889-R(G)B-WhtBal.png'],
                               'Context':{'NUVfile':'2021-09-23-0453_6-Jupiter_380NUV-Wavelets-Derot-Stack.png',
                                          'RGBfile':'2021-09-23-0442_4-Jupiter-Wavelets-RGB-RAllGB-WhtBal-ClrSmth-Smth-Wavelets.png'}},
 
                 '20210927UT':{'NH3file':'2021-09-27-0312_1-Jupiter-DR-ST-R656G647B632-WhtBal.png',
                               'CH4file':['2021-09-27-0326_7-Jupiter_NoWV-R656(G)B889-R(G)B-WhtBal.png'],
                               'Context':{'NUVfile':'2021-09-27-0329_0-Jupiter_380NUV-WV-DR-ST.png',
                                          'RGBfile':'2021-09-27-0319_8-Jupiter-Wavelets-R(ALL)GB-WV-RGB-WhtBal-ClrSmth-Smth-Wavelets.png'}},
 
                 '20211019UT':{'NH3file':'2021-10-19-0405_6-Jupiter-NoWV-R656G647B632-WhtBal.png',
                               'CH4file':['2021-10-19-0423_8-Jupiter_NoWV-R656(G)B889-R(G)B-WhtBal.png'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2021-10-19-0414_7-Jupiter-Wavelets-R685GB-WV-RGB-WhtBal-ClrSmth-Smth-Wavelets.png'}},
 
                 '20211202UT':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'X/10','Transparency':'X/10'},
                               'NH3file':'2021-12-02-0029_4-Jupiter-noWV-R656G647B632-WhtBal.png',
                               'CH4file':['2021-12-02-0031_0-Jupiter_NoWV-R656G730B632-RGB-WhtBal.png'],
                               'CH4labels':['Synth. Continuum @ 730nm','730nm (CH4)','730/Cont. (CH4)'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2021-12-02-0037_8-Jupiter-WV-R685GB-RGB-WhtBal-ClrSmth-Smth-Wavelets.png'}},
 
                 '20211203UT':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'X/10','Transparency':'X/10'},
                               'NH3file':'2021-12-03-0030_9-Jupiter-NoWV-R656G647B632-RGB-WhtBal.png',
                               'CH4file':['2021-12-03-0032_8-Jupiter_NoWV-R656G730B632-RGB-WhtBal.png'],
                               'CH4labels':['Synth. Continuum @ 730nm','730nm (CH4)','730/Cont. (CH4)'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2021-12-03-0040_5-Jupiter-WV-R685GB-RGB-WhtBal-ClrSmth-Smth-Wavelets.png'}},
                  #Special 889 CH4 Verson
                  '20220810UTa':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                             'Seeing':'7/10','Transparency':'7/10'},
                               'NH3file':'2022-08-10-1013_0-Jupiter_NoWV-R656G647B632-RGB-WhtBal.png',
                               'CH4file':['2022-08-10-1031_8-Jupiter_NoWV-R656(G)B889-R(G)B-WhtBal.png'],
                               'CH4labels':['656nm (Cont.)','889nm (CH4)','889/Cont. (CH4)'],
                               'Context':{'NUVfile':'2022-08-10-1052_4-Jupiter_889CH4-WV-DR-ST.png',
                                          'RGBfile':'2022-08-10-1030_0-Jupiter_WV-R(AllRED)GB-RGB-WhtBal-Wavelets-Str.png'},
                               'Contextlabels':['889nm','RGB']},
                  
                  '20220810UT_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                             'Seeing':'7/10','Transparency':'7/10'},
                               'NH3file':'2022-08-10-1013_0-Jupiter_NoWV-R656G647B632-RGB-WhtBal_CM2_L360_MAP-BARE.png',
                               'CH4file':['2022-08-10-1014_8-Jupiter_NoWV-R656G620B632-RGB-WhtBal_CM2_L360_MAP-BARE.png'],
                               'CH4labels':['656nm (Cont.)','889nm (CH4)','889/Cont. (CH4)'],
                               'Context':{'NUVfile':'2022-08-10-1052_4-Jupiter_889CH4-WV-DR-ST.png',
                                          'RGBfile':'2022-08-10-1030_0-Jupiter_WV-R(AllRED)GB-RGB-WhtBal-Wavelets-Str_CM2_L360_MAP-BARE.png'},
                               'Contextlabels':['889nm','RGB']},

                  '20220810UT':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'7/10','Transparency':'7/10'},
                               'NH3file':'2022-08-10-1013_0-Jupiter_NoWV-R656G647B632-RGB-WhtBal.png',
                               'CH4file':['2022-08-10-1014_8-Jupiter_NoWV-R656G620B632-RGB-WhtBal.png'],
                               'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4)'],
                               'Context':{'NUVfile':'2022-08-10-1052_4-Jupiter_889CH4-WV-DR-ST.png',
                                          'RGBfile':'2022-08-10-1030_0-Jupiter_WV-R(AllRED)GB-RGB-WhtBal-Wavelets-Str.png'},
                               'Contextlabels':['889nm','RGB']},

                  '20220812UT_Img':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'7/10','Transparency':'3/10'}, 
                               'CH4file':['2022-08-12-1025_5-Jupiter_NoWV-R656G620B632-RGB-WhtBal-CorrectedLabel.png'],
                               'NH3file':'2022-08-12-1025_5-Jupiter_NoWV-R656G647B632-RGB-WhtBalHP-CorrectedLabel.png',
                               'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'2022-08-12-1046_3-Jupiter_889CH4-WV-DR-ST.png',
                                          'RGBfile':'2022-08-12-1033_8-Jupiter_WV-R(AllRed)GB-RGB-WhtBal-Wavelets.png'},
                               'Contextlabels':['889nm','RGB']},
                  '20220812UT_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'7/10','Transparency':'3/10'}, 
                               'CH4file':['2022-08-12-1025_5-Jupiter_NoWV-R656G620B632-RGB-WhtBal-CorrectedLabel_CM2_L360_MAP-BARE.png'],
                               'NH3file':'2022-08-12-1025_5-Jupiter_NoWV-R656G647B632-RGB-WhtBalHP-CorrectedLabel_CM2_L360_MAP-BARE.png',
                               'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'2022-08-12-1046_3-Jupiter_889CH4-WV-DR-ST.png',
                                          'RGBfile':'2022-08-12-1033_8-Jupiter_WV-R(AllRed)GB-RGB-WhtBal-Wavelets_CM2_L360_MAP-BARE.png'},
                               'Contextlabels':['889nm','RGB']},

                  '20220818UT':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'7/10','Transparency':'7/10'}, 
                               'CH4file':['2022-08-18-0734_3-Jupiter_NoWV-R656G620B632-RGB-WhtBal.png'],
                               'NH3file':'2022-08-18-0733_4-Jupiter_NoWV-R656G647B632-RGB-WhtBal.png',
                               'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'2022-08-18-0801_5-Jupiter_889CH4-WV-ST-DR.png',
                                          'RGBfile':'2022-08-18-0745_4-Jupiter_AllRED-WV-RGB-WhtBal-Wavelets.png'},
                               'Contextlabels':['889nm','RGB']},

                  '20220818UT_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'7/10','Transparency':'7/10'}, 
                               'CH4file':['2022-08-18-0734_3-Jupiter_NoWV-R656G620B632-RGB-WhtBal_CM2_L360_MAP-BARE.png'],
                               'NH3file':'2022-08-18-0733_4-Jupiter_NoWV-R656G647B632-RGB-WhtBal_CM2_L360_MAP-BARE.png',
                               'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'2022-08-18-0745_4-Jupiter_AllRED-WV-RGB-WhtBal-Wavelets_CM2_L360_MAP-BARE.png',
                                          'RGBfile':'2022-08-18-0745_4-Jupiter_AllRED-WV-RGB-WhtBal-Wavelets_CM2_L360_MAP-BARE.png'},
                               'Contextlabels':['889nm','RGB']},

                  '20220828UT':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'6/10','Transparency':'7/10'}, 
                               'CH4file':['2022-08-28-0608_2-Jupiter_NoWV-R656G620B632-RGB-WhtBal.png',
                                          '2022-08-28-0608_0-Jupiter_NoWV-R656G730B632-RGB-WhtBal.png',
                                          '2022-08-28-0615_6-Jupiter_NoWV-R940G889B940-RGB-WhtBal.png'],
                               'NH3file':'2022-08-28-0608_2-Jupiter_NoWV-R656G647B632-RGB-WhtBal.png',
                               'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'2022-08-28-0611_6-Jupiter_889CH4-WV-DR-ST.png',
                                          'RGBfile':'2022-08-28-0601_3-Jupiter-RGB-JamesWillinghan-j220828a1-HalfSize.jpg'},
                               'Contextlabels':['889nm','']},
                  '20220828UT_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'6/10','Transparency':'7/10'}, 
                               'CH4file':['2022-08-28-0608_2-Jupiter_NoWV-R656G620B632-RGB-WhtBal_CM2_L360_MAP-BARE.png',
                                          '2022-08-28-0608_0-Jupiter_NoWV-R656G730B632-RGB-WhtBal_CM2_L360_MAP-BARE.png',
                                          '2022-08-28-0615_6-Jupiter_NoWV-R940G889B940-RGB-WhtBal_CM2_L360_MAP-BARE.png'],
                               'NH3file':'2022-08-28-0608_2-Jupiter_NoWV-R656G647B632-RGB-WhtBal_CM2_L360_MAP-BARE.png',
                               'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'2022-08-28-0611_6-Jupiter_889CH4-WV-DR-ST.png',
                                          'RGBfile':'2022-08-28-0601_3-Jupiter-RGB-JamesWillinghan-j220828a1_CM2_L360_MAP-BARE.jpg'},
                               'Contextlabels':['889nm','']},

                  '20220830UT_Img':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'7/10','Transparency':'7/10'}, 
                               'CH4file':['2022-08-30-0559_2-Jupiter_NoWV-R656G620B632-RGB-WhtBal.png',
                                          '2022-08-30-0559_1-Jupiter_NoWV-R656G730B632-RGB-WhtBal.png',
                                          '2022-08-30-0605_6-Jupiter_NoWV-R940G889B940-RGB-WhtBal.png'],
                               'NH3file':'2022-08-30-0559_1-Jupiter_NoWV-R656G647B632-RGB-WhtBal.png',
                               'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'2022-08-30-0604_2-Jupiter_889CH4-WV-DR-ST.png',
                                          'RGBfile':'2022-08-28-1429_7-Jupiter-RGB-Arakawa-j220828e1.bmp'},
                               'Contextlabels':['889nm','']},
                  '20220830UT_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'7/10','Transparency':'7/10'}, 
                               'CH4file':['2022-08-30-0559_2-Jupiter_NoWV-R656G620B632-RGB-WhtBal_CM2_L360_MAP-BARE.png',
                                          '2022-08-30-0559_1-Jupiter_NoWV-R656G730B632-RGB-WhtBal_CM2_L360_MAP-BARE.png',
                                          '2022-08-30-0605_6-Jupiter_NoWV-R940G889B940-RGB-WhtBal_CM2_L360_MAP-BARE.png'],
                               'NH3file':'2022-08-30-0559_1-Jupiter_NoWV-R656G647B632-RGB-WhtBal_CM2_L360_MAP-BARE.png',
                               'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'2022-08-30-0604_2-Jupiter_889CH4-WV-DR-ST.png',
                                          'RGBfile':'2022-08-28-1429_7-Jupiter-RGB-Arakawa-j220828e1_CM2_L360_MAP-BARE.bmp'},
                               'Contextlabels':['889nm','']},
                                    
                  '20220901UT_Img':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'6/10','Transparency':'7/10'}, 
                               'CH4file':['2022-09-01-0604_9-Jupiter_NoWV-R656G620B632-RGB-WhtBal.png',
                                          '2022-09-01-0604_8-Jupiter_NoWV-R656G730B632-RGB-WhtBal.png',
                                          '2022-09-01-0601_7-Jupiter_NoWV-R940G889B940-RGB-WhtBal.png'],
                               'NH3file':'2022-09-01-0604_9-Jupiter_NoWV-R656G647B632-RGB-WhtBal.png',
                               'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'2022-09-01-0603_5-Jupiter_889CH4-WV-DR-ST.png',
                                          'RGBfile':'2022-09-01-0618_0-Jupiter-Rivera-j220901l1.png'},
                               'Contextlabels':['889nm','']},
                  '20220901UT_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'6/10','Transparency':'7/10'}, 
                               'CH4file':['2022-09-01-0604_9-Jupiter_NoWV-R656G620B632-RGB-WhtBal_CM2_L360_MAP-BARE.png',
                                          '2022-09-01-0604_8-Jupiter_NoWV-R656G730B632-RGB-WhtBal_CM2_L360_MAP-BARE.png',
                                          '2022-09-01-0601_7-Jupiter_NoWV-R940G889B940-RGB-WhtBal_CM2_L360_MAP-BARE.png'],
                               'NH3file':'2022-09-01-0604_9-Jupiter_NoWV-R656G647B632-RGB-WhtBal_CM2_L360_MAP-BARE.png',
                               'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'2022-09-01-0603_5-Jupiter_889CH4-WV-DR-ST.png',
                                          'RGBfile':'2022-09-01-0618_0-Jupiter-Rivera-j220901l1_CM2_L360_MAP-BARE.png'},
                               'Contextlabels':['889nm','']},

                  '20220904UT_Img':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'7/10','Transparency':'7/10'}, 
                               'CH4file':['2022-09-04-0638_9-Jupiter_NoWV-R656G620B632-RGB-WhtBal.png',
                                          '2022-09-04-0638_6-Jupiter_NoWV-R656G730B632-RGB-WhtBal.png',
                                          '2022-09-04-0638_6-Jupiter_NoWV-R940G889B940-RGB-WhtBal.png'],
                               'NH3file':'2022-09-04-0638_2-Jupiter_NoWV-R656G647B632-RGB-WhtBal.png',
                               'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'2022-09-04-0639_1-Jupiter_889CH4-WV-DR-ST.png',
                                          'RGBfile':'NA'},
                               'Contextlabels':['889nm','']},
                  '20220904UT_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'7/10','Transparency':'7/10'}, 
                               'CH4file':['2022-09-04-0638_9-Jupiter_NoWV-R656G620B632-RGB-WhtBal_CM2_L360_MAP-BARE.png',
                                          '2022-09-04-0638_6-Jupiter_NoWV-R656G730B632-RGB-WhtBal_CM2_L360_MAP-BARE.png',
                                          '2022-09-04-0638_6-Jupiter_NoWV-R940G889B940-RGB-WhtBal_CM2_L360_MAP-BARE.png'],
                               'NH3file':'2022-09-04-0638_2-Jupiter_NoWV-R656G647B632-RGB-WhtBal_CM2_L360_MAP-BARE.png',
                               'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'2022-09-04-0639_1-Jupiter_889CH4-WV-DR-ST.png',
                                          'RGBfile':'2022-09-04-1644_7-Jupiter-Yamane-j220904j4_CM2_L360_MAP-BARE.jpg'},
                               'Contextlabels':['889nm','']},

                  '20220905UT_Img':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'7/10','Transparency':'7/10'}, 
                               'CH4file':['2022-09-05-0559_2-Jupiter_NoWV-R656G620B632-RGB-WhtBal.png',
                                          '2022-09-05-0559_2-Jupiter_NoWV-R656G730B632-RGB-WhtBal.png',
                                          '2022-09-05-0559_0-Jupiter_NoWV-R940G889B940-RGB-WhtBal.png'],
                               'NH3file':'2022-09-05-0559_1-Jupiter_NoWV-R656G647B632-RGB-WhtBal.png',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'2022-09-05-0558_9-Jupiter_889CH4-WV-DR-ST.png',
                                          'RGBfile':'2022-09-05-1559_0-Jupiter_j220905l1_Michael_Wong.png'},
                               'Contextlabels':['889nm','']},
                  '20220905UT_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'7/10','Transparency':'7/10'}, 
                               'CH4file':['2022-09-05-0559_2-Jupiter_NoWV-R656G620B632-RGB-WhtBal_CM2_L360_MAP-BARE.png',
                                          '2022-09-05-0559_2-Jupiter_NoWV-R656G730B632-RGB-WhtBal_CM2_L360_MAP-BARE.png',
                                          '2022-09-05-0559_0-Jupiter_NoWV-R940G889B940-RGB-WhtBal_CM2_L360_MAP-BARE.png'],
                               'NH3file':'2022-09-05-0559_1-Jupiter_NoWV-R656G647B632-ImgShiftS10-RGB-WhtBal_CM2_L360_MAP-BARE.png',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'2022-09-05-0558_9-Jupiter_889CH4-WV-DR-ST.png',
                                          'RGBfile':'2022-09-05-1559_0-Jupiter_j220905l1_Michael_Wong_CM2_L360_MAP-BARE.png'},
                               'Contextlabels':['889nm','']},

                  '20220905UTSmth_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'7/10','Transparency':'7/10'}, 
                               'CH4file':['2022-09-05-0559_2-Jupiter_NoWV-R656G620B632-RGB-WhtBal-Smth_CM2_L360_MAP-BARE.png',
                                          '2022-09-05-0559_2-Jupiter_NoWV-R656G730B632-RGB-WhtBal_CM2_L360_MAP-BARE.png',
                                          '2022-09-05-0559_0-Jupiter_NoWV-R940G889B940-RGB-WhtBal_CM2_L360_MAP-BARE.png'],
                               'NH3file':'2022-09-05-0559_1-Jupiter_NoWV-R656G647B632-RGB-WhtBal-Smth_CM2_L360_MAP-BARE.png',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'2022-09-05-0558_9-Jupiter_889CH4-WV-DR-ST.png',
                                          'RGBfile':'2022-09-05-1559_0-Jupiter_j220905l1_Michael_Wong_CM2_L360_MAP-BARE.png'},
                               'Contextlabels':['889nm','']},
                  
                  '20220905UTSa':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'7/10','Transparency':'7/10'}, 
                               'CH4file':['2022-09-05-0454_5-Saturn_R656G620B632-RGB-WhtBal.png',
                                          '2022-09-05-0454_4-Saturn_R656G730B632-RGB-WhtBal.png',
                                          '2022-09-05-0453_7-Saturn_R940G889B940-RGB-WhtBal.png'],
                               'NH3file':'2022-09-05-0454_5-Saturn_R656G647B632-RGB-WhtBal.png',
                               'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'2022-09-05-0453_6-Saturn_889CH4-WV-DR-ST.png',
                                          'RGBfile':'NA'},
                               'Contextlabels':['889nm','']},
                  
                  '20220912UT':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'8/10','Transparency':'7/10'}, 
                               'CH4file':['2022-09-12-0533_4-Jupiter_WV-R656G620B632-RGB-WhtBal.png'],
                               'NH3file':'2022-09-12-0533_4-Jupiter_WV-R656G647B632-RGB-WhtBal.png',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2022-09-12-0532_3-Jupiter_WV-R685G550B450-RGB-WhtBal-ClrSmth-Wavelets.png'},
                               'Contextlabels':['889nm','']},
                  '20220912UT_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'8/10','Transparency':'7/10'}, 
                               'CH4file':['2022-09-12-0533_4-Jupiter_WV-R656G620B632-RGB-WhtBal_CM2_L360_MAP-BARE.png'],
                               'NH3file':'2022-09-12-0533_4-Jupiter_WV-R656G647B632-RGB-WhtBal_CM2_L360_MAP-BARE.png',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2022-09-12-0532_3-Jupiter_WV-R685G550B450-RGB-WhtBal-ClrSmth-Wavelets_CM2_L360_MAP-BARE.png'},
                               'Contextlabels':['889nm','']},

                  '20220912UTSa':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'8/10','Transparency':'7/10'}, 
                               'CH4file':['2022-09-12-0432_5-Saturn_NoWV-R656G620B632-RGB-WhtBal.png'],
                               'NH3file':'2022-09-12-0432_5-Saturn_NoWV-R656G647B632-RGB-WhtBal.png',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2022-09-12-0413_1-Saturn_WV-RGB-WhtBal.png'},
                               'Contextlabels':['889nm','']},
                  
                  '20220913UT_Img':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'6/10','Transparency':'7/10'}, 
                               'CH4file':['2022-09-13-0457_4-Jupiter_NoWV-R656G620B632-RGB-WhtBal.png'],
                               'NH3file':'2022-09-13-0457_4-Jupiter_NoWV-R656G647B632-RGB-WhtBal.png',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2022-09-13-0455_6-Jupiter_WV-R685G550B450-RGB-WhtBal-ClrSmth-Wavelets.png'},
                               'Contextlabels':['889nm','']},
                  '20220913UT_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'6/10','Transparency':'7/10'}, 
                               'CH4file':['2022-09-13-0457_4-Jupiter_NoWV-R656G620B632-RGB-WhtBal_CM2_L360_MAP-BARE.png'],
                               'NH3file':'2022-09-13-0457_4-Jupiter_NoWV-R656G647B632-RGB-WhtBal_CM2_L360_MAP-BARE.png',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2022-09-13-0455_6-Jupiter_WV-R685G550B450-RGB-WhtBal-ClrSmth-Wavelets_CM2_L360_MAP-BARE.png'},
                               'Contextlabels':['889nm','']},

                  '20220919UTb_Img':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'6/10','Transparency':'7/10'}, 
                               'CH4file':['2022-09-19-0452_7-Jupiter_NoWV-R656G620B632-RGB-WhtBal.png'],
                               'NH3file':'2022-09-19-0453_4-Jupiter_NoWV-R656G647B632-RGB-WhtBal.png',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2022-09-19-0518_7-Jupiter_WV-R685G550B450-RGB-WhtBal-Wavelets.png'},
                               'Contextlabels':['889nm','']},
                  '20220919UTb_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'6/10','Transparency':'7/10'}, 
                               'CH4file':['2022-09-19-0452_7-Jupiter_NoWV-R656G620B632-RGB-WhtBal_CM2_L360_MAP-BARE.png'],
                               'NH3file':'2022-09-19-0453_4-Jupiter_NoWV-R656G647B632-RGB-WhtBal_CM2_L360_MAP-BARE.png',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2022-09-19-0518_7-Jupiter_WV-R685G550B450-RGB-WhtBal-Wavelets_CM2_L360_MAP-BARE.png'},
                               'Contextlabels':['889nm','']},
                  
                  '20220919UTa_Img':{'Metadata':{'Telescope':'VLT','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'6/10','Transparency':'7/10'}, 
                               'CH4file':['2022-09-19-0352_3-Jupiter-VLT-R656G620B632-RGB-WhtBal-MedFilt.png'],
                               'NH3file':'2022-09-19-0352_3-Jupiter-VLT-R656G647B632-RGB-WhtBal-MedFilt.png',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2022-09-19-0352_3-Jupiter-MUSE-R650G550B450-RGB-WhtBal-Filtered.png'},
                               'Contextlabels':['889nm','']},
                  '20220919UTa_Map':{'Metadata':{'Telescope':'VLT','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'6/10','Transparency':'7/10'}, 
                               'CH4file':['2022-09-19-0352_3-Jupiter-VLT-R656G620B632-RGB-WhtBal-MedFilt_CM2_L360_MAP-BARE.png'],
                               'NH3file':'2022-09-19-0352_3-Jupiter-VLT-R656G647B632-RGB-WhtBal-MedFilt_CM2_L360_MAP-BARE.png',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2022-09-19-0352_3-Jupiter-MUSE-R650G550B450-RGB-WhtBal-Filtered_CM2_L360_MAP-BARE.png'},
                               'Contextlabels':['889nm','']},

                  '20220919UTSa':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'6/10','Transparency':'7/10'}, 
                               'CH4file':['2022-09-19-0412_5-Saturn_NoWV-R656G620B632-RGB-WhtBal.png'],
                               'NH3file':'2022-09-19-0412_5-Saturn_NoWV-R656G647B632-RGB-WhtBal.png',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2022-09-19-0354_0-Saturn_WV-R685G550B450-RGB-WhtBal-Wavelets.png'},
                               'Contextlabels':['889nm','']},

                  '20220925UT_Img':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'5/10','Transparency':'7/10'}, 
                               'CH4file':['2022-09-25-0615_4-Jupiter_NoWV-R656G620B632-RGB-WhtBal.png'],
                               'NH3file':'2022-09-25-0615_6-Jupiter_NoWV-R656G647B632-RGB-WhtBal.png',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2022-09-25-0546_6-Jupiter_WV-R685G550B450-RGB-WhtBal-Wavelets.png'},
                               'Contextlabels':['889nm','RGB']},
                  '20220925UT_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'5/10','Transparency':'7/10'}, 
                               'CH4file':['2022-09-25-0615_4-Jupiter_NoWV-R656G620B632-RGB-WhtBal_CM2_L360_MAP-BARE.png'],
                               'NH3file':'2022-09-25-0615_6-Jupiter_NoWV-R656G647B632-RGB-WhtBal_CM2_L360_MAP-BARE.png',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2022-09-25-0546_6-Jupiter_WV-R685G550B450-RGB-WhtBal-Wavelets_CM2_L360_MAP-BARE.png'},
                               'Contextlabels':['889nm','RGB']},
                  '20220925UTS_Img':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'5/10','Transparency':'7/10'}, 
                               'CH4file':['2022-09-25-0500_4-Saturn_NoWV-R656G620B632-RGB-WhtBal.png'],
                               'NH3file':'2022-09-25-0500_4-Saturn_NoWV-R656G647B632-RGB-WhtBal.png',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2022-09-25-0441_4-Saturn_WV-R685G550B450-RGB-WhtBal-Wavelets.png'},
                               'Contextlabels':['889nm','RGB']},

                  '20221009UTa_Img':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'7/10','Transparency':'7/10'}, 
                               'CH4file':['2022-10-09-0401_5-Jupiter_NoWV-R656G620B632-RGB-WhtBal.png'],
                               'NH3file':'2022-10-09-0401_5-Jupiter_NoWV-R656G647B632-RGB-WhtBal.png',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2022-10-09-0339_0-Jupiter_WV2-R685G550B450-RGB-WhtBal-Wavelets.png'},
                               'Contextlabels':['889nm','RGB']},

                  '20221009UTa_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'7/10','Transparency':'7/10'}, 
                               'CH4file':['2022-10-09-0401_5-Jupiter_NoWV-R656G620B632-RGB-WhtBal_CM2_L360_MAP-BARE.png'],
                               'NH3file':'2022-10-09-0401_5-Jupiter_NoWV-R656G647B632-RGB-WhtBal_CM2_L360_MAP-BARE.png',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2022-10-09-0339_0-Jupiter_WV-R685G550B450-RGB-WhtBal-Wavelets_CM2_L360_MAP-BARE.png'},
                               'Contextlabels':['889nm','RGB']},

                  '20221009UTb_Img':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'7/10','Transparency':'7/10'}, 
                               'CH4file':['2022-10-09-0524_5-Jupiter_NoWV-R656G620B632-RGB-WhtBal.png'],
                               'NH3file':'2022-10-09-0524_5-Jupiter_NoWV-R656G647B632-RGB-WhtBal.png',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2022-10-09-0542_8-Jupiter_NoWV-R685G550B450-RGB-WhtBal-Str0to160-Wavelets.png'},
                               'Contextlabels':['889nm','RGB']},

                  '20221009UTb_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'7/10','Transparency':'7/10'}, 
                               'CH4file':['2022-10-09-0524_5-Jupiter_NoWV-R656G620B632-RGB-WhtBal_CM2_L360_MAP-BARE.png'],
                               'NH3file':'2022-10-09-0524_5-Jupiter_NoWV-R656G647B632-RGB-WhtBal_CM2_L360_MAP-BARE.png',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2022-10-09-0542_8-Jupiter_NoWV-R685G550B450-RGB-WhtBal-Str0to160-Wavelets_CM2_L360_MAP-BARE.png'},
                               'Contextlabels':['889nm','RGB']},
                  
                  '20221013UT_Img':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'5/10','Transparency':'8/10'}, 
                               'CH4file':['2022-10-13-0345_5-Jupiter-NoWV-R656G620B632-RGB-WhtBal.png'],
                               'NH3file':'2022-10-13-0345_5-Jupiter-NoWV-R656G647B632-RGB-WhtBal.png',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2022-10-13-0402_0-Jupiter-WV-R685G550B450-RGB-WhtBal-Wavelets.png'},
                               'Contextlabels':['889nm','RGB']},

                  '20221013UT_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'5/10','Transparency':'8/10'}, 
                               'CH4file':['2022-10-13-0345_5-Jupiter-NoWV-R656G620B632-RGB-WhtBal_CM2_L360_MAP-BARE.png'],
                               'NH3file':'2022-10-13-0345_5-Jupiter-NoWV-R656G647B632-RGB-WhtBal_CM2_L360_MAP-BARE.png',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2022-10-13-0402_0-Jupiter-WV-R685G550B450-RGB-WhtBal-Wavelets_CM2_L360_MAP-BARE.png'},
                               'Contextlabels':['889nm','RGB']},

                  '20221019UT_Img':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'6/10','Transparency':'7/10'}, 
                               'CH4file':['2022-10-19-0342_4-Jupiter-NoWV-R656G620B632-RGB-WhtBal.png'],
                               'NH3file':'2022-10-19-0342_4-Jupiter-NoWV-R656G647B632-RGB-WhtBal.png',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2022-10-19-0358_2-Jupiter-WV-R685G550B450-RGB-WhtBal-Wavelets.png'},
                               'Contextlabels':['889nm','RGB']},

                  '20221019UT_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'6/10','Transparency':'7/10'}, 
                               'CH4file':['2022-10-19-0342_4-Jupiter-NoWV-R656G620B632-RGB-WhtBal_CM2_L360_MAP-BARE.png'],
                               'NH3file':'2022-10-19-0342_4-Jupiter-NoWV-R656G647B632-RGB-WhtBal_CM2_L360_MAP-BARE.png',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2022-10-19-0358_2-Jupiter-WV-R685G550B450-RGB-WhtBal-Wavelets_CM2_L360_MAP-BARE.png'},
                               'Contextlabels':['889nm','RGB']},
                  '20221020UT_Img':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'6/10','Transparency':'7/10'}, 
                               'CH4file':['2022-10-20-0440_4-Jupiter-NoWV-R656G620B632-RGB-WhtBal.png'],
                               'NH3file':'2022-10-20-0440_4-Jupiter-NoWV-R656G647B632-RGB-WhtBal.png',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2022-10-20-0422_6-Jupiter-WV-R685G550B450-RGB-WhtBal-Wavelets.png'},
                               'Contextlabels':['889nm','RGB']},
                  '20221020UT_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'6/10','Transparency':'7/10'}, 
                               'CH4file':['2022-10-20-0440_4-Jupiter-NoWV-R656G620B632-RGB-WhtBal_CM2_L360_MAP-BARE.png'],
                               'NH3file':'2022-10-20-0440_4-Jupiter-NoWV-R656G647B632-RGB-WhtBal_CM2_L360_MAP-BARE.png',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2022-10-20-0422_6-Jupiter-WV-R685G550B450-RGB-WhtBal-Wavelets_CM2_L360_MAP-BARE.png'},
                               'Contextlabels':['889nm','RGB']},
                  '20221021UT_Img':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'6/10','Transparency':'7/10'}, 
                               'CH4file':['2022-10-21-0358_6-Jupiter-NoWV-R656G620B632-RGB-WhtBal.png'],
                               'NH3file':'2022-10-21-0358_6-Jupiter-NoWV-R656G647B632-RGB-WhtBal.png',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2022-10-21-0342_1-Jupiter-WV-R685G550B450-RGB-WhtBal-Wavelets.png'},
                               'Contextlabels':['889nm','RGB']},
                  '20221021UT_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'6/10','Transparency':'7/10'}, 
                               'CH4file':['2022-10-21-0358_6-Jupiter-NoWV-R656G620B632-RGB-WhtBal_CM2_L360_MAP-BARE.png'],
                               'NH3file':'2022-10-21-0358_6-Jupiter-NoWV-R656G647B632-RGB-WhtBal_CM2_L360_MAP-BARE.png',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2022-10-21-0342_1-Jupiter-WV-R685G550B450-RGB-WhtBal-Wavelets_CM2_L360_MAP-BARE.png'},
                               'Contextlabels':['889nm','RGB']},
                  '20230113UT_Img':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'8/10','Transparency':'7/10'}, 
                               'CH4file':['2023-01-13-0046_2-Jupiter-NoWV-R656G620B632-RGB-WhtBal.png'],
                               'NH3file':'2023-01-13-0046_2-Jupiter-NoWV-R656G647B632-RGB-WhtBal.png',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2023-01-13-0046_0-Jupiter-WV-R685G550B450-RGB-WhtBal-Wavelets.png'},
                               'Contextlabels':['889nm','RGB']},
                  '20230113UT_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'6/10','Transparency':'7/10'}, 
                               'CH4file':['2023-01-13-0046_2-Jupiter-NoWV-R656G620B632-RGB-WhtBal_CM2_L360_MAP-BARE.png'],
                               'NH3file':'2023-01-13-0046_2-Jupiter-NoWV-R656G647B632-RGB-WhtBal_CM2_L360_MAP-BARE.png',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2023-01-13-0046_0-Jupiter-WV-R685G550B450-RGB-WhtBal-Wavelets_CM2_L360_MAP-BARE.png'},
                               'Contextlabels':['889nm','RGB']},
                  '20230815UT_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'9/10','Transparency':'9/10'}, 
                               'CH4file':['2023-08-15-1112_8-Jupiter_R656G620B632-NoWV-RGB-Smth-WhtBal_CM2_L360_MAP-BARE.png'],
                               'NH3file':'2023-08-15-1112_8-Jupiter_R656G647B632-NoWV-RGB-Smth-WhtBal_CM2_L360_MAP-BARE.png',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2023-08-15-1129_3-Jupiter_R685G550B450-WV-RGB-WhtBal-ClrSmth-WV_CM2_L360_MAP-BARE.png'},
                               'Contextlabels':['889nm','RGB']},
                  '20230815UTSwap_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'9/10','Transparency':'9/10'}, 
                               'CH4file':['2023-08-15-1112_8-Jupiter_R656G620B632-NoWV-620swap632-RGB-Smth-WhtBal_CM2_L360_MAP-BARE.png'],
                               'NH3file':'2023-08-15-1112_9-Jupiter_R656G647B632-NoWV-620swap632-RGB-Smth-WhtBal_CM2_L360_MAP-BARE.png',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2023-08-15-1129_3-Jupiter_R685G550B450-WV-RGB-WhtBal-ClrSmth-WV_CM2_L360_MAP-BARE.png'},
                               'Contextlabels':['889nm','RGB']},
                  '20230815UTGauss632_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'9/10','Transparency':'9/10'}, 
                               'CH4file':['2023-08-15-1112_8-Jupiter_R656G620B632Gauss3pix-NoWV-RGB-Smth-WhtBal_CM2_L360_MAP-BARE.png'],
                               'NH3file':'2023-08-15-1112_8-Jupiter_R656G647B632-NoWV-RGB-Smth-WhtBal-Gauss632NH3_CM2_L360_MAP-BARE.png',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2023-08-15-1129_3-Jupiter_R685G550B450-WV-RGB-WhtBal-ClrSmth-WV_CM2_L360_MAP-BARE.png'},
                               'Contextlabels':['889nm','RGB']},
                  '20230816UT_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'7/10','Transparency':'7/10'}, 
                               'CH4file':['2023-08-16-1130_1-Jupiter_R656G620B632-NoWV-RGB-Smth-WhtBal_CM2_L360_MAP-BARE.png'],
                               'NH3file':'2023-08-16-1130_1-Jupiter_R656G647B632-NoWV-RGB-Smth-WhtBal_CM2_L360_MAP-BARE.png',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2023-08-16-1145_3-Jupiter_R685G550B450-WV-RGB-WhtBal-ClrSmth-WV_CM2_L360_MAP-BARE.png'},
                               'Contextlabels':['889nm','RGB']},
                  '20230817UT_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'8/10','Transparency':'7/10'}, 
                               'CH4file':['2023-08-17-1121_4-Jupiter_R656G620B632-RGB-Smth-WhtBal_CM2_L360_MAP-BARE.png'],
                               'NH3file':'2023-08-17-1121_4-Jupiter_R656G647B632-RGB-Smth-WhtBal_CM2_L360_MAP-BARE.png',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2023-08-17-1136_7-Jupiter_R685B550B450-WV-RGB-ClrSmth-WhtBal-WV_CM2_L360_MAP-BARE.png'},
                               'Contextlabels':['889nm','RGB']},
                  '20230818UT_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'7/10','Transparency':'4/10'}, 
                               'CH4file':['2023-08-18-1137_1-Jupiter_R656G620B632-NoWV-RGB-ClrSmth-WhtBal_CM2_L360_MAP-BARE.png'],
                               'NH3file':'2023-08-18-1137_1-Jupiter_R656G647B632-NoWV-RGB-ClrSmth-WhtBal_CM2_L360_MAP-BARE.png',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2023-08-18-1151_5-Jupiter_R685B550B450-WV-RGB-ClrSmth-WhtBal-WV_CM2_L360_MAP-BARE.png'},
                               'Contextlabels':['889nm','RGB']},
                  '20230827UTa_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'10/10','Transparency':'7/10'}, 
                               'CH4file':['2023-08-27-1130_6-Jupiter_R656G620B632-RGB-Smth-WhtBal_CM2_L360_MAP-BARE.png'],
                               'NH3file':'2023-08-27-1130_6-Jupiter_R656G647B632-RGB-Smth-WhtBal_CM2_L360_MAP-BARE.png',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2023-08-27-1144_1-Jupiter_R656B647B632-WV-RGB-WhtBal-ClrSmth-WV-Str0to192_CM2_L360_MAP-BARE.png'},
                               'Contextlabels':['889nm','RGB']},
                  '20230827UTb_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'10/10','Transparency':'7/10'}, 
                               'CH4file':['2023-08-27-1152_6-Jupiter_R656G620B632-RGB-Smth-WhtBal_CM2_L360_MAP-BARE.png'],
                               'NH3file':'2023-08-27-1153_8-Jupiter_R656G647B632-RGB-Smth-WhtBal_CM2_L360_MAP-BARE.png',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2023-08-27-1144_1-Jupiter_R656B647B632-WV-RGB-WhtBal-ClrSmth-WV-Str0to192_CM2_L360_MAP-BARE.png'},
                               'Contextlabels':['889nm','RGB']},
                  '20230830UT_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'7/10','Transparency':'8/10'}, 
                               'CH4file':['2023-08-30-1144_1-Jupiter_R656G620B632-RGB-Smth-WhtBal_CM2_L360_MAP-BARE.png'],
                               'NH3file':'2023-08-30-1144_1-Jupiter_R656G647B632-RGB-Smth-WhtBal_CM2_L360_MAP-BARE.png',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2023-08-30-1158_5-Jupiter_R685G550B450-RGB-WhtBal-ClrSmth-WV_CM2_L360_MAP-BARE.png'},
                               'Contextlabels':['889nm','RGB']},
                  '20230831UTa_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'6/10','Transparency':'8/10'}, 
                               'CH4file':['2023-08-31-1130_8-Jupiter_R656G620B632-NoWV-RGB-Smth-WhtBal_CM2_L360_MAP-BARE.png'],
                               'NH3file':'2023-08-31-1130_7-Jupiter_R656G647B632-NoWV-RGB-Smth-WhtBal_CM2_L360_MAP-BARE.png',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2023-08-31-1146_4-Jupiter_R685G550B450-RGB-ClrSmth-WhtBal-WV_CM2_L360_MAP-BARE.png'},
                               'Contextlabels':['889nm','RGB']},
                  '20230831UTb_Map':{'Metadata':{'Telescope':'C11','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'6/10','Transparency':'8/10'}, 
                               'CH4file':['2023-08-31-1153_5-Jupiter_R656G620B632-NoWV-RGB-Smth-WhtBal_CM2_L360_MAP-BARE.png'],
                               'NH3file':'2023-08-31-1153_7-Jupiter_R656G647B632-NoWV-RGB-Smth-WhtBal_CM2_L360_MAP-BARE.png',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2023-08-31-1146_4-Jupiter_R685G550B450-RGB-ClrSmth-WhtBal-WV_CM2_L360_MAP-BARE.png'},
                               'Contextlabels':['889nm','RGB']}}

    ###########################################################################
    # OBTAIN IMAGES TO DISPLAY AND DETERMINE IMAGE ARRAY SIZE
    ###########################################################################             
    path='c:/Astronomy/Projects/Planets/'+target+'/Imaging Data/'+obsdate[0:10]+'/'
    NH3file=sourcefiles[sourcedata]['NH3file']
    CH4file=sourcefiles[sourcedata]['CH4file']
    NUVfile=sourcefiles[sourcedata]['Context']['NUVfile']
    RGBfile=sourcefiles[sourcedata]['Context']['RGBfile']

    print(NUVfile,RGBfile)
    nny=1
    #print(nny)
    if NH3file != 'NA':
        NH3_RGB=load_png(path+NH3file)
        nny=nny+1
    if CH4file != ['NA']:
        nny=nny+len(CH4file)
        print("in CH4")
    if NUVfile != 'NA':
        NUV_RGB=load_png(path+NUVfile)
    else:
        NUV_RGB=[]
    if RGBfile != 'NA':
        RGB_RGB=imread(path+RGBfile)
    else:
        RGB_RGB=[]
    #print(nny)
    indices=(NH3_RGB>0.05*NH3_RGB.max())
    mask=np.zeros(NH3_RGB.shape)
    mask[indices]=1.
    
    if NUVfile !='NA' or RGBfile != 'NA':
        nny=nny+1
        
    print(nny)
    if nny==3:
        ny,nx,dx,dy=3,3,6.0,6.0
    elif nny==4:
        ny,nx,dx,dy=4,3,5.0,10.0
    elif nny==5:
        ny,nx,dx,dy=5,3,5.0,13.0
    elif nny==6:
        ny,nx,dx,dy=6,3,5.0,10.0

    ###########################################################################
    # SET UP PLOT/CANVAS AND CREATE FIGURE HEADER
    ###########################################################################

    fig,ax=pl.subplots(ny,nx,figsize=(dx,dy), dpi=150, facecolor="black")#,
                          #sharey=True,sharex=True)
    fig.suptitle(NH3file[0:10],x=0.5,ha='center',color='w',fontsize=10)
    #fig.suptitle("Methane and Ammonia Experiment "+NH3file[0:10],
    #             x=0.5,ha='center',color='w',fontsize=12)
    ax[0,1].annotate("Steven Hill",[0.5,1.2],ha='center',xycoords='axes fraction',color='w',fontsize=8)
    ax[0,1].annotate("Denver, Colorado",[0.5,1.1],ha='center',xycoords='axes fraction',color='w',fontsize=8)
    
    ax[0,0].annotate(sourcefiles[sourcedata]['Metadata']['Telescope'],[0.07,1.3],
                 ha='left',xycoords='axes fraction',color='w',fontsize=8)
    ax[0,0].annotate("5.60m FL",[0.07,1.2],
                 ha='left',xycoords='axes fraction',color='w',fontsize=8)
    ax[0,0].annotate("ASI120MM",[0.07,1.1],
                 ha='left',xycoords='axes fraction',color='w',fontsize=8)
    ax[0,2].annotate("Seeing="+sourcefiles[sourcedata]['Metadata']['Seeing'],[0.95,1.3],
                 ha='right',xycoords='axes fraction',color='w',fontsize=8)
    ax[0,2].annotate("Transparency="+sourcefiles[sourcedata]['Metadata']['Transparency'],[0.95,1.2],
                 ha='right',xycoords='axes fraction',color='w',fontsize=8)
    
    ###########################################################################
    # COMPUTE AMMONIA IMAGE PANELS ALONG WITH COLOR SLOPE (VALID AT NH3 EPOCH)
    ###########################################################################
    clrslp,nh3abs=Ammonia_Panels(NH3file,NH3_RGB,mask,ax) #INCLUDE "ROW" HERE BASED ON nx
    #WRITE PNG FILES WITH MANUAL SCALING
    #  !!!!NEED BETTER FILE NAMING CONVENTION HERE
    #  !!!!OR, FOR SEPARATION OF FUNCTION, CREATE A STANDALONG FILE-WRITING FUNCTION
    pathClrSlp='C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Analysis Data/Color Slope/'
    if imagetype=="Map":
        fnout=pathClrSlp+NH3file[0:26]+'ClrSlpMap'
    elif imagetype=="Img":
        fnout=pathClrSlp+NH3file[0:26]+'ClrSlpAbsImg'
    print('#############NH3file[0:26]=',NH3file[0:26])    
    print('#############CH4file[0:26]=',CH4file[0][0:26])    
    hdu = fits.PrimaryHDU(clrslp.astype(np.float32))
    hdul = fits.HDUList([hdu])
    hdul[0].header['BITPIX']=-32
    print(hdul[0].header[:])
    #fnout=path+'/'+NH3file[0:26]+'NH3Abs647.fits'
    try:
        os.remove(fnout+'.fits')
    except: 
        print("file doesn't exist")
    hdul.writeto(fnout+'.fits')
    hdul.close()
        
    clrslp16bit = np.nan_to_num((32000+500*(clrslp)).astype(np.uint16))
    imwrite(fnout+'.png', clrslp16bit.astype(np.uint16))
    #  !!!!OR, FOR SEPARATION OF FUNCTION, CREATE A STANDALONG FILE-WRITING FUNCTION   
    pathL2FITS='C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Analysis Data/L2 FITS/'
    if imagetype=="Map":
        fnout=pathL2FITS+NH3file[0:26]+'647NH3AbsMap'
    elif imagetype=="Img":
        fnout=pathL2FITS+NH3file[0:26]+'647NH3AbsImg'
    maskednh3=nh3abs*mask[:,:,1]
    normnh3=maskednh3/maskednh3[maskednh3>0].mean()
    hdu = fits.PrimaryHDU(normnh3.astype(np.float32))
    hdul = fits.HDUList([hdu])
    hdul[0].header['BITPIX']=-32
    print(hdul[0].header[:])
    #fnout=path+'/'+NH3file[0:26]+'NH3Abs647.fits'
    try:
        os.remove(fnout+'.fits')
    except: 
        print("file doesn't exist")
    hdul.writeto(fnout+'.fits')
    hdul.close()
    #nh3abs16bit = np.nan_to_num(((5.*65535.*(normnh3 - 0.9))*mask[:,:,1]).astype(np.uint16))
    normnh3scaled=np.nan_to_num(((5.*65535.*(normnh3*mask[:,:,1] - 0.9))))
    normnh3scaled[normnh3scaled<=0.]=0.0
    nh3abs16bit = normnh3scaled.astype(np.uint16)
    imwrite(fnout+'.png', nh3abs16bit)#.astype(np.uint16))

    ###########################################################################
    # COMPUTE METHANE IMAGE PANELS ALONG WITH COLOR SLOPE (VALID AT NH3 EPOCH)
    #   UP TO THREE SETS OF METHANE PANELS ARE COMPUTED (620NM, 730NM, 889NM)
    #   COULD ADD >1000NM, BUT THAT WOULD BE FOR A LATER DATE
    ###########################################################################
    #
    if CH4file != ['NA']: #NEED LOOPED CALLS HERE BASED ON NUMBER OF FILES AND LOOK UP FOR LABELS
        for ifile in range(0,len(CH4file)):
            CH4_RGB=load_png(path+CH4file[ifile])
            print(CH4file[ifile])
            if search("G620", CH4file[ifile]):
                wvstr='620'
                print("***********620")
                CH4labels=['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))']
            elif search("G730", CH4file[ifile]):
                wvstr='730'
                print("***********730")
                CH4labels=['Synth. Continuum @ 730nm','730nm (CH4)','730/Cont. (CH4))']
            elif search("G889", CH4file[ifile]):
                wvstr='889'
                print("***********889")
                CH4labels=['Continuum @ 940nm','889nm (CH4)','889/940 (CH4))']
    
            ch4abs=Methane_Panels(CH4file[ifile],CH4_RGB,mask,ax,ifile,CH4labels=CH4labels) #INCLUDE "ROW" HERE BASED ON nx
            ch4abs[ch4abs == inf] = 0
            ch4abs[ch4abs == -inf] = 0
            ch4abs[ch4abs == np.nan] = 0
            
            #nh3abs16bit = np.nan_to_num(((0.6667*65535*(ch4abs-0.5))*mask[:,:,0]).astype(np.uint16))
            #ch4abs16bit = np.nan_to_num(((5.*65535.*(ch4abs)-0.1)*mask[:,:,1]).astype(np.uint16))
            #imwrite(path+'/'+CH4file[0][0:26]+'CH4'+wvstr+'AbsPython.png', ch4abs16bit.astype(np.uint16))
            ###!!!!FIX CH4 FILE TIME IN NAMING!!!!

            if search("G620", CH4file[ifile]):
                if imagetype=="Map":
                    fnout=pathL2FITS+CH4file[0][0:26]+'620CH4AbsMap'
                elif imagetype=="Img":
                    fnout=pathL2FITS+CH4file[0][0:26]+'620CH4AbsImg'
                maskedch4=ch4abs*mask[:,:,1]
                normch4=maskedch4/maskednh3[maskedch4>0].mean()
                hdu = fits.PrimaryHDU(normch4.astype(np.float32))
                hdul = fits.HDUList([hdu])
                hdul[0].header['BITPIX']=-32
                print(hdul[0].header[:])
                #print(CH4file[ifile][0:26])
                try:
                    os.remove(fnout+'.fits')
                except: 
                    print("file doesn't exist")
                hdul.writeto(fnout+'.fits')
                hdul.close()

                normch4scaled=np.nan_to_num(((5.*65535.*(normch4*mask[:,:,1] - 0.9))))
                normch4scaled[normch4scaled<=0.]=0.0
                ch4abs16bit = normch4scaled.astype(np.uint16)
                imwrite(fnout+'.png', ch4abs16bit)#.astype(np.uint16))

    # CONTEXT PANELS (3)
    if NUVfile != 'NA' or RGBfile != 'NA':
        try:
            Contextlabels = sourcefiles[sourcedata]['Contextlabels']
        except KeyError:
            Contextlabels=['380nm','RGB']

        tmp=Context_Panels(ax,ny,NUVFile=NUVfile,NUV_RGB=NUV_RGB,RGBFile=RGBfile,RGB_RGB=RGB_RGB,mask=[],
                           Contextlabels=Contextlabels) #INCLUDE "ROW" HERE BASED ON nx


    for i in range(0,3):
        for j in range(0,3):
            print(i,j)
            ax[i,j].axis('off')

    fig.subplots_adjust(left=0.00, bottom=0.05, right=1.0, top=0.90,
                wspace=0.001, hspace=0.10)
    pathmaparraybyfilter=\
        'C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/Analysis Data/Map Array by Filter/'

    fig.savefig(pathmaparraybyfilter+obsdate+'-'+target+"_"+
                imagetype+"_SMHill.png",dpi=150,bbox_inches = 'tight')
    fig.savefig(pathmaparraybyfilter+obsdate+'-'+target+"_"+
                imagetype+"_SMHill.jpg",dpi=150,bbox_inches = 'tight')
    #imageio.imwrite(path+"PublicImageArray16bit"+obsdate+".png",fig.astype(np.uint16))
    
    #figx,axx=pl.subplots(1,1,figsize=(4,4), dpi=150, facecolor="black",
    #                      sharex=False)
    #ch4abs_scaled=(np.array(ch4abs)-1.0)/10.0+1.0
    #axx.imshow((np.array(nh3abs)/np.array(ch4abs_scaled))*np.array(mask[:,:,0]),'gist_gray')

###############################################################################
###############################################################################
    
def Ammonia_Panels(NH3file,NH3_RGB,mask,ax,
                   NH3labels=['656nm (Red Cont.)','647nm (NH3)','632nm (Blue Cont.)']):
    ###########################################################################
    # PURPOSE: COMPUTE AND PLOT PANELS FOR NH3 ANALYSIS
    # INPUTS:  NH3 FILE NAME, NH3_RGB PNG FILE, mask, and axis system to plot 
    #          on
    # OUTPUTS: FLOATING POINT ARRAYS OF THE COLOR SLOPE AND THE CONTINUUM
    #          DIVIDED 647NM AMMONIA OBSERVATIONS
    # ASSUMES: INPUT IS A NORMALIZED (WHITE BALANCED) RGB ARRAY WITH
    #          R(656), G(647), AND B(632)
    # IMPROVEMENTS:
    #          !!!!SPLIT COMPUTATIONS FROM PLOTTING INTO SEPARATE FUNCTIONS
    #          !!!!COMBINE FILE WRITING WITH COMPUTATIONS
    #          !!!!BETTER PLOT SCALING
    import numpy as np
    #from astropy.io import fits
    #import pylab as pl
    #import sharpen as sharp
    ###########################################################################
    # PARSE WINJUPOS TIME AND GET EPHEMERIS
    ###########################################################################
    sec=str(int(str(NH3file[16:17]))*6) #COMPUTE FROM FRACTIONAL WINJUPOS MINUTE
    NH3time=(NH3file[0:10]+"_"+NH3file[11:13]+":"+NH3file[13:15]+":"+sec.zfill(2))
    eph=get_WINJupos_ephem(NH3time)
    AmmoniaHeader=NH3time[11:19]+" UT; CM1"+eph[0]+"; CM2"+eph[1]+"; CM3"+eph[2]+"; Alt"+eph[3]   
    ###########################################################################
    # COMPUTE COLOR SLOPE BY PIXEL ASSUMING CHANNEL 0 IS R(656) AND CHANNEL 2
    #   IS B(632) AND THERE IS NO BULK COLOR SLOPE. THEN COMPUTE THE EFFECTIVE
    #   CONTINUUM AT 647NM (NH3) AND THE CONTINUUM DIVIDED NH3 IMAGE.
    ###########################################################################
    clrslp=(np.array(NH3_RGB[:,:,0]).astype(float)-np.array(NH3_RGB[:,:,2]).astype(float))/24.0 
    CNT647=15.0*clrslp+np.array(NH3_RGB[:,:,2])
    nh3abs=np.array(NH3_RGB[:,:,1])/CNT647
    ###########################################################################
    # PLOT 1ST ROW: CONTINUUM R(656), B(632), COLOR SLOPE AND ANNOTATE
    #   !!!!SCALING IS MANUAL HERE, PERHAPS IT SHOULD BE AUTOMATED?
    #   !!!!NOT SURE HOW SHARPENING IS FUNCTIONING HERE, IT DOESN'T SEEM TO 
    #       BE IMPORTED. MAYBE IT IS NATIVE TO imshow?
    ###########################################################################
    ax[0,0].imshow(sharpen(NH3_RGB[:,:,0]),'gist_gray',vmin=0.0,vmax=np.nanmax(NH3_RGB[:,:,0])*1.4)
    ax[0,0].annotate(NH3labels[0],[0.5,0.03],color='white',ha='center',
                    xycoords='axes fraction',fontsize=8)
    ax[0,1].imshow(sharpen(NH3_RGB[:,:,2])*mask[:,:,1],'gist_gray',vmin=0.0,vmax=np.nanmax(NH3_RGB[:,:,0])*1.4)
    ax[0,1].annotate(NH3labels[2],[0.5,0.03],color='white',ha='center',
                    xycoords='axes fraction',fontsize=8)
    # CUSTOM MASKING FOR COLOR SLOPE SINCE IT CAN BE POSITIVE AND NEGATIVE
    ax[0,2].imshow(clrslp*mask[:,:,1]+50*(mask[:,:,1]-1.0),'gist_gray',vmin=-50.,vmax=50.)
    ax[0,2].annotate('656/632 (Color Slope)',[0.5,0.03],color='white',ha='center',
                    xycoords='axes fraction',fontsize=8)
    ###########################################################################
    # PLOT 2ND ROW: 647NM SYNTH CONTINUUM, 647NM OBSERVATION, CONTINUUM DIVIDED
    #   647NM NH3 OBSERVATION AND ANNOTATE
    #   !!!!SCALING IS MANUAL HERE, PERHAPS IT SHOULD BE AUTOMATED?
    #   !!!!NOT SURE HOW SHARPENING IS FUNCTIONING HERE, IT DOESN'T SEEM TO 
    #       BE IMPORTED. MAYBE IT IS NATIVE TO imshow?
    ###########################################################################       
    ax[1,0].imshow(sharpen(CNT647)*mask[:,:,1],'gist_gray',vmin=0.0,vmax=np.nanmax(NH3_RGB[:,:,0])*1.4)
    ax[1,0].annotate('Synth. Continuum @ 647nm',[0.5,0.03],color='white',ha='center',
                    xycoords='axes fraction',fontsize=8)
    ax[1,1].imshow(sharpen(NH3_RGB[:,:,1]),'gist_gray',vmin=0.0,vmax=np.nanmax(NH3_RGB[:,:,0])*1.4)
    ax[1,1].annotate(NH3labels[1],[0.5,0.03],color='white',ha='center',
                    xycoords='axes fraction',fontsize=8)
    ax[1,2].imshow(nh3abs*mask[:,:,1],'gist_gray',vmin=0.9,vmax=1.1)
    ax[1,2].annotate('647/Cont. (NH3)',[0.51,0.03],color='white',ha='center',
                    xycoords='axes fraction',fontsize=8)
    ###########################################################################
    # ANNOTATE CONTINUUM ROW AND AMMONIA ROW OF PANELS
    ###########################################################################
    ax[0,1].annotate("CONTINUUM",[-1.2,0.93],color='white',ha='left',
                    xycoords='axes fraction',fontsize=9)
    ax[0,1].annotate(AmmoniaHeader,[2.2,0.93],color='white',ha='right',
                    xycoords='axes fraction',fontsize=8)
    ax[1,1].annotate(r"AMMONIA (NH$_3$)",[-1.2,0.93],color='white',ha='left',
                    xycoords='axes fraction',fontsize=9)
    ax[1,1].annotate(AmmoniaHeader,[2.2,0.93],color='white',ha='right',
                    xycoords='axes fraction',fontsize=8)
    ax[0,1].set_zorder(1) #required so not blocked by the ax[0,2] image
    ax[1,1].set_zorder(1) #required so not blocked by the ax[0,2] image
    #see https://stackoverflow.com/questions/29735743/getting-text-to-display-in-front-of-subplot-images
    
    return clrslp,nh3abs
    
def Methane_Panels(CH4file,CH4_RGB,mask,ax,ifile,
                   CH4labels=['656nm (Cont.)','889nm (CH4)','889/Cont. (CH4))']):
    import numpy as np
    #import sharpen as sharp
    ###########################################################################
    # PARSE WINJUPOS TIME AND GET EPHEMERIS
    #   !!!!COULD ROLL A TIME FORMAT OPTION INTO THE EPHEMERIS ROUTINE AND
    #       ALSO HAVE IT RETURN THE HEADER AS A GENERIC HEADER FOR ANY FILE
    ###########################################################################
    sec=str(int(str(CH4file[16:17]))*6)
    CH4time=(CH4file[0:10]+"_"+CH4file[11:13]+":"+CH4file[13:15]+":"+sec.zfill(2))
    eph=get_WINJupos_ephem(CH4time)
    MethaneHeader=CH4time[11:19]+" UT; CM1"+eph[0]+"; CM2"+eph[1]+"; CM3"+eph[2]+"; Alt"+eph[3]
    ###########################################################################
    # COMPUTE COLOR SLOPE BY PIXEL ASSUMING IF CH4 CHANNEL IS EITHER 620NM OR
    #   730NM. 620NM IS CLEARLY ADJACENT TO 632NM AND 656NM CONTINUA AND
    #   *SHOULD* BE REASONABLE TO EXTRAPOLATE COLOR SLOPE TO CREATE A SYNTHETIC
    #   CONTINUUM THERE. HOWEVER, TESTS SHOULD BE DONE TO SEE IF THIS APPROACH
    #   REMAINS VALID FOR 730NM. FOR 889NM, THERE IS NO COLOR SLOPE AVAILABLE
    #   SO THE CONTINUUM REFERENCE IS A NORMALIZED 940NM IMAGE. THE CONTINUUM
    #   DIVIDED METHANE IMAGE IS THEN COMPUTED.
    #   !!!!COLOR SLOPE, WHERE COMPUTED SHOULD BE RETURNED AND ALSO WRITTEN
    #       AS FITS AND PNG FILES, EVEN IF JUST FOR VALIDATION.
    ###########################################################################
    if CH4labels[1]=='889nm (CH4)':
        ch4abs=(np.array(CH4_RGB[:,:,1])+0.0001)/(np.array(CH4_RGB[:,:,0])+0.0001)
        ax[2+ifile,0].imshow(sharpen(CH4_RGB[:,:,0]),'gist_gray',vmin=0.0,vmax=np.nanmax(CH4_RGB[:,:,0])*1.2)
        ax[2+ifile,1].imshow(CH4_RGB[:,:,1],'gist_gray',vmin=0.0,vmax=np.nanmax(CH4_RGB[:,:,0])*1.2)
        maxmin=[0.5,2.0]
    elif CH4labels[1]=='730nm (CH4)':# or CH4labels[1]=='620nm (CH4)':
        clrslp=(np.array(CH4_RGB[:,:,0]).astype(float)-np.array(CH4_RGB[:,:,2]).astype(float))/24.0 
        meta={'620nm (CH4)':{'dwave':-12,'maxmin':[0.9,1.1]},
              '730nm (CH4)':{'dwave':98,'maxmin':[0.8,1.2]}}
        CNTSynth=meta[CH4labels[1]]['dwave']*clrslp+np.array(CH4_RGB[:,:,2])
        ch4abs=(np.array(CH4_RGB[:,:,1])+0.0001)/(CNTSynth+0.0001)
        ax[2+ifile,0].imshow(sharpen(CNTSynth)*mask[:,:,1],'gist_gray',vmin=0.0,vmax=np.nanmax(CNTSynth)*1.2)
        ax[2+ifile,1].imshow(sharpen(CH4_RGB[:,:,1])*mask[:,:,1],'gist_gray',vmin=0.0,vmax=np.nanmax(CH4_RGB[:,:,0])*1.2)
        maxmin=meta[CH4labels[1]]['maxmin']
    elif CH4labels[1]=='620nm (CH4)': ###ADDITION FOR PIECEWISE CONTINUOUS RATHER THAN EXTRAPOLATED
        meta={'620nm (CH4)':{'dwave':-12,'maxmin':[0.9,1.1]},
              '730nm (CH4)':{'dwave':98,'maxmin':[0.8,1.2]}}
        CNTSynth=np.array(CH4_RGB[:,:,2])
        ch4abs=(np.array(CH4_RGB[:,:,1])+0.0001)/(CNTSynth+0.0001)
        ax[2+ifile,0].imshow(sharpen(CNTSynth)*mask[:,:,1],'gist_gray',vmin=0.0,vmax=np.nanmax(CNTSynth)*1.2)
        ax[2+ifile,1].imshow(sharpen(CH4_RGB[:,:,1])*mask[:,:,1],'gist_gray',vmin=0.0,vmax=np.nanmax(CH4_RGB[:,:,0])*1.2)
        maxmin=meta[CH4labels[1]]['maxmin']

    ax[2+ifile,2].imshow(ch4abs*mask[:,:,1],'gist_gray',vmin=maxmin[0],vmax=maxmin[1])
    
    for icol in range(0,3):
        ax[2+ifile,icol].annotate(CH4labels[icol],[0.5,0.03],color='white',ha='center',
                        xycoords='axes fraction',fontsize=8)
    
    ax[2+ifile,1].annotate(r"METHANE (CH$_4$)",[-1.2,0.93],color='white',ha='left',
                    xycoords='axes fraction',fontsize=9)
    ax[2+ifile,1].annotate(MethaneHeader,[2.2,0.93],color='white',ha='right',
                    xycoords='axes fraction',fontsize=8)
    ax[2+ifile,1].set_zorder(1) #required so not blocked by the ax[2,2] image
    #see https://stackoverflow.com/questions/29735743/getting-text-to-display-in-front-of-subplot-images
    return ch4abs

def Context_Panels(ax,ny,NUVFile='NA',NUV_RGB=[],RGBFile='NA',RGB_RGB=[],mask=[],
                   Contextlabels=['380nm','RGB']):
    import numpy as np
   
    if NUVFile != 'NA': 
        sec=str(int(str(NUVFile[16:17]))*6)
        NUVtime=(NUVFile[0:10]+"_"+NUVFile[11:13]+":"+NUVFile[13:15]+":"+sec.zfill(2))
        eph=get_WINJupos_ephem(NUVtime)
        if np.array(mask).size != 0:
            ax[ny-1,0].imshow(np.array(NUV_RGB[:,:,0])*np.array(mask[:,:,0]),'gist_gray')
        else:
            ax[ny-1,0].imshow(NUV_RGB[:,:,0],'gist_gray')
        ax[ny-1,0].annotate(Contextlabels[0]+'; '+NUVtime[11:19]+" UT; Alt"+eph[3],[0.5,0.03],color='white',ha='center',
                        xycoords='axes fraction',fontsize=8)
        ax[ny-1,0].annotate("CM1"+eph[0]+"; CM2"+eph[1]+"; CM3"+eph[2],
                         [0.5,-0.05],color='white',ha='center',
                        xycoords='axes fraction',fontsize=7)
        ax[ny-1,1].imshow(np.zeros(NUV_RGB.shape))
        ax[ny-1,2].imshow(np.zeros(NUV_RGB.shape))
        
    print("RGBFile=",RGBFile)
    if RGBFile != 'NA': 
        sec=str(int(str(RGBFile[16:17]))*6)
        RGBtime=(RGBFile[0:10]+"_"+RGBFile[11:13]+":"+RGBFile[13:15]+":"+sec.zfill(2))
        eph=get_WINJupos_ephem(RGBtime)
        ax[ny-1,1].imshow(RGB_RGB)
        ax[ny-1,1].annotate(Contextlabels[1]+'; '+RGBtime[11:19]+" UT; Alt"+eph[3],[0.5,0.03],color='white',ha='center',
                        xycoords='axes fraction',fontsize=8)
        ax[ny-1,1].annotate("CM1"+eph[0]+"; CM2"+eph[1]+"; CM3"+eph[2],
                         [0.5,-0.05],color='white',ha='center',
                        xycoords='axes fraction',fontsize=7)
    
        ax[ny-1,0].imshow(np.zeros(RGB_RGB.shape))
        ax[ny-1,2].imshow(np.zeros(RGB_RGB.shape))
    
    ContextHeader="CONTEXT"
    ax[ny-1,1].annotate(ContextHeader,[-1.2,0.93],color='white',ha='left',
                    xycoords='axes fraction',fontsize=9)

    ax[ny-1,2].annotate("NOTES",[0.5,0.93],color='white',ha='center',
                     xycoords='axes fraction',fontsize=8)
    ax[ny-1,2].annotate("Ammonia and methane \n"+
                     "absorption are explored \n"+
                     "using in-band images \n"+
                     "divided by continuum \n"+
                     "images. Computed ratio \n"
                     "images are not sharpened \n"
                     "to minimize artifacts.",
                     [0.0,0.85],color='white',ha='left',va='top',
                     xycoords='axes fraction',fontsize=7)

    ax[ny-1,1].set_zorder(1) #required so not blocked by the ax[2,2] image
    #see https://stackoverflow.com/questions/29735743/getting-text-to-display-in-front-of-subplot-images 
    
def load_png(file_path):
    """
    Purpose: Properly load a 48-bit PNG file
    Read from KITTI .png file
    Args:
        file_path string: file path(absolute)
    Returns:
        data (numpy.array): data of image in (Height, Width, 3) layout
    
    FROM: https://www.programcreek.com/python/example/98900/png.Reader
    """
    import png
    import numpy as np

    flow_object = png.Reader(filename=file_path)
    flow_direct = flow_object.asDirect()
    flow_data = list(flow_direct[2])
    (w, h) = flow_direct[3]['size']

    flow = np.zeros((h, w, 3), dtype=np.float64)
    for i in range(len(flow_data)):
        flow[i, :, 0] = flow_data[i][0::3]
        flow[i, :, 1] = flow_data[i][1::3]
        flow[i, :, 2] = flow_data[i][2::3]

    return flow.astype(np.uint16) 

def get_WINJupos_ephem(dateobs):
    # Example call: get_WINJupos_ephem('2021-09-05_04:09:00')
    import win32com.shell.shell as shell
    import time
    #shell.ShellExecuteEx(lpVerb='runas', lpFile='cmd.exe', lpParameters='/c '+commands) #run as admin
    ###########################################################################
    # WRITE *.BAT FILE WITH COMMAND SCRIPT FOR WINJUPOS
    ###########################################################################   
    batfile=open("WINJupos_CM.bat",'w')
    Line1='cd "\Program Files\WinJUPOS 12.1.1"\r\n'
    #Line2='WinJUPOS.x64.exe Jupiter /GetCM:2021-09-05_04:09:00 /GeoLong:-104.9 /GeoLat:39.7 /GetAlt >"c:\Astronomy\Projects\SAS 2021 Ammonia\Jupiter_NH3_Analysis_P3\cm.txt"\r\n'
    Line2='WinJUPOS.x64.exe Jupiter /GetCM:'+dateobs+' /GeoLong:-104.9 /GeoLat:39.7 /GetAlt >"c:\Astronomy\Projects\SAS 2021 Ammonia\Jupiter_NH3_Analysis_P3\cm.txt"\r\n'
    #batfile.writelines(["Line1\r\n","Line2\r\n"])
    batfile.writelines([Line1,Line2])
    batfile.close()
    ###########################################################################
    # EXECUTE *.BAT COMMAND FILE FOR WINJUPOS AND WAIT TO READ RESULT FILE
    ###########################################################################
    commands = "WINJupos_CM.bat"  
    shell.ShellExecuteEx(lpFile='cmd.exe', lpParameters='/c '+commands)
    time.sleep(1)
    ephemfile=open("C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/cm.txt",'r')
    Lines=ephemfile.readlines()
    ephemfile.close()
    LineString=str(Lines[0])
    print(LineString)
    ###########################################################################
    # PARSE OUTPUT FILE AND CREATE STRING ARRAY EPH FOR CM1, CM2, CM3, AND ALT
    ###########################################################################
    start=[i for i, letter in enumerate(LineString) if letter == "="]
    end=[i for i, letter in enumerate(LineString) if letter == "°"]
    eph=[]
    for i in range(0,3):
        temp=LineString[int(start[i])+1:int(end[i])]
        #print("CM"+str(i+1)+" = "+LineString[int(start[i])+1:int(end[i])])#,Linestring[start[1]:end[1]],Linestring[start[2]:end[2]])
        print("CM"+str(i+1)+" = "+temp)#,Linestring[start[1]:end[1]],Linestring[start[2]:end[2]])
        eph.extend([temp])        
    temp=LineString[int(start[3])+1:int(end[3])]
    print("Alt =  "+temp)
    eph.extend([temp])        
    return eph

def sharpen(image):
    from astropy.convolution import Gaussian2DKernel
    from astropy.convolution import convolve
    kernel = Gaussian2DKernel(x_stddev=5)
    blurred=convolve(image,kernel)
    tst=image+0.99*(image-blurred) #Need to revalidate that this is correct
    return tst


    
