# -*- coding: utf-8 -*-
"""
Created on Mon Sep  4 19:33:21 2023

@author: smhil
"""

def get_variations_list():
    """
    Created on Tue Aug 22 10:37:40 2023
    
    @author: smhil
    """
    import json
    sourcefiles={'20220919UTc':{'Telescope':'VLT','FL':'5600mm','Camera':'ASI120MM',
                                     'Seeing':'6/10','Transparency':'7/10','Variation':'620Gauss2pix', 
                               'CH4file':['2022-09-19-0352_3-Jupiter-VLT-R656G620B632-RGB-WhtBal-MedFilt-620Gauss2pix'],
                               'NH3file':'2022-09-19-0352_3-Jupiter-VLT-R656G647B632-RGB-WhtBal-MedFilt',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2022-09-19-0352_3-Jupiter-MUSE-R650G550B450-RGB-WhtBal-Filtered'},
                               'Contextlabels':['889nm','']},
                  '20220919UTd':{'Telescope':'VLT','FL':'5600mm','Camera':'ASI120MM',
                                     'Seeing':'6/10','Transparency':'7/10','Variation':'647Gauss2pix', 
                               'CH4file':['2022-09-19-0352_3-Jupiter-VLT-R656G620B632-RGB-WhtBal-MedFilt'],
                               'NH3file':'2022-09-19-0352_3-Jupiter-VLT-R656G647B632-RGB-WhtBal-MedFilt-647Gauss2pix',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2022-09-19-0352_3-Jupiter-MUSE-R650G550B450-RGB-WhtBal-Filtered'},
                               'Contextlabels':['889nm','']},
                  '20220919UTe':{'Telescope':'VLT','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'6/10','Transparency':'7/10','Variation':'632Gauss2pix', 
                               'CH4file':['2022-09-19-0352_3-Jupiter-VLT-R656G620B632-RGB-WhtBal-MedFilt-632Gauss2pix'],
                               'NH3file':'2022-09-19-0352_3-Jupiter-VLT-R656G647B632-RGB-WhtBal-MedFilt-632Gauss2pix',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2022-09-19-0352_3-Jupiter-MUSE-R650G550B450-RGB-WhtBal-Filtered'},
                               'Contextlabels':['889nm','']},
                  '20220919UTf':{'Telescope':'VLT','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'6/10','Transparency':'7/10','Variation':'656Gauss2pix', 
                               'CH4file':['2022-09-19-0352_3-Jupiter-VLT-R656G620B632-RGB-WhtBal-MedFilt-656Gauss2pix'],
                               'NH3file':'2022-09-19-0352_3-Jupiter-VLT-R656G647B632-RGB-WhtBal-MedFilt-656Gauss2pix',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2022-09-19-0352_3-Jupiter-MUSE-R650G550B450-RGB-WhtBal-Filtered'},
                               'Contextlabels':['889nm','']},
                  '20220919UTg':{'Telescope':'VLT','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'6/10','Transparency':'7/10','Variation':'656Gauss2pix', 
                               'CH4file':['2022-09-19-0352_3-Jupiter-VLT-R656G620B632-RGB-WhtBal-MedFilt-656Gauss2pix'],
                               'NH3file':'2022-09-19-0352_3-Jupiter-VLT-R656G647B632-RGB-WhtBal-MedFilt-656Gauss2pix',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2022-09-19-0352_3-Jupiter-MUSE-R650G550B450-RGB-WhtBal-Filtered'},
                               'Contextlabels':['889nm','']},
                  '20220919UTh':{'Telescope':'VLT','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'6/10','Transparency':'7/10','Variation':'656Gauss2pix', 
                               'CH4file':['2022-09-19-0352_3-Jupiter-VLT-R656G620B632-RGB-WhtBal-MedFilt-656Gauss2pix'],
                               'NH3file':'2022-09-19-0352_3-Jupiter-VLT-R656G647B632-RGB-WhtBal-MedFilt-656Gauss2pix',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2022-09-19-0352_3-Jupiter-MUSE-R650G550B450-RGB-WhtBal-Filtered'},
                               'Contextlabels':['889nm','']},
                  '20220919UTi':{'Telescope':'VLT','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'6/10','Transparency':'7/10','Variation':'656Gauss2pix', 
                               'CH4file':['2022-09-19-0352_3-Jupiter-VLT-R656G620B632-RGB-WhtBal-MedFilt-656Gauss2pix'],
                               'NH3file':'2022-09-19-0352_3-Jupiter-VLT-R656G647B632-RGB-WhtBal-MedFilt-656Gauss2pix',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2022-09-19-0352_3-Jupiter-MUSE-R650G550B450-RGB-WhtBal-Filtered'},
                               'Contextlabels':['889nm','']},
                  '20220919UTj':{'Telescope':'VLT','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'6/10','Transparency':'7/10','Variation':'656Gauss2pix', 
                               'CH4file':['2022-09-19-0352_3-Jupiter-VLT-R656G620B632-RGB-WhtBal-MedFilt-656Gauss2pix'],
                               'NH3file':'2022-09-19-0352_3-Jupiter-VLT-R656G647B632-RGB-WhtBal-MedFilt-656Gauss2pix',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2022-09-19-0352_3-Jupiter-MUSE-R650G550B450-RGB-WhtBal-Filtered'},
                               'Contextlabels':['889nm','']},
                  '20220919UTk':{'Telescope':'VLT','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'6/10','Transparency':'7/10','Variation':'656Gauss2pix', 
                               'CH4file':['2022-09-19-0352_3-Jupiter-VLT-R656G620B632-RGB-WhtBal-MedFilt-656Gauss2pix'],
                               'NH3file':'2022-09-19-0352_3-Jupiter-VLT-R656G647B632-RGB-WhtBal-MedFilt-656Gauss2pix',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2022-09-19-0352_3-Jupiter-MUSE-R650G550B450-RGB-WhtBal-Filtered'},
                               'Contextlabels':['889nm','']},
                  '20220919UTl':{'Telescope':'VLT','FL':'5600mm','Camera':'ASI120MM',
                                            'Seeing':'6/10','Transparency':'7/10','Variation':'656Gauss2pix', 
                               'CH4file':['2022-09-19-0352_3-Jupiter-VLT-R656G620B632-RGB-WhtBal-MedFilt-656Gauss2pix'],
                               'NH3file':'2022-09-19-0352_3-Jupiter-VLT-R656G647B632-RGB-WhtBal-MedFilt-656Gauss2pix',
                               #'CH4labels':['Synth. Continuum @ 620nm','620nm (CH4)','620/Cont. (CH4))'],
                               'Context':{'NUVfile':'NA',
                                          'RGBfile':'2022-09-19-0352_3-Jupiter-MUSE-R650G550B450-RGB-WhtBal-Filtered'},
                               'Contextlabels':['889nm','']},
                  
                  }

    filename="C:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis_P3/L1PNG.json"
    with open(filename, "w") as fp:
        json.dump(sourcefiles , fp) 

    return(sourcefiles)