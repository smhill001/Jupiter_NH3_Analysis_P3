def get_spice_ephem(dateobs,planet="Jupiter"):
    """
    Created on Fri Sep 15 07:57:17 2023

    @author: smhil
    """
    #Example call: get_spice_ephem('2021-09-05_04:09:00')
    import planetmapper as pm
    import convert_system3_to_I_II_spice as convert
    
    new_dateobs=dateobs.replace("_"," ")
    
    body = pm.BodyXY('Jupiter', new_dateobs)
    CM3=body.subpoint_lon
    ephem=convert.convert_system3_to_I_II_spice(new_dateobs, CM3)

    ###########################################################################
    # PARSE OUTPUT FILE AND CREATE STRING ARRAY EPH FOR CM1, CM2, CM3, AND ALT
    ###########################################################################
    eph=[str(ephem["System I"]),str(ephem["System II"]),str(ephem["System III"])]      
    #print("Alt =  "+temp)
    return eph