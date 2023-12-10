def makeQR(fn,URL):
    """
    Created on Thu Dec  7 11:25:49 2023
    
    @author: smhil
    """
    import qrcode
    img = qrcode.make(URL)
    type(img)
    img.save("C:/Astronomy/Projects/SAS 2021 Ammonia/Publications/"+fn+".png")