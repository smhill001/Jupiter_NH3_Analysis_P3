# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 12:30:42 2023

@author: smhil
"""
def make_perlin_noise(dim=100,octv=10,seed=1):
    import matplotlib.pyplot as plt
    from perlin_noise import PerlinNoise
    
    noise = PerlinNoise(octaves=octv, seed=20)
    xpix, ypix = dim, dim
    pic = [[noise([i/xpix, j/ypix]) for j in range(xpix)] for i in range(ypix)]
    
    plt.imshow(pic, cmap='gray')
    plt.show()
    return(pic)