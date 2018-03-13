#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 14:42:50 2018

@author: peterfreeman
"""

import numpy as np

def genCStatistic(img):
    r,c = np.where(img!=0)
    intensity = np.sort(img[r,c].flatten())[::-1]
    
    totIntensity = np.sum(intensity)
    
    sumIntensity  = 0
    for ii in range(len(intensity)):
        sumIntensity += intensity[ii]
        if sumIntensity/totIntensity > 0.20:
            a20 = ii+1
            break
        
    for ii in range(a20,len(intensity)):
        sumIntensity += intensity[ii]
        if sumIntensity/totIntensity > 0.80:
            a80 = ii+1
            break
             
    return 5*np.log10(np.sqrt(a80/a20))
    
        
        