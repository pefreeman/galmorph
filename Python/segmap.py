#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  9 16:29:51 2018

@author: peterfreeman
"""

import numpy as np
from scipy import ndimage

def segmap(img,numSigma=5):
    nx, ny = np.shape(img)
    xcen = int(nx/2)
    ycen = int(ny/2)
    
    labeled, nr_objects = ndimage.label(img<0)
    if np.max(labeled) == 0:
        r,c = np.shape(img)
        return np.zeros((r,c)),np.zeros((r,c)),-9.0
    
    r,c = np.where(labeled>0)
    sigmaHat = np.sqrt(np.sum(np.square(img[r,c]))/len(r))
        
    labeled, nr_objects = ndimage.label(img>numSigma*sigmaHat)
    if labeled[xcen,ycen] == 0:
        r,c = np.shape(img)
        return np.zeros((r,c)),np.zeros((r,c)),-9.0
        
    allc = np.copy(labeled)
    r,c = np.where(allc>0)
    allc[r,c] = 1
    
    r,c = np.where(labeled != labeled[xcen,ycen])
    labeled[r,c] = 0
    labeled = labeled/np.max(labeled)
    istop = 0
    while istop == 0:
        fill = np.zeros((nx,ny))
        for jj in range(1,nx-1):
            for kk in range(1,ny-1):
                if labeled[jj,kk] == 0:
                    if np.sum(labeled[range(jj-1,jj+2),range(kk-1,kk+2)]) > 4:
                        fill[jj,kk] = 1
        if np.sum(fill) == 0:
            istop = 1
        labeled = labeled+fill
    return labeled,allc,sigmaHat