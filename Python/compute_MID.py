#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 10:12:16 2018

@author: peterfreeman
"""

import numpy as np
from scipy import ndimage
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve

def defineClump(img): # it really seems like there should be a python module code to do this...
    r,c = np.shape(img)
    clump = np.full((r,c),-1,dtype=int)
    xPeak = []
    yPeak = []
    for jj in range(r):
        for kk in range(c):
            if img[jj,kk] != 0:
                jjcl = jj
                kkcl = kk
                istop = 0
                while istop == 0:
                    jjmax = jjcl
                    kkmax = kkcl
                    imgmax = img[jjcl,kkcl]
                    for mm in np.arange(-1,2):
                        if jjcl+mm >= 0 and jjcl+mm < r:
                            for nn in np.arange(-1,2):
                                if kkcl+nn >=0 and kkcl+nn < c:
                                    if img[jjcl+mm,kkcl+nn] > imgmax:
                                        imgmax = img[jjcl+mm,kkcl+nn]
                                        jjmax = jjcl+mm
                                        kkmax = kkcl+nn
                    if jjmax == jjcl and kkmax == kkcl:
                        n = len(xPeak)
                        if n > 0:
                            cltmp = -1
                            for pp in range(n):
                                if xPeak[pp] == jjmax and yPeak[pp] == kkmax:
                                    cltmp = pp
                            if cltmp >= 0:
                                clump[jj,kk] = cltmp
                            else:
                                xPeak.append(jjmax)
                                yPeak.append(kkmax)
                                clump[jj,kk] = len(xPeak)-1
                        else:
                            xPeak.append(jjmax)
                            yPeak.append(kkmax)
                            clump[jj,kk] = len(xPeak)-1
                        istop = 1
                    else:
                        jjcl = jjmax
                        kkcl = kkmax
    return clump,xPeak,yPeak
              
def MStatistic(img,levels=np.arange(0,0.99,0.025)):
    r,c = np.shape(img)
    rw,_ = np.where(img==0)
    nPix = r*c-len(rw)
    nImg = img/np.max(img)
    nLev = len(levels)
    areaRatio = np.zeros(nLev)
    
    r,c = np.where(nImg!=0)
    v = np.sort(nImg[r,c].flatten())
    
    for ii in range(nLev):
        thr = round((nPix-1)*levels[ii]).astype(int)
        if thr == 0:
            continue
        labeled, nr_objects = ndimage.label(nImg>=v[thr])
        if np.max(labeled) == 0:
            continue
        areas = np.zeros(nr_objects)
        if nr_objects > 1:
            for jj in range(nr_objects):
                r,_ = np.where(labeled==jj+1)
                areas[jj] = len(r)
            areas = np.sort(areas)[::-1]
            areaRatio[ii] = (areas[1]/areas[0])*(areas[1]/nPix)
    return np.max(areaRatio)
    
def IStatistic(img,scale=0):
    if scale > 0:
        kernel = Gaussian2DKernel(scale)
        smoothImg = convolve(img,kernel)
    else:
        smoothImg = np.copy(img)
    clump,xPeak,yPeak = defineClump(smoothImg) # check for gradient ascent algorithms
    n = len(xPeak)
    if n == 1:
        return 0.0,xPeak[0],yPeak[0]
    else:
        clumpIntensity = np.zeros(n)
        r,c = np.shape(img)
        for jj in range(r):
            for kk in range(c):
                if clump[jj,kk] >= 0:    # beware indexing here!!
                    clumpIntensity[clump[jj,kk]] += img[jj,kk]
        w = np.argmax(clumpIntensity)
        xPeak = xPeak[w]
        yPeak = yPeak[w]
        s = np.sort(clumpIntensity)[::-1]
        return s[1]/s[0],xPeak,yPeak
        
def DStatistic(img,xPeak,yPeak):
    if xPeak == None or yPeak == None:
        return -9
    tot = xCen = yCen = 0
    r,c = np.shape(img)
    for jj in range(r):
        for kk in range(c):
            if img[jj,kk] > 0:
                xCen += jj*img[jj,kk]
                yCen += kk*img[jj,kk]
                tot  += img[jj,kk]
    xCen /= tot
    yCen /= tot
    r,_ = np.where(img!=0)
    area = len(r)
    deviation = np.sqrt(np.square(xPeak-xCen)+np.square(yPeak-yCen))/np.sqrt(area/np.pi)
    return deviation,xCen,yCen