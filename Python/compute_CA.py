#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 14:57:55 2018

@author: peterfreeman
"""

import math
import numpy as np

def distCircle(image,xc,yc):
    r,c = np.shape(image)
    aimage = np.full((r,c),-9.0)
    for ii in range(r):
        dx = ii-xc
        for jj in range(c):
            dy = jj-yc
            aimage[ii,jj] = np.sqrt(np.square(dx)+np.square(dy))
    return aimage

def petRadCir2(image,xc,yc,eta=0.2):
    r,c = np.shape(image)
    n = np.ceil(np.sqrt(np.square(r/2)+np.square(c/2))).astype(int)
    cir = distCircle(image,xc,yc)
    rPet = np.zeros(n)
    
    # compute curve of growth
    for r in range(n-1):
        if r == 0:
            annx,anny = np.where(cir<r+0.5)
        else:
            annx,anny = np.where((cir>=r-0.5)&(cir<r+0.5))
        if len(annx) == 0:
            rPet[r] = 1
            continue
        mu = np.mean(image[annx,anny])
        if mu < 0:
            mu = 0
        totx,toty = np.where(cir<r+0.5)
        avgMu = np.mean(image[totx,toty])
        if avgMu <= 0:
            avgMu = 0
            if r > 0:
                rPet[r] = rPet[r-1]
            else:
                rPet[r] = 1
        else:
            rPet[r] = mu/avgMu
        if rPet[r] < eta:
            break
    
    rCir = np.interp(eta,rPet[::-1],np.arange(n,dtype='float')[::-1])
    return cir,rCir

def fluxCir2(image,imgCir,rMax):
    rhi = np.ceil(rMax).astype(int)
    sumInt = np.zeros(rhi)
    aver   = np.zeros(rhi)
    
    for r in range(rhi):
        if r == 0:
            annx,anny = np.where(imgCir<r+0.5)
            sumInt[r] = 0
            if np.sum(image[annx,anny]) > 0:
                sumInt[r] = np.sum(image[annx,anny])
        else:
            annx,anny = np.where((imgCir>=r-0.5)&(imgCir<r+0.5))
            tmp = 0
            if np.sum(image[annx,anny]) > 0:
                tmp = np.sum(image[annx,anny])
            sumInt[r] = sumInt[r-1] + tmp
        if len(annx) > 0:
            aver[r] = np.mean(imgCir[annx,anny])
        else:
            aver[r] = 0
    
    # 3/25/14: add zero-point to avoid NA output for very concentrated galaxies
    if sumInt[rhi-1] == 0:  # may not be necessary given python indexing
        return -9.0,-9.0
    
    sumInt = sumInt/sumInt[rhi-1]  # used to add zero at beginning...
    #aver = c(0,aver)            # ditto...again, python indexing may make moot
    
    r20 = np.interp(0.20,sumInt,aver)
    r80 = np.interp(0.80,sumInt,aver)
    return r20,r80
    
def CStatistic(oimg,cimg,rPet):
    rMax = 1.5*rPet
    if rMax > np.floor(np.min(np.shape(oimg))*np.sqrt(2)/2).astype(int):
        rMax = np.floor(np.min(np.shape(oimg))*np.sqrt(2)/2).astype(int)
    r20,r80 = fluxCir2(oimg,cimg,rMax)
    if math.isnan(r20) == True or math.isnan(r80) == True or r20 <= 0.0 or r80 <= 0.0: 
        return -9.0
    return 5*np.log10(r80/r20)

def AStatistic(oimg,cimg,smapAll,rPet):
    rMax = 1.5*rPet
    if rMax > np.floor(np.min(np.shape(oimg))*np.sqrt(2)/2).astype(int):
        rMax = np.floor(np.min(np.shape(oimg))*np.sqrt(2)/2).astype(int)
    # Rotate 180 degrees
    rimg = np.rot90(np.rot90(oimg))
    wx,wy = np.where(cimg>rMax)
    cimg[wx,wy] = 0
    wx,wy = np.where(cimg>0)
    cimg[wx,wy] = 1
    ngal = len(wx)
    ogal = cimg*oimg
    rgal = cimg*rimg
    
    Aden = np.sum(np.abs(ogal))
    Agal = np.sum(np.abs(ogal-rgal))
    cimg = cimg-1
    cimg = np.abs(cimg)
    oann = cimg*oimg*np.abs(smapAll-1)
    wx,wy = np.where(oann!=0)
    if len(wx) > 0:
        rann = np.rot90(np.rot90(oann))
        wx,wy = np.where((oann!=0)&(rann!=0))
        if len(wx) > 0:
            nbkg = len(wx)
            Abkg = np.sum(np.abs(oann[wx,wy]-rann[wx,wy]))
        else:
            nbkg = 1
            Abkg = 0
    else:
        nbkg = 1
        Abkg = 0
    return (Agal-(ngal/nbkg)*Abkg)/(2*Aden) # factor of 2 from Conselice et al. 2000
