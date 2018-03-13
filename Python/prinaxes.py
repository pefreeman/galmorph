#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 15:38:12 2018

@author: peterfreeman
"""

import numpy as np

#
# assumes image is *normalized to 1*
#
def prinaxes(img):
    x = np.arange(np.shape(img)[0])
    y = np.arange(np.shape(img)[1])
    marx = np.sum(img,axis=1)
    mary = np.sum(img,axis=0)
    w = np.where(marx<0)
    marx[w] = 0
    w = np.where(mary<0)
    mary[w] = 0
    xcen = np.sum(x*marx)
    ycen = np.sum(y*mary)
    xvar = np.sum(x*x*marx)-np.square(xcen)
    yvar = np.sum(y*y*mary)-np.square(ycen)
    xycov = np.sum(np.outer(x,y)*img)-xcen*ycen
    
    b = -xvar-yvar
    c = xvar*yvar-np.square(xycov)
    sd = np.sqrt(np.square(b)-4*1*c)
    axmin = np.sqrt(0.5*(-b-sd))
    axmax = np.sqrt(0.5*(-b+sd))
    rad2deg = 180.0/np.pi
    a = xvar-np.square(axmax)
    b = xycov
    c = yvar-np.square(axmax)
    
    angle = np.arctan((a-b)/(c-b))*rad2deg
    return axmax,axmin,angle
