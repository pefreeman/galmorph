#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 09:15:15 2018

@author: peterfreeman
"""

import numpy as np

# based on equation 3 of Lotz, Primack, & Madau (2004)
def GStatistic(img):
    r,c = np.where(img!=0)
    v = np.sort(img[r,c].flatten())
    n = len(v)
    coeff = 2*np.arange(1,n+1)-n-1
    return np.sum(coeff*v)/np.mean(v)/n/(n-1)

# based on equations 7 and 8 of Lotz, Primack, & Madau (2004)
# this can probably be made faster/more efficient via use of "apply"
def M20Statistic(img):
    # Technically: "the center is computed by finding x_c, y_c such that M_tot is minimized."
    # We do not do this.
    nx,ny = np.shape(img)
    xi = range(nx)
    yi = range(ny)
    xcnum = 0
    xcden = 0
    for ii in xi:
        xcnum += np.sum(img[ii,:]*ii)
        xcden += np.sum(img[ii,:])
    xc = xcnum/xcden
    ycnum = 0
    ycden = 0
    for ii in yi:
        ycnum += np.sum(img[:,ii]*ii)
        ycden += np.sum(img[:,ii])
    yc = ycnum/ycden
    mtot = 0
    for ii in xi:
        for jj in yi:
            mtot += img[ii,jj]*(np.square(ii-xc)+np.square(jj-yc))
    r,c = np.where(img!=0)
    v = np.sort(img[r,c].flatten())[::-1]
    ftot = np.sum(v)
    fsum = 0
    ii = 0
    m20 = 0
    while fsum < 0.2*ftot:
        r,c = np.where(img==v[ii])
        for jj in range(len(r)):
            m20 += img[r[jj],c[jj]]*(np.square(r[jj]-xc)+np.square(c[jj]-yc))
        fsum += img[r[jj],c[jj]]
        ii += 1
    return np.log10(m20/mtot)
