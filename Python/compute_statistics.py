#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 16:49:58 2018

@author: peterfreeman
"""

import numpy as np
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
import segmap as sm
import prinaxes as pa
import compute_MID as mid
import compute_GM20 as gm20
import compute_CA as ca
import compute_genCA as gca

try:
    import progressbar as pb
except:
    pass

def compute_statistics(filename,trim=False,trimAmount=0,scalePSF=0,scalePix=0.06,
                       scaleSmooth=1,numSigma=5,incCounter=100):
    
    M      = np.full(len(filename),-9.0)
    I      = np.full(len(filename),-9.0)
    D      = np.full(len(filename),-9.0)
    G      = np.full(len(filename),-9.0)
    M20    = np.full(len(filename),-9.0)
    C      = np.full(len(filename),-9.0)
    A      = np.full(len(filename),-9.0)
    GC     = np.full(len(filename),-9.0)
    angle  = np.full(len(filename),-9.0)
    axmax  = np.full(len(filename),-9.0)
    axmin  = np.full(len(filename),-9.0)
    R      = np.full(len(filename),-9.0)
    SN     = np.full(len(filename),-9.0)
    sigma  = np.full(len(filename),-9.0)
    segpix = np.full(len(filename),-9.0)
    
    errFlag = np.zeros(len(filename))
        
    try:
        widgets = [pb.Percentage(),pb.Bar()]
        bar = pb.ProgressBar(widgets=widgets,max_value=len(filename))
    except:
        pass
    
    for ii in range(len(filename)):
        
        try:
                    
            hdu = fits.open(filename[ii])
            img = hdu[0].data
            hdu.close()
                
            r,c = np.shape(img)
            if trim == True:
                img = img[trimAmount:(r+2-trimAmount),
                          trimAmount:(c+2-trimAmount)]
                    
            if scalePSF > 0:
                kernel = Gaussian2DKernel(3*scalePSF)
                img = convolve(img,kernel)
          
            smap,smapAll,sigma[ii] = sm.segmap(img,numSigma=numSigma)
            if np.max(smap) == 0:
                errFlag[ii] = 1
                continue
                    
            oimg = np.copy(img)
            img = img*smap
            
            r,c = np.where(smap>0)
            if len(np.unique(r)) == 1 or len(np.unique(c)) == 1:
                errFlag[ii] = 2
                continue
            segpix[ii] = len(r)
   
            # Simple estimator of galaxy "size" (effective "radius")
            r,_ = np.where(smap>0)
            R[ii] = np.sqrt(len(r)/np.pi)
        
            # Simple estimator of signal-to-noise (S/N)
            sigmaImg = img/sigma[ii]
            r,c = np.where(sigmaImg>numSigma)
            SN[ii] = np.median(sigmaImg[r,c])
        
            # Convert size from pixels to arc-seconds
            pixelSize = R[ii]
            R[ii] = scalePix*R[ii]
            
            r,c = np.where(img<0)
            img[r,c] = 0.0
            nimg = img/np.sum(img)
            axmax[ii],axmin[ii],angle[ii] = pa.prinaxes(nimg)
    
            # Compute corrected M statistic of Freeman et al.
            M[ii] = mid.MStatistic(img,levels=np.arange(0,0.99,0.025))
        
            # Compute I and D statistics of Freeman et al.,
            # with pre-smoothing of data to mitigate noise (kernel size ~ 1 pixel)
            I[ii],XPeak,YPeak = mid.IStatistic(img,scale=scaleSmooth)
            D[ii],DXCen,DYCen = mid.DStatistic(img,XPeak,YPeak)

            # Compute Gini and M_20 based on prescription of Lotz et al. (2004)
            G[ii]   = gm20.GStatistic(img)
            M20[ii] = gm20.M20Statistic(img)
        
            # Comppute generalized C statistic (new)
            GC[ii] = gca.genCStatistic(img)
          
            # Image "shift" algorithm to deal with mis-centroided galaxies
            # Primarily intended to reduce A statistic value in such galaxies
            xc,yc = np.shape(oimg)/np.array([2,2])
            rX = np.round(DXCen-xc).astype(int)
            if rX > 0:
                rmX = np.arange(2*rX)
                DXCen -= 2*rX
            elif rX < 0:
                rmX = np.arange(np.shape(oimg)[0]+2*rX,np.shape(oimg)[0])    
            if rX != 0:
                oimg = np.delete(oimg,rmX,0)
                smapAll = np.delete(smapAll,rmX,0)
            rY = np.round(DYCen-yc).astype(int)
            if rY > 0:
                rmY = np.arange(2*rY)
                DYCen -= 2*rY
            elif rY < 0:
                rmY = np.arange(np.shape(oimg)[1]+2*rY,np.shape(oimg)[1])       
            if rY != 0:
                oimg = np.delete(oimg,rmY,1)
                smapAll = np.delete(smapAll,rmY,1)
            # end new algorithm

            # Compute C and A (Conselice 2003)
            cimg,rPet = ca.petRadCir2(oimg,DXCen,DYCen)
            if rPet > 2*pixelSize:
                rPet = pixelSize
            C[ii] = ca.CStatistic(oimg,cimg,rPet)
            A[ii] = ca.AStatistic(oimg,cimg,smapAll,rPet)
            
            try:
                bar.update(ii+1)
            except:
                if ii > 0 and ii%incCounter == 0:
                    print('Finished file index ',ii)
                
        except:
            errFlag[ii] = 3
            continue
    
    try:
        bar.finish()
    except:
        pass
        
    return M,I,D,G,M20,C,A,GC,R,SN,sigma,segpix,axmax,axmin,angle,errFlag


