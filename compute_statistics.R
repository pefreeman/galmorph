library(FITSio)
source("./prinaxes.R")
source("./segmap.R")
source("./compute_CA.R")
source("./compute_GM20.R")
source("./compute_MID.R")

compute_statistics = function(filename.input=c(),load=FALSE,increm.save=200,
                              scale.psf=0,scale.pix=0.06,scale.smooth=1,
                              eta=0.2,thrlev=10,
                              filename.output="compute_statistics.Rdata",
                              filename.increm="compute_statistics_increm.Rdata",
                              verbose=FALSE,delta.file=100,id.extract=FALSE,
                              id.extract.function=function(d){return(d)},
                              file.start=-9,image.hdu=1,seed.user=101,
                              noise.add=FALSE,noise.boost=NULL,subset=NULL)
{
  nfile = length(filename.input)

  if ( length(filename.input) == 0 && load == FALSE ) {
    return(list(finish=FALSE,file.number=0))
  }

  if ( load == FALSE ) {
    M       = M_o       = M_p       = rep(-9,nfile)
    M_level = M_level_o = M_level_p = rep(-9,nfile)
    I     = D     = rep(-9,nfile)
    Gini  = M20   = rep(-9,nfile)
    C     = A     = rep(-9,nfile)
    axmax = axmin = angle = rep(-9,nfile)
    sn    = rep(-9,nfile)
    size  = rep(-9,nfile)
    noise = sigma.hat.old = sigma.hat.new = rep(-9,nfile)
    id    = 1:nfile
    iilo  = 1
    if ( file.start > 1 ) iilo = file.start  # for debugging
  } else {
    load(filename.increm)
    nfile = length(M)
    iilo  = ii+1
  }

  set.seed(seed.user)

  for ( ii in iilo:nfile ) {

    if ( iilo > nfile ) break;

    if ( is.null(subset) == FALSE ) {
      if ( subset[ii] == 0 ) next;
    }

    result = tryCatch(
      {
        # Input galaxy image from postage stamp
        if ( id.extract == TRUE ) {
          id[ii] = id.extract.function(filename.input[ii])
        }

        x = readFITS(filename.input[ii],hdu=image.hdu)
        img = x$imDat
        rm(x)

        # Check for zeroes in the postage stamp image, indicating
        # possible edge-of-field issues.
        w = which(img==0,arr.ind=TRUE)
        if ( length(w[,1])/length(img) > 0.05 ) next;  # ad hoc

        # Add noise to the image
        if ( noise.add == TRUE ) {
          out = segmap(img,eta=eta,thrlev=thrlev)
          smap = out$c
          smap.all = out$allc
          if ( is.null(smap.all) == T ) next;
          w = which(smap.all==0,arr.ind=T)
          sigma.hat.old[ii] = sd(img[w])
          sigma.boost = sqrt(noise.boost^2-1)*sigma.hat.old[ii]
          img = img + rnorm(length(img),mean=0,sd=sigma.boost)
          sigma.hat.new[ii] = sd(img[w])
        }

        # Smooth the image if necessary to match PSFs
        if ( scale.psf > 0 ) img = conv_image(img,scale.psf)

        # Create segmentation map (mask) based on algorithm of Freeman et al.
        out = segmap(img,eta=eta,thrlev=thrlev)
        smap = out$c
        smap.all = out$allc  # for A computation (and adding noise)

        w = which(smap>0,arr.ind=T)
        r_smap = sqrt(length(w)/2/pi)
        if ( is.null(smap) == T ) next;

        # Estimator of galaxy "size" (effective "radius")
        w = which(smap>0,arr.ind=T)
        size[ii] = sqrt(length(w)/2/pi)

        # My estimator of signal-to-noise (S/N)
        w = which(smap.all==0,arr.ind=T)
        muhat = mean(img[w],trim=0.05)
        sdhat = sd(img[w])
        z = c()
        for ( jj in 1:dim(img)[1] ) {
          for ( kk in 1:dim(img)[2] ) {
            if ( smap[jj,kk] > 0 ) {
              if ( sqrt((jj-(dim(img)[1])/2)^2+(kk-(dim(img)[2])/2)^2) < 
                   0.5*size[ii] ) {
                z = append(z,(img[jj,kk]-muhat)/sdhat)
              }
            }
          }
        }
        if ( length(z) == 0 ) next;
        sn[ii] = median(z)

        # Convert size from pixels to arc-seconds
        size[ii] = scale.pix*size[ii]

        # Zero out all pixels outside the mask
        oimg   = img
        img    = img*smap
        w      = which(img<0,arr.ind=T)
        img[w] = 0

        # Normalize the image
        nimg = img/sum(img)

        # Compute principal axes (semi-major/semi-minor axes, rotation angle)
        pa = prinaxes(nimg)
        axmax[ii] = pa$axmax
        axmin[ii] = pa$axmin
        angle[ii] = pa$angle

        # Compute M statistic of Freeman et al. (M_o) plus two variants
        out = M_statistic(img,levels=seq(0.00,0.99,by=0.025))
        M[ii] = out$M
        M_o[ii] = out$M_o
        M_p[ii] = out$M_p
        M_level[ii] = out$level
        M_level_o[ii] = out$level_o
        M_level_p[ii] = out$level_p

        # Compute I and D statistics of Freeman et al., 
        # with pre-smoothing of data to mitigate noise (kernel size ~ 1 pixel)
        out = I_statistic(img,scale=scale.smooth)
        I[ii] = out$intensity_ratio
        D[ii] = D_statistic(img,out$x_peak,out$y_peak)

        # Compute Gini and M_20 based on prescription of Lotz et al. (2004)
        Gini[ii] = Gini_statistic(img)
        M20[ii]  = M20_statistic(img)

        # Compute CA (Conselice 2003)
        xc = dim(oimg)[1]/2
        yc = dim(oimg)[2]/2
        if ( xc%%2 == 0 ) xc=xc+0.5
        if ( yc%%2 == 0 ) yc=yc+0.5
        cir   = pet_rad_cir2(oimg,xc,yc)
        cimg  = cir$img_cir
        r_pet = cir$r_cir
        if ( r_pet > 2*r_smap ) r_pet = r_smap
        C[ii] = C_statistic(oimg,cimg,r_pet)
        A[ii] = A_statistic(oimg,cimg,smap.all,r_pet)
      },
      error = function(c) { return(FALSE) }
    )

    if ( result == FALSE ) {
      return(list(finish=FALSE,file.number=ii))
    }

    if ( ii%%increm.save == 0 ) {
      save(ii,scale.smooth,scale.pix,scale.psf,eta,thrlev,id,filename.input,
           M_level,M,M_level_o,M_o,M_level_p,M_p,I,D,axmax,axmin,angle,sn,size,
           Gini,M20,C,A,noise,sigma.hat.old,sigma.hat.new,subset,
           file=filename.increm)
    }
    if ( verbose == TRUE ) {
      if ( ii%%delta.file == 0 ) cat("Completed ",ii," of ",nfile,"...\n")
    }
  }

  save(ii,scale.smooth,scale.pix,scale.psf,eta,thrlev,id,filename.input,
       M_level,M,M_level_o,M_o,M_level_p,M_p,I,D,axmax,axmin,angle,sn,size,
       Gini,M20,C,A,noise,sigma.hat.old,sigma.hat.new,subset,
       file=filename.output)

  return(list(finish=TRUE,file.number=ii))
}

