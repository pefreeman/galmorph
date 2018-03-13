########################################################
# R Functions for Computing Galaxy Morphology Statistics
########################################################
#

Use:

  - create a file that calls the function compute_statistics
    (see compute_driver for an example)
  - run the driver file (e.g., "R CMD BATCH compute_driver.R")

Input to compute_statistics:

  - filename.input:  character array of names of files containing 
                     individual galaxy images (in FITS format)
  - id.extract:      if TRUE, use the function specified by id.extract.function
                     to extract the galaxy ID number from its filename;
                     if FALSE, assume the ID number for the first analyzed
                     galaxy is 1, etc.
  - id.extract.function: function for extracting the galaxy ID number from
                         its filename
  - filename.output: name of Rdata file containing all galaxy statistics
                     (which are recovered via the load command)
  - filename.increm: name of Rdata file containing an incremental dump
                     of galaxy statistics (see increm.save)
  - increm.save:     output statistics to disk every increm.save galaxies
                     (see filename.increm)
  - load:            if TRUE, load results from the file specified by
                     the argument filename.increm and re-commence 
                     analysis where that file left off (e.g., if the
                     suite crashes at file number 307, and delta.file = 100,
                     the program will re-commence analysis at file 201
                     (assuming increm.save has its default value, 200),
                     as the output statistics for files 1 -> 200 would have
                     been saved); if FALSE, begin analysis with the first file
  - scale.psf:       convolve all images with Gaussian kernel with 
                     sigma = scale.psf (units: pixels)
  - scale.pix:       plate scale (units: arcsec/pixel)
  - scale.smooth:    additionally convolve all images with Gaussian kernel
                     with sigma = scale.smooth, for I and D statistic
                     computation (units: pixels)
  - eta,thrlev:      parameters for segmentation map construction;
                     for experts only!
  - verbose:         if TRUE, compute_statistics will inform you of progress
                     every delta.file files
  - delta.file:      how many files to procress between messages to STDOUT
                     indicating the codes progress (moot if verbose=FALSE)
  - file.start:      start analysis at the file.start'th file in the 
                     filename.input array; primarily for debugging
  - image.hdu:       the HDU of the image in each input FITS file 
  - seed.user:       random number seed; used when adding noise to images
  - noise.add:       if TRUE, add noise to data:
                       Y -> Y + epsilon
                     where epsilon ~ N(0,sigma.boost^2) and
                       sigma.boost = sqrt(noise.boost^2-1)*sigma.hat.old
  - noise.boost:     see noise.add

Output from compute_statistics (in Rdata file; see filename.output above):

  - scale.smooth,scale.pix,scale.psf,eta,thrlev,filename.input: see above
  - id:            array of galaxy image ID numbers
  - M,M_level:     array of M statistics and level-set thresholds, 
                   current algorithm
  - M_o,M_level_o: array of M statistics and level-set thresholds, 
                   old algorithm (Freeman et al. 2013)
  - M_p,M_level_p: array of M statistics and level-set thresholds, 
                   Peth et al. (2015) algorithm
  - I:             array of I statistics
  - D:             array of D statistics
  - Gini:          array of Gini statistics
  - M20:           array of M20 statistics
  - C:             array of C statistics
  - A:             array of A statistics
  - axmax,axmin:   length of principal axes (units: pixels)
  - sn:            S/N estimate
  - size:          galaxy size estimate, based on segmentation map 
                   (units: arcseconds)
  - sigma.hat.old, array of noise estimates (applicable if noise.add=TRUE)
    sigma.hat.new:

