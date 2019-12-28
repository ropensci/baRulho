#' Measure blur ratio in the frequency domain 
#' 
#' \code{blur_ratio_spectrum} measures blur ratio of frequency spectra from signals referenced in an extended selection table.
#' @usage blur_ratio_spectrum(X, parallel = 1, pb = TRUE, method = 1, ssmooth = 50, 
#' output = "est", img = FALSE, res = 150,  wl = 512, 
#' ovlp = 70, pal = reverse.gray.colors.2, collevels = seq(-60, 0, 5), dest.path = NULL)
#' @param X object of class 'selection_table', 'extended_selection_table' created by the function \code{\link[warbleR]{selection_table}} from the warbleR package.
#' @param parallel Numeric vector of length 1. Controls whether parallel computing is applied by specifying the number of cores to be used. Default is 1 (i.e. no parallel computing).
#' @param pb Logical argument to control if progress bar is shown. Default is \code{TRUE}.
#' @param method Numeric vector of length 1 to indicate the 'experimental design' for measuring spectrum correlation. Two methods are available:
#' \itemize{
#' \item \code{1}: compare all signals with their counterpart that was recorded at the closest distance to source (e.g. compare a signal recorded at 5m, 10m and 15m with its counterpart recorded at 1m). This is the default method. 
#' \item \code{2}: compare all signals with their counterpart recorded at the distance immediately before (e.g. a signal recorded at 10m compared with the same signal recorded at 5m, then signal recorded at 15m compared with same signal recorded at 10m and so on).
#' }
#' @param ssmooth Numeric vector of length 1 determining the length of the sliding window used for a sum smooth for power spectrum calculation (in kHz). Default is 100.
#' @param output Character vector of length 1 to determine if an extended selection table ('est') or a data frame ('data.frame') is returned.
#' @param img Logical argument to control if image files in 'jpeg' format containing the images being compared and the corresponding spectra are produced. Default is no images ( \code{FALSE}).
#' @param res Numeric argument of length 1. Controls image resolution. Default is 150 (faster) although 300 - 400 is recommended for publication/presentation quality.
#' @param wl A numeric vector of length 1 specifying the window length of the spectrogram, default 
#'   is 512.
#' @param ovlp Numeric vector of length 1 specifying the percent overlap between two 
#'   consecutive windows, as in \code{\link[seewave]{spectro}}. Default is 70.
#' @param pal A color palette function to be used to assign colors in the 
#'   plot, as in \code{\link[seewave]{spectro}}. Default is reverse.gray.colors.2. 
#' @param collevels	Numeric vector indicating a set of levels which are used to partition the amplitude range of the spectrogram (in dB) as in \code{\link[seewave]{spectro}}. Default is \code{seq(-60, 0, 5)}. 
#' @param dest.path Character string containing the directory path where the image files will be saved. If NULL (default) then the folder containing the sound files will be used instead.
#' @return Data frame similar to input data, but also includes a new column (blur.ratio.spectrum)
#' with the blur ratio values. If \code{img = TRUE} it also returns 1 image file (in 'jpeg' format) for each comparison showing spectrograms of both signals and the overlaid power spectrum (as probability mass functions (PMF)).
#' @export
#' @name blur_ratio_spectrum
#' @details Blur ratio measures the degradation of sound as a function of the change in signal energy in the frequency domain, analogous to the blur ratio proposed by Dabelsteen et al (1993) for the time domain (and implemented in \code{\link{blur_ratio}}). Low values indicate low degradation of signals. The function measures the blur ratio of spectra from signals in which a reference playback has been re-recorded at different distances. Blur ratio is measured as the mismatch between power spectra (expressed as probability density functions) of the reference signal and the re-recorded signal. The function compare each signal type to the corresponding reference signal. The 'signal.id' column must be used to tell the function to only compare signals belonging to the same category (e.g. song-types). Two methods for settting the experimental design are provided. All wave objects in the extended selection table must have the same sampling rate so the length of spectra is comparable.   
#' @examples
#' {
#' # load example data
#' data("playback_est")
#' 
#' # remove noise selections
#' playback_est <- playback_est[playback_est$signal.id != "noise", ]
#' 
#' # using method 1
#' blur_ratio_spectrum(X = playback_est)
#' 
#' # using method 2
#' blur_ratio_spectrum(X = playback_est, method = 2)
#' }
#' 
#' @author Marcelo Araya-Salas (\email{marceloa27@@gmail.com}) #' @references {
#' Dabelsteen, T., Larsen, O. N., & Pedersen, S. B. (1993). Habitat-induced degradation of sound signals: Quantifying the effects of communication sounds and bird location on blur ratio, excess attenuation, and signal-to-noise ratio in blackbird song. The Journal of the Acoustical Society of America, 93(4), 2206.
#' 
#' Araya-Salas, M. (2019), baRulho: a R package to quantify habitat-induced degradation of (animal) acoustic signals. R package version 1.0.0
#' }
#last modification on dec-27-2019 (MAS)

blur_ratio_spectrum <- function(X, parallel = 1, pb = TRUE, method = 1,
                       ssmooth = 50, output = "est", 
                       img = FALSE, res = 150,  wl = 512, ovlp = 70, 
                       pal = reverse.gray.colors.2, collevels = seq(-60, 0, 5), 
                       dest.path = NULL){
  
  # set pb options 
  on.exit(pbapply::pboptions(type = .Options$pboptions$type), add = TRUE)
  
  # set dest.path if not provided
  if (is.null(dest.path)) 
    dest.path <- getwd() else 
      if (!dir.exists(dest.path)) 
    stop("'dest.path' provided does not exist")
  
  # If parallel is not numeric
  if (!is.numeric(parallel)) stop("'parallel' must be a numeric vector of length 1") 
  if (any(!(parallel %% 1 == 0),parallel < 1)) stop("'parallel' should be a positive integer")
  
  # If method is not numeric
  if (!is.numeric(method)) stop("'method' must be a numeric vector of length 1") 
  if (!any(method %in% 1:2)) stop("'method' must be either 1 or 2")
  
  # check signal.id column 
  if (is.null(X$signal.id)) stop("'X' must containe a 'signal.id' column")
  #check output
  if (!any(output %in% c("est", "data.frame"))) stop("'output' must be either 'est' or 'data.frame'")  
  
  # must have the same sampling rate
  if (length(unique(attr(X, "check.results")$sample.rate)) > 1) 
    stop("all wave objects in the extended selection table must have the same sampling rate (they can be homogenized using warbleR::resample_est())")
  
  # set pb options 
  pbapply::pboptions(type = ifelse(as.logical(pb), "timer", "none"))
  
  # set clusters for windows OS
  if (Sys.info()[1] == "Windows" & parallel > 1)
    cl <- parallel::makePSOCKcluster(getOption("cl.cores", parallel)) else cl <- parallel
  
  # print message
  if (pb) write(file = "", x = "calculating power spectra (step 1 of 2):")
  
  # calculate all spectra apply function
  specs <- pbapply::pblapply(X = 1:nrow(X), cl = cl, FUN = function(y, wl)   {
    clp <- warbleR::read_wave(X = X, index = y)
    seewave::spec(wave = clp, f = clp@samp.rate, plot = FALSE, wl = wl)
  }) 
  
  # add sound file selec column and names to spectra (weird column name so it does not overwrite user columns)
  X$TEMP....sgnl <- names(specs) <- paste(X$sound.files, X$selec, sep = "-")
  
  # add second column with names of the reference signals to be compare against
  X$TEMP....rfrnc <- sapply(1:nrow(X), function(x, meth = method){
    
    # extract for single signal and order by distance
    Y <- as.data.frame(X[X$signal.id == X$signal.id[X$TEMP....sgnl == X$TEMP....sgnl[x]], , drop = FALSE])
    Y <- Y[order(Y$distance), ]
    
    # method 1 compare to closest distance to source
    if (meth == 1) z <- Y$TEMP....sgnl[which.min(Y$distance)] else # if method 2
      # if not the first row then the previous row
      if (Y$TEMP....sgnl[1] != X$TEMP....sgnl[x]) z <- X$TEMP....sgnl[x - 1] else # else the first row
        z <- Y$TEMP....sgnl[1] 
    
    return(z)
  })
  
  ## function to measure blur ratio
  # y and z are the sound.files+selec names of the signals and reference signal (model)
  # spectrum mismatch ratio 
  blur_sp_FUN <- function(x, res, ovlp, wl, collevels, pal, ...){
    
    # get names of signal and reference
    sgnl <-  X$TEMP....sgnl[x]
    rfrnc <- X$TEMP....rfrnc[x]
    
    # if signals are the same or the selection is noise return NA
    if (sgnl == rfrnc | any(c(X$signal.id[X$TEMP....sgnl == sgnl], X$signal.id[X$TEMP....rfrnc == rfrnc]) == "noise")) out <- NA else {
      
      # extract spectrum for signal and model 
      sgnl.spc <- specs[[which(names(specs) == sgnl)]]
      rfrnc.spc <- specs[[which(names(specs) == rfrnc)]]
      
      # make them the same frequency range as reference
      bp <- c(X$bottom.freq[X$TEMP....sgnl == rfrnc], X$top.freq[X$TEMP....sgnl == rfrnc])
      bp <- bp + c(-0.2, 0.2)       # add 0.2 kHz buffer
      if (bp[1] < 0)  # force 0 if negative
        bp[1] <- 0
      if (bp[2] > ceiling(attr(X, "check.results")$sample.rate[1] / 2) - 1) 
        bp[2] <- ceiling(attr(X, "check.results")$sample.rate[1] / 2) - 1 # force lower than nyquist freq if higher
      
      # homogenize freq range and remove freq column
      sgnl.spc <- sgnl.spc[sgnl.spc[, 1] > bp[1] & sgnl.spc[, 2] < bp[2], 2]
      rfrnc.spc <- rfrnc.spc[rfrnc.spc[, 1] > bp[1] & rfrnc.spc[, 2] < bp[2], 2]
        
      # applied ssmooth
      if (!is.null(ssmooth)) {
        sgnl.spc <- as.matrix(seewave::sumsmooth(sgnl.spc, wl = ssmooth))
        rfrnc.spc <- as.matrix(seewave::sumsmooth(rfrnc.spc, wl = ssmooth))
        }

      # convert envelopes to PMF (probability mass function)
      rfrnc.spc <- rfrnc.spc / sum(rfrnc.spc)
      sgnl.spc <- sgnl.spc / sum(sgnl.spc)
      
      # get blur ratio as half the sum of absolute differences between spectra PMFs
      bl.rt <- sum(abs(rfrnc.spc - sgnl.spc)) / 2
      
      # plot
      if (img)
      {
        warbleR:::img_wrlbr_int(filename = paste0("blur_ratio_", X$signal.id[x], "-", rfrnc, "-", sgnl, ".jpeg"), path = dest.path, width = 10.16 * 1.5, 
                                height = 10.16 , units = "cm", res = res)
        
        # create time values for area calculation
        f.vals <- seq(bp[1], bp[2], length.out = length(rfrnc.spc))
        
        # difference between spectra
        spc.diff <- rfrnc.spc - sgnl.spc
        
        # matrix for layout
        ly.mat <- matrix(
          c(0, 0.3, 0, 0.5, # bottom left spectrogram
            0, 0.3, 0.5, 1, # top left spectrogram
            0.2, 1, 0, 1),  # right pannel spectra
          nrow = 3, byrow = TRUE)
        
        # close if open any screen
        invisible(close.screen(all.screens = TRUE))
        
        split.screen(ly.mat)
        
        # plot spectra
        screen(3)
        
        # set image margins
        par(mar = rep(4, 0, 4, 4))
        
        # reference spectrum first        
        plot(x = rfrnc.spc, y = f.vals,  type = "l", xlab = "", ylab = "", col = "#E37222", xlim = c(min(rfrnc.spc, sgnl.spc), max(rfrnc.spc, sgnl.spc) * 1.1), cex.main = 0.8, lwd = 1.2, yaxt = "n")
        
        # add x axis label
        mtext(text = "Power spectrum (PMF)", side = 1, line = 2.5)
        
        # add title
        mtext(text = paste("Signal type:", X$signal.id[x]), side = 3, line = 3, cex = 0.7)
        mtext(text = paste("Reference:", rfrnc), side = 3, line = 2, col = "#E37222", cex = 0.7)
        mtext(text = paste("Signal:", sgnl), side = 3, line = 1, col = "#07889B", cex = 0.7)
        
        # add y axis
        axis(side = 4)
        mtext(text = "Frequency (kHz)", side = 4, line = 2.5)
        
        # add signal spectrum
        lines(sgnl.spc, f.vals, col= "#07889B", lwd = 1.2)
        
        # signal spectrum on top
        polygon(y = c(f.vals, rev(f.vals)), x = c(sgnl.spc, rev(rfrnc.spc)), col =  "#07889B33", border = NA)
        
        # get plotting area limits
        usr <- par("usr")
        
        # and blu ratio value
        text(x = ((usr[1] + usr[2]) / 2) + usr[1], y = usr[4] * 0.9, paste("Blur ratio of spectrum:", round(bl.rt, 2)), cex = 0.8)
        
        # spectrogram of reference
        screen(1)
        
        # calculate margin for spectrogram
        mar.rf <- attr(X, "check.results")$duration[which(paste(X$sound.files, X$selec, sep = "-") == rfrnc)] / 4
        
        # index of reference
        rf.indx <- which(paste(X$sound.files, X$selec, sep = "-") == rfrnc)
        
        # freq limit of reference
        flim <- c(X$bottom.freq[rf.indx], X$top.freq[rf.indx])
        
        # extract clip reference and signal
        clp.sgnl <- warbleR::read_wave(X = X, index = x, from = X$start[x] - mar.rf, X$end[x] + mar.rf)
        clp.rfnc <- warbleR::read_wave(X = X, index = rf.indx, from = X$start[rf.indx] - mar.rf, X$end[rf.indx] + mar.rf)
        
        
        ## plot spectros
        # signal at bottom left
        screen(1)
        par(mar = c(0.3, 0.3, 0.15, 0.3))
        
        warbleR:::spectro_wrblr_int2(wave = clp.sgnl, f = clp.sgnl@samp.rate, 
                                     flim = flim,  axisX = FALSE, axisY = FALSE, 
                                     tlab = NULL, flab = NULL, main = NULL, grid = FALSE, rm.zero = TRUE, cexaxis = 1.2, add = TRUE, ovlp = ovlp, wl = wl, collevels = collevels, palette = pal)
        
        # add box with signal color
        box(col = "#07889B", lwd = 3)
        
        # reference at top left
        screen(2)
        par(mar = c(0.15, 0.3, 0.3, 0.3))
        
        warbleR:::spectro_wrblr_int2(wave = clp.sgnl, f = clp.sgnl@samp.rate, 
                                     flim = flim,  axisX = FALSE, axisY = FALSE, 
                                     tlab = NULL, flab = NULL, main = NULL, grid = FALSE, rm.zero = TRUE, cexaxis = 1.2, add = TRUE, ovlp = ovlp, wl = wl, collevels = collevels, palette = pal)
        
        # add box with reference color
        box(col = "#E37222", lwd = 3)
        
        # close graph    
        dev.off()        
      }
      
      # return maximum correlation
      return(bl.rt)
    }
    return(out)
  } 
  
  if (pb & !img) write(file = "", x = "calculating blur ratio (step 2 of 2):")
  if (pb & img) write(file = "", x = "calculating blur ratio and producing images (step 2 of 2):")
    
  # get blur ratio
  # calculate all spectra apply function
  X$blur.ratio.spectrum <- pbapply::pbsapply(X = 1:nrow(X), cl = cl, FUN = function(x, rs = res, wle = wl, colvs = collevels, pl = pal, ovp = ovlp, ...)   {
    blur_sp_FUN(x, res = rs, ovlp = ovp, wl = wle, collevels = colvs, pal = pl, ...)
  }) 
  
  
  # remove temporal columns
  X$TEMP....sgnl <- X$TEMP....rfrnc <- NULL
  
  # convert to data frame instead of extended selection table
  if (output == "data.frame") 
    X <- as.data.frame(X)
  
  return(X)
}
