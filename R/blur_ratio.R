#' Measure blur ratio in the time domain
#' 
#' \code{blur_ratio} measures blur ratio in sounds referenced in an extended selection table.
#' @usage blur_ratio(X, parallel = 1, cores = getOption("mc.cores", 1), 
#' pb = getOption("pb", TRUE), method = getOption("method", 1), ssmooth = 200, 
#' msmooth = NULL, output = "est", img = FALSE, res = 150, 
#' hop.size = getOption("hop.size", 11.6), wl = getOption("wl", NULL), 
#' ovlp = getOption("ovlp", 70), pal = viridis, collevels = seq(-120, 0, 5), 
#' dest.path = NULL, path = getOption("sound.files.path", "."))
#' @param X Object of class 'data.frame', 'selection_table' or 'extended_selection_table' (the last 2 classes are created by the function \code{\link[warbleR]{selection_table}} from the warbleR package) with the reference to the sounds in the master sound file. Must contain the following columns: 1) "sound.files": name of the .wav files, 2) "selec": unique selection identifier (within a sound file), 3) "start": start time and 4) "end": end time of selections, 5)  "bottom.freq": low frequency for bandpass, 6) "top.freq": high frequency for bandpass and 7) "sound.id": ID of sounds used to identify counterparts across distances. Each sound must have a unique ID within a distance.
#' @param parallel DEPRECATED. Use 'cores' instead.
#' @param cores Numeric vector of length 1. Controls whether parallel computing is applied by specifying the number of cores to be used. Default is 1 (i.e. no parallel computing).
#' @param pb Logical argument to control if progress bar is shown. Default is \code{TRUE}.
#' @param method Numeric vector of length 1 to indicate the 'experimental design' for measuring envelope correlation. Two methods are available:
#' \itemize{
#' \item \code{1}: compare all sounds with their counterpart that was recorded at the closest distance to source (e.g. compare a sound recorded at 5m, 10m and 15m with its counterpart recorded at 1m). This is the default method. 
#' \item \code{2}: compare all sounds with their counterpart recorded at the distance immediately before (e.g. a sound recorded at 10m compared with the same sound recorded at 5m, then sound recorded at 15m compared with same sound recorded at 10m and so on).
#' }
#' @param ssmooth Numeric vector of length 1 determining the length of the sliding window (in amplitude samples) used for a sum smooth for amplitude envelope calculation (used internally by \code{\link[seewave]{env}}). Default is 200.
#' @param msmooth Numeric vector of length 2 to smooth the amplitude envelope with a mean sliding window for amplitude envelope calculation. The first element is the window length (in number of amplitude values) and the second one the window overlap (used internally by \code{\link[seewave]{env}}). 
#' @param output Character vector of length 1 to determine if an extended selection table ('est', default), a data frame ('data.frame') or a list ("list") containing the extended selection table (first object in the list) and all (smoothed) wave envelopes (second object in the list) is returned. The envelope data can be used for plotting. 'est' format only available if 'X' is itself an extended selection table.
#' @param img Logical argument to control if image files in 'jpeg' format containing the images being compared and the corresponding envelopes are produced. Default is no images ( \code{FALSE}).
#' @param res Numeric argument of length 1. Controls image resolution. Default is 150 (faster) although 300 - 400 is recommended for publication/presentation quality.
#' @param hop.size A numeric vector of length 1 specifying the time window duration (in ms). Default is 11.6 ms, which is equivalent to 512 wl for a 44.1 kHz sampling rate. Ignored if 'wl' is supplied.
#' @param wl A numeric vector of length 1 specifying the window length of the spectrogram, default 
#' is NULL. If supplied, 'hop.size' is ignored.
#' @param ovlp Numeric vector of length 1 specifying the percent overlap between two 
#'   consecutive windows, as in \code{\link[seewave]{spectro}}. Only used when plotting. Default is 70. Applied to both spectra and spectrograms on image files.
#' @param pal A color palette function to be used to assign colors in the 
#'   plot, as in \code{\link[seewave]{spectro}}. Default is \code{\link[viridis]{viridis}}. 
#' @param collevels	Numeric vector indicating a set of levels which are used to partition the amplitude range of the spectrogram (in dB) as in \code{\link[seewave]{spectro}}. Default is \code{seq(-120, 0, 5)}. 
#' @param dest.path Character string containing the directory path where the image files will be saved. If NULL (default) then the folder containing the sound files will be used instead.
#' @return Data frame similar to input data, but also includes two new columns ('reference' and 'blur.ratio')
#' with the reference sound and blur ratio values. If \code{img = TRUE} it also returns 1 image file (in 'jpeg' format) for each comparison showing spectrograms of both sounds and the overlaid amplitude envelopes (as probability mass functions (PMF)). Spectrograms are shown within the frequency range of the reference sound and also show vertical lines with the start and end of sounds to allow users to visually check alignment. If \code{output = 'list'} the output would be a list including the data frame just described and a data frame with envelopes (amplitude values) for all sounds.
#' @param path Character string containing the directory path where the sound files are found. Only needed when 'X' is not an extended selection table.
#' @export
#' @name blur_ratio
#' @details Blur ratio measures the degradation of sound as a change in sound energy in the time domain as described by Dabelsteen et al (1993). Low values indicate low degradation of sounds. The function measures the blur ratio on sounds in which a reference playback has been re-recorded at different distances. Blur ratio is measured as the mismatch between amplitude envelopes (expressed as probability mass functions) of the reference sound and the re-recorded sound. By converting envelopes to probability mass functions the effect of energy attenuation is removed, focusing the analysis on the modification of the envelope shape. The function compares each sound to the corresponding reference sound within the supplied frequency range (e.g. bandpass) of the reference sound ('bottom.freq' and 'top.freq' columns in 'X'). The 'sound.id' column must be used to tell the function to only compare sounds belonging to the same category (e.g. song-types). Two methods for setting the experimental design are provided. All wave objects in the extended selection table must have the same sampling rate so the length of envelopes is comparable.
#' @seealso \code{\link{envelope_correlation}}, \code{\link{spectrum_blur_ratio}}
#' @examples
#' {
#' # load example data
#' data("playback_est")
#' 
#' # remove ambient selections
#' playback_est <- playback_est[playback_est$sound.id!= "ambient", ]
#' 
#' # using method 1
#'blur_ratio(X = playback_est)
#' 
#' # using method 2
#' blur_ratio(X = playback_est, method = 2)
#' }
#' 
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr}) 
#' @references {
#' Dabelsteen, T., Larsen, O. N., & Pedersen, S. B. (1993). Habitat-induced degradation of sound signals: Quantifying the effects of communication sounds and bird location on blur ratio, excess attenuation, and signal-to-noise ratio in blackbird song. The Journal of the Acoustical Society of America, 93(4), 2206.
#' 
#' Araya-Salas, M. (2020). baRulho: baRulho: quantifying habitat-induced degradation of (animal) acoustic signals in R. R package version 1.0.2
#' }
#last modification on dec-26-2019 (MAS)

blur_ratio <- function(X, parallel = 1, cores = getOption("mc.cores", 1), pb = getOption("pb", TRUE), method = getOption("method", 1),
                       ssmooth = 200, msmooth = NULL, output = "est", 
                       img = FALSE, res = 150, hop.size = getOption("hop.size", 11.6), wl = getOption("wl", NULL), ovlp = getOption("ovlp", 70), pal = viridis, collevels = seq(-120, 0, 5), dest.path = NULL, path = getOption("sound.files.path", ".")){
  
  # deprecated message
  if (parallel > 1) 
    stop2("'parallel' has been deprecated, Use 'cores' instead")
  
  # set path if not provided
  if (is.null(path)) 
    path <- getwd() else 
      if (!dir.exists(path)) 
        stop2("'path' provided does not exist")
  
  # set dest.path if not provided
  if (is.null(dest.path)) 
    dest.path <- getwd() else 
      if (!dir.exists(dest.path)) 
    stop2("'dest.path' provided does not exist")
  
  # make path null if extendeed selection table
  if (is_extended_selection_table(X))
    path <- NULL
  
  # hopsize  
  if (!is.numeric(hop.size) | hop.size < 0) stop2("'hop.size' must be a positive number") 
  
  # adjust wl based on hope.size
  if (is.null(wl))
    wl <- round(read_sound_file(X, index = 1, header = TRUE, path = path)$sample.rate * hop.size  / 1000, 0)
  
    # make wl even if odd
  if (!(wl %% 2) == 0) wl <- wl + 1
  
  # If cores is not numeric
  if (!is.numeric(cores)) stop2("'cores' must be a numeric vector of length 1") 
  if (any(!(cores %% 1 == 0), cores < 1)) stop2("'cores' should be a positive integer")
  
  # If method is not numeric
  if (!is.numeric(method)) stop2("'method' must be a numeric vector of length 1") 
  if (!any(method %in% 1:2)) stop2("'method' must be either 1 or 2")
  
  # check sound.idcolumn 
  if (is.null(X$sound.id)) stop2("'X' must contain a 'sound.id' column")
  
  #check output
  if (!any(output %in% c("est", "data.frame", "list"))) stop2("'output' must be 'est', 'data.frame' or 'list'")  
  
  # must have the same sampling rate
  if (length(unique(attr(X, "check.results")$sample.rate)) > 1) 
    stop2("all wave objects in the extended selection table must have the same sampling rate (they can be homogenized using warbleR::resample_est())")
  
  # set clusters for windows OS
  if (Sys.info()[1] == "Windows" & cores > 1)
    cl <- parallel::makePSOCKcluster(getOption("cl.cores", cores)) else cl <- cores
  
  # ingnore ssmooth if msmooth is supplied
  if (!is.null(msmooth)) 
    ssmooth <- NULL
  
  # add sound file selec column and names to X (weird column name so it does not overwrite user columns)
  if (pb) 
    write(file = "", x = paste0("Preparing data for analysis (step 1 out of 3):"))
  
  X <- prep_X_bRlo_int(X, method = method, cores = cores, pb = pb)
    
  # print message
  if (pb) write(file = "", x = "Calculating amplitude envelopes (step 2 out of 3):")
  
  # calculate all envelops apply function
  envs <- warbleR:::pblapply_wrblr_int(pbar = pb, X = 1:nrow(X), cl = cl, FUN = function(y, ssmth = ssmooth, msmth = msmooth, ov = ovlp)   {
    
    # load clip
    clp <- warbleR::read_sound_file(X = X, index = y, path = path)
    
    # define bandpass based on reference
    bp <- c(X$bottom.freq[X$TEMP....sgnl == X$reference[y]], X$top.freq[X$TEMP....sgnl == X$reference[y]])
    
    # bandpass filter
    clp <- seewave::ffilter(clp, from = bp[1] * 1000, 
                            ovlp = ov, to = bp[2] * 1000, bandpass = TRUE, 
                            wl = wl, output = "Wave")
    
    # calculate envelope
    nv <- seewave::env(wave = clp, f = clp@samp.rate, ssmooth = ssmth, msmooth = msmth, plot = FALSE, envt = "hil")[, 1]
    
    return(nv)
  }) 
  
  # add sound file selec column and names to envelopes
  names(envs) <- X$TEMP....sgnl
  
  ## function to measure blur ratio
  # y and z are the sound.files+selec names of the sounds and reference sound (model)
  # envelope mismatch ratio 
  blur_FUN <- function(x, res, ovlp, wl, collevels, pal, ...){
    
      # get names of sound and reference
      sgnl <-  X$TEMP....sgnl[x]
      rfrnc <- X$reference[x]
      
    # if sounds are the same or the selection is noise return NA
    if (sgnl == rfrnc | any(c(X$sound.id[X$TEMP....sgnl == sgnl], X$sound.id[X$reference == rfrnc]) == "ambient")) out <- NA else {
      
      # extract envelope for sound and model 
      sgnl.env <- envs[[which(names(envs) == sgnl)]]
      rfrnc.env <- envs[[which(names(envs) == rfrnc)]]
      
      # make them the same length as the shortest one
      if(length(sgnl.env) > length(rfrnc.env)) sgnl.env <- sgnl.env[1:length(rfrnc.env)]
      if(length(rfrnc.env) > length(sgnl.env)) rfrnc.env <- rfrnc.env[1:length(sgnl.env)]
      
      # duration (any works as they all must have the same sampling rate)
      # dur <- length(sgnl.env) / (attr(X, "check.results")$sample.rate[1] * 1000)
      dur <- length(sgnl.env) / warbleR::read_sound_file(X = X, index = x, path = path, header = TRUE)$sample.rate
      
      # convert envelopes to PMF (probability mass function)
      rfrnc.pmf <- rfrnc.env / sum(rfrnc.env)
      sgn.pmf <- sgnl.env / sum(sgnl.env)
      
      # get blur ratio as half the sum of absolute differences between envelope PMFs
      bl.rt <- sum(abs(rfrnc.pmf - sgn.pmf)) / 2
      
      # plot
      if (img)
      {
        img_bRlo_int(filename = paste0("blur_ratio_", X$sound.id[x], "-", rfrnc, "-", sgnl, ".jpeg"), path = dest.path, width = 10.16 * 1.5, 
                      height = 10.16 , units = "cm", res = res)
        
        # time values for plots
        time.vals <- seq(0, dur, length.out = length(sgnl.env))

        # difference between envelopes
        env.diff <- rfrnc.pmf - sgn.pmf
        
        # matrix for layout
        ly.mat <- matrix(
          c(0, 0.3, 0, 0.5, # bottom left spectrogram
            0, 0.3, 0.5, 1, # top left spectrogram
            0.2, 1, 0, 1),  # right pannel envelopes
          nrow = 3, byrow = TRUE)
        
        # save par settings
        oldpar <- par(no.readonly = TRUE)   
        on.exit(par(oldpar)) 
        
        # close if open any screen
        invisible(close.screen(all.screens = TRUE))
        
        # split screen
        split.screen(ly.mat)
        
        # plot envelopes
        screen(3)
        
        # set image margins
        par(mar = rep(4, 0, 4, 4))
        
        # reference envelope first        
        plot(time.vals, rfrnc.pmf, type = "l", xlab = "", ylab = "", col = "#31688E", ylim = c(min(rfrnc.pmf, sgn.pmf), max(rfrnc.pmf, sgn.pmf) * 1.1), cex.main = 0.8, lwd = 1.2, yaxt = "n")
      
        # add x axis label
        mtext(text = "Time (s)", side = 1, line = 2.5)
        
        # add title
        mtext(text = paste("Sound ID:", X$sound.id[x]), side = 3, line = 3, cex = 1)
        mtext(text = paste("Reference:", rfrnc), side = 3, line = 1.75, col = "#31688E", cex = 1)
        mtext(text = paste("Sound:", sgnl), side = 3, line = 0.5, col = "#B4DE2C", cex = 1)
            
        # add y axis
        axis(side = 4)
        mtext(text = "Amplitude (PMF)", side = 4, line = 2.5)
        
        # add sound envelope
        lines(time.vals, sgn.pmf, col= "#B4DE2CFF", lwd = 1.2)
        
        # sound envelope on top
        polygon(x = c(time.vals, rev(time.vals)), y = c(sgn.pmf, rev(rfrnc.pmf)), col = "#FDE72533", border = NA)
        
        # get plotting area limits
        usr <- par("usr")
        
        # and blu ratio value
        text(x = ((usr[1] + usr[2]) / 2) + usr[1], y = usr[4] * 0.9, paste("Blur ratio:", round(bl.rt, 2)), cex = 1)
        
        # index of reference
        rf.indx <- which(paste(X$sound.files, X$selec, sep = "-") == rfrnc)

        # freq limit of reference
        flim <- c(X$bottom.freq[rf.indx], X$top.freq[rf.indx])
        
        #####
        # end for sound and reference
        rf.info <- warbleR::read_sound_file(X = X, index = rf.indx, header = TRUE, path = path)
        rf.dur <- rf.info$samples / rf.info$sample.rate
        
        sgnl.info <- warbleR::read_sound_file(X = X, index = x, header = TRUE, path = path)
        sgnl.dur <- sgnl.info$samples / sgnl.info$sample.rate
        
        # calculate margin for spectrogram, before and after
        mar.rf.af <- mar.rf.bf <- (X$end[rf.indx] - X$start[rf.indx]) / 4
        
        # start for sound and reference
        strt.sgnl <- X$start[x] - mar.rf.bf
        if (strt.sgnl < 0) strt.sgnl <- 0
        strt.rf <- X$start[rf.indx] - mar.rf.bf
        if (strt.rf < 0) strt.rf <- 0
        
        end.sgnl <- X$end[x] + mar.rf.af
        if (end.sgnl > sgnl.dur) end.sgnl <- sgnl.dur
        end.rf <- X$end[rf.indx] + mar.rf.af
        if (end.rf > rf.dur) end.rf <- rf.dur
        
        # extract clip reference and sound
        clp.sgnl <- warbleR::read_sound_file(X = X, index = x, from = strt.sgnl, to = end.sgnl, path = path)
        clp.rfnc <- warbleR::read_sound_file(X = X, index = rf.indx, from = strt.rf, to = end.rf, path = path)
        
        ## plot spectros
        # sound at bottom left
        screen(1)
        par(mar = c(0.3, 0.3, 0.15, 0.3))
        
        warbleR:::spectro_wrblr_int2(wave = clp.sgnl, f = clp.sgnl@samp.rate, 
                                     flim = flim,  axisX = FALSE, axisY = FALSE, 
                                     tlab = NULL, flab = NULL, main = NULL, grid = FALSE, rm.zero = TRUE, cexaxis = 1.2, add = TRUE, ovlp = ovlp, wl = wl, collevels = collevels, palette = pal)
        
        # lines showing position of sound
        abline(v = c(mar.rf.bf, X$end[x] - X$start[x] + mar.rf.bf), col = "#B4DE2CFF", lty = 2)
                
        # add box with sound color
        box(col = "#B4DE2CFF", lwd = 3)
        
        # reference at top left
        screen(2)
        par(mar = c(0.15, 0.3, 0.3, 0.3))
        
        warbleR:::spectro_wrblr_int2(wave = clp.rfnc, f = clp.rfnc@samp.rate, 
           flim = flim, axisX = FALSE, axisY = FALSE, 
           tlab = NULL, flab = NULL, main = NULL, grid = FALSE, rm.zero = TRUE, cexaxis = 1.2, add = TRUE, ovlp = ovlp, wl = wl, collevels = collevels, palette = pal)

        # lines showing position of sound
        abline(v = c(mar.rf.bf, X$end[rf.indx] - X$start[rf.indx] + mar.rf.bf), col = "#31688ECC", lty = 2)
        
        # add box with reference color
        box(col = "#31688E", lwd = 3)
                
        # close graph    
        dev.off()        
      }
      
      # return maximum correlation
      return(bl.rt)
    }
    return(out)
  } 
  
  if (pb & !img) write(file = "", x = "Calculating blur ratio (step 3 out of 3):")
  if (pb & img) write(file = "", x = "Calculating blur ratio and producing images (step 3 out of 3):")
    
  # get blur ratio
  # calculate all envelops apply function
  X$blur.ratio <- pbapply::pbsapply(X = 1:nrow(X), cl = cl, FUN = function(x, rs = res, wle = wl, colvs = collevels, pl = pal, ovp = ovlp)   {
    blur_FUN(x, res = rs, ovlp = ovp, wl = wle, collevels = colvs, pal = pl)
  }) 
  
  # remove temporal column
  X$TEMP....sgnl <- NULL
  
  # convert to list instead of extended selection table, add envelopes
  if (output == "list") 
  {
    
    env.dfs <- warbleR:::pblapply_wrblr_int(pbar = pb, 1:length(envs), cl = cl, function(y){
      
      # extract 1 envelope
      x <- envs[[y]]
      
      # convert envelopes to PMF (probability mass function)
      x <- x / sum(x)       
      
      # put in data framme
      out <- data.frame(sound = names(envs)[y], sound.id= X$sound.id[paste(X$sound.files, X$selec, sep = "-") == names(envs)[y]], distance  = X$distance[paste(X$sound.files, X$selec, sep = "-") == names(envs)[y]], time = seq(from = 0, to = length(x) / (attr(X, "check.results")$sample.rate[1] * 1000), along.with =  x), amp = x)
      
      return(out)
    })
    
    # put together in a single data frame
    env.df <- do.call(rbind, env.dfs)
    
    # put est and envelopes in a list
    X <- list(est = X, envelopes = env.df)
    }
  
  # return data frame
  if (output == "data.frame") X <- as.data.frame(X) else
    attributes(X)$call <- base::match.call() # fix call attribute
  
  return(X)
}
