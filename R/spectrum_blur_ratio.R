#' Measure blur ratio in the frequency domain 
#' 
#' \code{spectrum_blur_ratio} measures blur ratio of frequency spectra from sounds referenced in an extended selection table.
#' @usage spectrum_blur_ratio(X, parallel = 1, cores = 1, pb = TRUE, method = 1, ssmooth = 50, 
#' output = "est", img = FALSE, res = 150, hop.size = 11.6, wl = NULL, 
#' ovlp = 70, pal = viridis, collevels = seq(-120, 0, 5), dest.path = NULL, path = NULL)
#' @param X Object of class 'data.frame', 'selection_table' or 'extended_selection_table' (the last 2 classes are created by the function \code{\link[warbleR]{selection_table}} from the warbleR package) with the reference to the sounds in the master sound file. Must contain the following columns: 1) "sound.files": name of the .wav files, 2) "selec": unique selection identifier (within a sound file), 3) "start": start time and 4) "end": end time of selections, 5)  "bottom.freq": low frequency for bandpass, 6) "top.freq": high frequency for bandpass and 7) "sound.id": ID of sounds used to identify counterparts across distances. Each sound must have a unique ID within a distance.
#' @param parallel DEPRECATED. Use 'cores' instead.
#' @param cores Numeric vector of length 1. Controls whether parallel computing is applied by specifying the number of cores to be used. Default is 1 (i.e. no parallel computing).
#' @param pb Logical argument to control if progress bar is shown. Default is \code{TRUE}.
#' @param method Numeric vector of length 1 to indicate the 'experimental design' for measuring spectrum correlation. Two methods are available:
#' \itemize{
#' \item \code{1}: compare all sounds with their counterpart that was recorded at the closest distance to source (e.g. compare a sound recorded at 5m, 10m and 15m with its counterpart recorded at 1m). This is the default method. 
#' \item \code{2}: compare all sounds with their counterpart recorded at the distance immediately before (e.g. a sound recorded at 10m compared with the same sound recorded at 5m, then sound recorded at 15m compared with same sound recorded at 10m and so on).
#' }
#' @param ssmooth Numeric vector of length 1 determining the length of the sliding window used for a sum smooth for power spectrum calculation (in kHz). Default is 100.
#' @param output Character vector of length 1 to determine if an extended selection table ('est', default), a data frame ('data.frame') or a list ("list") containing the extended selection table (first object in the list) and all (smoothed) wave envelopes (second object in the list) is returned. The envelope data can be used for plotting.
#' @param img Logical argument to control if image files in 'jpeg' format containing the images being compared and the corresponding spectra are produced. Default is no images ( \code{FALSE}).
#' @param res Numeric argument of length 1. Controls image resolution. Default is 150 (faster) although 300 - 400 is recommended for publication/presentation quality.
#' @param hop.size A numeric vector of length 1 specifying the time window duration (in ms). Default is 11.6 ms, which is equivalent to 512 wl for a 44.1 kHz sampling rate. Ignored if 'wl' is supplied.
#' @param wl A numeric vector of length 1 specifying the window length of the spectrogram, default 
#' is NULL. If supplied, 'hop.size' is ignored. Applied to both spectra and spectrograms on image files.
#' @param ovlp Numeric vector of length 1 specifying the percent overlap between two 
#'   consecutive windows, as in \code{\link[seewave]{spectro}}. Default is 70. Applied to both spectra and spectrograms on image files.
#' @param pal A color palette function to be used to assign colors in the 
#'   plot, as in \code{\link[seewave]{spectro}}. Default is \code{\link[viridis]{viridis}}. 
#' @param collevels	Numeric vector indicating a set of levels which are used to partition the amplitude range of the spectrogram (in dB) as in \code{\link[seewave]{spectro}}. Default is \code{seq(-120, 0, 5)}. 
#' @param dest.path Character string containing the directory path where the image files will be saved. If NULL (default) then the folder containing the sound files will be used instead.
#' @param path Character string containing the directory path where the sound files are found. Only needed when 'X' is not an extended selection table.
#' @return Data frame similar to input data, but also includes a new column (spectral.blur.ratio)
#' with the blur ratio values. If \code{img = TRUE} it also returns 1 image file (in 'jpeg' format) for each comparison showing spectrograms of both sounds and the overlaid power spectrum (as probability mass functions (PMF)). Spectrograms are shown within the frequency range of the reference sound and also show vertical lines with the start and end of sounds to allow users to visually check alignment. If \code{output = 'list'} the output would a list including the data frame just described and a data frame with spectra (amplitude values) for all sounds. 
#' @export
#' @name spectrum_blur_ratio
#' @details Spectral blur ratio measures the degradation of sound as a function of the change in sound energy in the frequency domain, analogous to the blur ratio proposed by Dabelsteen et al (1993) for the time domain (and implemented in \code{\link{blur_ratio}}). Low values indicate low degradation of sounds. The function measures the blur ratio of spectra from sounds in which a reference playback has been re-recorded at different distances. Spectral blur ratio is measured as the mismatch between power spectra (expressed as probability density functions) of the reference sound and the re-recorded sound. The function compares each sound type to the corresponding reference sound. The 'sound.id' column must be used to tell the function to only compare sounds belonging to the same category (e.g. song-types). Two methods for setting the experimental design are provided. All wave objects in the extended selection table must have the same sampling rate so the length of spectra is comparable.   
#' @examples
#' {
#' # load example data
#' data("playback_est")
#' 
#' # remove ambient selections
#' playback_est <- playback_est[playback_est$sound.id != "ambient", ]
#' 
#' # using method 1
#' spectrum_blur_ratio(X = playback_est)
#' 
#' # using method 2
#' spectrum_blur_ratio(X = playback_est, method = 2)
#' }
#' 
#' @seealso \code{\link{blur_ratio}}
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr}) 
#' @references {
#' Dabelsteen, T., Larsen, O. N., & Pedersen, S. B. (1993). Habitat-induced degradation of sound signals: Quantifying the effects of communication sounds and bird location on blur ratio, excess attenuation, and signal-to-noise ratio in blackbird song. The Journal of the Acoustical Society of America, 93(4), 2206.
#' 
#' Araya-Salas, M. (2020). baRulho: baRulho: quantifying habitat-induced degradation of (animal) acoustic signals in R. R package version 1.0.2
#' }

spectrum_blur_ratio <-
  function(X,
           parallel = 1, 
           cores = 1,
           pb = TRUE,
           method = 1,
           ssmooth = 50,
           output = "est",
           img = FALSE,
           res = 150,
           hop.size = 11.6,
           wl = NULL,
           ovlp = 70,
           pal = viridis,
           collevels = seq(-120, 0, 5),
           dest.path = NULL,
           path = NULL) {
    
    # set dest.path if not provided
    if (is.null(dest.path))
      dest.path <- getwd() else
      if (!dir.exists(dest.path))
        stop2("'dest.path' provided does not exist")
    
    # set path if not provided
    if (is.null(path))
      path <- getwd() else
      if (!dir.exists(path))
        stop2("'path' provided does not exist")
    
    # must have the same sampling rate
    if (is_extended_selection_table(X)) {
      if (length(unique(attr(X, "check.results")$sample.rate)) > 1)
        stop2(
          "all wave objects in the extended selection table must have the same sampling rate (they can be homogenized using warbleR::resample_est())"
        )
    } else
      print("assuming all sound files have the same sampling rate")
    
    # get sampling rate
    sampling_rate <-
      warbleR::read_sound_file(
        X = X,
        index = 1,
        path = path,
        header = TRUE
      )$sample.rate
    
    # If cores is not numeric
    if (!is.numeric(cores))
      stop2("'cores' must be a numeric vector of length 1")
    if (any(!(cores %% 1 == 0), cores < 1))
      stop2("'cores' should be a positive integer")
    
    # hopsize
    if (!is.numeric(hop.size) |
        hop.size < 0)
      stop2("'hop.size' must be a positive number")
    
    # adjust wl based on hope.size
    if (is.null(wl))
      wl <-
      round(read_sound_file(
        X,
        index = 1,
        header = TRUE,
        path = path
      )$sample.rate * hop.size  / 1000,
      0)
    
    # make wl even if odd
    if (!(wl %% 2) == 0)
      wl <- wl + 1
    
    # If method is not numeric
    if (!is.numeric(method))
      stop2("'method' must be a numeric vector of length 1")
    if (!any(method %in% 1:2))
      stop2("'method' must be either 1 or 2")
    
    # check sound.id column
    if (is.null(X$sound.id))
      stop2("'X' must contain a 'sound.id' column")
    
    #check output
    if (!any(output %in% c("est", "data.frame", "list")))
      stop2("'output' must be 'est', 'data.frame' or 'list'")
    
    # add sound file selec column and names to X (weird column name so it does not overwrite user columns)
    if (pb)
      write(file = "",
            x = paste0("Preparing data for analysis (step 1 out of 3):"))
    
    X <-
      prep_X_bRlo_int(X,
                      method = method,
                      cores = cores,
                      pb = pb)
    
    # set clusters for windows OS
    if (Sys.info()[1] == "Windows" & cores > 1)
      cl <-
      parallel::makePSOCKcluster(getOption("cl.cores", cores)) else
      cl <- cores
    
    # print message
    if (pb)
      write(file = "", x = "Calculating power spectra (step 2 out of 3):")
    
    # calculate all spectra apply function
    specs.list <-
      warbleR:::pblapply_wrblr_int(
        pbar = pb,
        X = 1:nrow(X),
        cl = cl,
        FUN = function(y, wl)   {
          # load clip
          clp <- warbleR::read_sound_file(X = X,
                                          index = y,
                                          path = path)
          
          # calculate spectrum
          clp.spc <-
            seewave::spec(
              wave = clp,
              f = clp@samp.rate,
              plot = FALSE,
              wl = wl
            )
          
          clp.spc[, 1] <-
            
            return(clp.spc)
        }
      )
    
    # add sound file selec names to spectra
    names(specs.list) <- X$TEMP....sgnl
    
    ## function to measure blur ratio
    # y and z are the sound.files+selec names of the sounds and reference sound (model)
    # spectrum mismatch ratio
    blur_sp_FUN <- function(x, res, ovlp, wl, collevels, pal, ...) {
      # get names of sound and reference
      sgnl <- X$TEMP....sgnl[x]
      rfrnc <- X$reference[x]
      
      # if sounds are the same or the selection is noise return NA
      if (sgnl == rfrnc |
          any(c(X$sound.id[X$TEMP....sgnl == sgnl], X$sound.id[X$reference == rfrnc]) == "ambient"))
        out <- NA else {
        # extract spectrum for sound and model
        sgnl.spc <- specs.list[[which(names(specs.list) == sgnl)]]
        rfrnc.spc <- specs.list[[which(names(specs.list) == rfrnc)]]
        
        # make them the same frequency range as reference
        bp <-
          c(X$bottom.freq[X$TEMP....sgnl == rfrnc], X$top.freq[X$TEMP....sgnl == rfrnc])
        
        bp <- bp + c(-0.2, 0.2)       # add 0.2 kHz buffer
        if (bp[1] < 0)
          # force 0 if negative
          bp[1] <- 0
        if (bp[2] > ceiling(sampling_rate / 2000) - 1)
          bp[2] <-
          ceiling(sampling_rate / 2000) - 1 # force lower than nyquist freq if higher
        
        # apply bandpass by shrinking freq range and remove freq column based on reference freq bins
        sgnl.spc <-
          sgnl.spc[rfrnc.spc[, 1] > bp[1] & rfrnc.spc[, 1] < bp[2], 2]
        rfrnc.spc <-
          rfrnc.spc[rfrnc.spc[, 1] > bp[1] & rfrnc.spc[, 1] < bp[2], 2]
        
        # applied ssmooth
        if (!is.null(ssmooth))
          if(length(rfrnc.spc) > ssmooth + 2){
          sgnl.spc <- as.matrix(seewave::sumsmooth(sgnl.spc, wl = ssmooth))
          rfrnc.spc <-
            as.matrix(seewave::sumsmooth(rfrnc.spc, wl = ssmooth))
        }
        
        # convert envelopes to PMF (probability mass function)
        rfrnc.pmf <- rfrnc.spc / sum(rfrnc.spc)
        sgnl.pmf <- sgnl.spc / sum(sgnl.spc)
        
        # get blur ratio as half the sum of absolute differences between spectra PMFs
        bl.rt <- sum(abs(rfrnc.pmf - sgnl.pmf)) / 2
        
        # plot
        if (img)
        {
          img_bRlo_int(
            filename = paste0(
              "blur_ratio_",
              X$sound.id[x],
              "-",
              rfrnc,
              "-",
              sgnl,
              ".jpeg"
            ),
            path = dest.path,
            width = 10.16 * 1.5,
            height = 10.16 ,
            units = "cm",
            res = res
          )
          
          # create time values for area calculation
          f.vals <- seq(bp[1], bp[2], length.out = length(rfrnc.pmf))
          
          # difference between spectra
          spc.diff <- rfrnc.pmf - sgnl.pmf
          
          # matrix for layout
          ly.mat <- matrix(c(0, 0.3, 0, 0.5, # bottom left spectrogram
                             0, 0.3, 0.5, 1, # top left spectrogram
                             0.2, 1, 0, 1),
                           # right pannel spectra
                           nrow = 3,
                           byrow = TRUE)
          
          # save par settings
          oldpar <- par(no.readonly = TRUE)
          on.exit(par(oldpar))
          
          # close if open any screen
          invisible(close.screen(all.screens = TRUE))
          
          split.screen(ly.mat)
          
          # plot spectra
          screen(3)
          
          # set image margins
          par(mar = rep(4, 0, 4, 4))
          
          # reference spectrum first
          plot(
            x = rfrnc.pmf,
            y = f.vals,
            type = "l",
            xlab = "",
            ylab = "",
            col = "#31688E",
            xlim = c(
              min(rfrnc.pmf, sgnl.pmf),
              max(rfrnc.pmf, sgnl.pmf) * 1.1
            ),
            cex.main = 0.8,
            lwd = 1.2,
            yaxt = "n"
          )
          
          # add x axis label
          mtext(text = "Power spectrum (PMF)",
                side = 1,
                line = 2.5)
          
          # add title
          mtext(
            text = paste("Sound ID:", X$sound.id[x]),
            side = 3,
            line = 3,
            cex = 1
          )
          mtext(
            text = paste("Reference:", rfrnc),
            side = 3,
            line = 1.75,
            col = "#31688E",
            cex = 1
          )
          mtext(
            text = paste("Sound ID:", sgnl),
            side = 3,
            line = 0.5,
            col = "#B4DE2C",
            cex = 1
          )
          
          # add y axis
          axis(side = 4)
          mtext(text = "Frequency (kHz)",
                side = 4,
                line = 2.5)
          
          # add sound spectrum
          lines(sgnl.pmf, f.vals, col = "#B4DE2C", lwd = 1.2)
          
          # sound spectrum on top
          polygon(
            y = c(f.vals, rev(f.vals)),
            x = c(sgnl.pmf, rev(rfrnc.pmf)),
            col =  "#FDE72533",
            border = NA
          )
          
          # get plotting area limits
          usr <- par("usr")
          
          # and blu ratio value
          text(
            x = ((usr[1] + usr[2]) / 2) + usr[1],
            y = usr[4] * 0.9,
            paste("Blur ratio of spectrum:", round(bl.rt, 2)),
            cex = 1
          )
          
          # index of reference
          rf.indx <-
            which(paste(X$sound.files, X$selec, sep = "-") == rfrnc)
          
          # freq limit of reference
          flim <- c(X$bottom.freq[rf.indx], X$top.freq[rf.indx])
          
          
          # end for sound and reference
          rf.info <-
            warbleR::read_sound_file(
              X = X,
              index = rf.indx,
              header = TRUE,
              path = path
            )
          rf.dur <- rf.info$samples / rf.info$sample.rate
          
          sgnl.info <-
            warbleR::read_sound_file(
              X = X,
              index = x,
              header = TRUE,
              path = path
            )
          sgnl.dur <- sgnl.info$samples / sgnl.info$sample.rate
          
          # calculate margin for spectrogram, before and after
          mar.rf.af <-
            mar.rf.bf <- (X$end[rf.indx] - X$start[rf.indx]) / 4
          
          # start for sound and reference
          strt.sgnl <- X$start[x] - mar.rf.bf
          if (strt.sgnl < 0)
            strt.sgnl <- 0
          strt.rf <- X$start[rf.indx] - mar.rf.bf
          if (strt.rf < 0)
            strt.rf <- 0
          
          
          end.sgnl <- X$end[x] + mar.rf.af
          if (end.sgnl > sgnl.dur)
            end.sgnl <- sgnl.dur
          end.rf <- X$end[rf.indx] + mar.rf.af
          if (end.rf > rf.dur)
            end.rf <- rf.dur
          
          # extract clip reference and sound
          clp.sgnl <-
            warbleR::read_sound_file(
              X = X,
              index = x,
              from = strt.sgnl,
              to = end.sgnl,
              path = path
            )
          clp.rfnc <-
            warbleR::read_sound_file(
              X = X,
              index = rf.indx,
              from = strt.rf,
              to = end.rf,
              path = path
            )
          
          ## plot spectros
          # sound at bottom left
          screen(1)
          par(mar = c(0.3, 0.3, 0.15, 0.3))
          
          warbleR:::spectro_wrblr_int2(
            wave = clp.sgnl,
            f = clp.sgnl@samp.rate,
            flim = flim,
            axisX = FALSE,
            axisY = FALSE,
            tlab = NULL,
            flab = NULL,
            main = NULL,
            grid = FALSE,
            rm.zero = TRUE,
            cexaxis = 1.2,
            add = TRUE,
            ovlp = ovlp,
            wl = wl,
            collevels = collevels,
            palette = pal
          )
          
          # lines showing position of sound
          abline(
            v = c(mar.rf.bf, X$end[x] - X$start[x] + mar.rf.bf),
            col = "#B4DE2CFF",
            lty = 2
          )
          
          # add box with sound color
          box(col = "#B4DE2C", lwd = 3)
          
          # reference at top left
          screen(2)
          par(mar = c(0.15, 0.3, 0.3, 0.3))
          
          warbleR:::spectro_wrblr_int2(
            wave = clp.rfnc,
            f = clp.rfnc@samp.rate,
            flim = flim,
            axisX = FALSE,
            axisY = FALSE,
            tlab = NULL,
            flab = NULL,
            main = NULL,
            grid = FALSE,
            rm.zero = TRUE,
            cexaxis = 1.2,
            add = TRUE,
            ovlp = ovlp,
            wl = wl,
            collevels = collevels,
            palette = pal
          )
          
          # lines showing position of sound
          abline(
            v = c(mar.rf.bf, X$end[rf.indx] - X$start[rf.indx] + mar.rf.bf),
            col = "#31688ECC",
            lty = 2
          )
          
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
    
    if (pb & !img)
      write(file = "", x = "Calculating spectrum blur ratio (step 3 out of 3):")
    if (pb & img)
      write(file = "", x = "Calculating blur ratio and producing images (step 3 out of 3):")
    
    # get blur ratio
    # calculate all spectra apply function
    X$spectral.blur.ratio <-
      pbapply::pbsapply(
        X = 1:nrow(X),
        cl = cl,
        FUN = function(x,
                       rs = res,
                       wle = wl,
                       colvs = collevels,
                       pl = pal,
                       ovp = ovlp,
                       ...)   {
          blur_sp_FUN(
            x,
            res = rs,
            ovlp = ovp,
            wl = wle,
            collevels = colvs,
            pal = pl,
            ...
          )
        }
      )
    
    # remove temporal columns
    X$TEMP....sgnl <- NULL
    
    # convert to list instead of extended selection table, add envelopes
    if (output == "list")
    {
      spec.dfs <-
        lapply(1:length(specs.list), function(y, ssmth = ssmooth) {
          # extract 1 envelope
          x <- specs.list[[y]]
          
          # convert envelopes to PMF (probability mass function)
          x[, 2] <- x[, 2] / sum(x[, 2])
          
          # applied ssmooth
          if (!is.null(ssmooth))
            x[, 2] <- as.matrix(seewave::sumsmooth(x[, 2], wl = ssmth))
          
          # put in data framme
          out <-
            data.frame(
              sound = names(specs.list)[y],
              sound.id = X$sound.id[paste(X$sound.files, X$selec, sep = "-") == names(specs.list)[y]],
              distance  = X$distance[paste(X$sound.files, X$selec, sep = "-") == names(specs.list)[y]],
              freq = x[, 1],
              amp = x[, 2]
            )
          
          return(out)
        })
      
      # put together in a single data frame
      spec.df <- do.call(rbind, spec.dfs)
      
      # put est and envelopes in a list
      X <- list(est = X, spectra = spec.df)
    }
    
    # return data frame
    if (output == "data.frame")
      X <- as.data.frame(X) else
        attributes(X)$call <- base::match.call() # fix call attribute
    
    
    return(X)
  }
