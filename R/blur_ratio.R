#' Measure blur ratio in the time domain
#'
#' \code{blur_ratio} measures blur ratio in sounds referenced in an extended selection table.
#' @usage blur_ratio(X, parallel = NULL, cores = getOption("mc.cores", 1),
#' pb = getOption("pb", TRUE), method = getOption("method", 1), ssmooth = 200,
#' msmooth = NULL, output = NULL, envelopes = FALSE, img = FALSE, res = 150,
#' hop.size = getOption("hop.size", 11.6), wl = getOption("wl", NULL),
#' ovlp = getOption("ovlp", 70), palette = viridis::viridis, collevels = seq(-120, 0, 5),
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
#' @param output DEPRECATED. Now the output format mirrors the class of the input 'X'. To obtain the amplitude envelopes use 'envelopes = TRUE'.
#' @param envelopes Logical to control if envelopes are returned (as attributes, 'attributes(X)$envelopes'). Default is \code{FALSE}.
#' @param img Logical argument to control if image files in 'jpeg' format containing the images being compared and the corresponding envelopes are produced. Default is no images (\code{FALSE}).
#' @param res Numeric argument of length 1. Controls image resolution. Default is 150 (faster) although 300 - 400 is recommended for publication/presentation quality.
#' @param hop.size A numeric vector of length 1 specifying the time window duration (in ms). Default is 11.6 ms, which is equivalent to 512 wl for a 44.1 kHz sampling rate. Ignored if 'wl' is supplied.
#' @param wl A numeric vector of length 1 specifying the window length of the spectrogram, default
#' is NULL. If supplied, 'hop.size' is ignored.
#' @param ovlp Numeric vector of length 1 specifying the percent overlap between two
#'   consecutive windows, as in \code{\link[seewave]{spectro}}. Only used when plotting. Default is 70. Applied to both spectra and spectrograms on image files.
#' @param palette A color palette function to be used to assign colors in the
#'   plot, as in \code{\link[seewave]{spectro}}. Default is \code{\link[viridis]{viridis}}.
#' @param collevels	Numeric vector indicating a set of levels which are used to partition the amplitude range of the spectrogram (in dB) as in \code{\link[seewave]{spectro}}. Default is \code{seq(-120, 0, 5)}.
#' @param dest.path Character string containing the directory path where the image files will be saved. If NULL (default) then the folder containing the sound files will be used instead.
#' @param path Character string containing the directory path where the sound files are found. Only needed when 'X' is not an extended selection table.
#' @return Object 'X' with two additional columns, 'reference' and 'blur.ratio', containing containing the id of the sound used as reference and the computed blur ratio values, respectively. If \code{img = TRUE} it also returns 1 image file (in 'jpeg' format) for each comparison showing spectrograms of both sounds and the overlaid amplitude envelopes (as probability mass functions (PMF)). Spectrograms are shown within the frequency range of the reference sound and also show vertical lines with the start and end of sounds to allow users to visually check alignment. If \code{envelopes = TRUE} the output would include amplitude envelopes for all sounds as attributes ('attributes(X)$envelopes').
#' @export
#' @name blur_ratio
#' @details Blur ratio measures the degradation of sound as a change in sound energy in the time domain as described by Dabelsteen et al (1993). Low values indicate low degradation of sounds. The function measures the blur ratio on sounds in which a reference playback has been re-recorded at different distances. Blur ratio is measured as the mismatch between amplitude envelopes (expressed as probability mass functions) of the reference sound and the re-recorded sound. By converting envelopes to probability mass functions the effect of energy attenuation is removed, focusing the analysis on the modification of the envelope shape. The function compares each sound to the corresponding reference sound within the supplied frequency range (e.g. bandpass) of the reference sound ('bottom.freq' and 'top.freq' columns in 'X'). The 'sound.id' column must be used to tell the function to only compare sounds belonging to the same category (e.g. song-types). Two methods for setting the experimental design are provided. All wave objects in the extended selection table must have the same sampling rate so the length of envelopes is comparable.
#' @seealso \code{\link{envelope_correlation}}, \code{\link{spectrum_blur_ratio}}
#' @examples
#' {
#'   # load example data
#'   data("degradation_est")
#'    
#'   # create subset of data with only re-recorded files
#'   rerecorded_est <- degradation_est[degradation_est$sound.files != "master.wav", ]
#'   
#'   # using method 1
#'   blur_ratio(X = rerecorded_est)
#'
#'   # using method 2
#'   # blur_ratio(X = rerecorded_est, method = 2)
#'
#'   # get envelopes
#'   br <- blur_ratio(X = rerecorded_est, envelopes = TRUE)
#'   envs <- attributes(br)$envelopes
#'
#'   # make distance a factor for plotting
#'   envs$distance <- as.factor(envs$distance)
#'
#'   \dontrun{
#'   # plot
#'   ggplot(envs, aes(x= time, y = amp, col = distance)) +
#'   geom_line() + facet_wrap(~ sound.id) +
#'   scale_color_manual(values = viridis::viridis(4)) +
#'   labs(x = "Time (s)", y = "Amplitude (PMF)") +
#'   theme_classic()
#'   }
#' }
#'
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})
#' @references {
#' Dabelsteen, T., Larsen, O. N., & Pedersen, S. B. (1993). Habitat-induced degradation of sound signals: Quantifying the effects of communication sounds and bird location on blur ratio, excess attenuation, and signal-to-noise ratio in blackbird song. The Journal of the Acoustical Society of America, 93(4), 2206.
#'
#' Araya-Salas, M. (2020). baRulho: baRulho: quantifying habitat-induced degradation of (animal) acoustic signals in R. R package version 1.0.2
#' }

blur_ratio <-
  function(X,
           parallel = NULL,
           cores = getOption("mc.cores", 1),
           pb = getOption("pb", TRUE),
           method = getOption("method", 1),
           ssmooth = 200,
           msmooth = NULL,
           output = NULL,
           envelopes = FALSE,
           img = FALSE,
           res = 150,
           hop.size = getOption("hop.size", 11.6),
           wl = getOption("wl", NULL),
           ovlp = getOption("ovlp", 70),
           palette = viridis::viridis,
           collevels = seq(-120, 0, 5),
           dest.path = NULL,
           path = getOption("sound.files.path", ".")) {
    # check arguments
    arguments <- as.list(base::match.call())
    
    # add objects to argument names
    for (i in names(arguments)[-1]) {
      arguments[[i]] <- get(i)
    }
    
    # check each arguments
    check_results <-
      check_arguments(fun = arguments[[1]], args = arguments)
    
    # report errors
    report_assertions2(check_results)
    
    # make path null if extendeed selection table
    if (is_extended_selection_table(X)) {
      path <- NULL
    }
    
    # adjust wl based on hope.size
    if (is.null(wl)) {
      wl <-
        round(
          read_sound_file(
            X,
            index = 1,
            header = TRUE,
            path = path
          )$sample.rate * hop.size / 1000,
          0
        )
    }
    
    # make wl even if odd
    if (!(wl %% 2) == 0)
      wl <- wl + 1
    
    # set clusters for windows OS
    if (Sys.info()[1] == "Windows" & cores > 1) {
      cl <- parallel::makePSOCKcluster(getOption("cl.cores", cores))
    } else {
      cl <- cores
    }
    
    # ingnore ssmooth if msmooth is supplied
    if (!is.null(msmooth)) {
      ssmooth <- NULL
    }
    
    # add sound file selec column and names to X (weird column name so it does not overwrite user columns)
    if (pb) {
      write(file = "",
            x = paste0("Preparing data for analysis (step 1 out of 3):"))
    }
    
    X <- prep_X_bRlo_int(X,
                         method = method,
                         cores = cores,
                         pb = pb)
    
    # print message
    if (pb)
      write(file = "", x = "Calculating amplitude envelopes (step 2 out of 3):")
    
    # calculate all envelops apply function
    envs <-
      warbleR:::pblapply_wrblr_int(
        pbar = pb,
        X = seq_len(nrow(X)),
        cl = cl,
        FUN = function(y,
                       ssmth = ssmooth,
                       msmth = msmooth,
                       ov = ovlp) {
          # load clip
          clp <- warbleR::read_sound_file(X = X,
                                          index = y,
                                          path = path)
          
          # define bandpass based on reference
          bp <-
            c(X$bottom.freq[X$.sgnl.temp == X$reference[y]], X$top.freq[X$.sgnl.temp == X$reference[y]])
          
          # bandpass filter
          clp <- seewave::ffilter(
            clp,
            from = bp[1] * 1000,
            ovlp = ov,
            to = bp[2] * 1000,
            bandpass = TRUE,
            wl = wl,
            output = "Wave"
          )
          
          # calculate envelope
          nv <-
            seewave::env(
              wave = clp,
              f = clp@samp.rate,
              ssmooth = ssmth,
              msmooth = msmth,
              plot = FALSE,
              envt = "hil"
            )[, 1]
          
          return(nv)
        }
      )
    
    # add sound file selec column and names to envelopes
    names(envs) <- X$.sgnl.temp
    
    # set options to run loops
    if (pb &
        !img)
      write(file = "", x = "Calculating blur ratio (step 3 out of 3):")
    if (pb &
        img)
      write(file = "", x = "Calculating blur ratio and producing images (step 3 out of 3):")
    
    # get blur ratio
    # calculate all envelops apply function
    X$blur.ratio <-
      pbapply::pbsapply(
        X = seq_len(nrow(X)),
        cl = cl,
        FUN = function(x,
                       rs = res,
                       wle = wl,
                       colvs = collevels,
                       pl = palette,
                       ovp = ovlp) {
          blur_FUN(
            X,
            envs,
            path,
            img,
            dest.path,
            x,
            res = rs,
            ovlp = ovp,
            wl = wle,
            collevels = colvs,
            palette = pl
          )
        }
      )
    
    # make NAs those sounds in which the reference is itself (only happens when method = 2)
    X$reference[X$reference ==  X$.sgnl.temp] <- NA
    
    # remove temporal column
    X$.sgnl.temp <- NULL
    
    # convert to list instead of extended selection table, add envelopes
    if (envelopes) {
      env.dfs <-
        warbleR:::pblapply_wrblr_int(pbar = pb, 1:length(envs), cl = cl, function(y) {
          # extract 1 envelope
          x <- envs[[y]]
          
          # convert envelopes to PMF (probability mass function)
          x <- x / sum(x)
          
          # put in data framme
          out <-
            data.frame(
              sound = names(envs)[y],
              sound.id = X$sound.id[paste(X$sound.files, X$selec, sep = "-") == names(envs)[y]],
              distance = X$distance[paste(X$sound.files, X$selec, sep = "-") == names(envs)[y]],
              time = seq(
                from = 0,
                to = length(x) / (attr(X, "check.results")$sample.rate[1] * 1000),
                along.with = x
              ),
              amp = x
            )
          
          return(out)
        })
      
      # put together in a single data frame
      env.df <- do.call(rbind, env.dfs)
      
      # add envelopes as attributes
      attributes(X)$envelopes <- env.df
    }
    
    # return data frame
    if (is_extended_selection_table(X) | is_selection_table(X)) {
      attributes(X)$call <- base::match.call()
    }
    return(X)
  }
