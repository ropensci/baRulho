#' Measure excess attenuation
#' 
#' \code{excess_attenuation} measures excess attenuation in sounds referenced in an extended selection table.
#' @usage excess_attenuation(X, parallel = 1, cores = getOption("mc.cores", 1), 
#' pb = getOption("pb", TRUE), method = getOption("method", 1), type = "Dabelsteen",
#'  output = "est", hop.size = getOption("hop.size", 1), wl = getOption("wl", NULL), 
#'  ovlp = getOption("ovlp", 50), gain = 0, bp = "freq.range", 
#'  path = getOption("sound.files.path", "."))
#' @param X Object of class 'data.frame', 'selection_table' or 'extended_selection_table' (the last 2 classes are created by the function \code{\link[warbleR]{selection_table}} from the warbleR package) with the reference to the sounds in the master sound file. Must contain the following columns: 1) "sound.files": name of the .wav files, 2) "selec": unique selection identifier (within a sound file), 3) "start": start time and 4) "end": end time of selections, 5)  "bottom.freq": low frequency for bandpass, 6) "top.freq": high frequency for bandpass and 7) "sound.id": ID of sounds used to identify counterparts across distances. Each sound must have a unique ID within a distance.
#' @param parallel DEPRECATED. Use 'cores' instead.
#' @param cores Numeric vector of length 1. Controls whether parallel computing is applied by specifying the number of cores to be used. Default is 1 (i.e. no parallel computing).
#' @param pb Logical argument to control if progress bar is shown. Default is \code{TRUE}.
#' @param method Numeric vector of length 1 to indicate the 'experimental design' for measuring excess attenuation. Two methods are available:
#' \itemize{
#' \item \code{1}: compare all sounds with their counterpart that was recorded at the closest distance to source (e.g. compare a sound recorded at 5m, 10m and 15m with its counterpart recorded at 1m). This is the default method. 
#' \item \code{2}: compare all sounds with their counterpart recorded at the distance immediately before (e.g. a sound recorded at 10m compared with the same sound recorded at 5m, then sound recorded at 15m compared with same sound recorded at 10m and so on).
#' }
#' @param type Character vector of length 1 to indicate the 'type' of excess attenuation to be used. Two types are available:
#' \itemize{
#' \item \code{Dabelsteen}: as described by Dabelsteen & Mathevon (2002): -20 x log(k) - 6/(2 x re-recorded distance) + gain. 'k' is the ratio of the mean amplitude envelopes of the re-recorded and reference sounds. This is the default method. 
#' \item \code{Darden}: as described by Darden et al. (2008): microphone_gain - 20 x log(reference distance / re-recorded distance) - 20 x log(k). 'k' is the ratio of the mean amplitude envelopes of the re-recorded and reference sounds. Microphone gain is the microphone gain of the reference and re-recorded sounds. 
#' }
#' If gain is not supplied (see 'gain' argument) gain is set as 0, which results in a relative measure of excess attenuation comparable only within the same experiment or between experiments with the same recording equipment and recording volume.
#' @param output Character vector of length 1 to determine if an extended selection table ('est', default) or a data frame ('data.frame'). 'est' format only available if 'X' is itself an extended selection table.
#' @param hop.size A numeric vector of length 1 specifying the time window duration (in ms). Default is 1 ms, which is equivalent to ~45 wl for a 44.1 kHz sampling rate. Ignored if 'wl' is supplied.
#' @param wl A numeric vector of length 1 specifying the window length of the spectrogram, default 
#' is \code{NULL}. If supplied, 'hop.size' is ignored.
#' Note that lower values will increase time resolution, which is more important for amplitude ratio calculations. 
#' @param ovlp Numeric vector of length 1 specifying the percent overlap between two 
#'   consecutive windows, as in \code{\link[seewave]{spectro}}. Default is 50. Only used for bandpass filtering.
#' @param gain Numeric vector of length 1 with the combined gain of the microphone and recorder (in dB). Default is 0, which results in a relative measure of excess attenuation comparable only within the same experiment, but not across experiments. Only used for \code{type = "Marten"}.  
#' @param bp Numeric vector of length 2 giving the lower and upper limits of a frequency bandpass filter (in kHz). Alternatively, when set to 'freq.range' (default) the function uses the 'bottom.freq' and 'top.freq' as the bandpass range.
#' @param path Character string containing the directory path where the sound files are found. Only needed when 'X' is not an extended selection table.
#' @return A data frame, or extended selection table similar to input data (depending on argument 'output'), but also includes a new column (excess.attenuation)
#' with the excess attenuation values.
#' @export
#' @name excess_attenuation
#' @details Excess attenuation is the amplitude loss of a sound in excess due to spherical spreading (observed attenuation - expected attenuation). With every doubling of distance, sounds attenuate with a 6 dB loss of amplitude (Morton, 1975; Marten & Marler, 1977). Any additional loss of amplitude results in energy loss in excess of that expected to occur with distance via spherical spreading. So it represents energy loss due to additional factors like vegetation or atmospheric conditions (Wiley & Richards, 1978). Low values indicate little additional attenuation. 
#' The goal of the function is to measure the excess attenuation on sounds in which a reference playback has been re-recorded at increasing distances. The 'sound.id' column must be used to indicate which sounds belonging to the same category (e.g. song-types). The function will then compare each sound type to the corresponding reference sound. Two methods for calculating excess attenuation are provided (see 'method' argument). 
#' @examples
#' {
#' # load example data
#' data("playback_est")
#' 
#' # using method 1
#'excess_attenuation(X = playback_est)
#' 
#' # using method 2
#' excess_attenuation(X = playback_est, method = 2)
#' }
#' 
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})
#' @seealso \code{\link{spcc_distortion}}; \code{\link{envelope_correlation}}
#' @references {
#' Araya-Salas, M. (2020). baRulho: baRulho: quantifying habitat-induced degradation of (animal) acoustic signals in R. R package version 1.0.2
#' 
#' Dabelsteen, T., Larsen, O. N., & Pedersen, S. B. (1993). Habitat-induced degradation of sound signals: Quantifying the effects of communication sounds and bird location on blur ratio, excess attenuation, and signal-to-noise ratio in blackbird song. The Journal of the Acoustical Society of America, 93(4), 2206.
#' 
#' Dabelsteen, T., & Mathevon, N. (2002). Why do songbirds sing intensively at dawn?. Acta ethologica, 4(2), 65-72.
#' 
#' Darden, SK, Pedersen SB, Larsen ON, & Dabelsteen T. (2008). Sound transmission at ground level in a short-grass prairie habitat and its implications for long-range communication in the swift fox *Vulpes velox*. The Journal of the Acoustical Society of America, 124(2), 758-766.
#' 
#' Marten K, & Marler P. (1977). Sound transmission and its significance for animal vocalization. Behavioral Ecology and Sociobiology, 2(3), 271-290.
#' 
#' Morton ES. (1975). Ecological sources of selection on avian sounds. The American Naturalist, 109(965), 17-34.
#' 
#' Wiley, R., & Richards, D. G. (1978). Physical constraints on acoustic communication in the atmosphere: implications for the evolution of animal vocalizations. Behavioral Ecology and Sociobiology, 3(1), 69-94.
#' }
#last modification on jul-19-2021 (MAS)

excess_attenuation <-
  function(X,
           parallel = 1, cores = getOption("mc.cores", 1),
           pb = getOption("pb", TRUE),
           method = getOption("method", 1),
           type = "Dabelsteen",
           output = "est",
           hop.size = getOption("hop.size", 1),
           wl = getOption("wl", NULL),
           ovlp = getOption("ovlp", 50),
           gain = 0,
           bp = "freq.range",
           path = getOption("sound.files.path", ".")) {
    
    # deprecated message
    if (parallel > 1) 
      stop2("'parallel' has been deprecated, Use 'cores' instead")
    
    # set path if not provided
    if (is.null(path))
      path <- getwd() else
      if (!dir.exists(path))
        stop2("'path' provided does not exist")
    
    # If cores is not numeric
    if (!is.numeric(cores))
      stop2("'cores' must be a numeric vector of length 1")
    if (any(!(cores %% 1 == 0), cores < 1))
      stop2("'cores' should be a positive integer")
    
    #check output
    if (!any(output %in% c("est", "data.frame")))
      stop2("'output' must be 'est' or 'data.frame'")
    # hopsize
    if (!is.numeric(hop.size) |
        hop.size < 0)
      stop2("'hop.size' must be a positive number")
    
    # must have the same sampling rate
    if (is_extended_selection_table(X)){
      if (length(unique(attr(X, "check.results")$sample.rate)) > 1)
        stop2(
          "all wave objects in the extended selection table must have the same sampling rate (they can be homogenized using warbleR::resample_est())"
        )
      } else 
          print("assuming all sound files have the same sampling rate")
    
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
    
    # add sound file selec column and names to X (weird column name so it does not overwrite user columns)
    if (pb)
      write(file = "",
            x = paste0("Preparing data for analysis (step 1 out of 3):"))
    X <-
      prep_X_bRlo_int(X,
                      method = method,
                      cores = cores,
                      pb = pb)
    
    # function to extract mean envelopes
    meanenv_FUN <- function(y, wl, ovlp) {
      # read sound clip
      clp <-
        warbleR::read_sound_file(
          X = X,
          index = y,
          from = 0,
          to = X$end[y],
          path = path
        )
      
      if (X$sound.id[y] != "ambient")
        noise_clp <-
          warbleR::read_sound_file(
            X = X,
            index = y,
            from = 0,
            to = X$start[y] - 0.001,
            path = path
          )
      
      # add band-pass frequency filter
      if (!is.null(bp)) {
        # filter to bottom and top freq range
        if (bp == "freq.range")
          bp <- c(X$bottom.freq[y], X$top.freq[y])
        
        clp <-
          seewave::ffilter(
            clp,
            f = clp@samp.rate,
            from = bp[1] * 1000,
            ovlp = ovlp,
            to = bp[2] * 1000,
            bandpass = TRUE,
            wl = wl,
            output = "Wave"
          )
        
        if (X$sound.id[y] != "ambient")
          noise_clp <-
          seewave::ffilter(
            noise_clp,
            f = noise_clp@samp.rate,
            from = bp[1] * 1000,
            ovlp = ovlp,
            to = bp[2] * 1000,
            bandpass = TRUE,
            wl = wl,
            output = "Wave"
          )
      }
      
      # get mean envelopes
      sig_env <-
        mean(seewave::env(
          clp,
          f = clp@samp.rate,
          envt = "hil",
          plot = FALSE
        ))
      
      return(data.frame((X[y, , drop = FALSE]), sig_env))
    }
    
    # set clusters for windows OS
    if (Sys.info()[1] == "Windows" & cores > 1)
      cl <-
      parallel::makePSOCKcluster(getOption("cl.cores", cores)) else
      cl <- cores
    
    if (pb)
      write(file = "",
            x = paste0("Calculating amplitude envelopes (step 2 out of 3):"))
    
    # run loop apply function
    mean_envs <-
      warbleR:::pblapply_wrblr_int(
        X = 1:nrow(X),
        pbar = pb,
        cl = cl,
        FUN = function(y)
          meanenv_FUN(y, wl, ovlp)
      )
    
    # put in a data frame
    X2 <- do.call(rbind, mean_envs)
    
    # split by sound ID
    sigtype_list <- split(X2, X2$sound.id)
    
    if (pb)
      write(file = "",
            x = paste0("Calculating excess attenuation (step 3 out of 3):"))
    
    # calculate excess attenuation
    X_list <-
      warbleR:::pblapply_wrblr_int(X = sigtype_list, pbar = pb, cl = cl, function(Y, meth = method, tp = type) {
        if (Y$sound.id[1] == "ambient")
          Y$excess.attenuation <- NA else {
          # method 1 compare to closest distance to source
          if (meth == 1) {
            # extract mean envelope of sounds
            sig_env_REF <- Y$sig_env[which.min(Y$distance)]
            dist_REF <- Y$distance[which.min(Y$distance)]
            
            ks <- Y$sig_env / sig_env_REF
            
            # type Dabelsteen
            if (tp == "Dabelsteen")
              ea <-
              (-20 * log(ks)) - (6 / (2 * (Y$distance - dist_REF))) + gain
            
            if (tp == "Darden")
              #EA = g - 20 log(d / 10) - 20 log(k)
              ea <-
              gain - 20 * log10(Y$distance / 10) - 20 * log(ks)
            
            Y$excess.attenuation <- ea
            Y$excess.attenuation[which.min(Y$distance)] <- NA
          }
          
          # compare to previous distance
          if (meth == 2) {
            # save original order
            Y$org....ord <- 1:nrow(Y)
            
            # sort by distance
            Y <- Y[order(Y$distance), ]
            
            ks <- Y$sig_env[-1] / Y$sig_env[-nrow(Y)]
            
            if (tp == "Dabelsteen")
              ea <-
              (-20 * log(ks)) - (6 / (2 * (Y$distance[-1] - Y$distance[1]))) + gain
            
            
            # type Darden
            if (tp == "Darden")
              ea <-
              gain - 20 * log10(Y$distance[-1] / 10) - 20 * log(ks)
            
            # add NA for first distance
            ea <- c(NA, ea)
            
            Y$excess.attenuation <- ea
            
            # reorder results
            Y <- Y[order(Y$org....ord), ]
            
            Y$org....ord <- NULL
          }
        }
        
        Y <- as.data.frame(Y)
        return(Y)
        
      })
    
    # put together in a data frame as X
    X2 <- do.call(rbind, X_list)
    
    # fix row names
    rownames(X2) <-  rownames(X)
    
    # remove temporal column
    X2$sigRMS <- X2$TEMP....sgnl <- X2$sig_env <- NULL
    
    # fix est
    if (output == "est" & is_extended_selection_table(X)){
      X2 <- warbleR::fix_extended_selection_table(X = X2, Y = X)
      
      # fix call attribute
      attributes(X2)$call <- base::match.call() 
      }
    return(X2)
  }
