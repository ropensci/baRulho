#' Measure spectrographic cross-correlation as a measure of sound distortion
#'
#' \code{detection_distance} measures spectrographic cross-correlation as a measure of sound distortion in sounds referenced in an extended selection table.
#' @usage detection_distance(X, cores = 1, pb = TRUE, output = "est",
#' hop.size = 11.6, wl = NULL, path = NULL, spl = NULL, spl.cutoff, temp = 20,
#' rh = 60, pa = 101325, hab.att.coef = 0.02, max.distance = 1000,
#'  resolution = 0.1, subtract.bgn = TRUE, envelope = "abs", mar = NULL)
#' @param X Object of class 'data.frame', 'selection_table' or 'extended_selection_table' (the last 2 classes are created by the function \code{\link[warbleR]{selection_table}} from the warbleR package) with the reference to the sounds in the master sound file. Must contain the following columns: 1) "sound.files": name of the .wav files, 2) "selec": unique selection identifier (within a sound file), 3) "start": start time and 4) "end": end time of selections, 5)  "bottom.freq": low frequency for bandpass, 6) "top.freq": high frequency for bandpass and 7) "sound.id": ID of sounds used to identify counterparts across distances. Each sound must have a unique ID within a distance.
#' @param cores Numeric vector of length 1. Controls whether parallel computing is applied by specifying the number of cores to be used. Default is 1 (i.e. no parallel computing).
#' @param pb Logical argument to control if progress bar is shown. Default is \code{TRUE}.
#' @param output Character vector of length 1 to determine if an extended selection table ('est', default) or a data frame ('data.frame').
#' @param hop.size A numeric vector of length 1 specifying the time window duration (in ms). Default is 11.6 ms, which is equivalent to 512 wl for a 44.1 kHz sampling rate. Ignored if 'wl' is supplied.
#' @param wl A numeric vector of length 1 specifying the window length of the spectrogram, default
#' is \code{NULL}. If supplied, 'hop.size' is ignored.
#' @param path Character string containing the directory path where the sound files are found. Only needed when 'X' is not an extended selection table.
#' @param spl A numeric vector of length 1 specifying the sound pressure level of sounds. If not supplied then it will be measured from the sounds  themselves.
#' @param spl.cutoff A numeric vector of length 1 specifying the sound pressure level cutoff to define if the sound is no longer detected. Ideally it should be estimated based on the sound detection threshold of the species.
#' @param temp Numeric vector of length 1 with frequency (in Celsius). Default is 20.
#' @param rh Numeric vector of length 1 with relative humidity (in percentage). Default is 60.
#' @param pa Numeric vector of length 1 with ambient pressure in Pa (standard: 101325, default). Used for Atmospheric attenuation.
#' @param hab.att.coef Attenuation coefficient of the habitat (in dB/kHz/m).
#' @return Returns the geometric, atmospheric and habitat attenuation (in dB) as well as the combined attenuation.
#' @param max.distance Numeric vector of length 1 with the maximum distance (in m) at which detection would be evaluated. Note that the function calculates the expected sound pressure level values along a vector of distances to find the distance at which the expected sound pressure level equates 'spl.cutoff'. Default is 1000 m.
#' @param resolution NUmeric vector of length 1 with the distance resolution (in m) for estimated detection distance. Higher resolutions take longer to estimate. Default is 0.1 m.
#' @param subtract.bgn Logical argument to control if SPL from background noise is excluded from the measured signal SPL. Default is \code{FALSE}.
#' @param envelope Character string vector with the method to calculate amplitude envelopes (in which SPL is measured, used required if 'spl' is not supplied), as in \code{\link[seewave]{env}}. Must be either 'abs' (absolute envelope, default) or 'hil' (Hilbert transformation).
#' @param mar numeric vector of length 1. Specifies the margins adjacent to
#'   the start and end points of selection over which to measure background noise. This is required to subtract background noise sound pressure level (so only needed when 'subtract.bgn = TRUE').
#' @return A data frame, or extended selection table similar to input data (depending on argument 'output'), but including two new columns (reference and detection.distance).
#' @export
#' @name detection_distance
#' @details The function computes the maximum distance at which a sound would be detected, which is calculated as the distance in which the sound pressure level (SPL) goes below the specified SPL cutoff ('spl.cutoff')). The function uses internally \code{\link{attenuation}} to estimate SPL at increasing values until it reaches the defined cutoff. The peak frequency (calculated on the power spectrum of the reference sound) of the reference sound for each sound ID is used as the carrier frequency for distance estimation. The sound recorded at the lowest distance is used as reference. \strong{This function assumes that all recordings have been made at the same recording volume}.
#' @examples
#' {
#' # load example data
#' data("playback_est")
#'
#'detection_distance(X = playback_est[playback_est$distance == 1,], spl.cutoff = 5, mar = 0.05)
#' }
#'
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})
#' @seealso \code{\link{attenuation}}
#' @references {
#' Araya-Salas, M. (2020). baRulho: baRulho: quantifying habitat-induced degradation of (animal) acoustic signals in R. R package version 1.0.2
#'
#' Clark, C.W., Marler, P. & Beeman K. (1987). Quantitative analysis of animal vocal phonology: an application to Swamp Sparrow song. Ethology. 76:101-115.
#' }
detection_distance <-
  function(X,
           cores = 1,
           pb = TRUE,
           output = "est",
           hop.size = 11.6,
           wl = NULL,
           path = NULL,
           spl = NULL,
           spl.cutoff = NULL,
           temp = 20,
           rh = 60,
           pa = 101325,
           hab.att.coef = 0.02,
           max.distance = 1000,
           resolution = 0.1,
           subtract.bgn = TRUE,
           envelope = "abs",
           mar = NULL) {
    
    # error if no spl.cutoff supplied
    if (is.null(spl.cutoff))
      stop2("'spl.cutoff' must be supplied")
    
    # error if no mar supplied when subtract.bgn
    if (is.null(mar) & subtract.bgn)
      stop2("'mar' must be supplied when 'subtract.bgn = TRUE'")
    
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
    
    # need mar if subtract.bgn TRUE
    if (subtract.bgn & is.null(mar))
      stop2("'mar' must be supplied if 'subtract.bgn = TRUE'")
    
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
                      method = 1,
                      cores = cores,
                      pb = pb)
    
    # set clusters for windows OS
    if (Sys.info()[1] == "Windows" & cores > 1)
      cl <-
      parallel::makePSOCKcluster(getOption("cl.cores", cores)) else
      cl <- cores
    
    # print message
    if (pb)
      if (is.null(spl))
        write(file = "", x = "Computing sound pressure level and peak frequency (step 2 out of 3):") else
      write(file = "", x = "Computing peak frequency (step 2 out of 3):")
    
    # calculate all spectra apply function
    peak_freq_list <-
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
          
          # get peak frequency
          peak_freq <-
            seewave::fpeaks(clp.spc, nmax = 1, plot = FALSE)
          
          if (is.null(spl)) {
            # get amplitude
            sigamp <-
              seewave::rms(seewave::env(clp, envt = envelope, plot = FALSE))
            
            # convert to dB
            signaldb <- 20 * log10(sigamp)
            
            # remove background SPL
            if (subtract.bgn) {
              bg.noise <-
                read_sound_file(
                  X,
                  index = y,
                  path = path,
                  from = if (X$start[y] - mar < 0)
                    0 else
                    X$start[y] - mar,
                  to = X$start[y]
                )
              
              noiseamp <-
                seewave::rms(seewave::env(
                  bg.noise,
                  f = bg.noise@samp.rate,
                  envt = envelope,
                  plot = FALSE
                ))
              noisedb <- 20 * log10(noiseamp)
              
              # remove noise SPL from signal SPL
              signaldb <-
                warbleR:::lessdB(signal.noise = signaldb, noise = noisedb)
            }
          } else
            signaldb <- spl
          # put output together
          output <- list(spl = signaldb, peakf = peak_freq[1, 1])
          
          return(output)
        }
      )
    
    # add sound file selec names to spectra
    names(peak_freq_list) <- X$TEMP....sgnl
    
    ## function to measure detection distance
    detection_dist_FUN <-
      function(x,
               spl.cutoff,
               temp,
               rh,
               pa,
               hab.att.coef,
               max.distance,
               resolution,
               spl,
               ...) {
        # get names of sound and reference
        sgnl <- X$TEMP....sgnl[x]
        rfrnc <- X$reference[x]
        
        # if sounds are the same or the selection is noise return NA
        if (any(c(X$sound.id[X$TEMP....sgnl == sgnl], X$sound.id[X$reference == rfrnc]) == "ambient"))
          detect_dist <- NA else {
          # extract spectrum for sound and model
          sgnl.spl <-
            peak_freq_list[[which(names(peak_freq_list) == sgnl)]]$spl
          rfrnc.pkf <-
            peak_freq_list[[which(names(peak_freq_list) == rfrnc)]]$peakf
          
          # get detection distance
          detect_dist <-
            detection_distance_ohn_int(
              spl.cutoff = spl.cutoff,
              spl = sgnl.spl,
              frequency = rfrnc.pkf * 1000,
              distance = X$distance[x],
              temp = temp,
              rh = rh,
              pa = pa,
              hab.att.coef = hab.att.coef,
              max.distance = max.distance,
              resolution = resolution
            )
          
        }
        return(detect_dist)
      }
    
    if (pb)
      write(file = "", x = "Computing detection distance (step 3 out of 3):")
    
    # get blur ratio
    # calculate all spectra apply function
    X$detection.distance <-
      pbapply::pbsapply(
        X = 1:nrow(X),
        cl = cl,
        FUN = function(x,
                       wle = wl,
                       spl.cutoffe = spl.cutoff,
                       tempe = temp,
                       rhe = rh,
                       pae = pa,
                       hab.att.coefe = hab.att.coef,
                       max.distancee = max.distance,
                       resolutione = resolution,
                       ...)   {
          detection_dist_FUN(
            x,
            wl = wle,
            spl.cutoff = spl.cutoffe,
            temp = tempe,
            rh = rhe,
            pa = pae,
            hab.att.coef = hab.att.coefe,
            max.distance = max.distancee,
            resolution = resolution,
            ...
          )
        }
      )
    
    # remove temporal columns
    X$TEMP....sgnl <- NULL
    
    # return data frame
    if (output == "data.frame")
      X <- as.data.frame(X) else
      attributes(X)$call <-
      base::match.call() # fix call attribute
    
    return(X)
  }

##############################################################################################################
# internal function to get detection distance from spl and cutoff
detection_distance_ohn_int <-
  function(spl.cutoff,
           spl,
           frequency,
           distance,
           temp = 20,
           rh = 60,
           pa = 101325,
           hab.att.coef = 0.02,
           max.distance = 1000,
           resolution = 0.1) {
    # initial SPL and distance
    L <- spl
    iter_dist <- 0
    
    # loop until SPL is equal or lower than background noise
    while (L - spl.cutoff > 0) {
      iter_dist <- iter_dist + resolution
      att <-
        attenuation(
          frequency = frequency,
          temp = temp,
          dist = iter_dist,
          dist0 = distance,
          rh = rh,
          pa = pa,
          hab.att.coef = hab.att.coef
        )
      L <- spl - att$combined.attenuation
      
      if (iter_dist >= max.distance) {
        iter_dist <- NA
        break
      }
    }
    
    return(iter_dist)
  }
