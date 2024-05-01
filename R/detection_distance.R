#' Measure detection distance of sound
#'
#' \code{detection_distance} detection distance of sounds.
#' @inheritParams template_params
#' @param spl A numeric vector of length 1 specifying the sound pressure level of sounds. If not supplied then it will be measured from the sounds  themselves.
#' @param spl.cutoff A numeric vector of length 1 specifying the sound pressure level cutoff to define if the sound is no longer detected. Ideally it should be estimated based on the sound detection threshold of the species.
#' @param temp Numeric vector of length 1 with frequency (in Celsius). Default is 20.
#' @param rh Numeric vector of length 1 with relative humidity (in percentage). Default is 60.
#' @param pa Numeric vector of length 1 with ambient pressure in Pa (standard: 101325, default). Used for Atmospheric attenuation.
#' @param hab.att.coef Attenuation coefficient of the habitat (in dB/kHz/m).
#' @param max.distance Numeric vector of length 1 with the maximum distance (in m) at which detection would be evaluated. Note that the function calculates the expected sound pressure level values along a vector of distances to find the distance at which the expected sound pressure level equates 'spl.cutoff'. Default is 1000 m.
#' @param resolution Numeric vector of length 1 with the distance resolution (in m) for estimated detection distance. Higher resolutions take longer to estimate. Default is 0.1 m.
#' @param subtract.bgn Logical argument to control if SPL from background noise is excluded from the measured signal SPL. Default is \code{FALSE}.
#' @param envelope Character string vector with the method to calculate amplitude envelopes (in which SPL is measured, used required if 'spl' is not supplied), as in \code{\link[seewave]{env}}. Must be either 'abs' (absolute envelope, default) or 'hil' (Hilbert transformation).
#' @param mar numeric vector of length 1. Specifies the margins adjacent to
#'   the start and end points of selection over which to measure background noise. This is required to subtract background noise sound pressure level (so only needed when 'subtract.bgn = TRUE').
#' @return Object 'X' with an additional column, 'detection.distance',
#' containing the computed detection distances (in m).
#' @export
#' @name detection_distance
#' @details The function computes the maximum distance at which a sound would be detected, which is calculated as the distance in which the sound pressure level (SPL) goes below the specified SPL cutoff ('spl.cutoff')). This is returned as an additional column 'detection.distance' (in m). The function uses internally \code{\link{attenuation}} to estimate SPL at increasing values until it reaches the defined cutoff. The peak frequency (calculated on the power spectrum of the reference sound) of the reference sound for each sound ID is used as the carrier frequency for distance estimation. The sound recorded at the lowest distance is used as reference. \strong{This function assumes that all recordings have been made at the same recording volume}.
#' @examples \dontrun{
#' # load example data
#' data("test_sounds_est")
#'
#' # add reference to X
#' X <- set_reference_sounds(X = test_sounds_est)
#'
#' detection_distance(X = X[X$distance %in% c(1, 10), ], spl.cutoff = 5, mar = 0.05)
#' }
#'
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})
#' @family quantify degradation
#' @seealso \code{\link{attenuation}}
#' @references {
#' Araya-Salas M., E. Grabarczyk, M. Quiroz-Oliva, A. Garcia-Rodriguez, A. Rico-Guevara. (2023), baRulho: an R package to quantify degradation in animal acoustic signals .bioRxiv 2023.11.22.568305.
#'
#' Clark, C.W., Marler, P. & Beeman K. (1987). Quantitative analysis of animal vocal phonology: an application to Swamp Sparrow song. Ethology. 76:101-115.
#' }
detection_distance <-
  function(X,
           cores = getOption("mc.cores", 1),
           pb = getOption("pb", TRUE),
           hop.size = getOption("hop.size", 11.6),
           wl = getOption("wl", NULL),
           path = getOption("sound.files.path", "."),
           spl = NULL,
           spl.cutoff = NULL,
           temp = 20,
           rh = 60,
           pa = 101325,
           hab.att.coef = 0.02,
           max.distance = 1000,
           resolution = 0.1,
           subtract.bgn = TRUE,
           envelope = c("abs", "hil"),
           mar = NULL) {
    
    # assign a value to type
    envelope <- rlang::arg_match(envelope)
             
    # check arguments
    arguments <- as.list(base::match.call())
    
    # add objects to argument names
    for (i in names(arguments)[-1]) {
      arguments[[i]] <- get(i)
    }
    
    # check each arguments
    check_results <-
      .check_arguments(fun = arguments[[1]], args = arguments)
    
    # report errors
    .report_assertions(check_results)
    
    # error if no mar supplied when subtract.bgn
    if (is.null(mar) & subtract.bgn) {
      .stop("'mar' must be supplied when 'subtract.bgn = TRUE'")
    }
    
    # adjust wl based on hop.size
    wl <- .adjust_wl(wl, X, hop.size, path)
    
    # set clusters for windows OS
    if (Sys.info()[1] == "Windows" & cores > 1) {
      cl <-
        parallel::makePSOCKcluster(getOption("cl.cores", cores))
    } else {
      cl <- cores
    }
    
    # print message
    if (pb) {
      if (is.null(spl)) {
        write(file = "", x = "Computing sound pressure level and peak frequency (step 1 out of 2):")
      } else {
        write(file = "", x = "Computing peak frequency (step 1 out of 2):")
      }
    }
    
    # add sound file selec colums to X (weird column name so it does not overwrite user columns)
    X$.sgnl.temp <- paste(X$sound.files, X$selec, sep = "-")
    
    # get names of envelopes involved (those as test with reference or as reference)
    target_sgnl_temp <-
      unique(c(X$.sgnl.temp[!is.na(X$reference)], X$reference[!is.na(X$reference)]))
    
    
    # calculate all spectra apply function
    peak_freq_list <-
      warbleR:::pblapply_wrblr_int(
        pbar = pb,
        X = target_sgnl_temp,
        cl = cl,
        FUN = function(y, wl) {
          # load clip
          clp <- warbleR::read_sound_file(X = X,
                                          index = which(X$.sgnl.temp == y),
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
                  index = which(X$.sgnl.temp == y),
                  path = path,
                  from = if (X$start[X$.sgnl.temp == y] - mar < 0) {
                    0
                  } else {
                    X$start[X$.sgnl.temp == y] - mar
                  },
                  to = X$start[X$.sgnl.temp == y]
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
          } else {
            signaldb <- spl
          }
          # put output together
          output <- list(spl = signaldb, peakf = peak_freq[1, 1])
          
          return(output)
        }
      )
    
    # add sound file selec names to spectra
    names(peak_freq_list) <- target_sgnl_temp
    
    if (pb) {
      write(file = "", x = "Computing detection distance (step 2 out of 2):")
    }
    
    # get detection distance
    # calculate all spectra apply function
    X$detection.distance <-
      unlist(warbleR:::pblapply_wrblr_int(
        X = seq_len(nrow(X)),
        pbar = pb,
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
                       Y = X,
                       pfl = peak_freq_list,
                       ...) {
          .detection_dist(
            x,
            wl = wle,
            spl.cutoff = spl.cutoffe,
            temp = tempe,
            rh = rhe,
            pa = pae,
            hab.att.coef = hab.att.coefe,
            max.distance = max.distancee,
            resolution = resolution,
            X = Y,
            peak_freq_list = pfl,
            ...
          )
        }
      ))
    
    # make NAs those sounds in which the reference is itself (only happens when method = 2) or is ambient noise
    X$reference[X$reference == X$.sgnl.temp |
                  X$sound.id == "ambient"] <- NA
    
    # remove temporal columns
    X$.sgnl.temp <- NULL
    
    # fix call if not a data frame
    if (!is.data.frame(X)) {
      attributes(X)$call <-
        base::match.call()
    } # fix call attribute
    
    return(X)
  }
