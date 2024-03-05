#' Measure excess attenuation
#'
#' \code{excess_attenuation} measures excess attenuation in sounds referenced in an extended selection table.
#' @inheritParams template_params
#' @param hop.size A numeric vector of length 1 specifying the time window duration (in ms). Default is 1 ms, which is equivalent to ~45 wl for a 44.1 kHz sampling rate. Ignored if 'wl' is supplied.
#' @param wl A numeric vector of length 1 specifying the window length of the spectrogram, default
#' is \code{NULL}. If supplied, 'hop.size' is ignored.
#' Note that lower values will increase time resolution, which is more important for amplitude calculations.
#' @param ovlp Numeric vector of length 1 specifying the percentage of overlap between two
#'   consecutive windows, as in \code{\link[seewave]{spectro}}. Default is 50. Only used for bandpass filtering.
#' @return Object 'X' with an additional column,  'excess.attenuation', containing the computed excess attenuation values (in dB).
#' @export
#' @name excess_attenuation
#' @details Excess attenuation is the amplitude loss of a sound in excess due to spherical spreading (observed attenuation - expected attenuation). With every doubling of distance, sounds attenuate with a 6 dB loss of amplitude (Morton, 1975; Marten & Marler, 1977). Any additional loss of amplitude results in energy loss in excess of that expected to occur with distance via spherical spreading. So it represents power loss due to additional factors like vegetation or atmospheric conditions (Wiley & Richards, 1978). Excess attenuation is computed as \code{-20 * log10(rms("test signal") / rms("reference signal"))) - (20 * log10(1 / "distance")} in which 'rms(..)' represents the root mean square of an amplitude envelope. Low values indicate little additional attenuation. The goal of the function is to measure the excess attenuation on sounds in which a reference playback has been re-recorded at increasing distances. The 'sound.id' column must be used to indicate which sounds belonging to the same category (e.g. song-types). The function will then compare each sound type to the corresponding reference sound. Two approaches for computing excess attenuation are provided (see 'type' argument). NAs will be returned if one of the envelopes is completely flat (e.g. no variation in amplitude).
#' @examples {
#'   # load example data
#'   data("test_sounds_est")
#'
#'   # using method 1
#'   # add reference to X
#'   X <- set_reference_sounds(X = test_sounds_est)
#'   excess_attenuation(X = X)
#'
#'   # using method 2
#'   X <- set_reference_sounds(X = test_sounds_est, method = 2)
#'   # excess_attenuation(X = X)
#' }
#'
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})
#' @seealso \code{\link{spcc}}; \code{\link{envelope_correlation}}
#' @references {
#' Araya-Salas, M. (2020). baRulho: baRulho: quantifying degradation of (animal) acoustic signals in R. R package version 1.0.2
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

excess_attenuation <-
  function(X,
           cores = getOption("mc.cores", 1),
           pb = getOption("pb", TRUE),
           hop.size = getOption("hop.size", 1),
           wl = getOption("wl", NULL),
           ovlp = getOption("ovlp", 50),
           bp = "freq.range",
           path = getOption("sound.files.path", ".")) {
    
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
    
    # adjust wl based on hop.size
    wl <- .adjust_wl(wl, X, hop.size, path)
    
    # set clusters for windows OS
    if (Sys.info()[1] == "Windows" & cores > 1) {
      cl <-
        parallel::makePSOCKcluster(getOption("cl.cores", cores))
    } else {
      cl <- cores
    }
    
    if (pb) {
      write(file = "",
            x = paste0("Computing amplitude envelopes (step 1 out of 2):"))
    }
    
    # add sound file selec colums to X (weird column name so it does not overwrite user columns)
    X$.sgnl.temp <- paste(X$sound.files, X$selec, sep = "-")
    
    # get names of envelopes involved (those as test with reference or as reference)
    target_sgnl_temp <-
      unique(c(X$.sgnl.temp[!is.na(X$reference)], X$reference[!is.na(X$reference)]))
    
    # run loop apply function
    mean_envs <-
      warbleR:::pblapply_wrblr_int(
        X = target_sgnl_temp,
        pbar = pb,
        cl = cl,
        FUN = function(y,
                       wln = wl,
                       ovl = ovlp,
                       Q = X,
                       pth = path,
                       bps = bp) {
          .mean.env(
            y,
            wl = wln,
            ovlp = ovl,
            X = Q,
            path = pth,
            bp = bps, 
            rms = TRUE
          )
        }
      )
    
    # add sound file selec column as names to envelopes
    names(mean_envs) <- target_sgnl_temp
    
    # put in a data frame
    X$sig_env <- vapply(seq_len(nrow(X)), function(x) {
      w <- if (any(names(mean_envs) == X$.sgnl.temp[x])) {
        mean_envs[[which(names(mean_envs) == X$.sgnl.temp[x])]]
      } else {
        NA
      }
      return(w)
    }, FUN.VALUE = numeric(1))
    
    # split by sound ID
    # sigtype_list <- split(X, X$sound.id)
    
    if (pb) {
      write(file = "",
            x = paste0("Computing excess attenuation (step 2 out of 2):"))
    }
    
    # calculate excess attenuation
    X$excess.attenuation <-
      unlist(warbleR:::pblapply_wrblr_int(
        X = seq_len(nrow(X)),
        pbar = pb,
        cl = cl,
        FUN = function(x) {
          .exc_att(y = x, X)
        }
      ))
    
    
    # remove temporal column
    X$.sgnl.temp <- X$sig_env <- NULL
    
    # fix call if not a data frame
    if (!is.data.frame(X)) {
      attributes(X)$call <-
        base::match.call()
    } # fix call attribute
    
    
    return(X)
  }
