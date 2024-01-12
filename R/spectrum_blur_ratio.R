#' Measure blur ratio in the frequency domain
#'
#' \code{spectrum_blur_ratio} measures blur ratio of frequency spectra from sounds referenced in an extended selection table.
#' @inheritParams template_params
#' @param spec.smooth Numeric vector of length 1 determining the length of the sliding window used for a sum smooth for power spectrum calculation (in kHz). Default is 5.
#' @param spectra Logical to control if power spectra are returned (as attributes). Default is \code{FALSE}.
#' @param res Numeric argument of length 1. Controls image resolution. Default is 150 (faster) although 300 - 400 is recommended for publication/presentation quality.
#' @param wl A numeric vector of length 1 specifying the window length of the spectrogram, default
#' is NULL. If supplied, 'hop.size' is ignored. Applied to both spectra and spectrograms on image files.
#' @param ovlp Numeric vector of length 1 specifying the percentage of overlap between two
#'   consecutive windows, as in \code{\link[seewave]{spectro}}. Default is 70. Applied to both spectra and spectrograms on image files.
#' @return Object 'X' with an additional column, 'spectrum.blur.ratio', containing the computed spectrum blur ratio values. If \code{spectra = TRUE} the output would include power spectra for all sounds as attributes ('attributes(X)$spectra').
#' @param n.bins Numeric vector of length 1 specifying the number of frequency bins to use for representing power spectra. Default is 100. If null the raw power spectrum is used (note that this can result in high RAM memory usage for large data sets). Power spectrum values are interpolated using \code{\link[stats]{approx}}.
#' @export
#' @name spectrum_blur_ratio
#' @details Spectral blur ratio measures the degradation of sound as a function of the change in sound power in the frequency domain, analogous to the blur ratio proposed by Dabelsteen et al (1993) for the time domain (and implemented in \code{\link{blur_ratio}}). Low values indicate low degradation of sounds. The function measures the blur ratio of spectra from sounds in which a reference playback has been re-recorded at different distances. Spectral blur ratio is measured as the mismatch between power spectra (expressed as probability density functions) of the reference sound and the re-recorded sound. The function compares each sound type to the corresponding reference sound. The 'sound.id' column must be used to tell the function to only compare sounds belonging to the same category (e.g. song-types). Two methods for setting the experimental design are provided. All wave objects in the extended selection table must have the same sampling rate so the length of spectra is comparable. The function uses \code{\link[seewave]{spec}} internally to compute power spectra. NA is returned if at least one the power spectra cannot be computed.
#' @examples {
#'   # load example data
#'   data("test_sounds_est")
#'
#' # add reference to X
#' X <- set_reference_sounds(X = test_sounds_est)
#'
#'   # get spetrum blur ratio
#'   spectrum_blur_ratio(X = X)
#'
#'   # using method 2
#'    X <- set_reference_sounds(X = test_sounds_est, method = 2)
#'   spectrum_blur_ratio(X = X)
#'
#'   # get power spectra
#'   sbr <- spectrum_blur_ratio(X = X, spectra = TRUE)
#'   spctr <- attributes(sbr)$spectra
#'
#'   # make distance a factor for plotting
#'   spctr$distance <- as.factor(spctr$distance)
#'
#'   
#'   # plot
#'   rlang::check_installed("ggplot2")
#'   library(ggplot2)
#'   
#'   ggplot(spctr[spctr$freq > 0.3, ], aes(y = amp, x = freq,
#'   col = distance)) +
#'   geom_line() +
#'   facet_wrap(~sound.id) +
#'   scale_color_viridis_d(alpha = 0.7) +
#'   labs(x = "Frequency (kHz)", y = "Amplitude (PMF)") +
#'   coord_flip() +
#'   theme_classic()
#' }
#'
#' @seealso \code{\link{blur_ratio}}
#' @family quantify degradation
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})
#' @references {
#' Dabelsteen, T., Larsen, O. N., & Pedersen, S. B. (1993). Habitat-induced degradation of sound signals: Quantifying the effects of communication sounds and bird location on blur ratio, excess attenuation, and signal-to-noise ratio in blackbird song. The Journal of the Acoustical Society of America, 93(4), 2206.
#'
#' Araya-Salas, M. (2020). baRulho: baRulho: quantifying degradation of (animal) acoustic signals in R. R package version 1.0.2
#' }

spectrum_blur_ratio <-
  function(X,
           cores = getOption("mc.cores", 1),
           pb = getOption("pb", TRUE),
           spec.smooth = getOption("spec.smooth", 5),
           spectra = FALSE,
           res = 150,
           hop.size = getOption("hop.size", 11.6),
           wl = getOption("wl", NULL),
           ovlp = getOption("ovlp", 70),
           path = getOption("sound.files.path", "."),
           n.bins = 100) {
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
    
    # total number of steps depending on whether envelopes are returned
    steps <- if (spectra)
      3
    else
      2
    
    # get sampling rate
    sampling_rate <-
      warbleR::read_sound_file(
        X = X,
        index = 1,
        path = path,
        header = TRUE
      )$sample.rate
    
    
    # adjust wl based on hop.size
    wl <- .adjust_wl(wl, X, hop.size, path)
    
    # add sound file selec colums to X (weird column name so it does not overwrite user columns)
    X$.sgnl.temp <- paste(X$sound.files, X$selec, sep = "-")
    
    # get names of envelopes involved (those as test with reference or as reference)
    target_sgnl_temp <-
      unique(c(X$.sgnl.temp[!is.na(X$reference)], X$reference[!is.na(X$reference)]))
    
    # set clusters for windows OS
    if (Sys.info()[1] == "Windows" & cores > 1) {
      cl <-
        parallel::makePSOCKcluster(getOption("cl.cores", cores))
    } else {
      cl <- cores
    }
    
    # print message
    if (pb)
      write(file = "",
            x = paste0("Computing power spectra (step 1 out of ", steps, "):"))
    
    # calculate all spectra apply function
    specs <-
      warbleR:::pblapply_wrblr_int(
        pbar = pb,
        X = target_sgnl_temp,
        cl = cl,
        FUN = function(y,
                       ssmth = spec.smooth,
                       wln = wl,
                       Q = X,
                       pth = path,
                       nb = n.bins) {
          .spctr(
            y,
            spec.smooth = ssmth,
            wl = wln,
            X = Q,
            path = pth,
            n.bins = nb
          )
        }
      )
    
    # add sound file selec names to spectra
    names(specs) <- target_sgnl_temp
    
    # print second message
    if (pb)
      write(
        file = "",
        x = paste0("Computing spectrum blur ratio (step 2 out of ", steps, "):")
      )
    
    # get blur ratio
    # calculate all spectra apply function
    X$spectrum.blur.ratio <-
      unlist(warbleR:::pblapply_wrblr_int(
        X = seq_len(nrow(X)),
        pbar = pb,
        cl = cl,
        FUN = function(x,
                       Q = X,
                       wle = wl,
                       ovp = ovlp,
                       spcs = specs,
                       sr = sampling_rate) {
          .blur_sp(
            x,
            X = Q,
            ovlp = ovp,
            wl = wle,
            specs = spcs,
            sampling_rate = sr
          )
        }
      ))
    
    # remove temporal columns
    X$.sgnl.temp <- NULL
    
    if (pb & spectra)
      write(file = "", x = "Saving spectra (step 3 out of 3):")
    
    
    # convert to list instead of extended selection table, add envelopes
    if (spectra) {
      spec.dfs <-
        warbleR:::pblapply_wrblr_int(
          X = seq_along(specs),
          cl = cl,
          pbar = pb,
          FUN = function(y) {
            # extract 1 envelope
            x <- specs[[y]]
            
            # convert envelopes to PMF (probability mass function)
            x[, 2] <- x[, 2] / sum(x[, 2])
            
            # put in data framme
            out <-
              data.frame(
                sound = names(specs)[y],
                sound.id = X$sound.id[paste(X$sound.files, X$selec, sep = "-") == names(specs)[y]],
                distance = X$distance[paste(X$sound.files, X$selec, sep = "-") == names(specs)[y]],
                freq = x[, 1],
                amp = x[, 2]
              )
            
            return(out)
          }
        )
      
      # put together in a single data frame
      spec.df <- do.call(rbind, spec.dfs)
      
      # add envelopes as attributes
      attributes(X)$spectra <- spec.df
    }
    
    # return data frame
    if (warbleR::is_extended_selection_table(X) | is_selection_table(X)) {
      attributes(X)$call <- base::match.call()
    }
    return(X)
  }
