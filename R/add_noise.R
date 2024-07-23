#' Add synthetic noise
#'
#' \code{add_noise} adds synthetic noise to annotations in extended selection tables
#' @inheritParams template_params
#' @param X Object of class 'extended_selection_table' (created by the function \code{\link[warbleR]{selection_table}} from the warbleR package), generated 'by element', with the reference to the test sounds (typically the output of \code{\link{align_test_files}}). Must contain the following columns: 1) "sound.files": name of the .wav files, 2) "selec": unique selection identifier (within a sound file), 3) "start": start time and 4) "end": end time of selections, 5) "bottom.freq": low frequency for bandpass, 6) "top.freq": high frequency for bandpass and 7) "sound.id": ID of sounds (needed to calculate signal to noise ratio internally using \code{\link{signal_to_noise_ratio}}).
#' @param mar numeric vector of length 1. Specifies the margins adjacent to
#'   the start point of the annotation over which to measure ambient noise.
#' @param target.snr numeric vector of length 1. Specifies the desired signal-to-noise ratio. Must be lower that the current signal-to-noise ratio. Annotations showing a signal-to-noise ratio higher than 'target.snr' will remain unchanged. Must be supplied.
#' @param precision  numeric vector of length 1. Specifies the precision of the adjusted signal-to-noise ratio (in dB).
#' @param max.iterations Numeric vector of length 1. Specifies the maximum number of iterations that the internal signal-to-noise adjusting routine will run before stopping. Note that in most cases the default maximum number of iterations (1000) is not reached.
#' @param kind Character vector of length 1 indicating the kind of noise, “white”, “pink”, “power”, "brown", or “red”. Noise is synthesized with a modified version of the function \code{\link[tuneR]{noise}}. Default is "pink" which is similar to background noise in natural environments.
#' @param alpha Numeric vector of length 1. The power for the power law noise (defaults are 1 for pink and 1.5 for red noise). Only used when \code{kind = "power"}.
#' @param ... Additional arguments to be passed internally to \code{\link{signal_to_noise_ratio}}.
#' @return Object 'X' in which the wave objects have been modified to match the target signal-to-noise ratio. It also includes an additional column, 'adjusted.snr', with the new signal-to-noise ratio values.
#' @export
#' @name add_noise
#' @details The function adds synthetic noise to sounds referenced in an extended selection table (class created by the function \code{\link[warbleR]{selection_table}} from the warbleR package) to decrease the signal-to-noise ratio. This can be useful, for instance, for evaluating the effect of background noise on signal structure. Note that the implementation is slow.
#' @examples \dontrun{
#' # load example data
#' data("test_sounds_est")
#'
#' # make it a 'by element' extended selection table
#' X <- warbleR::by_element_est(X = test_sounds_est)
#'
#' # add noise to the first five rows
#' X_noise <- add_noise(X = X[1:5, ], mar = 0.2, target.snr = 3)
#' }
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})
#' @family miscellaneous
#' @seealso \code{\link{signal_to_noise_ratio}}
#' @references {
#' Araya-Salas M., E. Grabarczyk, M. Quiroz-Oliva, A. Garcia-Rodriguez, A. Rico-Guevara. (2023), baRulho: an R package to quantify degradation in animal acoustic signals .bioRxiv 2023.11.22.568305.
#' Timmer. J and M. König (1995): On generating power law noise. Astron. Astrophys. 300, 707-710.
#' }

add_noise <-
  function(X,
           mar = NULL,
           target.snr = 2,
           precision = 0.1,
           cores = getOption("mc.cores", 1),
           pb = getOption("pb", TRUE),
           max.iterations = 1000,
           kind = c("pink", "white", "brown", "red", "power"),
           alpha = 1,
           ...) {
    
    # assign a value to kind
    kind <- rlang::arg_match(kind)
    
    # check arguments
    arguments <- as.list(base::match.call())
    
    # add objects to argument names
    for (i in names(arguments)[-1]) {
      # use try to avoid errors with argumets from dots (...)
      try(arguments[[i]] <- get(i), silent = TRUE)
    }
    
    # check each arguments
    check_results <-
      .check_arguments(fun = arguments[[1]], args = arguments)
    
    # report errors
    .report_assertions(check_results)
    
    # get index number
    target_rows <- which(!X$sound.id %in% c("ambient", "marker"))
    
    wav_snr_list <-
      warbleR:::.pblapply(
        pbar = pb,
        X = target_rows,
        cl = cores,
        FUN = .add_noise,
        mar = mar,
        target.snr = target.snr,
        precision = precision,
        max.iterations = max.iterations,
        Y = X,
        kind = kind,
        alpha = alpha
      )
    
    names(wav_snr_list) <- target_rows
    
    # add empty column to fill with adjusted SNR values
    X$adjusted.snr <- NA
    
    # add modified waves to original extended selection table
    for (i in target_rows) {
      attributes(X)$wave.objects[[i]] <-
        wav_snr_list[[which(target_rows == i)]]$wave
      X$adjusted.snr[i] <-
        wav_snr_list[[which(target_rows == i)]]$snr
    }
    
    # check if some were not modified
    if (any(!sapply(wav_snr_list, "[[", 3))) {
      n_wavs <- sum(!sapply(wav_snr_list, "[[", 3))
      
      cli::cli_alert_warning(text = .colortext(
        paste0(
          "Warning: {n_wavs} wave{?s} already have a signal-to-noise ratio lower than ",
          target.snr,
          " (the target SNR) and w{?as/ere} not modified (run signal_to_noise_ratio(X, mar = ",
          mar,
          ") to spot {?it/them})"
        ),
        as = "magenta"
      ),
      wrap = TRUE)
    }
    return(X)
  }
