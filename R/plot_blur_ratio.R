#' Plot blur ratio 
#'
#' \code{plot_blur_ratio} plots time and frequency blur ratio in sounds referenced in an extended selection table.
#' @inheritParams template_params
#' @param type Character vector of length 1 indicating the type of blur ratio to plot. The two options are 'envelope' (for regular blur ratio as in \code{\link{blur_ratio}}, default) and 'spectrum' (for spectrum blur ratio as in \code{\link{spectrum_blur_ratio}}).
#' @param env.smooth Numeric vector of length 1 determining the length of the sliding window (in amplitude samples) used for a sum smooth for amplitude envelope calculation (used internally by \code{\link[seewave]{env}}). Default is 200.
#' @param spec.smooth Numeric vector of length 1 determining the length of the sliding window used for a sum smooth for power spectrum calculation (in kHz). Default is 5.
#' @param res Numeric argument of length 1. Controls image resolution. Default is 150 (faster) although 300 - 400 is recommended for publication/presentation quality.
#' @param flim A numeric vector of length 2 indicating the highest and lowest frequency limits (kHz) of the spectrograms, as in \code{\link[seewave]{spectro}}. Default is \code{NULL}. Alternatively, a character vector similar to \code{c("-1", "1")} in which the first number is the value to be added to the minimum bottom frequency in 'X' and the second the value to be added to the maximum top frequency in 'X'. This is computed independently for each sound id so the frequency limit better fits the frequency range of the annotated signals. This is useful when test sounds show marked differences in their frequency ranges.
#' @param ovlp Numeric vector of length 1 specifying the percentage of overlap between two
#'   consecutive windows, as in \code{\link[seewave]{spectro}}. Only used when plotting. Default is 70. Applied to both spectra and spectrograms on image files.
#' @param palette A color palette function to be used to assign colors in the
#'   plot, as in \code{\link[seewave]{spectro}}. Default is \code{\link[viridis]{viridis}}.
#' @param collevels	Numeric vector indicating a set of levels which are used to partition the amplitude range of the spectrogram (in dB) as in \code{\link[seewave]{spectro}}. Default is \code{seq(-120, 0, 5)}.
#' @param colors Character vector of length 4 containing the colors to be used for the color to identify the reference sound  (element 1),  the color to identify the test sound (element 2) and the color of blurred region (element 3).
#' @param n.samples Numeric vector of length 1 specifying the number of amplitude samples (or frequency bins if \code{spectrum = TRUE}) to use for representing power distributions. Default is 100. If null the raw power distribution is used (note that this can result in high RAM memory usage for large data sets).
#' @return It returns 1 image file (in 'jpeg' format) for each blur ratio estimation, showing spectrograms of both sounds and the overlaid amplitude envelopes (or power spectra if \code{spectrum = TRUE}) as probability mass functions (PMF). Spectrograms are shown within the frequency range of the reference sound. It also returns the file path of the images invisibly.
#' @export
#' @name plot_blur_ratio
#' @details The function generates image files (in 'jpeg' format) for each possible blur ratio estimation in 'X'. The image files show the spectrograms of both sounds and the overlaid power distribution (either amplitude envelopes or power spectrum, see argument 'type') as probability mass functions (PMF). The output graphs highlight the mismatch between the compared distribution which represent the estimated blur ratio returned by either \code{\link{blur_ratio}} or \code{\link{spectrum_blur_ratio}}. Spectrograms are shown within the frequency range of the reference sound and also show dotted lines with the time (type = "envelope") or frequency range (type = "spectrum") in which energy distributions where computed.
#' @family quantify degradation
#' @seealso \code{\link{envelope_correlation}}, \code{\link{spectrum_blur_ratio}}, \code{\link{blur_ratio}}
#' @examples {
#'   # load example data
#'   data("test_sounds_est")
#'
#'   # add reference to X
#'   X <- set_reference_sounds(X = test_sounds_est)
#'
#'   # create plots
#'   plot_blur_ratio(X = X, dest.path = tempdir())
#' }
#'
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})
#' @references {
#' Dabelsteen, T., Larsen, O. N., & Pedersen, S. B. (1993). Habitat-induced degradation of sound signals: Quantifying the effects of communication sounds and bird location on blur ratio, excess attenuation, and signal-to-noise ratio in blackbird song. The Journal of the Acoustical Society of America, 93(4), 2206.
#'
#' Araya-Salas M., E. Grabarczyk, M. Quiroz-Oliva, A. Garcia-Rodriguez, A. Rico-Guevara. (2023), baRulho: an R package to quantify degradation in animal acoustic signals .bioRxiv 2023.11.22.568305.
#' }

plot_blur_ratio <-
  function(X,
           type = c("envelope", "spectrum"),
           cores = getOption("mc.cores", 1),
           pb = getOption("pb", TRUE),
           env.smooth = getOption("env.smooth", 200),
           spec.smooth = getOption("spec.smooth", 5),
           res = 150,
           flim = c("-1", "+1"),
           hop.size = getOption("hop.size", 11.6),
           wl = getOption("wl", NULL),
           ovlp = getOption("ovlp", 70),
           palette = viridis::viridis,
           collevels = seq(-120, 0, 5),
           dest.path = getOption("dest.path", "."),
           path = getOption("sound.files.path", "."),
           colors = viridis::viridis(3),
           n.samples = 100) {
    
    # assign a value to type
    type <- rlang::arg_match(type)
    
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
    wl <- .adjust_wl(wl, X, hop.size)
    
    # set clusters for windows OS
    if (Sys.info()[1] == "Windows" & cores > 1) {
      cl <- parallel::makePSOCKcluster(getOption("cl.cores", cores))
    } else {
      cl <- cores
    }
    
    # add sound file selec colums to X (weird column name so it does not overwrite user columns)
    X$.sgnl.temp <- paste(X$sound.files, X$selec, sep = "-")
    
    # get names of envelopes involved (those as test with reference or as reference)
    target_sgnl_temp <-
      unique(c(X$.sgnl.temp[!is.na(X$reference)], X$reference[!is.na(X$reference)]))
    
    # print message
    if (pb) {
      if (type == "envelope") {
        write(file = "", x = "Computing amplitude envelopes (step 1 out of 2):")
      } else {
        write(file = "", x = "Computing power spectra (step 1 out of 2):")
      }
    }
    
    # calculate all envelops apply function
    if (type == "envelope") {
      energy_vectors <-
        warbleR:::pblapply_wrblr_int(
          pbar = pb,
          X = target_sgnl_temp,
          cl = cl,
          FUN = function(x,
                         ssmth = env.smooth,
                         ovl = ovlp,
                         Q = X,
                         wln = wl,
                         pth = path,
                         n.samp = n.samples) {
            .env(
              X = Q,
              y = x,
              env.smooth = ssmth,
              ovlp = ovl,
              wl = wln,
              path = pth,
              n.samples = n.samp
            )
          }
        )
    }
    
    if (type == "spectrum") {
      # calculate all spectra apply function
      energy_vectors <-
        warbleR:::pblapply_wrblr_int(
          pbar = pb,
          X = target_sgnl_temp,
          cl = cl,
          FUN = function(y,
                         ssmth = spec.smooth,
                         wln = wl,
                         Q = X,
                         pth = path,
                         n.b = n.samples) {
            .spctr(
              y,
              spec.smooth = ssmth,
              wl = wln,
              X = Q,
              path = pth,
              n.bins = n.b
            )
          }
        )
    }
    
    
    
    # add sound file selec column as names to envelopes
    names(energy_vectors) <- target_sgnl_temp
    
    # set options to run loops
    if (pb) {
      write(file = "", x = "Producing images (step 2 out of 2):")
    }
    
    # plot blur ratio
    file_paths <-
      warbleR:::pblapply_wrblr_int(
        pbar = pb,
        X = which(!is.na(X$.sgnl.temp)),
        cl = cl,
        FUN = function(x,
                       rs = res,
                       en.vctr = energy_vectors,
                       spct = if (type == "envelope") {
                         FALSE
                       } else {
                         TRUE
                       },
                       wle = wl,
                       colvs = collevels,
                       pl = palette,
                       ovp = ovlp,
                       dest.pth = dest.path,
                       fl = flim,
                       cols = colors,
                       Q = X) {
          .plot_blur(
            X = Q,
            energy_vectors = en.vctr,
            spectr = spct,
            path,
            dest.path = dest.pth,
            x,
            res = rs,
            ovlp = ovp,
            wl = wle,
            collevels = colvs,
            palette = pl,
            flim = fl,
            colors = cols
          )
        }
      )
    # message to let know users where the files have been saved
    .message(
      paste0(
        "The image files have been saved in the directory path '",
        normalizePath(dest.path),
        "'"
      )
    )
    
    # return file names without printing them
    invisible(unlist(file_paths))  
  }
