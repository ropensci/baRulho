#' Measure blur ratio in the time domain
#'
#' \code{plot_blur_ratio} measures blur ratio in sounds referenced in an extended selection table.
#' @usage plot_blur_ratio(X, type = "envelope", cores = getOption("mc.cores", 1),
#' pb = getOption("pb", TRUE), env.smooth = getOption("env.smooth", 200),
#' spec.smooth = getOption("spec.smooth", 5), res = 150, flim = c("-1", "+1"),
#' hop.size = getOption("hop.size", 11.6), wl = getOption("wl", NULL),
#' ovlp = getOption("ovlp", 70), palette = viridis::viridis,
#' collevels = seq(-120, 0, 5), dest.path = getOption("dest.path", "."),
#' path = getOption("sound.files.path", "."),
#' colors = viridis::viridis(3), n.samples = 100)
#' @param X The output of \code{\link{set_reference_sounds}} which is an object of class 'data.frame', 'selection_table' or 'extended_selection_table' (the last 2 classes are created by the function \code{\link[warbleR]{selection_table}} from the warbleR package) with the reference to the test sounds . Must contain the following columns: 1) "sound.files": name of the .wav files, 2) "selec": unique selection identifier (within a sound file), 3) "start": start time and 4) "end": end time of selections, 5)  "bottom.freq": low frequency for bandpass, 6) "top.freq": high frequency for bandpass, 7) "sound.id": ID of sounds used to identify counterparts across distances and 8) "reference": identity of sounds to be used as reference for each test sound (row). See \code{\link{set_reference_sounds}} for more details on the structure of 'X'.
#' @param type Character vector of length 1 indicating the type of blur ratio to plot. The two options are 'envelope' (for regular blur ratio as in \code{\link{blur_ratio}}, default) and 'spectrum' (for spectrum blur ratio as in \code{\link{spectrum_blur_ratio}}).
#' @param cores Numeric vector of length 1. Controls whether parallel computing is applied by specifying the number of cores to be used. Default is 1 (i.e. no parallel computing).
#' @param pb Logical argument to control if progress bar is shown. Default is \code{TRUE}.
#' @param env.smooth Numeric vector of length 1 determining the length of the sliding window (in amplitude samples) used for a sum smooth for amplitude envelope calculation (used internally by \code{\link[seewave]{env}}). Default is 200.
#' @param spec.smooth Numeric vector of length 1 determining the length of the sliding window used for a sum smooth for power spectrum calculation (in kHz). Default is 5.
#' @param res Numeric argument of length 1. Controls image resolution. Default is 150 (faster) although 300 - 400 is recommended for publication/presentation quality.
#' @param flim A numeric vector of length 2 indicating the highest and lowest frequency limits (kHz) of the spectrograms, as in \code{\link[seewave]{spectro}}. Default is \code{NULL}. Alternatively, a character vector similar to \code{c("-1", "1")} in which the first number is the value to be added to the minimum bottom frequency in 'X' and the second the value to be added to the maximum top frequency in 'X'. This is computed independently for each sound id so the frequency limit better fits the frequency range of the annotated signals. This is useful when test sounds show marked differences in their frequency ranges.
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
#' @param colors Character vector of length 4 containing the colors to be used for the background of column and row title panels (element 1), the color of amplitude envelopes (element 2), the color of power spectra (element 3), and the background color of envelopes and spectra (element 4).
#' @param n.samples Numeric vector of length 1 specifying the number of amplitude samples (or frequency bins if \code{spectrum = TRUE}) to use for representing power distributions. Default is 100. If null the raw power distribution is used (note that this can result in high RAM memory usage for large data sets).
#' @return It returns 1 image file (in 'jpeg' format) for each blur ratio estimation, showing spectrograms of both sounds and the overlaid amplitude envelopes (or power spectra if \code{spectrum = TRUE}) as probability mass functions (PMF). Spectrograms are shown within the frequency range of the reference sound.
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
#' Araya-Salas, M. (2020). baRulho: baRulho: quantifying degradation of (animal) acoustic signals in R. R package version 1.0.2
#' }

plot_blur_ratio <-
  function(X,
           type = "envelope",
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

    # adjust wl based on hop.size
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
    if (!(wl %% 2) == 0) {
      wl <- wl + 1
    }

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
            env_FUN(
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
            spctr_FUN(
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
    catch_out <-
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
          plot_blur_FUN(
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
  }
