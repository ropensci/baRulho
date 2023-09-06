#' Measure frequency spectrum correlation
#'
#' \code{spectrum_correlation} measures frequency spectrum correlation of sounds referenced in an extended selection table.
#' @usage spectrum_correlation(X, parallel = NULL, cores = getOption("mc.cores", 1),
#' pb = getOption("pb", TRUE), cor.method = "pearson",
#' spec.smooth = getOption("spec.smooth", 5), output = NULL,
#' hop.size = getOption("hop.size", 11.6), wl = getOption("wl", NULL),
#' ovlp = getOption("ovlp", 70), path = getOption("sound.files.path", "."))
#' @param X The output of \code{\link{set_reference_sounds}} which is an object of class 'data.frame', 'selection_table' or 'extended_selection_table' (the last 2 classes are created by the function \code{\link[warbleR]{selection_table}} from the warbleR package) with the reference to the sounds in the master sound file. Must contain the following columns: 1) "sound.files": name of the .wav files, 2) "selec": unique selection identifier (within a sound file), 3) "start": start time and 4) "end": end time of selections, 5)  "bottom.freq": low frequency for bandpass, 6) "top.freq": high frequency for bandpass, 7) "sound.id": ID of sounds used to identify counterparts across distances and 8) "reference": identity of sounds to be used as reference for each test sound (row). See \code{\link{set_reference_sounds}} for more details on the structure of 'X'.
#' @param parallel DEPRECATED. Use 'cores' instead.
#' @param cores Numeric vector of length 1. Controls whether parallel computing is applied by specifying the number of cores to be used. Default is 1 (i.e. no parallel computing).
#' If \code{NULL} (default) then the current working directory is used.
#' @param pb Logical argument to control if progress bar is shown. Default is \code{TRUE}.
#' @param cor.method Character string indicating the correlation coefficient to be applied ("pearson", "spearman", or "kendall", see \code{\link[stats]{cor}}).
#' @param spec.smooth Numeric vector of length 1 determining the length of the sliding window used for a sum smooth for power spectrum calculation (in kHz). Default is 5.
#' @param output DEPRECATED. Now the output format mirrors the class of the input 'X'.
#' @param hop.size A numeric vector of length 1 specifying the time window duration (in ms). Default is 11.6 ms, which is equivalent to 512 wl for a 44.1 kHz sampling rate. Ignored if 'wl' is supplied.
#' @param wl A numeric vector of length 1 specifying the window length of the spectrogram, default
#' is NULL. If supplied, 'hop.size' is ignored.
#' @param ovlp Numeric vector of length 1 specifying the percent overlap between two
#'   consecutive windows, as in \code{\link[seewave]{spectro}}. Default is 70.
#' @param path Character string containing the directory path where the sound files are found. Only needed when 'X' is not an extended selection table.
#' @return Object 'X' with two additional columns, 'reference' and 'spectrum.correlation', containing the id of the sound used as reference and the computed frequency spectrum correlation coefficients, respectively.
#' @export
#' @name spectrum_correlation
#' @details spectral correlation measures the similarity of two sounds in the frequency domain. The function measures the spectral correlation coefficients of sounds in which a reference playback has been re-recorded at increasing distances. Values range from 1 (identical frequency spectrum, i.e. no degradation) to 0. The 'sound.id' column must be used to indicate the function to only compare sounds belonging to the same category (e.g. song-types). The function will then compare each sound to the corresponding reference sound. Two methods for computing spectral correlation are provided (see 'method' argument). The function uses \code{\link[seewave]{meanspec}} internally to compute power spectra. Use \code{\link{spectrum_blur_ratio}} to extract raw spectra values.
#' @examples {
#'   # load example data
#'   data("test_sounds_est")
#'
#'   # method 1
#'   # add reference column
#'   Y <- set_reference_sounds(X = test_sounds_est)
#'
#'   # run spectrum correlation
#'   spectrum_correlation(X = Y)
#'
#'   # method 2
#'   Y <- set_reference_sounds(X = test_sounds_est, method = 2)
#'   # spectrum_correlation(X = Y)
#' }
#'
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})
#' @seealso \code{\link{envelope_correlation}}, \code{\link{spectrum_blur_ratio}}
#' @references {
#' Araya-Salas, M. (2020). baRulho: baRulho: quantifying degradation of (animal) acoustic signals in R. R package version 1.0.2
#'
#' Apol, C.A., Sturdy, C.B. & Proppe, D.S. (2017). Seasonal variability in habitat structure may have shaped acoustic signals and repertoires in the black-capped and boreal chickadees. Evol Ecol. 32:57-74.
#' }

spectrum_correlation <-
  function(X,
           parallel = NULL,
           cores = getOption("mc.cores", 1),
           pb = getOption("pb", TRUE),
           cor.method = "pearson",
           spec.smooth = getOption("spec.smooth", 5),
           output = NULL,
           hop.size = getOption("hop.size", 11.6),
           wl = getOption("wl", NULL),
           ovlp = getOption("ovlp", 70),
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
    if (!(wl %% 2) == 0) {
      wl <- wl + 1
    }

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
    if (pb) {
      write(file = "", x = "Computing power spectra (step 1 out of 2):")
    }

    # calculate all spectra apply function
    specs <-
      warbleR:::pblapply_wrblr_int(
        pbar = pb,
        X = target_sgnl_temp,
        cl = cl,
        FUN = function(x,
                       ssmth = spec.smooth,
                       wln = wl,
                       Q = X,
                       pth = path,
                       ovl = ovlp) {
          spctr_FUN(
            y = x,
            spec.smooth = ssmth,
            wl = wln,
            X = Q,
            meanspc = TRUE,
            path = pth,
            ovlp = ovl
          )
        }
      )

    # add sound file selec names to envelopes (weird column name so it does not overwrite user columns)
    names(specs) <- target_sgnl_temp

    if (pb) {
      write(file = "", x = "Computing spectrum correlations (step 2 out of 2):")
    }

    # calculate all envelops apply function
    X$spectrum.correlation <-
      unlist(warbleR:::pblapply_wrblr_int(
        X = seq_len(nrow(X)),
        pbar = pb,
        cl = cl,
        FUN =
          function(x,
                   spcs = specs,
                   cm = cor.method,
                   Q = X) {
            spctr_cor_FUN(
              y = x,
              specs = spcs,
              X = Q,
              cor.method = cm
            )
          }
      ))

    # remove temporal columns
    X$.sgnl.temp <- NULL

    # fix call if not a data frame
    if (!is.data.frame(X)) {
      attributes(X)$call <-
        base::match.call()
    } # fix call attribute


    return(X)
  }
