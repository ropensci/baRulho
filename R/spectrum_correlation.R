#' Measure frequency spectrum correlation
#'
#' \code{spectrum_correlation} measures frequency spectrum correlation of sounds referenced in an extended selection table.
#' @usage spectrum_correlation(X, parallel = NULL, cores = getOption("mc.cores", 1),
#' pb = getOption("pb", TRUE), method = getOption("method", 1), cor.method = "pearson",
#' output = NULL, hop.size = getOption("hop.size", 11.6),
#' wl = getOption("wl", NULL), ovlp = getOption("ovlp", 70),
#' path = getOption("sound.files.path", "."))
#' @param X Object of class 'data.frame', 'selection_table' or 'extended_selection_table' (the last 2 classes are created by the function \code{\link[warbleR]{selection_table}} from the warbleR package) with the reference to the sounds in the master sound file. Must contain the following columns: 1) "sound.files": name of the .wav files, 2) "selec": unique selection identifier (within a sound file), 3) "start": start time and 4) "end": end time of selections, 5)  "bottom.freq": low frequency for bandpass, 6) "top.freq": high frequency for bandpass, 7) "sound.id": ID of sounds used to identify counterparts across distances and 8) "distance": distance at which each test sound was re-recorded. Each sound must have a unique ID within a distance. An additional 'transect' column labeling those sounds recorded in the same transect is required if 'method = 2'.
#' @param parallel DEPRECATED. Use 'cores' instead.
#' @param cores Numeric vector of length 1. Controls whether parallel computing is applied by specifying the number of cores to be used. Default is 1 (i.e. no parallel computing).
#' If \code{NULL} (default) then the current working directory is used.
#' @param pb Logical argument to control if progress bar is shown. Default is \code{TRUE}.
#' @param method Numeric vector of length 1 to indicate the 'experimental design' to measure frequency spectrum correlation. Two methods are available:
#' \itemize{
#' \item \code{1}: compare all sounds with their counterpart that was recorded at the closest distance to source (e.g. compare a sound recorded at 5m, 10m and 15m with its counterpart recorded at 1m). This is the default method.
#' \item \code{2}: compare all sounds with their counterpart recorded at the distance immediately before (e.g. a sound recorded at 10m compared with the same sound recorded at 5m, then sound recorded at 15m compared with same sound recorded at 10m and so on).
#' }
#' @param cor.method Character string indicating the correlation coefficient to be applied ("pearson", "spearman", or "kendall", see \code{\link[stats]{cor}}).
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
#' @details spectral correlation measures the similarity of two sounds in the frequency domain. The function measures the spectral correlation coefficients of sounds in which a reference playback has been re-recorded at increasing distances. Values range from 1 (identical frequency spectrum, i.e. no degradation) to 0. The 'sound.id' column must be used to indicate the function to only compare sounds belonging to the same category (e.g. song-types). The function will then compare each sound to the corresponding reference sound. Two methods for calculating spectral correlation are provided (see 'method' argument). Use \code{\link{spectrum_blur_ratio}} to get spectra for plotting.
#' @examples
#' {
#'   # load example data
#'   data("degradation_est")
#'
#'   # create subset of data with only re-recorded files
#'   rerecorded_est <- degradation_est[degradation_est$sound.files != "master.wav", ]
#'
#'   # method 1
#'   spectrum_correlation(X = rerecorded_est)
#'
#'   # method 2
#'   # spectrum_correlation(X = rerecorded_est, method = 2)
#' }
#'
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})
#' @seealso \code{\link{envelope_correlation}}, \code{\link{spectrum_blur_ratio}}
#' @references {
#' Araya-Salas, M. (2020). baRulho: baRulho: quantifying habitat-induced degradation of (animal) acoustic signals in R. R package version 1.0.2
#'
#' Apol, C.A., Sturdy, C.B. & Proppe, D.S. (2017). Seasonal variability in habitat structure may have shaped acoustic signals and repertoires in the black-capped and boreal chickadees. Evol Ecol. 32:57-74.
#' }

spectrum_correlation <-
  function(X,
           parallel = NULL,
           cores = getOption("mc.cores", 1),
           pb = getOption("pb", TRUE),
           method = getOption("method", 1),
           cor.method = "pearson",
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
    checkmate::reportAssertions(check_results)

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

    # add sound file selec column and names to X (weird column name so it does not overwrite user columns)
    if (pb) {
      write(
        file = "",
        x = paste0("Preparing data for analysis (step 1 out of 3):")
      )
    }

    X <-
      prep_X_bRlo_int(X,
        method = method,
        cores = cores,
        pb = pb
      )

    # set clusters for windows OS
    if (Sys.info()[1] == "Windows" & cores > 1) {
      cl <-
        parallel::makePSOCKcluster(getOption("cl.cores", cores))
    } else {
      cl <- cores
    }

    if (pb) {
      write(file = "", x = "Calculating frequency spectrums (step 2 out of 3):")
    }

    # calculate all spectra apply function
    spcs <-
      warbleR:::pblapply_wrblr_int(
        pbar = pb,
        X = seq_len(nrow(X)),
        cl = cl,
        FUN = function(y, wle = wl, ovl = ovlp) {
          # load clip
          clp <- warbleR::read_sound_file(
            X = X,
            index = y,
            path = path
          )

          # mean spec
          mspc <-
            meanspec(
              wave = clp,
              f = clp@samp.rate,
              plot = FALSE,
              wl = wle,
              ovlp = ovl
            )

          return(mspc)
        }
      )

    # add sound file selec names to envelopes (weird column name so it does not overwrite user columns)
    names(spcs) <- X$.sgnl.temp

    # make a data frame with 2 columns with names of the sounds to be compare
    X$reference <- sapply(seq_len(nrow(X)), function(x, meth = method) {
      # extract for single sound and order by distance
      Y <-
        as.data.frame(X[X$sound.id == X$sound.id[X$.sgnl.temp == X$.sgnl.temp[x]], , drop = FALSE])
      Y <- Y[order(Y$distance), ]

      # method 1 compare to closest distance to source
      if (meth == 1) {
        z <- Y$.sgnl.temp[which.min(Y$distance)]
      } else
      # if method 2
      # if not the first row then the previous row
      if (Y$.sgnl.temp[1] != X$.sgnl.temp[x]) {
        z <- X$.sgnl.temp[x - 1]
      } else {
        # else the first row
        z <- Y$.sgnl.temp[1]
      }

      return(z)
    })

    if (pb) {
      write(file = "", x = "Calculating spectrum correlations (step 3 out of 3):")
    }

    # calculate all envelops apply function
    X$spectrum.correlation <-
      pbapply::pbsapply(
        X = seq_len(nrow(X)),
        cl = cl,
        FUN = function(x) {
          spctr_cor_FUN(y = X$.sgnl.temp[x], z = X$reference[x], spcs, X, cor.method)
        }
      )

    # make NAs those sounds in which the reference is itself (only happens when method = 2) or is ambient noise
    X$reference[X$reference == X$.sgnl.temp | X$sound.id == "ambient"] <- NA

    # remove temporal columns
    X$.sgnl.temp <- NULL

    # fix call if not a data frame
    if (!is.data.frame(X)) {
      attributes(X)$call <-
        base::match.call()
    } # fix call attribute


    return(X)
  }
