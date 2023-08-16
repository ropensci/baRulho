#' Measure spectrographic cross-correlation as a measure of sound distortion
#'
#' \code{spcc} measures spectrographic cross-correlation as a measure of sound distortion in sounds referenced in an extended selection table.
#' @usage spcc(X, parallel = NULL, cores = getOption("mc.cores", 1), pb = getOption("pb", TRUE),
#' method = getOption("method", 1), cor.method = "pearson", output = NULL,
#' hop.size = getOption("hop.size", 11.6), wl = getOption("wl", NULL),
#' ovlp = getOption("ovlp", 90), wn = 'hanning', path = getOption("sound.files.path", "."))
#' @param X Object of class 'data.frame', 'selection_table' or 'extended_selection_table' (the last 2 classes are created by the function \code{\link[warbleR]{selection_table}} from the warbleR package) with the reference to the sounds in the master sound file. Must contain the following columns: 1) "sound.files": name of the .wav files, 2) "selec": unique selection identifier (within a sound file), 3) "start": start time and 4) "end": end time of selections, 5)  "bottom.freq": low frequency for bandpass, 6) "top.freq": high frequency for bandpass and 7) "sound.id": ID of sounds used to identify counterparts across distances. Each sound must have a unique ID within a distance.
#' @param parallel DEPRECATED. Use 'cores' instead.
#' @param cores Numeric vector of length 1. Controls whether parallel computing is applied by specifying the number of cores to be used. Default is 1 (i.e. no parallel computing).
#' @param pb Logical argument to control if progress bar is shown. Default is \code{TRUE}.
#' @param method Numeric vector of length 1 to indicate the 'experimental design' for measuring envelope correlation. Two methods are available:
#' \itemize{
#' \item \code{1}: compare all sounds with their counterpart that was recorded at the closest distance to source (e.g. compare a sound recorded at 5m, 10m and 15m with its counterpart recorded at 1m). This is the default method.
#' \item \code{2}: compare all sounds with their counterpart recorded at the distance immediately before (e.g. a sound recorded at 10m compared with the same sound recorded at 5m, then sound recorded at 15m compared with same sound recorded at 10m and so on).
#' }
#' @param cor.method Character string indicating the correlation coefficient to be applied ("pearson", "spearman", or "kendall", see \code{\link[stats]{cor}}).
#' @param output DEPRECATED. Now the output format mirrors the class of the input 'X'.
#' @param hop.size A numeric vector of length 1 specifying the time window duration (in ms). Default is 11.6 ms, which is equivalent to 512 wl for a 44.1 kHz sampling rate. Ignored if 'wl' is supplied.
#' @param wl A numeric vector of length 1 specifying the window length of the spectrogram, default
#' is \code{NULL}. If supplied, 'hop.size' is ignored.
#' @param ovlp Numeric vector of length 1 specifying \% of overlap between two
#' consecutive windows, as in \code{\link[seewave]{spectro}}. Default is 90. High values of ovlp
#' slow down the function but produce more accurate results.
#' @param wn A character vector of length 1 specifying the window name as in \code{\link[seewave]{ftwindow}}.
#' @param path Character string containing the directory path where the sound files are found. Only needed when 'X' is not an extended selection table.
#' @return Object 'X' with two additional columns, 'reference' and 'cross.correlation', containing the id of the sound used as reference and the computed spectrogram cross-correlation coefficients, respectively.
#' @export
#' @name spcc
#' @details Spectrographic cross-correlation measures frequency distortion of sounds as a similarity metric. Values close to 1 means very similar spectrograms (i.e. little sound distortion has occurred). Cross-correlation is measured of sounds in which a reference playback has been re-recorded at increasing distances. The 'sound.id' column must be used to indicate the function to only compare sounds belonging to the same category (e.g. song-types). The function compares each sound to the corresponding reference sound within the supplied frequency range (e.g. bandpass) of the reference sound ('bottom.freq' and 'top.freq' columns in 'X'). Two methods for calculating cross-correlation are provided (see 'method' argument). The function is a wrapper on warbleR's \code{\link[warbleR]{cross_correlation}} function.
#' @examples
#' {
#'   # load example data
#'   data("degradation_est")
#'
#'   # create subset of data with only re-recorded files
#'   rerecorded_est <- degradation_est[degradation_est$sound.files != "master.wav", ]
#'
#'   # method 1
#'   spcc(X = rerecorded_est, method = getOption("method", 1))
#'
#'   # method 2
#'   # spcc(X = rerecorded_est, method = 2)
#' }
#'
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})
#' @seealso \code{\link{blur_ratio}}, \code{\link{realign_test_sounds}}, \code{\link[warbleR]{cross_correlation}}
#' @references {
#' Araya-Salas, M. (2020). baRulho: baRulho: quantifying habitat-induced degradation of (animal) acoustic signals in R. R package version 1.0.2
#'
#' Clark, C.W., Marler, P. & Beeman K. (1987). Quantitative analysis of animal vocal phonology: an application to Swamp Sparrow song. Ethology. 76:101-115.
#' }
# last modification on jan-06-2020 (MAS)

spcc <-
  function(X,
           parallel = NULL,
           cores = getOption("mc.cores", 1),
           pb = getOption("pb", TRUE),
           method = getOption("method", 1),
           cor.method = "pearson",
           output = NULL,
           hop.size = getOption("hop.size", 11.6),
           wl = getOption("wl", NULL),
           ovlp = getOption("ovlp", 90),
           wn = "hanning",
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


    # remove ambient if any from sound types
    sound.ids <- setdiff(unique(X$sound.id), "ambient")

    # create matrix containing pairwise comparisons of selections (2 columns)
    comp_mats <- lapply(sound.ids, function(x) {
      # extract for single sound and order by distance
      Y <- as.data.frame(X[X$sound.id == x, ])

      # create selec ID column (unique ID for each selection (row))
      Y$sf.selec <- paste(Y$sound.files, Y$selec, sep = "-")

      # create matrix with 2 columns of the selections to be compare
      if (method == 1) {
        # comparing to closest distance to source
        cmp.mt <-
          cbind(Y$sf.selec[which.min(Y$distance)], Y$sf.selec[-which.min(Y$distance)])
      } else {
        # comparing to previous distance
        cmp.mt <- cbind(Y$sf.selec[-nrow(Y)], Y$sf.selec[-1])
      }

      return(cmp.mt)
    })

    # put together in a single
    comp_mat <- do.call(rbind, comp_mats)

    # save previous warbleR options
    prev_wl <- .Options$warbleR

    on.exit(
      warbleR_options(
        wl = prev_wl$wl,
        ovlp = prev_wl$ovlp,
        wn = prev_wl$wn,
        parallel = prev_wl$parallel,
        pb = prev_wl$pb
      )
    )

    # steps for warbleR message
    options("int_warbleR_steps" = c(current = 0, total = 1))

    on.exit(options("int_warbleR_steps" = c(current = 0, total = 0)), add = TRUE)

    warbleR_options(
      wl = wl,
      ovlp = ovlp,
      wn = wn,
      parallel = cores,
      pb = pb,
      compare.matrix = comp_mat
    )

    # run spcc
    xcorrs <-
      warbleR::cross_correlation(
        X = X,
        cor.method = "pearson",
        path = path
      )$max.xcorr.matrix

    # put results back into X
    X$reference <- NA
    X$cross.correlation <- NA

    # add correlation and reference only for calculated correlations
    X$reference[match(xcorrs$X2, paste(X$sound.files, X$selec, sep = "-"))] <-
      as.character(xcorrs$X1)

    X$cross.correlation[match(xcorrs$X2, paste(X$sound.files, X$selec, sep = "-"))] <-
      xcorrs$score

    # fix call if not a data frame
    if (!is.data.frame(X)) {
      attributes(X)$call <-
        base::match.call()
    } # fix call attribute

    return(X)
  }


##############################################################################################################
#' alternative name for \code{\link{spcc}}
#'
#' @keywords internal
#' @details see \code{\link{spcc}} for documentation. \code{\link{spcc_distortion}} will be deprecated in future versions.
#' @export

spcc_distortion <- spcc
