#' Fix misaligned start and end of test sounds
#'
#' \code{realign_test_sounds} fixes misaligned start and end of test sounds in an extended selection table using spectrographic cross-correlation
#' @usage realign_test_sounds(X, parallel = NULL, cores = getOption("mc.cores", 1),
#' pb = getOption("pb", TRUE), hop.size = getOption("hop.size", 11.6),
#' wl = getOption("wl", NULL), ovlp = getOption("ovlp", 90), wn = 'hanning')
#' @param X object of class 'extended_selection_table' created by the function \code{\link[warbleR]{selection_table}} from the warbleR package. The object must include the following additional columns: 'sound.id', 'bottom.freq' and 'top.freq'.
#' @param parallel DEPRECATED. Use 'cores' instead.
#' @param cores Numeric vector of length 1. Controls whether parallel computing is applied by specifying the number of cores to be used. Default is 1 (i.e. no parallel computing).
#' @param pb Logical argument to control if progress bar is shown. Default is \code{TRUE}.
#' @param hop.size A numeric vector of length 1 specifying the time window duration (in ms). Default is 11.6 ms, which is equivalent to 512 wl for a 44.1 kHz sampling rate. Ignored if 'wl' is supplied.
#' @param wl A numeric vector of length 1 specifying the window length of the spectrogram, default
#' is NULL. If supplied, 'hop.size' is ignored.
#' @param ovlp Numeric vector of length 1 specifying \% of overlap between two
#' consecutive windows, as in \code{\link[seewave]{spectro}}. Default is 90. High values of ovlp
#' slow down the function but produce more accurate results.
#' @param wn A character vector of length 1 specifying the window name as in \code{\link[seewave]{ftwindow}}.
#' @return Object 'X' in which time parameters (columns 'start' and 'end') have been tailored to more closely match the start and end of the reference sound.
#' @export
#' @name realign_test_sounds
#' @details This function uses spectrographic cross-correlation to align the position in time of sounds with regard to a reference sound. The sound recorded at the closest distance to the source is used as reference. Precise alignment is crucial for downstream measures of sound degradation. The function calls warbleR's \code{\link[warbleR]{cross_correlation}} internally to align sounds using cross-correlation. The output extended selection table contains the new start and end values after alignment.
#' @examples
#' {
#'   # load example data
#'   data("degradation_est")
#'
#'   # create subset of data with only re-recorded files
#'   rerecorded_est <- degradation_est[degradation_est$sound.files != "master.wav", ]
#'
#'   # create "unaligned_rerecorded_est" by adding noise to "rerecorded_est" start and end
#'   unaligned_rerecorded_est <- rerecorded_est
#'   set.seed(123)
#'   noise_time <- sample(c(0.005, -0.005, 0.006, -0.006, 0, 0.002, -0.002),
#'     nrow(unaligned_rerecorded_est),
#'     replace = TRUE
#'   )
#'   attr(unaligned_rerecorded_est, "check.res")$start <-
#'     unaligned_rerecorded_est$start <- unaligned_rerecorded_est$start + noise_time
#'   attr(unaligned_rerecorded_est, "check.res")$end <- unaligned_rerecorded_est$end <-
#'     unaligned_rerecorded_est$end + noise_time
#'
#'   # re align
#'   realign_test_sounds(X = unaligned_rerecorded_est)
#' }
#'
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})
#' @seealso \code{\link{blur_ratio}}, \code{\link[warbleR]{cross_correlation}}
#' @references {
#' Araya-Salas, M. (2020). baRulho: baRulho: quantifying habitat-induced degradation of (animal) acoustic signals in R. R package version 1.0.2
#'
#' Clark, C.W., Marler, P. & Beeman K. (1987). Quantitative analysis of animal vocal phonology: an application to Swamp Sparrow song. Ethology. 76:101-115.
#' }

realign_test_sounds <-
  function(X,
           parallel = NULL,
           cores = getOption("mc.cores", 1),
           pb = getOption("pb", TRUE),
           hop.size = getOption("hop.size", 11.6),
           wl = getOption("wl", NULL),
           ovlp = getOption("ovlp", 90),
           wn = "hanning") {
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


    # set number of processes for printing message
    if (pb) {
      on.exit(options(warbleR.steps = 3), add = TRUE)

      options(warbleR.steps = 3)
    }

    # adjust wl based on hope.size
    if (is.null(wl)) {
      wl <- round(attr(X, "check.results")$sample.rate[1] * hop.size, 0)
    }

    # make wl even if odd
    if (!(wl %% 2) == 0) {
      wl <- wl + 1
    }

    # remove ambient if any from sound types
    sig.types <- setdiff(unique(X$sound.id), "ambient")

    # create matrix containing pairwise comparisons of selections (2 columns)
    comp_mats <- lapply(sig.types, function(x) {
      # extract for single sound and order by distance
      Y <- as.data.frame(X[X$sound.id == x, ])

      # create selec ID column (unique ID for each selection (row))
      Y$sf.selec <- paste(Y$sound.files, Y$selec, sep = "-")

      # create matrix with 2 columns of the selections to be compare
      # comparing to closest distance to source
      cmp.mt <-
        cbind(Y$sf.selec[-which.min(Y$distance)], Y$sf.selec[which.min(Y$distance)])

      return(cmp.mt)
    })

    # put together in a single
    comp_mat <- do.call(rbind, comp_mats)

    # resave X
    Y <- X

    # get index of sounds that would be align
    indx.algn <-
      which(paste(Y$sound.files, Y$selec, sep = "-") %in% comp_mat[, 1])

    # fix end to include half the duration of the selection at both sides
    attr(Y, "check.results")$end[indx.algn] <- Y$end[indx.algn] <-
      sapply(indx.algn, function(x) {
        wv.info <- warbleR::read_wave(X, index = x, header = TRUE)
        mxdur <- wv.info$samples / wv.info$sample.rate
        new.end <- X$end[x] + (X$end[x] - X$start[x]) * 0.7

        if (mxdur < new.end) {
          return(mxdur)
        } else {
          return(new.end)
        }
      })

    # fix start in the same way
    attr(Y, "check.results")$start[indx.algn] <-
      Y$start[indx.algn] <-
      X$start[indx.algn] - (X$end[indx.algn] - X$start[indx.algn]) * 0.7

    # make 0 any negatives
    attr(Y, "check.results")$start[Y$start < 0] <-
      Y$start[Y$start < 0] <- 0

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
    options("int_warbleR_steps" = c(current = 0, total = 2))

    on.exit(options("int_warbleR_steps" = c(current = 0, total = 0)), add = TRUE)

    warbleR_options(
      wl = wl,
      ovlp = ovlp,
      wn = wn,
      parallel = cores,
      pb = pb,
      compare.matrix = comp_mat,
      X = Y
    )

    # run spcc
    xcorrs <- warbleR::cross_correlation(output = "list")

    # message
    if (pb) {
      write(file = "", x = "finding peaks and aligning (step 2 out of 2)")
    }

    # find peaks and lags
    peaks <-
      find_peaks_bRlh_int(
        xc.output = xcorrs,
        cores = cores,
        max.peak = TRUE
      )

    # fix start and end in original data set and its attributes
    # start
    attr(X, "check.results")$start[paste(X$sound.files, X$selec, sep = "-") %in% comp_mat[, 1]] <-
      X$start[paste(X$sound.files, X$selec, sep = "-") %in% comp_mat[, 1]] <-
      peaks$start

    # end
    attr(X, "check.results")$end[paste(X$sound.files, X$selec, sep = "-") %in% comp_mat[, 1]] <-
      X$end[paste(X$sound.files, X$selec, sep = "-") %in% comp_mat[, 1]] <-
      X$start[paste(X$sound.files, X$selec, sep = "-") %in% comp_mat[, 1]] +
      (peaks$end - peaks$start)

    attr(X, "check.results")$duration <-
      attr(X, "check.results")$end - attr(X, "check.results")$start

    # fix call attribute
    attributes(X)$call <- base::match.call()

    return(X)
  }

##############################################################################################################
#' alternative name for \code{\link{realign_test_sounds}}
#'
#' @keywords internal
#' @details see \code{\link{realign_test_sounds}} for documentation. \code{\link{spcc_align}} will be deprecated in future versions.
#' @export

spcc_align <- realign_test_sounds
