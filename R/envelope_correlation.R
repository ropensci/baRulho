#' Measure amplitude envelope correlation
#'
#' \code{envelope_correlation} measures amplitude envelope correlation of sounds referenced in an extended selection table.
#' @usage envelope_correlation(X, parallel = NULL, cores = getOption("mc.cores", 1),
#' pb = getOption("pb", TRUE), method = getOption("method", 1),
#' cor.method = "pearson", ssmooth = NULL, msmooth = NULL, output = NULL,
#' hop.size = getOption("hop.size", 11.6), wl = getOption("wl", NULL),
#' ovlp = getOption("ovlp", 70), path = getOption("sound.files.path", "."))
#' @param X Object of class 'data.frame', 'selection_table' or 'extended_selection_table' (the last 2 classes are created by the function \code{\link[warbleR]{selection_table}} from the warbleR package) with the reference to the sounds in the master sound file. Must contain the following columns: 1) "sound.files": name of the .wav files, 2) "selec": unique selection identifier (within a sound file), 3) "start": start time and 4) "end": end time of selections, 5)  "bottom.freq": low frequency for bandpass, 6) "top.freq": high frequency for bandpass and 7) "sound.id": ID of sounds used to identify counterparts across distances. Each sound must have a unique ID within a distance.
#' @param parallel DEPRECATED. Use 'cores' instead.
#' @param cores Numeric vector of length 1. Controls whether parallel computing is applied by specifying the number of cores to be used. Default is 1 (i.e. no parallel computing).
#' If \code{NULL} (default) then the current working directory is used.
#' @param pb Logical argument to control if progress bar is shown. Default is \code{TRUE}.
#' @param method Numeric vector of length 1 to indicate the 'experimental design' to measure amplitude envelope correlation. Two methods are available:
#' \itemize{
#' \item \code{1}: compare all sounds with their counterpart that was recorded at the closest distance to source (e.g. compare a sound recorded at 5m, 10m and 15m with its counterpart recorded at 1m). This is the default method.
#' \item \code{2}: compare all sounds with their counterpart recorded at the distance immediately before (e.g. a sound recorded at 10m compared with the same sound recorded at 5m, then sound recorded at 15m compared with same sound recorded at 10m and so on).
#' }
#' @param cor.method Character string indicating the correlation coefficient to be applied ("pearson", "spearman", or "kendall", see \code{\link[stats]{cor}}).
#' @param ssmooth Numeric vector of length 1 to determine the length of the sliding window used for a sum smooth for amplitude envelope calculation (used internally by \code{\link[seewave]{env}}).
#' @param msmooth Numeric vector of length 2 to smooth the amplitude envelope with a mean sliding window for amplitude envelope calculation. The first element is the window length (in number of amplitude values) and the second one the window overlap (used internally by \code{\link[seewave]{env}}).
#' @param output DEPRECATED. Now the output format mirrors the class of the input 'X'.
#' @param hop.size A numeric vector of length 1 specifying the time window duration (in ms). Default is 11.6 ms, which is equivalent to 512 wl for a 44.1 kHz sampling rate. Ignored if 'wl' is supplied.
#' @param wl A numeric vector of length 1 specifying the window length of the spectrogram, default
#' is NULL. If supplied, 'hop.size' is ignored.
#' @param ovlp Numeric vector of length 1 specifying the percent overlap between two
#'   consecutive windows, as in \code{\link[seewave]{spectro}}. Default is 70.
#' @param path Character string containing the directory path where the sound files are found. Only needed when 'X' is not an extended selection table.
#' @return Object 'X' with two additional columns, 'reference' and 'envelope.correlation', containing the id of the sound used as reference and the computed envelope correlation coefficients, respectively.
#' @export
#' @name envelope_correlation
#' @details Amplitude envelope correlation measures the similarity of two sounds in the time domain. The function measures the envelope correlation coefficients of sounds in which a reference playback has been re-recorded at increasing distances. Values close to 1 means very similar amplitude envelopes (i.e. little degradation has occurred). If envelopes have different lengths (which means sounds have different lengths) cross-correlation is used and the maximum correlation coefficient is returned. Cross-correlation is achieved by sliding the shortest sound along the largest one and calculating the correlation at each step. The 'sound.id' column must be used to indicate the function to only compare sounds belonging to the same category (e.g. song-types). The function compares each sound to the corresponding reference sound within the supplied frequency range (e.g. bandpass) of the reference sound ('bottom.freq' and 'top.freq' columns in 'X'). Two methods for calculating envelope correlation are provided (see 'method' argument). Use \code{\link{blur_ratio}} to create envelopes graphs.
#' @seealso \code{\link{blur_ratio}}, \code{\link{spectrum_blur_ratio}}
#' @examples
#' {
#'   # load example data
#'   data("degradation_est")
#'
#'   # create subset of data with only re-recorded files
#'   rerecorded_est <- degradation_est[degradation_est$sound.files != "master.wav", ]
#'
#'   # remove ambient selections
#'   rerecorded_est <- rerecorded_est[rerecorded_est$sound.id != "ambient", ]
#'
#'   # method 1
#'   envelope_correlation(X = rerecorded_est)
#'
#'   # method 2
#'   # envelope_correlation(X = rerecorded_est, method = 2)
#' }
#'
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})
#' @references {
#' Araya-Salas, M. (2020). baRulho: baRulho: quantifying habitat-induced degradation of (animal) acoustic signals in R. R package version 1.0.2
#'
#' Apol, C.A., Sturdy, C.B. & Proppe, D.S. (2017). Seasonal variability in habitat structure may have shaped acoustic signals and repertoires in the black-capped and boreal chickadees. Evol Ecol. 32:57-74.
#' }
# last modification on nov-01-2019 (MAS)

envelope_correlation <- function(X, parallel = NULL, cores = getOption("mc.cores", 1), pb = getOption("pb", TRUE), method = getOption("method", 1), cor.method = "pearson", ssmooth = NULL, msmooth = NULL, output = NULL, hop.size = getOption("hop.size", 11.6), wl = getOption("wl", NULL), ovlp = getOption("ovlp", 70), path = getOption("sound.files.path", ".")) {
  # check arguments
  arguments <- as.list(base::match.call())

  # add objects to argument names
  for (i in names(arguments)[-1]) {
    arguments[[i]] <- get(i)
  }

  # check each arguments
  check_results <- check_arguments(fun = arguments[[1]], args = arguments)

  # report errors
  checkmate::reportAssertions(check_results)

  # adjust wl based on hope.size
  if (is.null(wl)) {
    wl <- round(read_sound_file(X, index = 1, header = TRUE, path = path)$sample.rate * hop.size / 1000, 0)
  }

  # make wl even if odd
  if (!(wl %% 2) == 0) wl <- wl + 1

  # add sound file selec column and names to X (weird column name so it does not overwrite user columns)
  if (pb) {
    write(file = "", x = paste0("Preparing data for analysis (step 1 out of 3):"))
  }

  X <- prep_X_bRlo_int(X, method = method, cores = cores, pb = pb)

  # set clusters for windows OS
  if (Sys.info()[1] == "Windows" & cores > 1) {
    cl <- parallel::makePSOCKcluster(getOption("cl.cores", cores))
  } else {
    cl <- cores
  }

  if (pb) write(file = "", x = "Calculating amplitude envelopes (step 2 out of 3):")

  # calculate all envelopes apply function
  envs <- warbleR:::pblapply_wrblr_int(pbar = pb, X = seq_len(nrow(X)), cl = cl, FUN = function(y) {
    # get clip
    clp <- warbleR::read_sound_file(X = X, index = y, path = path)

    # define bandpass based on reference
    bp <- c(X$bottom.freq[X$.sgnl.temp == X$reference[y]], X$top.freq[X$.sgnl.temp == X$reference[y]])

    # bandpass filter
    clp <- seewave::ffilter(clp,
      from = bp[1] * 1000,
      ovlp = ovlp, to = bp[2] * 1000, bandpass = TRUE,
      wl = wl, output = "Wave"
    )

    # calculate envelope
    nv <- env(wave = clp, f = clp@samp.rate, ssmooth = ssmooth, plot = FALSE, msmooth = msmooth)[, 1]

    return(nv)
  })

  # add sound file selec column and names to envelopes
  names(envs) <- X$.sgnl.temp

  # set options for loop
  if (pb) write(file = "", x = "Calculating envelope correlations (step 3 out of 3):")

  # calculate all envelops apply function
  X_list <- warbleR:::pblapply_wrblr_int(X = seq_len(nrow(X)), pbar = pb, cl = cl, FUN = function(x) {
    Y <- X[x, , drop = FALSE]
    Y$envelope.correlation <- env_cor_FUN(y = Y$.sgnl.temp, z = Y$reference, envs, cor.method)
    return(as.data.frame(Y))
  })

  X2 <- do.call(rbind, X_list)

  # make NAs those sounds in which the reference is itself (only happens when method = 2)
  X2$reference[X2$reference == X2$.sgnl.temp] <- NA

  # remove temporal columns
  X2$.sgnl.temp <- NULL

  # fix est
  if (is_extended_selection_table(X)) {
    X2 <- warbleR::fix_extended_selection_table(X = X2, Y = X)
  }

  # fix call if not a data frame
  if (!is.data.frame(X)) {
    attributes(X)$call <-
      base::match.call()
  } # fix call attribute

  return(X2)
}
