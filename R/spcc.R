#' Measure spectrographic cross-correlation as a measure of sound distortion
#'
#' \code{spcc} measures spectrographic cross-correlation as a measure of sound distortion in sounds referenced in an extended selection table.
#' @usage spcc(X, parallel = NULL, cores = getOption("mc.cores", 1), pb = getOption("pb", TRUE),
#'  cor.method = "pearson", output = NULL,
#' hop.size = getOption("hop.size", 11.6), wl = getOption("wl", NULL),
#' ovlp = getOption("ovlp", 90), wn = 'hanning', path = getOption("sound.files.path", "."))
#' @param X The output of \code{\link{set_reference_sounds}} which is an object of class 'data.frame', 'selection_table' or 'extended_selection_table' (the last 2 classes are created by the function \code{\link[warbleR]{selection_table}} from the warbleR package) with the reference to the sounds in the master sound file. Must contain the following columns: 1) "sound.files": name of the .wav files, 2) "selec": unique selection identifier (within a sound file), 3) "start": start time and 4) "end": end time of selections, 5)  "bottom.freq": low frequency for bandpass, 6) "top.freq": high frequency for bandpass, 7) "sound.id": ID of sounds used to identify counterparts across distances and 8) "reference": identity of sounds to be used as reference for each test sound (row). See \code{\link{set_reference_sounds}} for more details on the structure of 'X'.
#' @param parallel DEPRECATED. Use 'cores' instead.
#' @param cores Numeric vector of length 1. Controls whether parallel computing is applied by specifying the number of cores to be used. Default is 1 (i.e. no parallel computing).
#' @param pb Logical argument to control if progress bar is shown. Default is \code{TRUE}.
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
#' @return Object 'X' with an additional column,  'cross.correlation', containing the id of the sound used as reference and the computed spectrogram cross-correlation coefficients, respectively.
#' @export
#' @name spcc
#' @details Spectrographic cross-correlation measures frequency distortion of sounds as a similarity metric. Values close to 1 means very similar spectrograms (i.e. little sound distortion has occurred). Cross-correlation is measured of sounds in which a reference playback has been re-recorded at increasing distances. The 'sound.id' column must be used to indicate the function to only compare sounds belonging to the same category (e.g. song-types). The function compares each sound to the corresponding reference sound within the supplied frequency range (e.g. bandpass) of the reference sound ('bottom.freq' and 'top.freq' columns in 'X'). Two methods for computing cross-correlation are provided (see 'method' argument). The function is a wrapper on warbleR's \code{\link[warbleR]{cross_correlation}} function.
#' @examples {
#'   # load example data
#'   data("test_sounds_est")
#'
#'   # add reference to X
#'   X <- set_reference_sounds(X = test_sounds_est)
#'
#'   # get spcc
#'   spcc(X = X)
#' }
#'
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})
#' @family quantify degradation
#' @seealso \code{\link{blur_ratio}}, \code{\link{manual_realign}}, \code{\link[warbleR]{cross_correlation}}
#' @references {
#' Araya-Salas, M. (2020). baRulho: baRulho: quantifying degradation of (animal) acoustic signals in R. R package version 1.0.2
#'
#' Clark, C.W., Marler, P. & Beeman K. (1987). Quantitative analysis of animal vocal phonology: an application to Swamp Sparrow song. Ethology. 76:101-115.
#' }

spcc <-
  function(X,
           parallel = NULL,
           cores = getOption("mc.cores", 1),
           pb = getOption("pb", TRUE),
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
    
    
    # # put together in a single
    comp_mat <- cbind(X$.sgnl.temp, X$reference)
    
    # remove NA rows
    comp_mat <- comp_mat[stats::complete.cases(comp_mat),]
    
    
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
      warbleR::cross_correlation(X = X,
                                 cor.method = "pearson",
                                 path = path)$max.xcorr.matrix
    
    # put results back into X
    X$cross.correlation <- NA
    
    # fill score values on X
    X$cross.correlation <- vapply(seq_len(nrow(X)), function(x) {
      # get score for each row
      xc <-
        xcorrs$score[xcorrs$X1 == X$.sgnl.temp[x] &
                       xcorrs$X2 == X$reference[x]]
      
      # if empty then NA
      if (length(xc) == 0)
        xc <- NA
      
      return(xc)
    }, FUN.VALUE = numeric(1L))
    
    # fix call if not a data frame
    if (!is.data.frame(X)) {
      attributes(X)$call <-
        base::match.call()
    } # fix call attribute
    
    # remove temporary colu8mn
    X$.sgnl.temp <- NULL
    
    return(X)
  }


##############################################################################################################
#' alternative name for \code{\link{spcc}}
#'
#' @keywords internal
#' @details see \code{\link{spcc}} for documentation. \code{\link{spcc_distortion}} will be deprecated in future versions.
#' @export

spcc_distortion <- spcc
