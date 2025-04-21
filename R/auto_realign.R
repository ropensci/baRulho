#' Fix small misalignments in the time position test sounds
#'
#' \code{auto_realign} fixes small misalignments in the time position of test sounds in an extended selection table using spectrographic cross-correlation
#' @inheritParams template_params
#' @param X object of class 'extended_selection_table' (created by the function \code{\link[warbleR]{selection_table}} from the warbleR package) with the test sound files' annotations to be aligned. Must contain the following columns: 1) "sound.files": name of the .wav files, 2) "selec": unique selection identifier (within a sound file), 3) "start": start time and 4) "end": end time of selections, 5)  "bottom.freq": low frequency for bandpass, 6) "top.freq": high frequency for bandpass and 7) "sound.id": ID of sounds used to identify counterparts across distances. Each sound must have a unique ID within a given distance.. The object must include the following additional columns: 'sound.id', 'bottom.freq' and 'top.freq'.
#' @param Y object of class 'extended_selection_table' (a class created by the function \code{\link[warbleR]{selection_table}} from the warbleR package) with the master sound file annotations. This should be the same data than that was used for finding the position of markers in \code{\link{find_markers}}. It should also contain a 'sound.id' column.
#' @param ovlp Numeric vector of length 1 specifying the percentage of overlap between two
#' consecutive windows, as in \code{\link[seewave]{spectro}}. Default is 90. High values slow down the function but produce more accurate results. Can be set globally for the current R session via the "ovlp" option (see \code{\link[base]{options}}).
#' @param bp Numeric vector of length 2 giving the lower and upper limits of a frequency bandpass filter (in kHz). Default is \code{NULL}.
#' @return Object 'X' in which time parameters (columns 'start' and 'end') have been tailored to more closely match the start and end of the reference sound.
#' @export
#' @name auto_realign
#' @details Precise alignment is crucial for downstream measures of sound degradation. This function uses spectrogram cross-correlation to improve the time position alignment of test sounds. The master sound file is used as reference. The function calls warbleR's \code{\link[warbleR]{cross_correlation}} internally to align sounds using cross-correlation. The output extended selection table contains the new start and end values after alignment.  \strong{Note that 1) this function only works to further improve alignments if the estimated position of the test sound is already close to the actual position and 2) both 'X' and 'Y' must be extended selection tables sensu \code{\link[warbleR]{selection_table}}}. The function might not work properly with annotations with a small frequency range (e.g. pure tones).
#' 
#' @examples {
#'   # load example data
#'   data("test_sounds_est")
#'   data("master_est")
#'   
#'   # create "unaligned_test_sounds_est" by
#'   # adding error to "test_sounds_est" start and end
#'   unaligned_test_sounds_est <- test_sounds_est
#'   set.seed(123)
#'   noise_time <- sample(c(0.009, -0.01, 0.03, -0.03, 0, 0.07, -0.007),
#'   nrow(unaligned_test_sounds_est),
#'   replace = TRUE)
#'   
#'   attr(unaligned_test_sounds_est, "check.res")$start <- 
#'   unaligned_test_sounds_est$start <- 
#'   unaligned_test_sounds_est$start + noise_time
#'   attr(unaligned_test_sounds_est, "check.res")$end <- 
#'   unaligned_test_sounds_est$end  <- 
#'   unaligned_test_sounds_est$end + noise_time
#'
#' # re align
#' realigned_est <- auto_realign(X = unaligned_test_sounds_est, Y = master_est)
#' }
#'
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})
#' @family test sound alignment
#' @seealso \code{\link{blur_ratio}}, \code{\link[warbleR]{cross_correlation}}
#' @references 
#' Araya-Salas, M., Grabarczyk, E. E., Quiroz-Oliva, M., Garcia-Rodriguez, A., & Rico-Guevara, A. (2025). Quantifying degradation in animal acoustic signals with the R package baRulho. Methods in Ecology and Evolution, 00, 1-12. https://doi.org/10.1111/2041-210X.14481
#'
#' Clark, C.W., Marler, P. & Beeman K. (1987). Quantitative analysis of animal vocal phonology: an application to Swamp Sparrow song. Ethology. 76:101-115.

auto_realign <-
  function(X,
           Y,
           cores = getOption("mc.cores", 1),
           pb = getOption("pb", TRUE),
           hop.size = getOption("hop.size", 11.6),
           wl = getOption("wl", NULL),
           ovlp = getOption("ovlp", 90),
           wn = c("hanning", "hamming", "bartlett", "blackman", "flattop", "rectangle"), 
           bp = NULL) {
    
    # assign a value to wn
    wn <- rlang::arg_match(wn)
    
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
    
    # add column to discriminate between reference and unaligned
    X$..type.. <- "unaligned"
    Y$..type.. <- "reference"
    
    # combine the two extended selection tables
    common_cols <- intersect(names(X), names(Y))  
    W <- rbind(X[, common_cols], Y[, common_cols])
    
    # adjust wl based on hop.size
    wl <- .adjust_wl(wl, W, hop.size)
    
    # remove ambient if any from sound types
    sig.types <- setdiff(unique(W$sound.id), c("ambient", "start_marker", "end_marker"))
    
    # create matrix containing pairwise comparisons of selections (2 columns)
    comp_mats <- lapply(sig.types, function(x) {
      # extract for single sound and order by distance
      Q <- as.data.frame(W[W$sound.id == x,])
      
      # create selec ID column (unique ID for each selection (row))
      Q$sf.selec <- paste(Q$sound.files, Q$selec, sep = "-")
      
      # create matrix with 2 columns of the selections to be compare
      # comparing to closest distance to source
      cmp.mt <-
        cbind(Q$sf.selec[Q$..type.. == "unaligned"], Q$sf.selec[Q$..type.. == "reference"])
      
      return(cmp.mt)
    })
    
    # put together in a single
    comp_mat <- do.call(rbind, comp_mats)
    
    # resave W
    Z <- W
    
    # get index of sounds that would be align
    indx.algn <-
      which(paste(Z$sound.files, Z$selec, sep = "-") %in% comp_mat[, 1])
    
    # fix end to include half the duration of the selection at both sides
    attr(Z, "check.results")$end[indx.algn] <- Z$end[indx.algn] <-
      vapply(indx.algn, function(x) {
        wv.info <- warbleR::read_sound_file(W, index = x, header = TRUE)
        mxdur <- wv.info$samples / wv.info$sample.rate
        new.end <- W$end[x] + (W$end[x] - W$start[x]) * 0.7
        
        if (mxdur < new.end) {
          return(mxdur)
        } else {
          return(new.end)
        }
      }, FUN.VALUE = numeric(1))
    
    # fix start in the same way
    attr(Z, "check.results")$start[indx.algn] <-
      Z$start[indx.algn] <-
      W$start[indx.algn] - (W$end[indx.algn] - W$start[indx.algn]) * 0.7
    
    # make 0 any negatives
    attr(Z, "check.results")$start[Z$start < 0] <-
      Z$start[Z$start < 0] <- 0
    
    # save previous warbleR options
    prev_wl <- .Options$warbleR
    
    on.exit(
      warbleR_options(
        wl = prev_wl$wl,
        ovlp = prev_wl$ovlp,
        wn = prev_wl$wn,
        parallel = prev_wl$parallel,
        pb = prev_wl$pb
      ),
      add = TRUE
    )
    
    # steps for warbleR message
      if (pb)  {
    warbleR:::.update_progress(total = 2)
        }

    warbleR_options(
      wl = wl,
      ovlp = ovlp,
      wn = wn,
      parallel = cores,
      pb = pb,
      compare.matrix = comp_mat,
      X = Z, 
      bp = bp
    )
    
    # run spcc
    xcorrs <- warbleR::cross_correlation(output = "list",  wl = 256, method = 1)
    
    # find peaks and lags
    peaks <-
      .find_peaks(xc.output = xcorrs,
                  cores = cores,
                  max.peak = TRUE, 
                  pb = pb)
    
    
    # add column with original sound file selec labels
    peaks$sound.files.selec <- comp_mat[, 1]
    
    # fix start and end in original data set and its attributes
    # start
    for (i in peaks$sound.files.selec){
      attr(X, "check.results")$start[paste(X$sound.files, X$selec, sep = "-") == i] <-
      X$start[paste(X$sound.files, X$selec, sep = "-") == i] <- peaks$start[peaks$sound.files.selec == i]
    
      attr(X, "check.results")$end[paste(X$sound.files, X$selec, sep = "-") == i] <-
        X$end[paste(X$sound.files, X$selec, sep = "-") == i] <- peaks$end[peaks$sound.files.selec == i]
    }
    
    # recalculate duration in attribute metadata
    attr(X, "check.results")$duration <-
      attr(X, "check.results")$end - attr(X, "check.results")$start
    
    # fix call attribute
    attributes(X)$call <- base::match.call()
    
    # remove temporary column
    X$..type.. <- NULL
    
    return(X)
  }
