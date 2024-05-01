#' Fix small misalignments in the time position test sounds
#'
#' \code{auto_realign} fixes small misalignments in the time position of test sounds in an extended selection table using spectrographic cross-correlation
#' @inheritParams template_params
#' @param X object of class 'extended_selection_table' created by the function \code{\link[warbleR]{selection_table}} from the warbleR package. The object must include the following additional columns: 'sound.id', 'bottom.freq' and 'top.freq'.
#' @param Y object of class 'extended_selection_table' (a class created by the function \code{\link[warbleR]{selection_table}} from the warbleR package) with the master sound file annotations. This should be the same data than that was used for finding the position of markers in \code{\link{find_markers}}. It should also contain a 'sound.id' column.
#' @param ovlp Numeric vector of length 1 specifying the percentage of overlap between two
#' consecutive windows, as in \code{\link[seewave]{spectro}}. Default is 90. High values slow down the function but produce more accurate results.
#' @return Object 'X' in which time parameters (columns 'start' and 'end') have been tailored to more closely match the start and end of the reference sound.
#' @export
#' @name auto_realign
#' @details Precise alignment is crucial for downstream measures of sound degradation. This function uses spectrographic cross-correlation to align the position in time of test sounds. The master sound file is used as reference. The function calls warbleR's \code{\link[warbleR]{cross_correlation}} internally to align sounds using cross-correlation. The output extended selection table contains the new start and end values after alignment. Note that this function only works to further improve alignments if the estimated position of the test sound is already close to the actual position. Note that both 'X' and 'Y' must be extended selection tables sensu \code{\link[warbleR]{selection_table}}.
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
#' @references {
#' Araya-Salas M., E. Grabarczyk, M. Quiroz-Oliva, A. Garcia-Rodriguez, A. Rico-Guevara. (2023), baRulho: an R package to quantify degradation in animal acoustic signals .bioRxiv 2023.11.22.568305.
#'
#' Clark, C.W., Marler, P. & Beeman K. (1987). Quantitative analysis of animal vocal phonology: an application to Swamp Sparrow song. Ethology. 76:101-115.
#' }

auto_realign <-
  function(X,
           Y,
           cores = getOption("mc.cores", 1),
           pb = getOption("pb", TRUE),
           hop.size = getOption("hop.size", 11.6),
           wl = getOption("wl", NULL),
           ovlp = getOption("ovlp", 90),
           wn = c("hanning", "hamming", "bartlett", "blackman", "flattop", "rectangle")) {
    
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
    
    
    # combine the two extended selection tables
    common_cols <- intersect(names(X), names(Y))  
    W <- rbind(X[, common_cols], Y[, common_cols])
    
    # set number of processes for printing message
    if (pb) {
      on.exit(options(warbleR.steps = 3))
      
      options(warbleR.steps = 3)
    }
    
    # adjust wl based on hop.size
    wl <- .adjust_wl(wl, W, hop.size)

    # remove ambient if any from sound types
    sig.types <- setdiff(unique(W$sound.id), c("ambient", "start_marker", "end_marker"))
    
    # create matrix containing pairwise comparisons of selections (2 columns)
    comp_mats <- lapply(sig.types, function(x) {
      # extract for single sound and order by distance
      Y <- as.data.frame(W[W$sound.id == x,])
      
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
    
    # resave W
    Z <- W
    
    # get index of sounds that would be align
    indx.algn <-
      which(paste(Z$sound.files, Z$selec, sep = "-") %in% comp_mat[, 1])
    
    # fix end to include half the duration of the selection at both sides
    attr(Z, "check.results")$end[indx.algn] <- Z$end[indx.algn] <-
      vapply(indx.algn, function(x) {
        wv.info <- warbleR::read_wave(W, index = x, header = TRUE)
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
    options("int_warbleR_steps" = c(current = 0, total = 2))
    
    on.exit(options("int_warbleR_steps" = c(current = 0, total = 0)), add = TRUE)
    
    warbleR_options(
      wl = wl,
      ovlp = ovlp,
      wn = wn,
      parallel = cores,
      pb = pb,
      compare.matrix = comp_mat,
      X = Z
    )
    
    # run spcc
    xcorrs <- warbleR::cross_correlation(output = "list")
    
    # message
    if (pb) {
      write(file = "", x = "finding peaks and aligning (step 2 out of 2)")
    }
    
    # find peaks and lags
    peaks <-
      .find_peaks(xc.output = xcorrs,
                  cores = cores,
                  max.peak = TRUE)
    
    
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
    
    return(X)
  }
