#' Align start and end of signal using spectrographic cross-correlation  
#' 
#' \code{spcc_align} aligns start and end of signal in an extended selection table using spectrographic cross-correlation  
#' @usage spcc_align(X = NULL, parallel = 1, pb = TRUE, 
#' hop.size = 11.6, wl = NULL, ovlp = 90, wn = 'hanning')
#' @param X object of class 'selection_table', 'extended_selection_table' created by the function \code{\link[warbleR]{selection_table}} from the warbleR package. The object must include the following additional columns: 'signal.type', 'bottom.freq' and 'top.freq'.
#' @param parallel Numeric vector of length 1. Controls whether parallel computing is applied by specifying the number of cores to be used. Default is 1 (i.e. no parallel computing).
#' @param pb Logical argument to control if progress bar is shown. Default is \code{TRUE}.
#' @param hop.size A numeric vector of length 1 specifying the time window duration (in ms). Default is 11.6 ms, which is equivalent to 512 wl for a 44.1 kHz sampling rate. Ignored if 'wl' is supplied.
#' @param wl A numeric vector of length 1 specifying the window length of the spectrogram, default 
#' is NULL. If supplied, 'hop.size' is ignored.
#' @param ovlp Numeric vector of length 1 specifying \% of overlap between two 
#' consecutive windows, as in \code{\link[seewave]{spectro}}. Default is 90. High values of ovlp 
#' slow down the function but produce more accurate results.
#' @param wn A character vector of length 1 specifying the window name as in \code{\link[seewave]{ftwindow}}. 
#' @return Extended selection table similar to input data in which time parameters (columns 'start' and 'end') have been tailored to more closely match the start and end of the reference signal. 
#' @export
#' @name spcc_align
#' @details This function uses spectrographic cross-correlation to align the position in time of signals with regard to a reference signal. The signal recorded at the closest distance to the source is used as reference. Precise alignment is crucial for downstream measures of signal degradation. The function calls warbleR's \code{\link[warbleR]{xcorr}} and \code{\link[warbleR]{find_peaks}} internally to align signals using cross-correlation. The output extended selection table contains the new start and end values after alignment. 
#' @examples
#' {
#' # load example data
#' data("playback_est_unaligned")
#' 
#' # method 1
#'spcc_align(X = playback_est_unaligned)
#' }
#' 
#' @author Marcelo Araya-Salas (\email{marceloa27@@gmail.com})
#' @seealso \code{\link{blur_ratio}}, \code{\link[warbleR]{xcorr}}
#' @references {
#' Araya-Salas, M. (2020), baRulho: baRulho: quantifying habitat-induced degradation of (animal) acoustic signals in R. R package version 1.0.0
#' 
#' Clark, C.W., Marler, P. & Beeman K. (1987). Quantitative analysis of animal vocal phonology: an application to Swamp Sparrow song. Ethology. 76:101-115. 
#' }
# last modification on nov-01-2019 (MAS)

spcc_align <- function(X = NULL, parallel = 1, pb = TRUE, hop.size = 11.6, wl = NULL, ovlp = 90, wn = 'hanning'){
  
  # set number of processes for printing message 
  if (pb){
    
    on.exit(options(warbleR.steps = 3), add = TRUE)
    
    options(warbleR.steps = 3)
  }
  
  # is extended sel tab
  if (!warbleR::is_extended_selection_table(X)) 
    stop("'X' must be and extended selection table")
  
  # If parallel is not numeric
  if (!is.numeric(parallel)) stop("'parallel' must be a numeric vector of length 1") 
  if (any(!(parallel %% 1 == 0),parallel < 1)) stop("'parallel' should be a positive integer")
  
  # must have the same sampling rate
  if (length(unique(attr(X, "check.results")$sample.rate)) > 1) 
    stop("all wave objects in the extended selection table must have the same sampling rate (they can be homogenized using warbleR::resample_est())")
  
  # hopsize  
  if (!is.numeric(hop.size) | hop.size < 0) stop("'parallel' must be a positive number") 
  
  # adjust wl based on hope.size
  if (is.null(wl))
    wl <- round(attr(X, "check.results")$sample.rate[1] * hop.size, 0)
  
  # check signal.type column 
  if (is.null(X$signal.type)) stop("'X' must containe a 'signal.type' column")
  
  #remove ambient if any from signal types
  sig.types <- setdiff(unique(X$signal.type), "ambient")
  
  # create matrix containing pairwise comparisons of selections (2 columns)
  comp_mats <- lapply(sig.types, function(x){
    
    # extract for single signal and order by distance
    Y <- as.data.frame(X[X$signal.type == x, ])
    
    # create selec ID column (unique ID for each selection (row)) 
    Y$sf.selec <- paste(Y$sound.files, Y$selec, sep = "-")
    
    # create matrix with 2 columns of the selections to be compare 
    # comparing to closest distance to source
    cmp.mt <- cbind(Y$sf.selec[-which.min(Y$distance)], Y$sf.selec[which.min(Y$distance)]) 
    
    return(cmp.mt)
  })
  
  # put together in a single
  comp_mat <- do.call(rbind, comp_mats)
  
  # resave X
  Y <- X
  
  # get index of signals that would be align
  indx.algn <- which(paste(Y$sound.files, Y$selec, sep = "-") %in% comp_mat[, 1])
  
  # fix end to include half the duration of the selection at both sides
  attr(Y, "check.results")$end[indx.algn] <- Y$end[indx.algn] <- 
    sapply(indx.algn, function(x)
    {
      wv.info <- warbleR::read_wave(X, index = x, header = TRUE)
      mxdur <- wv.info$samples / wv.info$sample.rate
      new.end <-  X$end[x] + (X$end[x] -  X$start[x]) * 0.7
      
      if(mxdur < new.end) return(mxdur) else return(new.end)
      })

  # fix start in the same way
  attr(Y, "check.results")$start[indx.algn] <- Y$start[indx.algn] <- X$start[indx.algn] - (X$end[indx.algn] -  X$start[indx.algn]) * 0.7
  attr(Y, "check.results")$start[Y$start < 0] <- Y$start[Y$start < 0] <- 0

  # run spcc 
  xcorrs <- warbleR::xcorr(X = Y, wl = wl, ovlp = ovlp, wn = wn, parallel = parallel, pb = pb, compare.matrix = comp_mat, output = "list")
  
  # message
  if (pb)  
    write(file = "", x = "finding peaks and aligning (step 3 of 3)")
  
  # find peaks and lags
  peaks <- warbleR::find_peaks(xc.output = xcorrs, parallel = parallel, max.peak = TRUE)
  
  # fix start and end in original data set and its attributes 
  # start
  attr(X, "check.results")$start[paste(X$sound.files, X$selec, sep = "-") %in% comp_mat[, 1]] <- X$start[paste(X$sound.files, X$selec, sep = "-") %in% comp_mat[, 1]] <- 
    peaks$start + Y$start[paste(Y$sound.files, X$selec, sep = "-") %in% comp_mat[, 1]]
  
  # end
  attr(X, "check.results")$end[paste(X$sound.files, X$selec, sep = "-") %in% comp_mat[, 1]] <- 
  X$end[paste(X$sound.files, X$selec, sep = "-") %in% comp_mat[, 1]] <- 
    peaks$end + Y$start[paste(Y$sound.files, X$selec, sep = "-") %in% comp_mat[, 1]]

  attr(X, "check.results")$duration <- attr(X, "check.results")$end - attr(X, "check.results")$start
  
  return(X)
}
