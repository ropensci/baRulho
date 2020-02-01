#' Measure spectrographic cross-correlation as a measure of signal distortion
#' 
#' \code{spcc} measures spectrographic cross-correlation as a measure of signal distortion in signals referenced in an extended selection table.
#' @usage spcc(X, parallel = 1, pb = TRUE,  method = 1, 
#' cor.method = "pearson", hop.size = 11.6, wl = NULL, ovlp = 90, wn = 'hanning')
#' @param X object of class 'extended_selection_table' created by the function \code{\link[warbleR]{selection_table}} from the warbleR package. The object must include the following additional columns: 'signal.type', 'bottom.freq' and 'top.freq'.
#' @param parallel Numeric vector of length 1. Controls whether parallel computing is applied by specifying the number of cores to be used. Default is 1 (i.e. no parallel computing).
#' @param pb Logical argument to control if progress bar is shown. Default is \code{TRUE}.
#' @param method Numeric vector of length 1 to indicate the 'experimental design' for measuring envelope correlation. Two methods are available:
#' \itemize{
#' \item \code{1}: compare all signals with their counterpart that was recorded at the closest distance to source (e.g. compare a signal recorded at 5m, 10m and 15m with its counterpart recorded at 1m). This is the default method. 
#' \item \code{2}: compare all signals with their counterpart recorded at the distance immediately before (e.g. a signal recorded at 10m compared with the same signal recorded at 5m, then signal recorded at 15m compared with same signal recorded at 10m and so on).
#' }
#' @param cor.method Character string indicating the correlation coefficient to be applied ("pearson", "spearman", or "kendall", see \code{\link[stats]{cor}}).
#' @param hop.size A numeric vector of length 1 specifying the time window duration (in ms). Default is 11.6 ms, which is equivalent to 512 wl for a 44.1 kHz sampling rate. Ignored if 'wl' is supplied.
#' @param wl A numeric vector of length 1 specifying the window length of the spectrogram, default 
#' is NULL. If supplied, 'hop.size' is ignored.
#' @param ovlp Numeric vector of length 1 specifying \% of overlap between two 
#' consecutive windows, as in \code{\link[seewave]{spectro}}. Default is 90. High values of ovlp 
#' slow down the function but produce more accurate results.
#' @param wn A character vector of length 1 specifying the window name as in \code{\link[seewave]{ftwindow}}. 
#' @return Extended selection table similar to input data, but includes a new column (cross.correlation)
#' with the spectrogram cross-correlation coefficients.
#' @export
#' @name spcc
#' @details Spectrographic cross-correlation measures frequency distortion of signals as a similarity metric. Values close to 1 means very similar spectrograms (i.e. little signal distortion has occurred). Cross-correlation is measured of signals in which a reference playback has been re-recorded at increasing distances. The 'signal.type' column must be used to indicate the function to only compare signals belonging to the same category (e.g. song-types). The function compares each signal type to the corresponding reference signal within the supplied frequency range (e.g. bandpass) of the reference signal ('bottom.freq' and 'top.freq' columns in 'X'). Two methods for calculating cross-correlation are provided (see 'method' argument). The function is a wrapper on warbleR's \code{\link[warbleR]{xcorr}} function.
#' @examples
#' {
#' # load example data
#' data("playback_est")
#' 
#' # method 1
#'spcc(X = playback_est, method = 1)
#' 
#' # method 2
#' spcc(X = playback_est, method = 2)
#' }
#' 
#' @author Marcelo Araya-Salas (\email{marceloa27@@gmail.com})
#' @seealso \code{\link{blur_ratio}}, \code{\link{spcc_align}}, \code{\link[warbleR]{xcorr}}
#' @references {
#' Araya-Salas, M. (2020). baRulho: baRulho: quantifying habitat-induced degradation of (animal) acoustic signals in R. R package version 1.0.0
#' 
#' Clark, C.W., Marler, P. & Beeman K. (1987). Quantitative analysis of animal vocal phonology: an application to Swamp Sparrow song. Ethology. 76:101-115. 
#' }
# last modification on jan-06-2020 (MAS)

spcc <- function(X, parallel = 1, pb = TRUE,  method = 1, cor.method = "pearson", hop.size = 11.6, wl = NULL, ovlp = 90, wn = 'hanning'){
  
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

  # If method is not numeric
  if (!is.character(cor.method)) stop("'cor.method' must be a character vector of length 1")
  
  if (!any(cor.method %in%  c("pearson", "kendall", "spearman"))) stop("'method' must be either  'pearson', 'kendall' or 'spearman'")
  
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
    if (method == 1) # comparing to closest distance to source
      cmp.mt <- cbind(Y$sf.selec[which.min(Y$distance)], Y$sf.selec[-which.min(Y$distance)]) else # comparing to previous distance
        cmp.mt <- cbind(Y$sf.selec[-nrow(Y)], Y$sf.selec[-1])
    
    return(cmp.mt)
  })
  
  # put together in a single
  comp_mat <- do.call(rbind, comp_mats)
  
  # run spcc 
  xcorrs <- warbleR::xcorr(X = X, wl = wl, ovlp = ovlp, wn = wn, cor.method = "pearson", parallel = parallel, pb = pb, compare.matrix = comp_mat)
  
  # put results back into X
  X$reference <- NA
  X$cross.correlation <- NA
  
  # add correlation and reference only for calculated correlations
  X$reference[match(xcorrs$X2, paste(X$sound.files, X$selec, sep = "-"))] <- as.character(xcorrs$X1)
  
  X$cross.correlation[match(xcorrs$X2, paste(X$sound.files, X$selec, sep = "-"))] <- xcorrs$score
  
  return(X)
}
