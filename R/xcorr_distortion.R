#' Measure spectrographic cross-correlation as a measure of signal distortion
#' 
#' \code{xcorr_distortion} Measures spectrographic cross-correlation as a measure of signal distortion in signals referenced in an extended selection table.
#' @usage xcorr_distortion(X = NULL, parallel = 1, pb = TRUE,  method = 1, cor.method = "pearson",
#' wl = 512, ovlp = 90, wn = 'hanning', output = "est")
#' @param X object of class 'selection_table', 'extended_selection_table' created by the function \code{\link[warbleR]{selection_table}} from the warbleR package.
#' @param parallel Numeric vector of length 1. Controls whether parallel computing is applied by specifying the number of cores to be used. Default is 1 (i.e. no parallel computing).
#' @param pb Logical argument to control if progress bar is shown. Default is \code{TRUE}.
#' @param method Numeric vector of length 1 to indicate the 'experimental design' for measuring envelope correlation. Two methods are available:
#' \itemize{
#' \item \code{1}: compare all signals with their counterpart that was recorded at the closest distance to source (e.g. compare a signal recorded at 5m, 10m and 15m with its counterpart recorded at 1m). This is the default method. 
#' \item \code{2}: compare all signals with their counterpart recorded at the distance immediately before (e.g. a signal recorded at 10m compared with the same signal recorded at 5m, then signal recorded at 15m compared with same signal recorded at 10m and so on).
#' }
#' @param cor.method Character string indicating the correlation coefficient to be applied ("pearson", "spearman", or "kendall", see \code{\link[stats]{cor}}).
#' @param wl A numeric vector of length 1 specifying the window length of the spectrogram, default 
#' is 512.
#' @param ovlp Numeric vector of length 1 specifying \% of overlap between two 
#' consecutive windows, as in \code{\link[seewave]{spectro}}. Default is 90. High values of ovlp 
#' slow down the function but produce more accurate results.
#' @param wn A character vector of length 1 specifying the window name as in \code{\link[seewave]{ftwindow}}. 
#' @param output Character vector of length 1 to determine if an extended selection table ('est', default) or a data frame ('data.frame') is returned.
#' @return Data frame or extended selection table (depending on 'output' argument) similar to input data, but includes a new column (cross.correlation)
#' with the spectrogram cross-correlation coefficients.
#' @export
#' @name xcorr_distortion
#' @details The spectrographic cross-correlation measures  frequency distortion of signals as a similarity metric where values range from 1 (completely equal, no distortion) and decays towards 0 (highly distorted). Cross-correlation is measured of signals in which a reference playback has been re-recorded at increasing distances. The 'signal.id' column must be used to indicate the function to only compare signals belonging to the same category (e.g. song-types). The function will then compare each signal type to the corresponding reference signal. Two methods for calculating cross-correlation are provided (see 'method' argument). The function is a wrapper on warbleR's \code{\link[warbleR]{xcorr}} function.
#' @examples
#' {
#' # load example data
#' data("playback_est")
#' 
#' # remove noise selections
#' playback_est <- playback_est[playback_est$signal.id != "noise", ]
#' 
#' # method 1
#'xcorr_distortion(X = playback_est, method = 1)
#' 
#' # method 2
#' xcorr_distortion(X = playback_est, method = 2)
#' }
#' 
#' @author Marcelo Araya-Salas (\email{marceloa27@@gmail.com})
#' @references {
#' Araya-Salas, M. (2019), baRulho: a R package to quantify habitat-induced degradation of (animal) acoustic signals. R package version 1.0.0
#' 
#' Clark, C.W., Marler, P. & Beeman K. (1987). Quantitative analysis of animal vocal phonology: an application to Swamp Sparrow song. Ethology. 76:101-115. 
#' }
# last modification on nov-01-2019 (MAS)

xcorr_distortion <- function(X = NULL, parallel = 1, pb = TRUE,  method = 1, cor.method = "pearson", wl = 512, ovlp = 90, wn = 'hanning', output = "est"){
  
  # If parallel is not numeric
  if (!is.numeric(parallel)) stop("'parallel' must be a numeric vector of length 1") 
  if (any(!(parallel %% 1 == 0),parallel < 1)) stop("'parallel' should be a positive integer")
  
  # If method is not numeric
  if (!is.character(cor.method)) stop("'cor.method' must be a character vector of length 1") 
  if (!any(cor.method %in%  c("pearson", "kendall", "spearman"))) stop("'method' must be either  'pearson', 'kendall' or 'spearman'")
  
  # check signal.id column 
  if (is.null(X$signal.id)) stop("'X' must containe a 'signal.id' column")
  
  #check output
  if (!any(output %in% c("est", "data.frame"))) stop("'output' must be either 'est' or 'data.frame'")  
  
  # create matrix containing pairwise comparisons of selections (2 columns)
  comp_mats <- lapply(unique(X$signal.id), function(x){
    
    # extract for single signal and order by distance
    Y <- as.data.frame(X[X$signal.id == x, ])
   
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
  
  # run xcorr 
  xcorrs <- warbleR::xcorr(X = X, wl = wl, ovlp = ovlp, wn = wn, cor.method = cor.method, parallel = parallel, pb = pb, compare.matrix = comp_mat)
  
  # put results back into X
  X$cross.correlation <- NA
  
  X$cross.correlation[match(xcorrs$X2, paste(X$sound.files, X$selec, sep = "-"))] <- xcorrs$score
  
  # convert to data frame instead of extended selection table
  if (output == "data.frame") 
    X <- as.data.frame(X)
  
  return(X)
}
