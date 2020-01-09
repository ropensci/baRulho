#' Measure excess attenuation
#' 
#' \code{excess_attenuation} measures excess attenuation in signals referenced in an extended selection table.
#' @usage excess_attenuation(X, parallel = 1, pb = TRUE, method = 1, 
#' bp = NULL, hop.size = 1, wl = NULL, ovlp = 70)
#' @param X object of class 'selection_table', 'extended_selection_table' created by the function \code{\link[warbleR]{selection_table}} from the warbleR package. The data frame must include the following additional columns: 'distance', 'signal.type', 'bottom.freq' and 'top.freq'.
#' @param parallel Numeric vector of length 1. Controls whether parallel computing is applied by specifying the number of cores to be used. Default is 1 (i.e. no parallel computing).
#' @param pb Logical argument to control if progress bar is shown. Default is \code{TRUE}.
#' @param method Numeric vector of length 1 to indicate the 'experimental design' for measuring excess attenuation. Two methods are available:
#' \itemize{
#' \item \code{1}: compare all signals with their counterpart that was recorded at the closest distance to source (e.g. compare a signal recorded at 5m, 10m and 15m with its counterpart recorded at 1m). This is the default method. 
#' \item \code{2}: compare all signals with their counterpart recorded at the distance immediately before (e.g. a signal recorded at 10m compared with the same signal recorded at 5m, then signal recorded at 15m compared with same signal recorded at 10m and so on).
#' }
#' @param bp Numeric vector of length 2 giving the lower and upper limits of a frequency bandpass filter (in kHz). Default is \code{NULL}.
#' @param hop.size A numeric vector of length 1 specifying the time window duration (in ms). Default is 1 ms, which is equivalent to ~45 wl for a 44.1 kHz sampling rate. Ignored if 'wl' is supplied.
#' @param wl A numeric vector of length 1 specifying the window length of the spectrogram, default 
#' is NULL. Ignored if \code{bp = NULL}. If supplied, 'hop.size' is ignored.
#' Note that lower values will increase time resolution, which is more important for amplitude ratio calculations. 
#' @param ovlp Numeric vector of length 1 specifying the percent overlap between two 
#'   consecutive windows, as in \code{\link[seewave]{spectro}}. Only used when plotting. Default is 70. Only used for bandpass filtering.
#' @return Extended selection table similar to input data, but also includes a new column (excess.attenuation)
#' with the excess attenuation values.
#' @export
#' @name excess_attenuation
#' @details Excess attenuation is the amplitude loss of a sound in excess of that due to spherical spreading. With every doubling of distance, sounds attenuate with a 6 dB loss of amplitude (Morton, 1975; Marten & Marler, 1977). Any additional loss of amplitude results in excess attenuation, or energy loss in excess of that expected to occur with distance via spherical spreading due to atmospheric conditions or habitat (Wiley & Richards, 1978). Low values indicate little signal attenuation. 
#' The goal of the function is to measure the excess attenuation on signals in which a reference playback has been re-recorded at increasing distances. The 'signal.type' column must be used to indicate which signals belonging to the same category (e.g. song-types). The function will then compare each signal type to the corresponding reference signal the supplied frequency range (e.g. bandpass) of the reference signal ('bottom.freq' and 'top.freq' columns in 'X'). Two methods for calculating excess attenuation are provided (see 'method' argument). 
#' @examples
#' {
#' # load example data
#' data("playback_est")
#' 
#' # using method 1
#'excess_attenuation(X = playback_est)
#' 
#' # using method 2
#' excess_attenuation(X = playback_est, method = 2)
#' }
#' 
#' @author Marcelo Araya-Salas (\email{marceloa27@@gmail.com})
#' @seealso \code{\link{spcc}}
#' @references {
#' Araya-Salas, M. (2020), baRulho: baRulho: quantifying habitat-induced degradation of (animal) acoustic signals in R. R package version 1.0.0
#' 
#' Marten, K., & Marler, P. (1977). Sound transmission and its significance for animal vocalization. Behavioral Ecology and Sociobiology, 2(3), 271-290.
#' 
#' Morton, E. S. (1975). Ecological sources of selection on avian sounds. The American Naturalist, 109(965), 17-34.
#' }
#last modification on nov-01-2019 (MAS)

excess_attenuation <- function(X, parallel = 1, pb = TRUE, method = 1, 
                               bp = NULL, hop.size = 1, wl = NULL, ovlp = 70){
  
  # set pb options 
  on.exit(pbapply::pboptions(type = .Options$pboptions$type), add = TRUE)
  
  # is extended sel tab
  if (!warbleR::is_extended_selection_table(X)) 
    stop("'X' must be and extended selection table")
  
  # If parallel is not numeric
  if (!is.numeric(parallel)) stop("'parallel' must be a numeric vector of length 1") 
  if (any(!(parallel %% 1 == 0),parallel < 1)) stop("'parallel' should be a positive integer")
  
  # hopsize  
  if (!is.numeric(hop.size) | hop.size < 0) stop("'parallel' must be a positive number") 
  
  # adjust wl based on hope.size
  if (is.null(wl))
    wl <- round(attr(X, "check.results")$sample.rate[1] * hop.size, 0)
  
  # If method is not numeric
  if (!is.numeric(method)) stop("'method' must be a numeric vector of length 1") 
  if (!any(method %in% 1:2)) stop("'method' must be either 1 or 2")
  
  # check signal.type column 
  if (is.null(X$signal.type)) stop("'X' must containe a 'signal.type' column")
  
  # add sound file selec column and names to X (weird column name so it does not overwrite user columns)
  X <- prep_X_bRlo_int(X, method = method)
  
  # function to measure RMS for signal and noise
  rms_FUN <- function(y, bp, wl, ovlp){
    
    # read signal clip
    clp <- warbleR::read_wave(X = X, index = y)
    
    # define bandpass based on reference
    bp <- c(X$bottom.freq[X$TEMP....sgnl == X$reference[y]], X$top.freq[X$TEMP....sgnl == X$reference[y]])
    
    # bandpass filter
    clp <- seewave::ffilter(clp, from = bp[1] * 1000, 
                            ovlp = ovlp, to = bp[2] * 1000, bandpass = TRUE, 
                            wl = wl, output = "Wave")
    
    # get RMS for signal
    sigRMS <- seewave::rms(seewave::env(clp, f = clp@samp.rate, envt = "abs", plot = FALSE))
    sigRMS <- 20 * log10(sigRMS)
    
    return(data.frame(X[y, , drop = FALSE], sigRMS))
  }
  
  # set pb options 
  pbapply::pboptions(type = ifelse(as.logical(pb), "timer", "none"))
  
  # set clusters for windows OS
  if (Sys.info()[1] == "Windows" & parallel > 1)
    cl <- parallel::makePSOCKcluster(getOption("cl.cores", parallel)) else cl <- parallel
  
  # run loop apply function
  RMS <- pbapply::pblapply(X = 1:nrow(X), cl = cl, FUN = function(y)  rms_FUN(y, bp, wl, ovlp)) 
  
  # put in a data frame
  RMS_df <- do.call(rbind, RMS)
  
  # split by signal ID
  RMS_list <- split(RMS_df, RMS_df$signal.type)
  
  # calculate excess attenuation
  X_list <- sapply(RMS_list, function(Y, meth = method){
    
    if (Y$signal.type[1] == "ambient") Y$excess.attenuation <- NA else {
    
    # method 1 compare to closest distance to source
    if (meth == 1){
      
      # extract RMS of signal and background references
      sig_RMS_REF <- Y$sigRMS[which.min(Y$distance)]
      dist_REF <- Y$distance[which.min(Y$distance)]
      
      # term 1: decrease in signal amplitude (RMS) of reference (Ref) vs re-recorded (RR)
      term1 <- sig_RMS_REF - Y$sigRMS 
      
      # lost due to spheric spreading
      term2 <- -20 * log10(1 / Y$distance)
      
      # distance traveled by sound
      term3 <- Y$distance - dist_REF
      
      # excess attenuation = (total attenuation - spheric spreading attenuation) / distance
      Y$excess.attenuation <- (term1 + term2) / term3
      Y$excess.attenuation[which.min(Y$distance)] <- NA
      }
    
    # compare to previous distance 
    if (meth == 2){
      
      # save original order
      Y$org....ord <- 1:nrow(Y)
      
      # sort by distance
      Y <- Y[order(Y$distance), ]
      
      # term 1: decrease in signal amplitude (RMS) of reference (Ref) vs re-recorded (RR)
      term1 <- Y$sigRMS[-nrow(Y)] - Y$sigRMS[-1] 
      
      # lost due to spheric spreading
      term2 <- -20 * log10(1 / (Y$distance[-1]))
      
      # distance traveled by sound
      term3 <- Y$distance[-1] - Y$distance[-nrow(Y)]
      
      # excess attenuation = (total attenuation - spheric spreading attenuation) / distance
      EA <- (term1 + term2) / term3
      
      # add NA for first distance
      Y$excess.attenuation <- c(NA, EA)      
      
      # reorder results
      Y <- Y[order(Y$org....ord), ]
      
      Y$org....ord <- NULL
      }
    }
    
    return(Y)
    
  }, simplify = FALSE)
  
  # put together in a data frame as X
  X2 <- do.call(rbind, X_list)
  
  # fix row names 
  rownames(X2) <-  rownames(X)
  
  # remove temporal column
  X2$sigRMS <- X2$TEMP....sgnl <- NULL
  
  # fix est
  X2 <- warbleR::fix_extended_selection_table(X = X2, Y = X)
  
  return(X2)
}
