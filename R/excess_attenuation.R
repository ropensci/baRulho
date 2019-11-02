#' Measure excess attenuation
#' 
#' \code{excess_attenuation} measures excess attenuation in signals referenced in a extended selection table.
#' @usage excess_attenuation(X, parallel = 1, pb = TRUE, method = 1, 
#' bp = NULL, wl = 10, output = "est")
#' @param X object of class 'selection_table', 'extended_selection_table' created by the function \code{\link[warbleR]{selection_table}} from the warbleR package.
#' @param parallel Numeric. Controls whether parallel computing is applied.
#' It specifies the number of cores to be used. Default is 1 (i.e. no parallel computing).
#' If \code{NULL} (default) then the current working directory is used.
#' @param pb Logical argument to control if progress bar is shown. Default is \code{TRUE}. It can also be
#' set globally using the 'pb' option (see \code{\link{warbleR_options}}).
#' @param method Numeric vector of length 1 to indicate the 'experimental design' for measuring excess attenuation. Two methods are available:
#' \itemize{
#' \item \code{1}: compare all signals with their counterpart that was recorded at the closest distance to source (e.g. compare a signal recorded at 5m, 10m and 15m with its counterpart recorded at 1m). This is the default method. 
#' \item \code{2}: compare all signals with their counterpart recorded at the distance immediately before (e.g. a signal recorded at 10m compared with the same signal recorded at 5m, then signal recorded at 15m compared with same signal recorded at 10m and so on).
#' }
#' @param bp Numeric vector of length 2 giving the lower and upper limits of a frequency bandpass filter (in kHz). Default is \code{NULL}.
#' @param wl A numeric vector of length 1 specifying the window length of the spectrogram for applying bandpass. Default 
#'   is 10. Ignored if \code{bp = NULL}.
#'  Note that lower values will increase time resolution, which is more important for amplitude ratio calculations. 
#' @param output Character vector of length 1 to determine if an extended selection table ('est', default) or a data frame ('data.frame') is returned.
#' @return Data frame or extended selection table (depending on 'output' argument) similar to input data, but also includes a new column (excess.attenuation)
#' with the excess attenuation values.
#' @export
#' @name excess_attenuation
#' @details Excess attenuation is the attenuation of a sound in excess of that due to spherical spreading as described by Dabelsteen et al (1993). The goal of the function is to measure the excess attenuation on signals in which a master playback has been re-recorded at different distances. The 'signal.id' column must be used to indicate which signals belonging to the same category (e.g. song-types). The function will then compared each signal type to its reference. Two methods for calculating excess attenuation are provided (see 'method' argument).   
#' @examples
#' {
#' # load example data
#' data("playback_est")
#' 
#' # remove noise selections
#' playback_est <- playback_est[playback_est$signal.id != "noise", ]
#'
#' # using method 1
#'excess_attenuation(X = playback_est)
#' 
#' # using method 2
#' excess_attenuation(X = playback_est, method = 2)
#' }
#' 
#' @author Marcelo Araya-Salas (\email{marceloa27@@gmail.com}) #' @references {
#' Dabelsteen, T., Larsen, O. N., & Pedersen, S. B. (1993). Habitat-induced degradation of sound signals: Quantifying the effects of communication sounds and bird location on blur ratio, excess attenuation, and signal-to-noise ratio in blackbird song. The Journal of the Acoustical Society of America, 93(4), 2206.
#' 
#' Araya-Salas, M. (2019), baRulho: a R package to evaluate habitat-induced degradation of (animal) acoustic signals. R package version 1.0.0
#' }
#last modification on nov-01-2019 (MAS)

excess_attenuation <- function(X, parallel = 1, pb = TRUE, method = 1, 
                       bp = NULL, wl = 10, output = "est"){
  
  # set pb options 
  on.exit(pbapply::pboptions(type = .Options$pboptions$type), add = TRUE)
  
  # If parallel is not numeric
  if (!is.numeric(parallel)) stop("'parallel' must be a numeric vector of length 1") 
  if (any(!(parallel %% 1 == 0),parallel < 1)) stop("'parallel' should be a positive integer")
  
  # If method is not numeric
  if (!is.numeric(method)) stop("'method' must be a numeric vector of length 1") 
  if (!any(method %in% 1:2)) stop("'method' must be either 1 or 2")
  
  # check signal.id column 
  if (is.null(X$signal.id)) stop("'X' must containe a 'signal.id' column")
  
  #check output
  if (!any(output %in% c("est", "data.frame"))) stop("'output' must be either 'est' or 'data.frame'")  
  
  # function to measure RMS for signal and noise
  rms_FUN <- function(y, bp, wl){
    
    # read signal clip
    signal <- warbleR::read_wave(X = X, index = y)
    
    # get RMS for signal
    sigRMS <- seewave::rms(seewave::env(signal, f = signal@samp.rate, envt = "abs", plot = FALSE))
    sigRMS <- 20*log10(sigRMS)
    
    return(data.frame(X[y, , drop = FALSE], sigRMS))
  }
  
  # set pb options 
  pbapply::pboptions(type = ifelse(as.logical(pb), "timer", "none"))
  
  # set clusters for windows OS
  if (Sys.info()[1] == "Windows" & parallel > 1)
    cl <- parallel::makePSOCKcluster(getOption("cl.cores", parallel)) else cl <- parallel
  
  # run loop apply function
  RMS <- pbapply::pblapply(X = 1:nrow(X), cl = cl, FUN = function(y)  rms_FUN(y, bp, wl)) 
  
  # put in a data frame
  RMS_df <- do.call(rbind, RMS)
  
  # split by signal ID
  RMS_list <- split(RMS_df, RMS_df$signal.id)
  
  # calculate excess attenuation
  X$excess.attenuation <- unlist(sapply(RMS_list, function(Y, meth = method){
    
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
      EA <- (term1 + term2) / term3
      EA[which.min(Y$distance)] <- NA
    }
    
    # compare to previous distance 
    if (meth == 2){
      
      # save original order
      Y$org.ord <- 1:nrow(Y)
      
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
      EA <- c(NA, EA)      
      
      # reorder results
      EA <- EA[order(Y$org.ord)]
    }
    
    return(EA)
    
  }, simplify = FALSE))
  
  # convert to data frame instead of extended selection table
  if (output == "data.frame") 
    X <- as.data.frame(X)
  
  return(X)
}
