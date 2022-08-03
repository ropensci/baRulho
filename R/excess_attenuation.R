#' Measure excess attenuation
#' 
#' \code{excess_attenuation} measures excess attenuation in signals referenced in an extended selection table.
#' @usage excess_attenuation(X, parallel = 1, pb = TRUE, method = 1, type = "Darden",
#'  output = "est", hop.size = 1, wl = NULL, ovlp = 70, gain = 0)
#' @param X object of class 'extended_selection_table' created by the function \code{\link[warbleR]{selection_table}} from the warbleR package. The data frame must include the following additional columns: 'distance', 'signal.type', 'bottom.freq' and 'top.freq'.
#' @param parallel Numeric vector of length 1. Controls whether parallel computing is applied by specifying the number of cores to be used. Default is 1 (i.e. no parallel computing).
#' @param pb Logical argument to control if progress bar is shown. Default is \code{TRUE}.
#' @param method Numeric vector of length 1 to indicate the 'experimental design' for measuring excess attenuation. Two methods are available:
#' \itemize{
#' \item \code{1}: compare all signals with their counterpart that was recorded at the closest distance to source (e.g. compare a signal recorded at 5m, 10m and 15m with its counterpart recorded at 1m). This is the default method. 
#' \item \code{2}: compare all signals with their counterpart recorded at the distance immediately before (e.g. a signal recorded at 10m compared with the same signal recorded at 5m, then signal recorded at 15m compared with same signal recorded at 10m and so on).
#' }
#' @param type Character vector of length 1 to indicate the 'type' of excess attenuation to be used. Two types are available:
#' \itemize{
#' \item \code{Darden}: as described by Darden et al. (2008): microphone_gain - 20 x log(reference distance / re-recorded distance) - 20 x log(envelope_correlation). The function \code{\link{envelope_correlation}} is used internally. Microphone gain is the microphone gain of the reference and re-recorded signals. This is the default method. If gain is not supplied (see 'gain' argument) gain is set as 0, which results in a relative measure of excess attenuation comparable only within the same experiment or between experiments with the same recording equipment and recording volume.
#' \item \code{Marten}: as described by Marten et al. (1977): Observed attenuation is calculated as the difference in sound pressure level between the reference and the re-recorded signal. Sound pressure level is measured as the root mean square of the amplitude envelope and as such it is not calibrated. Background noise amplitude (measured on the noise right before the signals) is subtracted from the sound pressure level estimates.
#' }
#' @param output Character vector of length 1 to determine if an extended selection table ('est', default) or a data frame ('data.frame').
#' @param hop.size A numeric vector of length 1 specifying the time window duration (in ms). Default is 1 ms, which is equivalent to ~45 wl for a 44.1 kHz sampling rate. Ignored if 'wl' is supplied.
#' @param wl A numeric vector of length 1 specifying the window length of the spectrogram, default 
#' is \code{NULL}. If supplied, 'hop.size' is ignored.
#' Note that lower values will increase time resolution, which is more important for amplitude ratio calculations. 
#' @param ovlp Numeric vector of length 1 specifying the percent overlap between two 
#'   consecutive windows, as in \code{\link[seewave]{spectro}}. Only used when plotting. Default is 70. Only used for bandpass filtering.
#' @param gain Numeric vector of length 1 with the combined gain of the microphone and recorder (in dB). Default is 0, which results in a relative measure of excess attenuation comparable only within the same experiment, but not across experiments. Only used for \code{type = "Marten"}.  
#' @return Extended selection table similar to input data, but also includes a new column (excess.attenuation)
#' with the excess attenuation values.
#' @export
#' @name excess_attenuation
#' @details Excess attenuation is the amplitude loss of a sound in excess due to spherical spreading (observed attenuation - expected attenuation). With every doubling of distance, sounds attenuate with a 6 dB loss of amplitude (Morton, 1975; Marten & Marler, 1977). Any additional loss of amplitude results in energy loss in excess of that expected to occur with distance via spherical spreading. So it represents energy loss due to additional factors like vegetation or atmospheric conditions (Wiley & Richards, 1978). Low values indicate little additional attenuation. 
#' The goal of the function is to measure the excess attenuation on signals in which a reference playback has been re-recorded at increasing distances. The 'signal.type' column must be used to indicate which signals belonging to the same category (e.g. song-types). The function will then compare each signal type to the corresponding reference signal. Two methods for calculating excess attenuation are provided (see 'method' argument). 
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
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})
#' @seealso \code{\link{spcc_distortion}}; \code{\link{envelope_correlation}}
#' @references {
#' Araya-Salas, M. (2020). baRulho: baRulho: quantifying habitat-induced degradation of (animal) acoustic signals in R. R package version 1.0.2
#' 
#' Dabelsteen, T., Larsen, O. N., & Pedersen, S. B. (1993). Habitat-induced degradation of sound signals: Quantifying the effects of communication sounds and bird location on blur ratio, excess attenuation, and signal-to-noise ratio in blackbird song. The Journal of the Acoustical Society of America, 93(4), 2206.
#' 
#' Darden, SK, Pedersen SB, Larsen ON, & Dabelsteen T. (2008). Sound transmission at ground level in a short-grass prairie habitat and its implications for long-range communication in the swift fox *Vulpes velox*. The Journal of the Acoustical Society of America, 124(2), 758-766.
#' 
#' Marten K, & Marler P. (1977). Sound transmission and its significance for animal vocalization. Behavioral Ecology and Sociobiology, 2(3), 271-290.
#' 
#' Morton ES. (1975). Ecological sources of selection on avian sounds. The American Naturalist, 109(965), 17-34.
#' 
#' Wiley, R., & Richards, D. G. (1978). Physical constraints on acoustic communication in the atmosphere: implications for the evolution of animal vocalizations. Behavioral Ecology and Sociobiology, 3(1), 69-94.
#' }
#last modification on jul-19-2021 (MAS)

excess_attenuation <- function(X, parallel = 1, pb = TRUE, method = 1, type = "Darden", 
                               output = "est", hop.size = 1, wl = NULL, ovlp = 70, gain = 0){
  
  # is extended sel tab
  if (!warbleR::is_extended_selection_table(X)) 
    stop("'X' must be and extended selection table")
  
  # If parallel is not numeric
  if (!is.numeric(parallel)) stop("'parallel' must be a numeric vector of length 1") 
  if (any(!(parallel %% 1 == 0),parallel < 1)) stop("'parallel' should be a positive integer")
  
  #check output
  if (!any(output %in% c("est", "data.frame"))) stop("'output' must be 'est' or 'data.frame'")  
  
  # hopsize  
  if (!is.numeric(hop.size) | hop.size < 0) stop("'hop.size' must be a positive number") 
  
  # adjust wl based on hope.size
  if (is.null(wl))
    wl <- round(attr(X, "check.results")$sample.rate[1] * hop.size, 0)
  
  # make wl even if odd
  if (!(wl %% 2) == 0) wl <- wl + 1
  
  # If method is not numeric
  if (!is.numeric(method)) stop("'method' must be a numeric vector of length 1") 
  if (!any(method %in% 1:2)) stop("'method' must be either 1 or 2")
  
  # check signal.type column 
  if (is.null(X$signal.type)) stop("'X' must containe a 'signal.type' column")
  
  # add sound file selec column and names to X (weird column name so it does not overwrite user columns)
  if (pb) 
    if (type == "Marten")
    write(file = "", x = paste0("Preparing data for analysis (step 1 out of 3):")) else
      write(file = "", x = paste0("Preparing data for analysis (step 1 out of 2):"))
  
  X <- prep_X_bRlo_int(X, method = method, parallel = parallel, pb = pb)
  
  # # function to measure RMS for signal and noise
  spl_FUN <- function(y, wl, ovlp){
    
    # read signal clip
    clp <- warbleR::read_wave(X = X, index = y, from = 0, to = X$end[y])
    
    # get RMS for signal
    sigRMS <- seewave::rms(seewave::env(clp, f = clp@samp.rate, envt = "abs", plot = FALSE))
    sigSPL <- 20 * log10(sigRMS)

    # measure ambient SPL
    # read noise before signal
    if (X$signal.type[y] != "ambient"){
      
      noise_clp <- warbleR::read_wave(X = X, index = y, from = 0, to = X$start[y]- 0.001)
      noiseRMS <- seewave::rms(seewave::env(noise_clp, f = noise_clp@samp.rate, envt = "abs", plot = FALSE))
      noiseSPL <- 20 * log10(noiseRMS)
      
      # remove noise SPL from signal SPL
      sigSPL <- lessdB(signal.noise = sigSPL, noise = noiseSPL)
      } 
    return(data.frame((X[y, , drop = FALSE]), sigSPL))
  }
  
  # set clusters for windows OS
  if (Sys.info()[1] == "Windows" & parallel > 1)
    cl <- parallel::makePSOCKcluster(getOption("cl.cores", parallel)) else cl <- parallel
  
  if (type == "Marten"){
  if (pb) 
    write(file = "", x = paste0("Measuring sound pressure level (step 2 out of 3):"))
  
  # run loop apply function
  SPLs <- warbleR:::pblapply_wrblr_int(X = 1:nrow(X), pbar = pb, cl = cl, FUN = function(y)  spl_FUN(y, wl, ovlp))
  
  # put in a data frame
  X2 <- do.call(rbind, SPLs)
}  else 
  X2 <- X
  
  # split by signal ID
  sigtype_list <- split(X2, X2$signal.type)
  
  if (pb)
    if (type == "Marten")
    write(file = "", x = paste0("Calculating excess attenuation (step 3 out of 3):")) else
      write(file = "", x = paste0("Calculating excess attenuation (step 2 out of 2):"))
  
  # calculate excess attenuation
  X_list <- warbleR:::pblapply_wrblr_int(X = sigtype_list, pbar = pb, cl = cl, function(Y, meth = method, tp = type){
    
    if (Y$signal.type[1] == "ambient") Y$excess.attenuation <- NA else {
      
      # method 1 compare to closest distance to source
      if (meth == 1){
        
        # extract SPL of signal and background references
        sig_SPL_REF <- Y$sigSPL[which.min(Y$distance)]
        dist_REF <- Y$distance[which.min(Y$distance)]
        
        # type Marten
        if (tp == "Marten"){
          # ea <- observed_attenuation - expected_attenuation
          
          # term 1: decrease in signal amplitude (RMS) of reference (Ref) vs re-recorded (RR)
          term1 <- sig_SPL_REF - Y$sigSPL
          
          # lost due to spherical spreading
          term2 <- 20 * log10(Y$distance / dist_REF)
          
          ea <- (term2 - term1)
  
          # segun chirras
          ## EA = -20logK - 6dB / dd + dB anadidos
           
        } 
        
        if (tp == "Darden"){
          
          #EA = g - 20 log(d / 10) - 20 log(k)
          # term1 = g (combined mic gain)
          term1 <- gain
          
          # term2 = -20 log(d / 10)
          # term2 <-  -20 * log(Y$distance / 10)
          term2 <- -20 * log10(dist_REF / Y$distance)
          
          # term3 = - 20 log(k)
          # get envelope correlation (k)
          k <- envelope_correlation(X[X$signal.type == Y$signal.type[1],], output = "data.frame", pb = FALSE)$envelope.correlation
          term3 <- -20 * log(k)
          
          # excess attenuation = (total attenuation - spheric spreading attenuation) / distance
          ea <- term1 + term2 + term3
        } 
        
        Y$excess.attenuation <- ea
        Y$excess.attenuation[which.min(Y$distance)] <- NA
      }
      
      # compare to previous distance 
      if (meth == 2){
        
        # save original order
        Y$org....ord <- 1:nrow(Y)
        
        # sort by distance
        Y <- Y[order(Y$distance), ]
        
        
        if (tp == "Marten"){
          # term 1: decrease in signal amplitude (RMS) of reference (Ref) vs re-recorded (RR)
          term1 <- Y$sigSPL[-nrow(Y)] - Y$sigSPL[-1] 
          
          # lost due to spheric spreading
          term2 <- 20 * log10(Y$distance[-1] / Y$distance[-nrow(Y)])
          
          # ea <- observed_attenuation - expected_attenuation
          ea <- term2 - term1
          
          # add NA for first distance
          ea <- c(NA, ea)      
        }
        
        # type Darden
        if (tp == "Darden"){
          
          #EA = g - 20 log(d / 10) - 20 log(k)
          # term1 = g (combined mic gain)
          term1 <- gain
          
          # term2 = -20 log(ref d / d)
          term2 <- -20 * log10(Y$distance[-nrow(Y)] / Y$distance[-1])

          # term3 = - 20 log(k)
          # get envelope correlation (k)
          k <- envelope_correlation(X[X$signal.type == Y$signal.type[1],], output = "data.frame", pb = FALSE, method = 2)$envelope.correlation
          
          # order by distance too
          k <- k[order(X$distance[X$signal.type == Y$signal.type[1]])]
          term3 <- -20 * log10(k[-1])
          
          # excess attenuation
          ea <- term1 + term2 + term3
          
          # add NA for first distance
          ea <- c(NA, ea)  
        } 
        
        Y$excess.attenuation <- ea
        # reorder results
        Y <- Y[order(Y$org....ord), ]
        
        Y$org....ord <- NULL
      }
    }
    
    Y <- as.data.frame(Y)
    return(Y)
    
  })
  
  # put together in a data frame as X
  X2 <- do.call(rbind, X_list)
  
  # fix row names 
  rownames(X2) <-  rownames(X)
  
  # remove temporal column
  X2$sigRMS <- X2$TEMP....sgnl <- NULL
  
  # fix est
  if (output == "est")
    X2 <- warbleR::fix_extended_selection_table(X = X2, Y = X)
  
  return(X2)
}

## internal function to subtract SPL from background noise
# signal = signal SPL
# noise = noise SPL
lessdB <- function(signal.noise, noise){
  
  puttative_SPLs <- seq(0.01, signal.noise, by = 0.01)
  
  sum_SPLs <-  20 * log10((10^(puttative_SPLs/20)) + (10^(noise/20)))
  
  signal_SPL <- puttative_SPLs[which.min(abs(sum_SPLs - signal.noise))]

  return(signal_SPL)
  }