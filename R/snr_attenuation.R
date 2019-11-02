#' Measure attenuation as signal-to-noise ratio 
#' 
#' \code{snr_attenuation} measures attenuation as signal-to-noise ratio of signals referenced in a extended selection table.
#' @usage snr_attenuation(X, mar, parallel = 1, pb = TRUE, eq.dur = FALSE,
#' noise.ref = "adjacent", type = 1, bp = NULL, wl = 10, output = "est")
#' @param X object of class 'selection_table', 'extended_selection_table' created by the function \code{\link[warbleR]{selection_table}} from the warbleR package.
#' @param mar numeric vector of length 1. Specifies the margins adjacent to
#'   the start and end points of selection over which to measure noise.
#' @param parallel Numeric. Controls whether parallel computing is applied.
#' It specifies the number of cores to be used. Default is 1 (i.e. no parallel computing).
#' If \code{NULL} (default) then the current working directory is used.
#' @param pb Logical argument to control if progress bar is shown. Default is \code{TRUE}. It can also be
#' set globally using the 'pb' option (see \code{\link{warbleR_options}}).
#' @param eq.dur Logical. Controls whether the noise segment that is measured has the same duration 
#' than the signal (if \code{TRUE}. Default is \code{FALSE}). If \code{TRUE} then 'mar' and 'noise.ref' arguments are ignored.
#' @param noise.ref Character vector of length 1 to determined if a noise segment to be used for measuring ambient noise. Two options are available: 
#' \itemize{
#' \item \code{adjacent}: measure ambient noise right before the signal (using argument 'mar' to define duration of noise segments). If several 'noise' selections by sound file are supplied, then the root mean square of the amplitude envelope will be averaged across those selections.
#' \item \code{custom}: measure noise segments referenced in the selection table (labeled as 'noise' in the 'signal.id' column). Those segments will be used to apply the same noise reference to all signals in a sound file. Therefore, at least one 'noise' selection for each sound file must be provided.
#' }
#' @param type Numeric. Determine the formula to be used to calculate the signal-to-noise ratio (S = signal
#' , N = background noise): 
#' \itemize{
#' \item \code{1}: ratio of S amplitude envelope quadratic mean to N amplitude envelope quadratic mean
#'  (\code{rms(env(S))/rms(env(N))})
#' \item \code{2}: ratio of the difference between S amplitude envelope quadratic mean and N amplitude envelope quadratic mean to N amplitude envelope quadratic mean (\code{(rms(env(S)) - rms(env(N)))/rms(env(N))}, as proposed by Dabelsteen et al (1993))
#' }
#' @param bp Numeric vector of length 2 giving the lower and upper limits of a frequency bandpass filter (in kHz). Default is \code{NULL}.
#' @param wl A numeric vector of length 1 specifying the window length of the spectrogram for applying bandpass. Default 
#'   is 10. Ignored if \code{bp = NULL}.
#'  Note that lower values will increase time resolution, which is more important for amplitude ratio calculations. 
#' @param output Character vector of length 1 to determine if an extended selection table ('est', default) or a data frame ('data.frame') is returned.
#' @return Data frame or extended selection table (depending on 'output' argument) similar to input data, but also includes a new column (snr.attenuation)
#' with the signal-to-noise ratio values.
#' @export
#' @name snr_attenuation
#' @details Signal-to-noise ratio is the attenuation is the attenuation of a sound in excess of that due to spherical spreading as described by Dabelsteen et al (1993). The goal of the function is to measure the excess attenuation on signals in which a master playback has been re-recorded at different distances. The 'signal.id' column must be used to indicate which signals belonging to the same category (e.g. song-types). The function will then compared each signal type to its reference. Two methods for calculating excess attenuation are provided (see 'method' argument).   
#' @examples
#' {
#' # load example data
#' data("playback_est")
#' 
#' # using noise reference selections 
#' snr_attenuation(X = playback_est, mar = 0.05, noise.ref = 'custom')
#' 
#' #' # remove noise selections
#' playback_est <- playback_est[playback_est$signal.id != "noise", ]
#' # using margin for noise of 0.05 and adjacent noise reference
#'snr_attenuation(X = playback_est, mar = 0.05, noise.ref = 'adjacent')
#' }
#' 
#' @author Marcelo Araya-Salas (\email{marceloa27@@gmail.com}) #' @references {
#' Dabelsteen, T., Larsen, O. N., & Pedersen, S. B. (1993). Habitat-induced degradation of sound signals: Quantifying the effects of communication sounds and bird location on blur ratio, excess attenuation, and signal-to-noise ratio in blackbird song. The Journal of the Acoustical Society of America, 93(4), 2206.
#' 
#' Araya-Salas, M. (2019), baRulho: a R package to evaluate habitat-induced degradation of (animal) acoustic signals. R package version 1.0.0
#' }
#last modification on nov-01-2019 (MAS)

snr_attenuation <- function(X, mar, parallel = 1, pb = TRUE, eq.dur = FALSE,
                       noise.ref = "adjacent", type = 1, bp = NULL, wl = 10, output = "est"){
  
  # set pb options 
  on.exit(pbapply::pboptions(type = .Options$pboptions$type), add = TRUE)
  
  # If parallel is not numeric
  if (!is.numeric(parallel)) stop("'parallel' must be a numeric vector of length 1") 
  if (any(!(parallel %% 1 == 0),parallel < 1)) stop("'parallel' should be a positive integer")
  
  # check signal.id column 
  if (is.null(X$signal.id)) stop("'X' must containe a 'signal.id' column")
  
  #check noise.ref
  if (!any(noise.ref %in% c("custom", "adjacent"))) stop("'noise.ref' must be either 'custom' or 'adjacent'")  
  
  # check if 'noise' is found in  signal.id column 
  if (!any(X$signal.id %in% 'noise') & noise.ref == "custom") stop("'noise' selections must be contained in 'X' (and label in 'signal.id' column) when 'noise.red == TRUE'")
  
  # check if 'noise' is found in  signal.id column 
  if (!any(X$signal.id %in% 'noise') & noise.ref == "custom") stop("'noise' selections must be contained in 'X' (and label in 'signal.id' column) when 'noise.red == TRUE'")
  
  if (noise.ref == "custom" & any(sapply(unique(X$sound.files), function(x) sum(X$sound.files == x & X$signal.id == "noise")) == 0)) stop("Each sound file referenced in 'X' must have at least 1 'noise' selection when 'noise.ref == custom'")
  
  #check output
  if (!any(output %in% c("est", "data.frame"))) stop("'output' must be either 'est' or 'data.frame'")  
  
  # set pb options 
  pbapply::pboptions(type = ifelse(as.logical(pb), "timer", "none"))
  
  # set clusters for windows OS
  if (Sys.info()[1] == "Windows" & parallel > 1)
    cl <- parallel::makePSOCKcluster(getOption("cl.cores", parallel)) else cl <- parallel
  
  # calculate all envelops with a apply function
  envs <- pbapply::pblapply(X = 1:nrow(X), cl = cl, FUN = function(y)   {
    if (noise.ref == "custom"){
      
      # read signal clip
      signal <- warbleR::read_wave(X = X, index = y)
      
      # add band-pass frequency filter
      if (!is.null(bp)) {
        
        signal <- seewave::ffilter(signal, f = signal@samp.rate, from = bp[1] * 1000, ovlp = 0,
                              to = bp[2] * 1000, bandpass = TRUE, wl = wl, 
                              output = "Wave")
      }
      
      # get RMS for signal
      sig.env <- seewave::env(signal, f = signal@samp.rate, envt = "abs", plot = FALSE)
      
      bg.env <- NA
    } 
    
    if (noise.ref == "adjacent"){
      
      
      # Read sound files to get sample rate and length
      r <- warbleR::read_wave(X = X, index = y, header = TRUE)
      
      # read sample rate
      f <- r$sample.rate
      
      # set margin to half of signal duration
      if (eq.dur) mar <- (X$end[y] - X$start[y])/2
      
      #reset time coordinates of signals if lower than 0 o higher than duration
      stn <- X$start[y] - mar
      enn <- X$end[y] + mar
      mar1 <- mar
      
      if (stn < 0) { 
        mar1 <- mar1  + stn
        stn <- 0
      }
      
      mar2 <- mar1 + X$end[y] - X$start[y]
      
      if (enn > r$samples/f) enn <- r$samples/f
      
      # read signal and margin
      r <- warbleR::read_wave(X = X, index = y, from = stn, to = enn)
      
      # read clip with signal  
      signal <- warbleR::read_wave(X = X, index = y)
      
      # get RMS for signal
      sig.env <- seewave::env(signal, f = signal@samp.rate, envt = "abs", plot = FALSE)
      
      # cut noise before signal
      noise1 <- seewave::cutw(r, from =  0, to = mar1, f = f)
      
      # get RMS for background noise
      bg.env <- seewave::env(noise1, f = f, envt = "abs", plot = FALSE)
    }
    return(list(sig.env = sig.env, bg.env = bg.env))
  }) 
  
  # add sound file selec column and names to envelopes (weird column name so it does not overwrite user columns)
  X$TEMP....y <- names(envs) <- paste(X$sound.files, X$selec, sep = "-")

  # calculate SNR 
  X$snr.attenuation <- sapply(1:nrow(X), function(y){
    
    if (X$signal.id[y] != "noise"){
      
    # signal RMS
    sig_RMS <- seewave::rms(envs[[X$TEMP....y[y]]]$sig.env)  
    
    # get referene noise RMS
    if (noise.ref == "adjacent"){
      bg_RMS <- seewave::rms(envs[[X$TEMP....y[y]]]$bg.env)  
    } else {
      # get envelopes from noise selections
      bg_envs <- sapply(envs[X$TEMP....y[X$sound.files == X$sound.files[y] & X$signal.id == "noise"]], "[", 'sig.env')

      # get mean RMS from combined envelopes
      bg_RMS <- seewave::rms(unlist(sapply(bg_envs, as.vector)))
      }
    
    # Calculate signal-to-noise ratio
   if (type == 1)
    snr <- sig_RMS / bg_RMS
   
   if (type == 2)
     snr <- (sig_RMS - bg_RMS) / bg_RMS
   
    return(20*log10(snr))  
    } else return(NA) # return NA if current row is noise
  })
  
  # remove temporary column
  X$TEMP....y <- NULL  

  # convert to data frame instead of extended selection table
  if (output == "data.frame") 
    X <- as.data.frame(X)
  
  return(X)
  }
