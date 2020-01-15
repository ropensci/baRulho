#' Measure attenuation as signal-to-noise ratio 
#' 
#' \code{snr} measures attenuation as signal-to-noise ratio of signals referenced in an extended selection table.
#' @usage snr(X, mar, parallel = 1, pb = TRUE, eq.dur = FALSE,
#' noise.ref = "adjacent", type = 1, bp = NULL, hop.size = 1, wl = NULL)
#' @param X object of class 'extended_selection_table' created by the function \code{\link[warbleR]{selection_table}} from the warbleR package.
#' @param mar numeric vector of length 1. Specifies the margins adjacent to
#'   the start and end points of selection over which to measure ambient noise.
#' @param parallel Numeric vector of length 1. Controls whether parallel computing is applied by specifying the number of cores to be used. Default is 1 (i.e. no parallel computing).
#' @param pb Logical argument to control if progress bar is shown. Default is \code{TRUE}.
#' @param eq.dur Logical. Controls whether the ambient noise segment that is measured has the same duration 
#' to that of the signal (if \code{TRUE}. Default is \code{FALSE}). If \code{TRUE} then 'mar' and 'noise.ref' arguments are ignored.
#' @param noise.ref Character vector of length 1 to determined if a measure ambient noise segment to be used for measuring ambient noise. Two options are available: 
#' \itemize{
#' \item \code{adjacent}: measure ambient noise right before the signal (using argument 'mar' to define duration of ambient noise segments). If several 'ambient' selections by sound file are supplied, then the root mean square of the amplitude envelope will be averaged across those selections.
#' \item \code{custom}: measure ambient noise segments referenced in the selection table (labeled as 'ambient' in the 'signal.type' column). Those segments will be used to apply the same ambient noise reference to all signals in a sound file. Therefore, at least one 'ambient' selection for each sound file must be provided.
#' }
#' @param type  Numeric vector of length 1. Selects the formula to be used to calculate the signal-to-noise ratio (S = signal
#' , N = background noise): 
#' \itemize{
#' \item \code{1}: ratio of S amplitude envelope root mean square to N amplitude envelope root mean square
#'  (\code{rms(env(S))/rms(env(N))})
#' \item \code{2}: ratio of the difference between S amplitude envelope root mean square and N amplitude envelope root mean square to N amplitude envelope root mean square (\code{(rms(env(S)) - rms(env(N)))/rms(env(N))}, as proposed by Dabelsteen et al (1993))
#' }
#' @param bp Numeric vector of length 2 giving the lower and upper limits of a frequency bandpass filter (in kHz). Default is \code{NULL}.
#' @param hop.size A numeric vector of length 1 specifying the time window duration (in ms). Default is 1 ms, which is equivalent to ~45 wl for a 44.1 kHz sampling rate. Ignored if 'wl' is supplied.
#' @param wl A numeric vector of length 1 specifying the window length of the spectrogram, default 
#' is NULL. Ignored if \code{bp = NULL}. If supplied, 'hop.size' is ignored.
#' Note that lower values will increase time resolution, which is more important for amplitude ratio calculations. 
#' @return Extended selection table similar to input data, but also includes a new column (snr.attenuation)
#' with the signal-to-noise ratio values.
#' @export
#' @name snr
#' @details Signal-to-noise ratio (SNR) measures signal amplitude level in relation to ambient noise. A general margin in which ambient noise will be measured must be specified. Alternatively, a selection of ambient noise can be used as reference (see 'noise.ref' argument). When margins overlap with another acoustic signal nearby, SNR will be inaccurate, so margin length should be carefully considered. Any SNR less than or equal to one suggests background noise is equal to or overpowering the acoustic signal. The 'signal.type' column must be used to indicate which signals belong to the same category (e.g. song-types). The function will measure signal-to-noise ratio within the supplied frequency range (e.g. bandpass) of the reference signal ('bottom.freq' and 'top.freq' columns in 'X'). Two methods for calculating signal-to-noise ratio are provided (see 'type' argument).   
#' @examples
#' {
#' # load example data
#' data("playback_est")
#' 
#' # using measure ambient noise reference selections 
#' snr(X = playback_est, mar = 0.05, noise.ref = 'custom')
#' 
#' # remove ambient selections
#' playback_est <- playback_est[playback_est$signal.type != "ambient", ]
#' # using margin for ambient noise of 0.05 and adjacent measure ambient noise reference
#'snr(X = playback_est, mar = 0.05, noise.ref = 'adjacent')
#' }
#' 
#' @author Marcelo Araya-Salas (\email{marceloa27@@gmail.com})
#' @seealso \code{\link{excess_attenuation}}
#' @references {
#' Dabelsteen, T., Larsen, O. N., & Pedersen, S. B. (1993). Habitat-induced degradation of sound signals: Quantifying the effects of communication sounds and bird location on blur ratio, excess attenuation, and signal-to-noise ratio in blackbird song. The Journal of the Acoustical Society of America, 93(4), 2206.
#' 
#' Araya-Salas, M. (2020). baRulho: baRulho: quantifying habitat-induced degradation of (animal) acoustic signals in R. R package version 1.0.0
#' }
#last modification on nov-01-2019 (MAS)

snr <- function(X, mar, parallel = 1, pb = TRUE, eq.dur = FALSE,
                       noise.ref = "adjacent", type = 1, bp = NULL, hop.size = 1, 
                       wl = NULL){
  
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
  
  # check signal.type column 
  if (is.null(X$signal.type)) stop("'X' must containe a 'signal.type' column")
  
  #check noise.ref
  if (!any(noise.ref %in% c("custom", "adjacent"))) stop("'noise.ref' must be either 'custom' or 'adjacent'")  
  
  # check if 'ambient' is found in  signal.type column 
  if (!any(X$signal.type %in% 'ambient') & noise.ref == "custom") stop("'ambient' selections must be contained in 'X' (and label in 'signal.type' column) when 'noise.ref == TRUE'")
  
  # check if 'ambient' is found in  signal.type column 
  if (!any(X$signal.type %in% 'ambient') & noise.ref == "custom") stop("'ambient' selections must be contained in 'X' (and label in 'signal.type' column) when 'noise.ref == TRUE'")
  
  if (noise.ref == "custom" & any(sapply(unique(X$sound.files), function(x) sum(X$sound.files == x & X$signal.type == "ambient")) == 0)) stop("Each sound file referenced in 'X' must have at least 1 'ambient' selection when 'noise.ref == custom'")
  
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
      
      # cut ambient noise before signal
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
    
    if (X$signal.type[y] != "ambient"){
      
    # signal RMS
    sig_RMS <- seewave::rms(envs[[X$TEMP....y[y]]]$sig.env)  
    
    # get reference ambient noise RMS
    if (noise.ref == "adjacent"){
      bg_RMS <- seewave::rms(envs[[X$TEMP....y[y]]]$bg.env)  
    } else {
      # get envelopes from ambient selections
      bg_envs <- sapply(envs[X$TEMP....y[X$sound.files == X$sound.files[y] & X$signal.type == "ambient"]], "[", 'sig.env')

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
  
  return(X)
  }
