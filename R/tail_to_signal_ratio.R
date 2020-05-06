#' Measure reverberations as tail-to-signal ratio 
#' 
#' \code{tail_to_signal_ratio} measures reverberations as tail-to-signal ratio of signals referenced in an extended selection table.
#' @usage tail_to_signal_ratio(X, mar, parallel = 1, pb = TRUE,  type = 1, 
#' bp = 'freq.range', output = "est", hop.size = 1, wl = NULL)
#' @param X object of class 'extended_selection_table' created by the function \code{\link[warbleR]{selection_table}} from the warbleR package.
#' @param mar numeric vector of length 1. Specifies the margins adjacent to
#'   the start and end points of selection over which to measure ambient noise.
#' @param parallel Numeric vector of length 1. Controls whether parallel computing is applied by specifying the number of cores to be used. Default is 1 (i.e. no parallel computing).
#' @param pb Logical argument to control if progress bar is shown. Default is \code{TRUE}.
#' @param bp Numeric vector of length 2 giving the lower and upper limits of a frequency bandpass filter (in kHz). Alternatively, when set to 'freq.range' (default), which will make the funcition use the 'bottom.freq' and 'top.freq' as the bandpass range.
#' @param output Character vector of length 1 to determine if an extended selection table ('est', default) or a data frame ('data.frame').
#' @param hop.size A numeric vector of length 1 specifying the time window duration (in ms). Default is 1 ms, which is equivalent to ~45 wl for a 44.1 kHz sampling rate. Ignored if 'wl' is supplied.
#' @param type Numeric. Determine the formula to be used to calculate the tail-to-signal ratio (S = signal, T = tail, N = background noise): 
#' \itemize{
#' \item \code{1}: ratio of T amplitude envelope quadratic mean to S amplitude envelope quadratic mean
#'  (\code{rms(env(T))/rms(env(S))}) as described by Dabelsteen et al. (1993).
#' \item \code{2}: ratio of T amplitude envelope quadratic mean to N amplitude envelope quadratic mean (\code{rms(env(T))/rms(env(N))}). N is measure in the margin right before the signal. So type 2 actually measures tail-to-noise ratio.
#' }
#' @param wl A numeric vector of length 1 specifying the window length of the spectrogram, default 
#' is NULL. Ignored if \code{bp = NULL}. If supplied, 'hop.size' is ignored.
#' Note that lower values will increase time resolution, which is more important for amplitude ratio calculations. 
#' @return Extended selection table similar to input data, but also includes a new column (tail.to.signal.ratio)
#' with the tail-to-signal ratio values.
#' @export
#' @name tail_to_signal_ratio
#' @details Tail-to-signal ratio (TSR) measures ratio of energy in the tail of reverberations to energy in the signal. A general margin in which reverberation tail will be measured must be specified. The function will measure TSR within the supplied frequency range (e.g. bandpass) of the reference signal ('bottom.freq' and 'top.freq' columns in 'X'). Two methods for calculating reverberations are provided (see 'type' argument). Note that 'type' 2 is not equivalent to the original description of TSR in Dabelsteen et al. (1993) and  is better referred to as tail-to-noise ratio.  
#' @examples
#' {
#' # load example data
#' data("playback_est")
#' 
#'  # remove noise selections
#'  pe <- playback_est[playback_est$signal.type != "ambient", ]
#'  
#'  # using margin for noise of 0.01
#'  tail_to_signal_ratio(X = pe, mar = 0.01, bp = NULL)
#'  
#'  # tail-to-noise ratio (type 2)
#'  tail_to_signal_ratio(X = playback_est, mar = 0.01, type = 2)
#' }
#' 
#' @author Marcelo Araya-Salas (\email{marceloa27@@gmail.com})
#' @seealso \code{\link{excess_attenuation}}
#' @references {
#' Araya-Salas, M. (2020). baRulho: baRulho: quantifying habitat-induced degradation of (animal) acoustic signals in R. R package version 1.0.2
#' 
#' Dabelsteen, T., Larsen, O. N., & Pedersen, S. B. (1993). Habitat-induced degradation of sound signals: Quantifying the effects of communication sounds and bird location on blur ratio, excess attenuation, and signal-to-noise ratio in blackbird song. The Journal of the Acoustical Society of America, 93(4), 2206.
#' 
#' Darden, SK, Pedersen SB, Larsen ON, & Dabelsteen T. (2008). Sound transmission at ground level in a short-grass prairie habitat and its implications for long-range communication in the swift fox *Vulpes velox*. The Journal of the Acoustical Society of America, 124(2), 758-766.
#' 
#' }
#last modification on nov-01-2019 (MAS)

tail_to_signal_ratio <- function(X, mar, parallel = 1, pb = TRUE,
                      type = 1, bp = 'freq.range', output = "est", hop.size = 1, 
                       wl = NULL){
  
  # get call argument names
  argus <- names(as.list(base::match.call()))
  
  # set pb options 
  on.exit(pbapply::pboptions(type = .Options$pboptions$type), add = TRUE)
  
  # is extended sel tab
  if (!warbleR::is_extended_selection_table(X)) 
    stop("'X' must be and extended selection table")
  
  # If parallel is not numeric
  if (!is.numeric(parallel)) stop("'parallel' must be a numeric vector of length 1") 
  if (any(!(parallel %% 1 == 0),parallel < 1)) stop("'parallel' should be a positive integer")
  
  #check output
  if (!any(output %in% c("est", "data.frame"))) stop("'output' must be 'est' or 'data.frame'")  
  
  # hopsize  
  if (!is.numeric(hop.size) | hop.size < 0) stop("'parallel' must be a positive number") 
  
  # adjust wl based on hope.size
  if (is.null(wl))
    wl <- round(attr(X, "check.results")$sample.rate[1] * hop.size, 0)
  
  # make wl even if odd
  if (!(wl %% 2) == 0) wl <- wl + 1
  
  # check signal.type column 
  if (is.null(X$signal.type)) stop("'X' must containe a 'signal.type' column")
  
  # set pb options 
  pbapply::pboptions(type = ifelse(as.logical(pb), "timer", "none"))
  
  # set clusters for windows OS
  if (Sys.info()[1] == "Windows" & parallel > 1)
    cl <- parallel::makePSOCKcluster(getOption("cl.cores", parallel)) else cl <- parallel
  

  # calculate STR 
  X$tail.to.signal.ratio <- pbapply::pbsapply(1:nrow(X), cl = cl, function(y){
    
    if (X$signal.type[y] != "ambient"){
     
      # Read sound files to get sample rate and length
      r <- warbleR::read_wave(X = X, index = y, from = 0, to = Inf, header = TRUE)
      
      # read sample rate
      f <- r$sample.rate
      
      #reset time coordinates of signals if higher than duration
      enn <- X$end[y] + mar
      if (enn > r$samples/f) enn <- r$samples/f
    
      # read signal and margin
      tail.wv <- warbleR::read_wave(X = X, index = y, from = X$end[y], to = enn)
      
      # read signal  
      if (type == 1)
      signal <- warbleR::read_wave(X = X, index = y)
  
      # read background noise right before the signal
      if (type == 2)
      signal <- warbleR::read_wave(X = X, index = y, from = X$start[y] - mar, to = X$start[y])
      
      # add band-pass frequency filter
      if (!is.null(bp)) {
        
        # filter to bottom and top freq range
        if (bp == "freq.range") 
          bp <- c(X$bottom.freq[y], X$top.freq[y])
        
        signal <- seewave::ffilter(signal, f = signal@samp.rate, from = bp[1] * 1000, ovlp = 0,
                                   to = bp[2] * 1000, bandpass = TRUE, wl = wl, 
                                   output = "Wave")
        
        tail.wv <- seewave::ffilter(tail.wv, f = tail.wv@samp.rate, from = bp[1] * 1000, ovlp = 0,
                                   to = bp[2] * 1000, bandpass = TRUE, wl = wl, 
                                   output = "Wave")  
      }
      
      
      # get RMS for signal (or noise if type 2)
      sig.env <- seewave::env(signal, f = signal@samp.rate, envt = "abs", plot = FALSE)
      
      # get RMS for background noise
      tail.env <- seewave::env(tail.wv, f = f, envt = "abs", plot = FALSE)

    # signal (or noise) RMS
    sig_RMS <- seewave::rms(sig.env)  
    
    # get reference ambient noise RMS
    tail_RMS <- seewave::rms(tail.env)  
    
    # Calculate tail.to.signal ratio
    str <- sig_RMS / tail_RMS
   
    return(20*log10(str))  
    } else return(NA) # return NA if current row is noise
  })
  
  if (output == "data.frame") X <- as.data.frame(X)
  
  return(X)
  }
