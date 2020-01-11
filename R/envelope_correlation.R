#' Measure amplitude envelope correlation
#' 
#' \code{envelope_correlation} measures amplitude envelope correlation of signals referenced in an extended selection table.
#' @usage envelope_correlation(X, parallel = 1, pb = TRUE, method = 1, cor.method = "pearson", 
#' ssmooth = NULL, msmooth = NULL, hop.size = 11.6, wl = NULL, ovlp = 70)
#' @param X object of class 'selection_table', 'extended_selection_table' created by the function \code{\link[warbleR]{selection_table}} from the warbleR package, or can be generated from Raven and Syrinx selection tables using the \code{\link[Rraven]{imp_raven}} or \code{\link[Rraven]{imp_syrinx}} functions from the Rraven package.
#' @param parallel Numeric vector of length 1. Controls whether parallel computing is applied by specifying the number of cores to be used. Default is 1 (i.e. no parallel computing).
#' If \code{NULL} (default) then the current working directory is used.
#' @param pb Logical argument to control if progress bar is shown. Default is \code{TRUE}.
#' @param method Numeric vector of length 1 to indicate the 'experimental design' to measure amplitude envelope correlation. Two methods are available:
#' \itemize{
#' \item \code{1}: compare all signals with their counterpart that was recorded at the closest distance to source (e.g. compare a signal recorded at 5m, 10m and 15m with its counterpart recorded at 1m). This is the default method. 
#' \item \code{2}: compare all signals with their counterpart recorded at the distance immediately before (e.g. a signal recorded at 10m compared with the same signal recorded at 5m, then signal recorded at 15m compared with same signal recorded at 10m and so on).
#' }
#' @param cor.method Character string indicating the correlation coefficient to be applied ("pearson", "spearman", or "kendall", see \code{\link[stats]{cor}}).
#' @param ssmooth Numeric vector of length 1 to determine the length of the sliding window used for a sum smooth for amplitude envelope calculation (used internally by \code{\link[seewave]{env}}).
#' @param msmooth Numeric vector of length 2 to smooth the amplitude envelope with a mean sliding window for amplitude envelope calculation. The first element is the window length (in number of amplitude values) and the second one the window overlap (used internally by \code{\link[seewave]{env}}). 
#' @param hop.size A numeric vector of length 1 specifying the time window duration (in ms). Default is 11.6 ms, which is equivalent to 512 wl for a 44.1 kHz sampling rate. Ignored if 'wl' is supplied.
#' @param wl A numeric vector of length 1 specifying the window length of the spectrogram, default 
#' is NULL. If supplied, 'hop.size' is ignored.
#' @param ovlp Numeric vector of length 1 specifying the percent overlap between two 
#'   consecutive windows, as in \code{\link[seewave]{spectro}}. Default is 70.
#' @return Extended selection table similar to input data, but also includes two new columns ('reference' and  'envelope.correlation')
#' with the reference signal and the amplitude envelope correlation coefficients.
#' @export
#' @name envelope_correlation
#' @details Amplitude envelope correlation measures the similarity of two signals in the time domain. The  function measures the envelope correlation coefficients of signals in which a reference playback has been re-recorded at increasing distances. Values close to 1 means very similar amplitude envelopes (i.e. little degradation has occurred). If envelopes have different lengths (which means signals have different lengths) cross-correlation is used and the maximum correlation coefficient is returned. Cross-correlation is achieved by sliding the shortest signal along the largest one and calculating the correlation at each step. The 'signal.type' column must be used to indicate the function to only compare signals belonging to the same category (e.g. song-types).The function compares each signal type to the corresponding reference signal within the supplied frequency range (e.g. bandpass) of the reference signal ('bottom.freq' and 'top.freq' columns in 'X'). Two methods for calculating envelope correlation are provided (see 'method' argument). Use \code{\link{blur_ratio}} to extract envelopes. 
#' @seealso \code{\link{blur_ratio}}, \code{\link{spectral_blur_ratio}}
#' @examples
#' {
#' # load example data
#' data("playback_est")
#' 
#' # remove ambient selections
#' playback_est <- playback_est[playback_est$signal.type != "ambient", ]
#' 
#' # method 1
#'envelope_correlation(X = playback_est)
#' 
#' # method 2
#' envelope_correlation(X = playback_est, method = 2)
#' }
#' 
#' @author Marcelo Araya-Salas (\email{marceloa27@@gmail.com})
#' @references {
#' Araya-Salas, M. (2020). baRulho: baRulho: quantifying habitat-induced degradation of (animal) acoustic signals in R. R package version 1.0.0
#' 
#' Apol, C.A., Sturdy, C.B. & Proppe, D.S. (2017). Seasonal variability in habitat structure may have shaped acoustic signals and repertoires in the black-capped and boreal chickadees. Evol Ecol. 32:57-74.
#' }
#last modification on nov-01-2019 (MAS)

envelope_correlation <- function(X, parallel = 1, pb = TRUE, method = 1,  cor.method = "pearson", ssmooth = NULL, msmooth = NULL, hop.size = 11.6, wl = NULL, ovlp = 70){
  
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
  if (!is.character(cor.method)) stop("'cor.method' must be a character vector of length 1") 
  if (!any(cor.method %in%  c("pearson", "kendall", "spearman"))) stop("'method' must be either  'pearson', 'kendall' or 'spearman'")
  
  # check signal.type column 
  if (is.null(X$signal.type)) stop("'X' must containe a 'signal.type' column")
  
  # add sound file selec column and names to X (weird column name so it does not overwrite user columns)
  X <- prep_X_bRlo_int(X, method = method)
  
  # set pb options 
  pbapply::pboptions(type = ifelse(as.logical(pb), "timer", "none"))
  
  # set clusters for windows OS
  if (Sys.info()[1] == "Windows" & parallel > 1)
    cl <- parallel::makePSOCKcluster(getOption("cl.cores", parallel)) else cl <- parallel
  
  if (pb) write(file = "", x = "calculating amplitude envelopes (step 1 of 2):")
  
  # calculate all envelopes apply function
  envs <- pbapply::pblapply(X = 1:nrow(X), cl = cl, FUN = function(y)   {
    
    # get clip
    clp <- warbleR::read_wave(X = X, index = y)
    
    # define bandpass based on reference
    bp <- c(X$bottom.freq[X$TEMP....sgnl == X$reference[y]], X$top.freq[X$TEMP....sgnl == X$reference[y]])
    
    # bandpass filter
    clp <- seewave::ffilter(clp, from = bp[1] * 1000, 
                            ovlp = ovlp, to = bp[2] * 1000, bandpass = TRUE, 
                            wl = wl, output = "Wave")
    
    # calculate envelope
    nv <- env(wave = clp, f = clp@samp.rate, ssmooth = ssmooth, plot = FALSE, msmooth = msmooth)[, 1]
  
    return(nv)
  }) 
  
  # add sound file selec column and names to envelopes
  names(envs) <- X$TEMP....sgnl
  
  # function to measure envelope correlation
  # y and z are the sound.files+selec names of the signals and reference signal (model)
  env_cor_FUN <- function(y, z){
    
    # if names are the same return NA
    if (y == z) out <- NA else {
      
      # extract envelope for signal and model 
      sgnl.env <- envs[[which(names(envs) == y)]]
      mdl.env <- envs[[which(names(envs) == z)]]
      
      # define short and long envelope for sliding one (short) over the other (long)
      if(length(mdl.env) > length(sgnl.env)) {
        lg.env <- mdl.env
        shrt.env <- sgnl.env
      } else {
        lg.env <- sgnl.env
        shrt.env <- mdl.env
      }
      
      # get length of shortest minus 1 (1 if same length so it runs a single correlation)
      shrt.lgth <- length(shrt.env) - 1
      
      # steps for sliding one signal over the other  
      stps <- length(lg.env) - shrt.lgth
      
      # calculate correlations at each step
      cors <- sapply(1:stps, function(x) {
        cor(lg.env[x:(x + shrt.lgth)], shrt.env, method = cor.method)
      })
    
    # return maximum correlation
    out <- max(cors, na.rm = TRUE)
    }
  
    return(out)
    }
  

  if (pb) write(file = "", x = "calculating envelope cross-correlations (step 2 of 2):")
  
  # calculate all envelops apply function
  X$env.cor <- pbapply::pbsapply(X = 1:nrow(X), cl = cl, FUN = function(x) {
    env_cor_FUN(y = X$TEMP....sgnl[x], z = X$reference[x])
  }) 
  
  # remove temporal columns
  X$TEMP....sgnl <- NULL
  
  return(X)
}
