#' Measure frequency spectrum correlation
#' 
#' \code{spectral_correlation} measures frequency spectrum correlation of signals referenced in an extended selection table.
#' @usage spectral_correlation(X, parallel = 1, pb = TRUE, method = 1, 
#' cor.method = "pearson", hop.size = 11.6, wl = NULL, ovlp = 70)
#' @param X object of class 'extended_selection_table' created by the function \code{\link[warbleR]{selection_table}} from the warbleR package.
#' @param parallel Numeric vector of length 1. Controls whether parallel computing is applied by specifying the number of cores to be used. Default is 1 (i.e. no parallel computing).
#' If \code{NULL} (default) then the current working directory is used.
#' @param pb Logical argument to control if progress bar is shown. Default is \code{TRUE}.
#' @param method Numeric vector of length 1 to indicate the 'experimental design' to measure frequency spectrum correlation. Two methods are available:
#' \itemize{
#' \item \code{1}: compare all signals with their counterpart that was recorded at the closest distance to source (e.g. compare a signal recorded at 5m, 10m and 15m with its counterpart recorded at 1m). This is the default method. 
#' \item \code{2}: compare all signals with their counterpart recorded at the distance immediately before (e.g. a signal recorded at 10m compared with the same signal recorded at 5m, then signal recorded at 15m compared with same signal recorded at 10m and so on).
#' }
#' @param cor.method Character string indicating the correlation coefficient to be applied ("pearson", "spearman", or "kendall", see \code{\link[stats]{cor}}).
#' @param hop.size A numeric vector of length 1 specifying the time window duration (in ms). Default is 11.6 ms, which is equivalent to 512 wl for a 44.1 kHz sampling rate. Ignored if 'wl' is supplied.
#' @param wl A numeric vector of length 1 specifying the window length of the spectrogram, default 
#' is NULL. If supplied, 'hop.size' is ignored.
#' @param ovlp Numeric vector of length 1 specifying the percent overlap between two 
#'   consecutive windows, as in \code{\link[seewave]{spectro}}. Default is 70.
#' @return Extended selection table similar to input data, but also includes a new column ('spectrum.correlation')
#' with the calculated frequency spectrum correlation coefficients.
#' @export
#' @name spectral_correlation
#' @details spectrum correlation measures the similarity of two signals in the frequency domain. The function measures the spectrum correlation coefficients of signals in which a reference playback has been re-recorded at increasing distances. Values range from 1 (identical frequency spectrum, i.e. no degradation) to 0. The 'signal.type' column must be used to indicate the function to only compare signals belonging to the same category (e.g. song-types). The function will then compare each signal type to the corresponding reference signal. Two methods for calculating spectrum correlation are provided (see 'method' argument). Use \code{\link{spectral_blur_ratio}} to get spectra for plotting. 
#' @examples
#' {
#' # load example data
#' data("playback_est")
#' 
#' # remove ambient selections
#' pe <- playback_est[playback_est$signal.type != "ambient", ]
#' 
#' # method 1
#'spectral_correlation(X = pe)
#' 
#' # method 2
#' spectral_correlation(X = pe, method = 2)
#' }
#' 
#' @author Marcelo Araya-Salas (\email{marceloa27@@gmail.com})
#' @seealso \code{\link{envelope_correlation}}, \code{\link{spectral_blur_ratio}} 
#' @references {
#' Araya-Salas, M. (2020). baRulho: baRulho: quantifying habitat-induced degradation of (animal) acoustic signals in R. R package version 1.0.0
#' 
#' Apol, C.A., Sturdy, C.B. & Proppe, D.S. (2017). Seasonal variability in habitat structure may have shaped acoustic signals and repertoires in the black-capped and boreal chickadees. Evol Ecol. 32:57-74.
#' }
#last modification on nov-01-2019 (MAS)

spectral_correlation <- function(X, parallel = 1, pb = TRUE, method = 1, cor.method = "pearson", hop.size = 11.6, wl = NULL, ovlp = 70){
  
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
  
  if (pb) write(file = "", x = "calculating frequency spectrums (step 1 of 2):")
  
  # calculate all spectra apply function
  spcs <- pbapply::pblapply(X = 1:nrow(X), cl = cl, FUN = function(y, wle = wl, ovl = ovlp){
   
    # load clip
    clp <- warbleR::read_wave(X = X, index = y)
    
    # mean spec
    mspc <- meanspec(wave = clp, f = clp@samp.rate, plot = FALSE, wl = wle, ovlp = ovl)  
  
    return(mspc)
    }) 
  
  # add sound file selec names to envelopes (weird column name so it does not overwrite user columns)
  names(spcs) <- X$TEMP....sgnl
  
  # function to measure envelope correlation
  # y and z are the sound.files+selec names of the signals and reference signal (model)
  spctr_cor_FUN <- function(y, z){
    
    # if names are the same return NA
    if (y == z) cor.spctr <- NA else {
      
      # extract envelope for signal and model 
      sgnl.spctr <- spcs[[which(names(spcs) == y)]]
      mdl.spctr <- spcs[[which(names(spcs) == z)]]
      
      
      ### filter to freq range of signals and remove freq column
      # get range as lowest bottom and highest top
      frng <- c(min(X$bottom.freq[X$TEMP....sgnl %in% c(y, z)]), max(X$top.freq[X$TEMP....sgnl %in% c(y, z)]))
      sgnl.spctr <- sgnl.spctr[sgnl.spctr[, 1] > frng[1] & sgnl.spctr[, 1] < frng[2], 2]
      mdl.spctr <- mdl.spctr[mdl.spctr[, 1] > frng[1] & mdl.spctr[, 1] < frng[2], 2]
      
      # get correlation assuming they have same length  
      cor.spctr <- cor(sgnl.spctr, mdl.spctr, method = cor.method)
    }
    
    return(cor.spctr)
  }
  
  # make a data frame with 2 columns with names of the signals to be compare
  X$reference <- sapply(1:nrow(X), function(x, meth = method){
    
    # extract for single signal and order by distance
    Y <- as.data.frame(X[X$signal.type == X$signal.type[X$TEMP....sgnl == X$TEMP....sgnl[x]], , drop = FALSE])
    Y <- Y[order(Y$distance), ]
    
    # method 1 compare to closest distance to source
    if (meth == 1) z <- Y$TEMP....sgnl[which.min(Y$distance)] else # if method 2
      # if not the first row then the previous row
      if (Y$TEMP....sgnl[1] != X$TEMP....sgnl[x]) z <- X$TEMP....sgnl[x - 1] else # else the first row
        z <- Y$TEMP....sgnl[1] 
    
    return(z)
  })
  
  if (pb) write(file = "", x = "calculating envelope cross-correlations (step 2 of 2):")
  
  # calculate all envelops apply function
  X$spectrum.correlation <- pbapply::pbsapply(X = 1:nrow(X), cl = cl, FUN = function(x) {
    spctr_cor_FUN(y = X$TEMP....sgnl[x], z = X$reference[x])
  }) 
  
  # remove temporal columns
  X$TEMP....sgnl <-NULL
  
  return(X)
}
