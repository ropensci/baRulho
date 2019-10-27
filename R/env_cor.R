#' Measure amplitude envelope correlation
#' 
#' \code{env_cor} Measures amplitude envelope correlation in signals referenced in a extended selection table.
#' @usage env_cor(X, parallel = 1, pb = TRUE, method = 1, cor.method = "pearson", 
#' ssmooth = NULL, msmooth = NULL, output = "est")
#' @param X object of class 'selection_table', 'extended_selection_table' created by the function \code{\link[warbleR]{selection_table}} from the warbleR package.
#' @param parallel Numeric. Controls whether parallel computing is applied.
#' It specifies the number of cores to be used. Default is 1 (i.e. no parallel computing).
#' If \code{NULL} (default) then the current working directory is used.
#' @param pb Logical argument to control if progress bar is shown. Default is \code{TRUE}. It can also be
#' set globally using the 'pb' option (see \code{\link{warbleR_options}}).
#' @param method Numeric vector of length 1 to indicate the 'experimental design' for measuring envelope correlation. Two methods are available:
#' \itemize{
#' \item \code{1}: compare all signals with their counterpart that was recorded at the closest distance to source (e.g. compare a signal recorded at 5m, 10m and 15m with its counterpart recorded at 1m). This is the default method. 
#' \item \code{2}: compare all signals with their counterpart recorded at the distance immediately before (e.g. a signal recorded at 10m compared with the same signal recorded at 5m, then signal recorded at 15m compared with same signal recorded at 10m and so on).
#' }
#' @param cor.method Character string indicating the correlation coefficient to be applied ("pearson", "spearman", or "kendall", see \code{\link[stats]{cor}}).
#' @param ssmooth Numeric vector of length 1 determining the length of the sliding window used for a sum smooth for amplitude envelope calculation (used internally by \code{\link[seewave]{env}}).
#' @param msmooth Numeric vector of length 2 to smooth the amplitude envelope with a mean sliding window for amplitude envelope calculation. The first element is the window length (in number of amplitude values) and the second one the window overlap (used internally by \code{\link[seewave]{env}}). 
#' @param output Character vector of length 1 to determine if an extended selection table ('est', default) or a data frame ('data.frame') is returned.
#' @return Data frame or extended selection table (depending on 'output' argument) similar to input data, but also includes a new column 
#' with the amplitude envelope correlation coefficients.
#' @export
#' @name env_cor
#' @details The correlation of amplitude envelopes is intended to measure the distortion of signals in the time domain. The goal of the function is to measure the envelope correlation on signals in which a master playback has been re-recorded at different distances. The 'signal.id' column must be used to indicate the function to only compare signals belonging to the same category (e.g. song-types). The function will then compared each signal type to its reference. Two methods for calculating envelope correlation are provided (see 'method' argument).
#' @examples
#' {
#' # First set temporary folder
#' # setwd(tempdir())
#' 
#' data("playback_est")
#' 
#' # method 1
#'env_cor(X = playback_est)
#' 
#' # method 2
#' env_cor(X = playback_est, method = 2)
#' }
#' 
#' @author Marcelo Araya-Salas (\email{marceloa27@@gmail.com}) 
#last modification on oct-22-2019 (MAS)

env_cor <- function(X, parallel = 1, pb = TRUE, method = 1,  cor.method = "pearson", ssmooth = NULL, msmooth = NULL, output = "est"){
  
  # set pb options 
  on.exit(pbapply::pboptions(type = .Options$pboptions$type), add = TRUE)
  
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

  # set pb options 
  pbapply::pboptions(type = ifelse(as.logical(pb), "timer", "none"))
  
  # set clusters for windows OS
  if (Sys.info()[1] == "Windows" & parallel > 1)
    cl <- parallel::makePSOCKcluster(getOption("cl.cores", parallel)) else cl <- parallel
  
  if (pb) write(file = "", x = "calculating amplitude envelopes (step 1 of 2):")
  
  
  # calculate all envelops apply function
  envs <- pbapply::pblapply(X = 1:nrow(X), cl = cl, FUN = function(y)   {
    clp <- warbleR::read_wave(X = X, index = y)
    env(wave = clp, f = clp@samp.rate, ssmooth = ssmooth, plot = FALSE, msmooth = msmooth)[, 1]
    }) 

  # add sound file selec column and names to envelopes (weird column name so it does not overwrite user columns)
  X$TEMP....y <- names(envs) <- paste(X$sound.files, X$selec, sep = "-")
  
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
    out <- max(cors)
    }
  
    return(out)
    }
  
  # make a data frame with 2 columns with names of the signals to be compare
  X$TEMP....z <- sapply(1:nrow(X), function(x, meth = method){
    
     # extract for single signal and order by distance
     Y <- as.data.frame(X[X$signal.id == X$signal.id[X$TEMP....y == X$TEMP....y[x]], , drop = FALSE])
     Y <- Y[order(Y$distance), ]
     
    # method 1 compare to closest distance to source
    if (meth == 1) z <- Y$TEMP....y[which.min(Y$distance)] else # if method 2
      # if not the first row then the previous row
    if (Y$TEMP....y[1] != X$TEMP....y[x]) z <- X$TEMP....y[x - 1] else # else the first row
    z <- Y$TEMP....y[1] 
    
    return(z)
  })
  
  if (pb) write(file = "", x = "calculating envelope cross-correlations (step 2 of 2):")
  
  # calculate all envelops apply function
  X$env.cor <- pbapply::pbsapply(X = 1:nrow(X), cl = cl, FUN = function(x) {
    env_cor_FUN(y = X$TEMP....y[x], z = X$TEMP....z[x])
  }) 
  
  # remove temporal columns
  X$TEMP....y <- X$TEMP....z <- NULL
  
  # convert to data frame instead of extended selection table
  if (output == "data.frame") 
    X <- as.data.frame(X)
  
  return(X)
}
