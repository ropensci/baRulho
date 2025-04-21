#' Measure amplitude envelope correlation
#'
#' \code{envelope_correlation} measures amplitude envelope correlation of sounds referenced in an extended selection table.
#' @inheritParams template_params
#' @param env.smooth Numeric vector of length 1 to determine the length of the sliding window used for a sum smooth for amplitude envelope calculation (used internally by \code{\link[seewave]{env}}). Can be set globally for the current R session via the "env.smooth" option (see \code{\link[base]{options}}).
#' @param ovlp Numeric vector of length 1 specifying the percentage of overlap between two
#'   consecutive windows, as in \code{\link[seewave]{spectro}}. Default is 70. Can be set globally for the current R session via the "ovlp" option (see \code{\link[base]{options}}).
#' @return Object 'X' with an additional column, 'envelope.correlation', containing the computed envelope correlation coefficients.
#' @export
#' @name envelope_correlation
#' @details Amplitude envelope correlation measures the similarity of two sounds in the time domain. The function measures the envelope correlation coefficients of sounds in which a reference playback has been re-recorded at increasing distances. Values close to 1 means very similar amplitude envelopes (i.e. little degradation has occurred). If envelopes have different lengths (which means sounds have different lengths) cross-correlation is used and the maximum correlation coefficient is returned. Cross-correlation is achieved by sliding the shortest sound along the largest one and computing the correlation at each step. The 'sound.id' column must be used to indicate the function to only compare sounds belonging to the same category (e.g. song-types). The function compares each sound to the corresponding reference sound within the supplied frequency range (e.g. bandpass) of the reference sound ('bottom.freq' and 'top.freq' columns in 'X'). Two methods for computing envelope correlation are provided (see 'method' argument). Use \code{\link{blur_ratio}} to create envelopes graphs.
#' @seealso \code{\link{blur_ratio}}, \code{\link{spectrum_blur_ratio}}
#' @family quantify degradation
#' @examples {
#'   # load example data
#'   data("test_sounds_est")
#'
#'   # add reference to X
#'   X <- set_reference_sounds(X = test_sounds_est)
#'
#'   envelope_correlation(X = X)
#'
#'   # method 2
#'   # add reference to X
#'   X <- set_reference_sounds(X = test_sounds_est, method = 2)
#'   envelope_correlation(X = X)
#' }
#'
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})
#' @references
#' Araya-Salas, M., Grabarczyk, E. E., Quiroz-Oliva, M., Garcia-Rodriguez, A., & Rico-Guevara, A. (2025). Quantifying degradation in animal acoustic signals with the R package baRulho. Methods in Ecology and Evolution, 00, 1-12. https://doi.org/10.1111/2041-210X.14481
#'
#' Apol, C.A., Sturdy, C.B. & Proppe, D.S. (2017). Seasonal variability in habitat structure may have shaped acoustic signals and repertoires in the black-capped and boreal chickadees. Evol Ecol. 32:57-74.

envelope_correlation <-
  function(X,
           cores = getOption("mc.cores", 1),
           pb = getOption("pb", TRUE),
           cor.method = c("pearson", "spearman", "kendall"),
           env.smooth = getOption("env.smooth", 200),
           hop.size = getOption("hop.size", 11.6),
           wl = getOption("wl", NULL),
           ovlp = getOption("ovlp", 70),
           path = getOption("sound.files.path", ".")) {
    
    # assign a value to cor.method
    cor.method <- rlang::arg_match(cor.method)
    
    # check arguments
    arguments <- as.list(base::match.call())
    
    # add objects to argument names
    for (i in names(arguments)[-1]) {
      arguments[[i]] <- get(i)
    }
    
    # check each arguments
    check_results <-
      .check_arguments(fun = arguments[[1]], args = arguments)
    
    # report errors
    .report_assertions(check_results)
    
    # adjust wl based on hop.size
    wl <- .adjust_wl(wl, X, hop.size, path)
    
    # add sound file selec colums to X (weird column name so it does not overwrite user columns)
    X$.sgnl.temp <- paste(X$sound.files, X$selec, sep = "-")
    
    # get names of envelopes involved (those as test with reference or as reference)
    target_sgnl_temp <-
      unique(c(X$.sgnl.temp[!is.na(X$reference)], X$reference[!is.na(X$reference)]))
    
    
    # set clusters for windows OS
    if (Sys.info()[1] == "Windows" & cores > 1) {
      cl <- parallel::makePSOCKcluster(cores)
    } else {
      cl <- cores
    }
    
    # calculate all envelopes apply function
    envs <-
      warbleR:::.pblapply(
        pbar = pb,
        X = target_sgnl_temp,
        cl = cl,
        message = "computing amplitude envelopes", 
        current = 1,
        total = 2,
        FUN = function(x,
                       ssmth = env.smooth,
                       ovl = ovlp,
                       Q = X,
                       wln = wl,
                       pth = path) {
          .env(
            X = Q,
            y = x,
            env.smooth = ssmth,
            ovlp = ovl,
            wl = wln,
            path = pth
          )
        }
      )
    
    # add sound file selec column as names to envelopes
    names(envs) <- target_sgnl_temp
    
    # calculate all envelops apply function
    envelope_correlation_list <- warbleR:::.pblapply(
      X = seq_len(nrow(X)),
      pbar = pb,
      cl = cl,
      message = "computing envelope correlations", 
      current = 2,
      total = 2,
      FUN =
        function(x,
                 nvs = envs,
                 cm = cor.method,
                 Q = X) {
          .env_cor(X = Q,
                   x,
                   envs = nvs,
                   cor.method = cm)
        }
    )
    
    # unlist
    X$envelope.correlation <- unlist(envelope_correlation_list)
    
    # # remove temporal columns
    X$.sgnl.temp <- NULL
    
    
    # fix call if not a data frame
    if (!is.data.frame(X)) {
      attributes(X)$call <-
        base::match.call()
    } # fix call attribute
    
    return(X)
  }
