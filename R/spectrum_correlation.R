#' Measure frequency spectrum correlation
#'
#' \code{spectrum_correlation} measures frequency spectrum correlation of sounds referenced in an extended selection table.
#' @inheritParams template_params
#' @param spec.smooth Numeric vector of length 1 determining the length of the sliding window used for a sum smooth for power spectrum calculation (in kHz). Default is 5.
#' @param ovlp Numeric vector of length 1 specifying the percentage of overlap between two
#'   consecutive windows, as in \code{\link[seewave]{spectro}}. Default is 70.
#' @param n.bins Numeric vector of length 1 specifying the number of frequency bins to use for representing power spectra. Default is 100. If null the raw power spectrum is used (note that this can result in high RAM memory usage for large data sets). Power spectrum values are interpolated using \code{\link[stats]{approx}}.
#' @return Object 'X' with an additional column,  'spectrum.correlation', containing the id of the sound used as reference and the computed frequency spectrum correlation coefficients, respectively.
#' @export
#' @name spectrum_correlation
#' @details spectral correlation measures the similarity of two sounds in the frequency domain. The function measures the spectral correlation coefficients of sounds in which a reference playback has been re-recorded at increasing distances. Values range from 1 (identical frequency spectrum, i.e. no degradation) to 0. The 'sound.id' column must be used to indicate the function to only compare sounds belonging to the same category (e.g. song-types). The function will then compare each sound to the corresponding reference sound. Two methods for computing spectral correlation are provided (see 'method' argument). The function uses \code{\link[seewave]{meanspec}} internally to compute power spectra. Use \code{\link{spectrum_blur_ratio}} to extract raw spectra values. NA is returned if at least one the power spectra cannot be computed.
#' @examples {
#'   # load example data
#'   data("test_sounds_est")
#'
#'   # method 1
#'   # add reference column
#'   Y <- set_reference_sounds(X = test_sounds_est)
#'
#'   # run spectrum correlation
#'   spectrum_correlation(X = Y)
#'
#'   # method 2
#'   Y <- set_reference_sounds(X = test_sounds_est, method = 2)
#'   # spectrum_correlation(X = Y)
#' }
#'
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})
#' @family quantify degradation
#' @seealso \code{\link{envelope_correlation}}, \code{\link{spectrum_blur_ratio}}
#' @references {
#' Araya-Salas, M. (2020). baRulho: baRulho: quantifying degradation of (animal) acoustic signals in R. R package version 1.0.2
#'
#' Apol, C.A., Sturdy, C.B. & Proppe, D.S. (2017). Seasonal variability in habitat structure may have shaped acoustic signals and repertoires in the black-capped and boreal chickadees. Evol Ecol. 32:57-74.
#' }

spectrum_correlation <-
  function(X,
           cores = getOption("mc.cores", 1),
           pb = getOption("pb", TRUE),
           cor.method = c("pearson", "spearman", "kendall"),
           spec.smooth = getOption("spec.smooth", 5),
           hop.size = getOption("hop.size", 11.6),
           wl = getOption("wl", NULL),
           ovlp = getOption("ovlp", 70),
           path = getOption("sound.files.path", "."),
           n.bins = 100) {
    
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
      cl <-
        parallel::makePSOCKcluster(getOption("cl.cores", cores))
    } else {
      cl <- cores
    }
    
    # print message
    if (pb) {
      write(file = "", x = "Computing power spectra (step 1 out of 2):")
    }
    
    # calculate all spectra apply function
    specs <-
      warbleR:::pblapply_wrblr_int(
        pbar = pb,
        X = target_sgnl_temp,
        cl = cl,
        FUN = function(x,
                       ssmth = spec.smooth,
                       wln = wl,
                       Q = X,
                       pth = path,
                       ovl = ovlp,
                       nb = n.bins) {
          .spctr(
            y = x,
            spec.smooth = ssmth,
            wl = wln,
            X = Q,
            meanspc = TRUE,
            path = pth,
            ovlp = ovl,
            n.bins = nb
          )
        }
      )
    
    # add sound file selec names to envelopes (weird column name so it does not overwrite user columns)
    names(specs) <- target_sgnl_temp
    
    if (pb) {
      write(file = "", x = "Computing spectrum correlations (step 2 out of 2):")
    }
    
    # calculate all envelops apply function
    X$spectrum.correlation <-
      unlist(warbleR:::pblapply_wrblr_int(
        X = seq_len(nrow(X)),
        pbar = pb,
        cl = cl,
        FUN =
          function(x,
                   spcs = specs,
                   cm = cor.method,
                   Q = X) {
            .spctr_cor(
              y = x,
              specs = spcs,
              X = Q,
              cor.method = cm
            )
          }
      ))
    
    # remove temporal columns
    X$.sgnl.temp <- NULL
    
    # fix call if not a data frame
    if (!is.data.frame(X)) {
      attributes(X)$call <-
        base::match.call()
    } # fix call attribute
    
    
    return(X)
  }
