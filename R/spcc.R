#' Measure spectrographic cross-correlation as a measure of sound distortion
#'
#' \code{spcc} measures spectrographic cross-correlation as a measure of sound distortion in sounds referenced in an extended selection table.
#' @inheritParams template_params
#' @param ovlp Numeric vector of length 1 specifying \% of overlap between two
#' consecutive windows, as in \code{\link[seewave]{spectro}}. Default is 90. High values of ovlp
#' slow down the function but produce more accurate results. Can be set globally for the current R session via the "ovlp" option (see \code{\link[base]{options}}).
#' @return Object 'X' with an additional column, 'cross.correlation', containing the computed spectrogram cross-correlation coefficients.
#' @export
#' @name spcc
#' @details Spectrographic cross-correlation measures frequency distortion of sounds as a similarity metric. Values close to 1 means very similar spectrograms (i.e. little sound distortion has occurred). Cross-correlation is measured of sounds in which a reference playback has been re-recorded at increasing distances. The 'sound.id' column must be used to indicate the function to only compare sounds belonging to the same category (e.g. song-types). The function compares each sound to the corresponding reference sound within the supplied frequency range (e.g. bandpass) of the reference sound ('bottom.freq' and 'top.freq' columns in 'X'). Two methods for computing cross-correlation are provided (see 'method' argument). The function is a wrapper on warbleR's \code{\link[warbleR]{cross_correlation}} function.
#' @examples {
#'   # load example data
#'   data("test_sounds_est")
#'
#'   # add reference to X
#'   X <- set_reference_sounds(X = test_sounds_est)
#'
#'   # get spcc
#'   spcc(X = X)
#' }
#'
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})
#' @family quantify degradation
#' @seealso \code{\link{blur_ratio}}, \code{\link{manual_realign}}, \code{\link[warbleR]{cross_correlation}}
#' @references {
#' Araya-Salas, M., Grabarczyk, E. E., Quiroz-Oliva, M., Garcia-Rodriguez, A., & Rico-Guevara, A. (2025). Quantifying degradation in animal acoustic signals with the R package baRulho. Methods in Ecology and Evolution, 00, 1-12. https://doi.org/10.1111/2041-210X.14481
#'
#' Clark, C.W., Marler, P. & Beeman K. (1987). Quantitative analysis of animal vocal phonology: an application to Swamp Sparrow song. Ethology. 76:101-115.
#' }

spcc <-
  function(X,
           cores = getOption("mc.cores", 1),
           pb = getOption("pb", TRUE),
           cor.method = c("pearson", "spearman", "kendall"),
           hop.size = getOption("hop.size", 11.6),
           wl = getOption("wl", NULL),
           ovlp = getOption("ovlp", 90),
           wn = "hanning",
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
    
    # # put together in a single
    comp_mat <- cbind(X$.sgnl.temp, X$reference)
    
    # remove NA rows
    comp_mat <- comp_mat[stats::complete.cases(comp_mat), , drop = FALSE]
    
    
    # save previous warbleR options
    prev_wl <- .Options$warbleR
    
    on.exit(
      warbleR_options(
        wl = prev_wl$wl,
        ovlp = prev_wl$ovlp,
        wn = prev_wl$wn,
        parallel = prev_wl$parallel,
        pb = prev_wl$pb
      )
    )
    
    # steps for warbleR message
    options("int_warbleR_steps" = list(current = 0, total = 1))
    
    on.exit(options("int_warbleR_steps" = list(current = 0, total = 0)), add = TRUE)
    
    warbleR_options(
      wl = wl,
      ovlp = ovlp,
      wn = wn,
      parallel = cores,
      pb = pb,
      compare.matrix = comp_mat
    )
    
    # run spcc
    xcorrs <-
      cross_correlation(X = X,
                                 cor.method = "pearson",
                                 path = path)$max.xcorr.matrix
    
    # put results back into X
    X$cross.correlation <- NA
    
    # fill score values on X
    X$cross.correlation <- vapply(seq_len(nrow(X)), function(x) {
      # get score for each row
      xc <-
        xcorrs$score[xcorrs$X1 == X$.sgnl.temp[x] &
                       xcorrs$X2 == X$reference[x]]
      
      # if empty then NA
      if (length(xc) == 0) {
        xc <- NA
      }
      
      return(xc)
    }, FUN.VALUE = numeric(1L))
    
    # fix call if not a data frame
    if (!is.data.frame(X)) {
      attributes(X)$call <-
        base::match.call()
    } # fix call attribute
    
    # remove temporary colu8mn
    X$.sgnl.temp <- NULL
    
    return(X)
  }
