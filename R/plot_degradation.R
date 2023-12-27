#' Save multipanel plots with reference and test sounds
#'
#' \code{plot_degradation} creates multipanel plots (as image files) with reference and test sounds by distance and transect.
#' @usage plot_degradation(X, nrow = 4, env.smooth = getOption("env.smooth", 200),
#' hop.size = getOption("hop.size", 11.6), wl = getOption("wl", NULL),
#' ovlp = getOption("ovlp", 70),  path = getOption("sound.files.path", "."),
#' dest.path = getOption("dest.path", "."), cores = getOption("mc.cores", 1),
#' pb = getOption("pb", TRUE), collevels = seq(-120, 0, 5),
#' palette = viridis::viridis, flim = c("-1", "+1"), envelope = TRUE, spectrum = TRUE,
#' heights = c(4, 1), widths = c(5, 1), margins = c(2, 1), row.height = 2, col.width = 2,
#' cols = viridis::mako(4, alpha = 0.3), res = 120, ...)
#' @param X The output of \code{\link{set_reference_sounds}} which is an object of class 'data.frame', 'selection_table' or 'extended_selection_table' (the last 2 classes are created by the function \code{\link[warbleR]{selection_table}} from the warbleR package) with the reference to the test sounds . Must contain the following columns: 1) "sound.files": name of the .wav files, 2) "selec": unique selection identifier (within a sound file), 3) "start": start time and 4) "end": end time of selections, 5)  "bottom.freq": low frequency for bandpass, 6) "top.freq": high frequency for bandpass, 7) "sound.id": ID of sounds used to identify counterparts across distances and 8) "reference": identity of sounds to be used as reference for each test sound (row). See \code{\link{set_reference_sounds}} for more details on the structure of 'X'.
#' @param nrow Numeric vector of length 1 with the number of rows per image file. Default is 4. This would be dynamically adjusted if more rows than needed are set.
#' @param env.smooth Numeric vector of length 1 determining the length of the sliding window (in amplitude samples) used for a sum smooth for amplitude envelope and power spectrum calculations (used internally by \code{\link[seewave]{env}}). Default is 200.
#' @param hop.size A numeric vector of length 1 specifying the time window duration (in ms). Default is 11.6 ms, which is equivalent to 512 wl for a 44.1 kHz sampling rate. Ignored if 'wl' is supplied.
#' @param wl A numeric vector of length 1 specifying the window length of the spectrogram, default
#' is NULL. If supplied, 'hop.size' is ignored.
#' @param ovlp Numeric vector of length 1 specifying the percent overlap between two
#'   consecutive windows, as in \code{\link[seewave]{spectro}}. Only used when plotting. Default is 70. Applied to both spectra and spectrograms on image files.
#' @param path Character string containing the directory path where the sound files are found. Only needed when 'X' is not an extended selection table.
#' @param dest.path Character string containing the directory path where the image files will be saved. If NULL (default) then the folder containing the sound files will be used instead.
#' @param cores Numeric vector of length 1. Controls whether parallel computing is applied by specifying the number of cores to be used. Default is 1 (i.e. no parallel computing).
#' @param pb Logical argument to control if progress bar is shown. Default is \code{TRUE}.
#' @param collevels	Numeric vector indicating a set of levels which are used to partition the amplitude range of the spectrogram (in dB) as in \code{\link[seewave]{spectro}}. Default is \code{seq(-120, 0, 5)}.
#' @param palette A color palette function to be used to assign colors in the
#'   plot, as in \code{\link[seewave]{spectro}}. Default is \code{\link[viridis]{viridis}}.
#' @param flim A numeric vector of length 2 indicating the highest and lowest frequency limits (kHz) of the spectrogram, as in \code{\link[seewave]{spectro}}. Default is \code{NULL}. Alternatively, a character vector similar to \code{c("-1", "1")} in which the first number is the value to be added to the minimum bottom frequency in 'X' and the second the value to be added to the maximum top frequency in 'X'. This is computed independently for each sound id so the frequency limit better fits the frequency range of the annotated signals. This is useful when test sounds show marked differences in their frequency ranges.
#' @param envelope Logical to control if envelopes are plotted. Default is \code{TRUE}.
#' @param spectrum Logical to control if power spectra are plotted. Default is \code{TRUE}.
#' @param heights Numeric vector of length 2 to control the relative heights of spectrogram (first number) and amplitude envelope (second number) when \code{envelope = TRUE}. Default is c(4, 1).
#' @param widths Numeric vector of length 2 to control the relative widths of spectrogram (first number) and power spectrum (second number) when \code{spectrum = TRUE}. Default is c(5, 1).
#' @param margins Numeric vector of length 2 to control the relative time of the test sound (first number) and adjacent margins (i.e. adjacent background noise, second number) to be included in the spectrogram \code{spectrum = TRUE}. Default is c(2, 1) which means that each margin next to the sound is half the duration of the sound. Note that all spectrograms will have the same time length so margins will be calculated to ensure all spectrograms match the duration of the spectrogram in the longest sound. As such, this argument controls the margin on the longest sound.
#' @param row.height Numeric vector of length 1 controlling the height (in inches) of sound panels in the output image file. Default is 2.
#' @param col.width Numeric vector of length 1 controlling the width (in inches) of sound panels in the output image file. Default is 2.
#' @param cols Character vector of length 4 containing the colors to be used for the background of column and row title panels (element 1), the color of amplitude envelopes (element 2), the color of power spectra (element 3), and the background color of envelopes and spectra (element 4).
#' @param res Numeric argument of length 1. Controls image resolution. Default is 120 (faster) although 300 - 400 is recommended for publication/presentation quality.
#' @param ... Additional arguments to be passed to the internal spectrogram
#' creating function for customizing graphical output. The function is a modified
#' version of \code{\link[seewave]{spectro}}, so it takes the same arguments.
#' @return One ore more image files with a multipanel figure of spectrograms of test sound by distance, sound id and transect.
#' @export
#' @name plot_degradation
#' @details The function aims to simplify the visual inspection of sound degradation by producing multipanel figures (saved in 'dest.path') containing visualizations of each test sound and its reference. Sounds are sorted by distance (columns) and transect (if more than 1). Visualizations include spectrograms, amplitude envelopes and power spectra (the last 2 are optional). Each row includes all the copies of a sound id for a given transect (the row label includes the sound id in the first line and transect in the second line), also including its reference if it comes from another transect. Ambient noise annotations (sound.id 'ambient') are excluded.
#' @examples {
#'   # load example data
#'   data("test_sounds_est")
#'
#'   # order so spectrograms from same sound id as close in the graph
#'   test_sounds_est <- test_sounds_est[order(test_sounds_est$sound.id), ]
#'
#'   # set directory to save image files
#'   options(dest.path = tempdir())
#'
#'   # method 1
#'   Y <- set_reference_sounds(X = test_sounds_est)
#'
#'   # plot degradation spectrograms
#'   plot_degradation(
#'     X = Y, nrow = 3, ovlp = 95
#'   )
#'
#'   # using other color palettes
#'   plot_degradation(
#'     X = Y, nrow = 3, ovlp = 95,
#'     cols = viridis::magma(4, alpha = 0.3),
#'     palette = viridis::magma
#'   )
#'
#'   # missing some data, 2 rows
#'   plot_degradation(
#'     X = Y[-3, ], nrow = 2, ovlp = 95,
#'     cols = viridis::mako(4, alpha = 0.4), palette = viridis::mako, wl = 200
#'   )
#'
#'   # changing marging and high overlap
#'   plot_degradation(X = Y, margins = c(5, 1), nrow = 6, ovlp = 95)
#'
#'   # more rows than needed (will adjust it automatically)
#'   plot_degradation(X = Y, nrow = 10, ovlp = 90)
#' }
#'
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})
#' @family quantify degradation
#' @seealso \code{\link{blur_ratio}}, \code{\link{plot_aligned_sounds}}, \code{\link{plot_degradation}}
#' @references {
#' Araya-Salas, M. (2020). baRulho: baRulho: quantifying degradation of (animal) acoustic signals in R. R package version 1.0.2
#' }
#'
plot_degradation <-
  function(X,
           nrow = 4,
           env.smooth = getOption("env.smooth", 200),
           hop.size = getOption("hop.size", 11.6),
           wl = getOption("wl", NULL),
           ovlp = getOption("ovlp", 70),
           path = getOption("sound.files.path", "."),
           dest.path = getOption("dest.path", "."),
           cores = getOption("mc.cores", 1),
           pb = getOption("pb", TRUE),
           collevels = seq(-120, 0, 5),
           palette = viridis::viridis,
           flim = c("-1", "+1"),
           envelope = TRUE,
           spectrum = TRUE,
           heights = c(4, 1),
           widths = c(5, 1),
           margins = c(2, 1),
           row.height = 2,
           col.width = 2,
           cols = viridis::mako(4, alpha = 0.3),
           res = 120,
           ...) {
    # check arguments
    arguments <- as.list(base::match.call())
    
    # add objects to argument names
    for (i in names(arguments)[-1]) {
      try(arguments[[i]] <- get(i), silent = TRUE)
    }
    
    # check each arguments
    check_results <-
      .check_arguments(fun = arguments[[1]], args = arguments)
    
    # report errors
    .report_assertions(check_results)
    
    # set colors
    bg_titles <- cols[1]
    env_fill <- cols[2]
    spc_fill <- cols[3]
    bg_sp_env <- cols[4]
    
    # adjust wl based on hop.size
    wl <- .adjust_wl(wl, X, hop.size, path)
    
    # set clusters for windows OS
    if (Sys.info()[1] == "Windows" & cores > 1) {
      cl <- parallel::makePSOCKcluster(getOption("cl.cores", cores))
    } else {
      cl <- cores
    }
    
    # remove ambient sounds
    X <- X[X$sound.id != "ambient",]
    
    # stop if more than 1 sample per distance but no transect info
    if (any(table(X$sound.id, X$distance)) > 1 &
        is.null(X$transect)) {
      .stop(
        "There are more than 1 test sound per sound.id/distance combination but no 'transect' column to group by transect"
      )
    }
    
    # stop if there are more than 1 sample for distance sound id and transect combination (only 1 samp)
    if (!is.null(X$transect)) {
      if (any(table(X$sound.id, X$distance, X$transect)) > 1) {
        .stop("There is more than 1 test sound per sound.id/distance/transect combination")
      }
    }
    
    if (is.null(X$transect)) {
      transects <- X$transect <- 1
    } else {
      transects <- unique(X$transect)
    }
    
    # add sound file selec colums to X (weird column name so it does not overwrite user columns)
    X$.sgnl.temp <- paste(X$sound.files, X$selec, sep = "-")
    
    # for each sound fix start and end based on margin
    X <- .fix_margins(X, i, path, margins)
    
    # get unique distances and maximum number of columns (distances) for any sound id
    distances <- sort(unique(X$distance))
    ncol <- length(distances)
    
    soundid_X <- .prep_data_plot_degrad(X, distances, nrow)
    
    # adjust nrow if less sound.id.seq than nrow
    if (nrow > length(unique(soundid_X$sound.id.seq))) {
      nrow <- length(unique(soundid_X$sound.id.seq))
    }
    
    # set image size fixing height for when only 1 row
    img_width <- ncol * col.width
    img_heigth <- nrow * row.height * if (nrow == 1) {
      1.4
    } else {
      1
    }
    
    # set basic layout for all pages
    page_layout <- .page_layout(ncol, nrow, spectrum, envelope, heights, widths)
    
    close.screen(all.screens = TRUE)
    
    # set clusters for windows OS
    if (Sys.info()[1] == "Windows" & cores > 1) {
      cl <- parallel::makePSOCKcluster(getOption("cl.cores", cores))
    } else {
      cl <- cores
    }
    
    out <-
      warbleR:::pblapply_wrblr_int(
        X = sort(unique(soundid_X$page)),
        pbar = pb,
        cl = cl,
        FUN = .plot_degrad,
        X,
        soundid_X,
        flim,
        path,
        dest.path,
        img_width,
        img_heigth,
        res,
        page_layout,
        nrow,
        ncol,
        envelope,
        spectrum,
        distances,
        wl,
        spc_fill,
        bg_sp_env,
        bg_titles,
        ovlp,
        collevels, 
        env.smooth,
        palette,
        ...
      )
  }
