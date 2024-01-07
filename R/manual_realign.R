#' Plot spectrograms to check test sound files alignment
#'
#' \code{manual_realign} plots spectrograms to visually inspect alignment precision on test sound files.
#' @param X Object of class 'data.frame', 'selection_table' or 'extended_selection_table' (the last 2 classes are created by the function \code{\link[warbleR]{selection_table}} from the warbleR package) with the reference to the test sounds (typically the output of \code{\link{align_test_files}}). Must contain the following columns: 1) "sound.files": name of the .wav files, 2) "selec": unique selection identifier (within a sound file), 3) "start": start time and 4) "end": end time of selections, 5)  "bottom.freq": low frequency for bandpass, 6) "top.freq": high frequency for bandpass and 7) "sound.id": ID of sounds used to identify counterparts across distances. Each sound must have a unique ID within a distance.
#' @param Y object of class 'data.frame', 'selection_table' or 'extended_selection_table' (the last 2 classes are created by the function \code{\link[warbleR]{selection_table}} from the warbleR package) with the master sound file annotations. This should be the same data than that was used for finding the position of markers in \code{\link{find_markers}}. It should also contain a 'sound.id' column.
#' @param hop.size A numeric vector of length 1 specifying the time window duration (in ms). Default is 1 1. 6ms, which is equivalent to 512 wl for a 44.1 kHz sampling rate. Ignored if 'wl' is supplied.
#' @param wl A numeric vector of length 1 specifying the window length of the spectrogram, default is NULL. Ignored if \code{bp = NULL}. If supplied, 'hop.size' is ignored.
#' @param ovlp Numeric vector of length 1 specifying the percent overlap between two consecutive windows, as in \code{\link[seewave]{spectro}}. Default is 0.
#' @param path Character string containing the directory path where the sound files are found. Only needed when 'X' is not an extended selection table.
#' @param collevels A numeric vector of length 3. Specifies levels to partition the amplitude range of the spectrogram (in dB). The more levels the higher the resolution of the spectrogram. Default is seq(-120, 0, 1).
#' @param palette Color palette function for spectrogram. Default is  \code{\link[viridis]{viridis}}. See \code{\link[seewave]{spectro}} for more palettes. Palettes as \code{\link[monitoR:specCols]{gray.2}} may work better when \code{fast.spec = TRUE}.
#' @param duration A numeric vector of length 1. Specifies the overall duration of the clip that will be plotted. Notice that only the initial part of the test files are plotted as this is enough to tell the precision of the alignment.
#' @param mar numeric vector of length 1. Specifies the minimum margins adjacent (before and after) to the start of the marker used for checking alignments (see 'marker' argument). Default is 0.2.
#' @param step.lengths Numeric vector of length 2 indicating the time length (in ms) of short (min(step.lengths)) and long steps (max(step.lengths)) for manually aligning spectrograms. Default is \code{c(5, 30)}.
#' @param flim A numeric vector of length 2 indicating the highest and lowest frequency limits (kHz) of the spectrogram, as in \code{\link[seewave]{spectro}}. Default is \code{NULL} which will plot spectrograms in the full frequency range (0 - nyquist frequency).
#' @param label.col Character string controlling the color of lines and sound ID labels.
#' @param ext.window Logical. If \code{TRUE} then and external graphic window is used.Dimensions can be set using the 'width' and 'height' arguments. Default is \code{TRUE}.
#' @param width Numeric vector of length 1. Single value (in inches) indicating the width of the output image files. Default is 10.
#' @param height Numeric vector of length 1. Single value (in inches) indicating the height of the output image files. Default is 5.
#' @param srt Numeric argument of length 1. The rotation (in degrees) of the sound id labels. Default is 0.
#' @param cex Numeric argument of length 1controlling the size of sound id text labels. Default is 1.
#' @param fast.spec Logical. If \code{TRUE} then image function is used internally to create spectrograms, which substantially increases performance (much faster), although some options become unavailable, as collevels (amplitude scale). Default is \code{FALSE}.
#' @param marker Character string with the name of the marker to be used as the main reference for checking/adjusting time alignments. Default is 'start_marker'. Note that this can take any of the sound IDs in 'Y$sound.id'.
#' @param grid Numeric vector of length 1 controlling the spacing between vertical lines on the spectrogram. Default is 0.2 s. Use 0 to remove grid.
#' @param ... Additional arguments to be passed to the internal spectrogram
#' creating function for customizing graphical output. The function is a modified
#' version of \code{\link[seewave]{spectro}}, so it takes the same arguments.
#' @return Creates a multipanel graph with spectrograms of master and test sound files in which users can interactively adjust their alignment in time. Return an object similar to the input object 'X' in which the start and end of the sounds have been adjusted.
#' @export
#' @name manual_realign
#' @details This function allows the interactive adjustment of the alignment of test sound files produced by \code{\link{align_test_files}}. The function generates a multipanel graph with the spectrogram of the master sound file in top of that from test sound files, highlighting the position of correspondent test sounds on both in order to simplify assessing and adjusting their alignment. Spectrograms include the first few seconds of the sound files (controlled by 'duration') which is usually enough to tell the precision of the alignment. The lower spectrogram shows a series of 'buttons' that users can click on to control if the test sound file spectrogram (low panel) needs to be moved to the left ("<") or right (">"). Users can also reset the spectrogram to its original position ('reset'), move on to the next sound file in 'X' (test sound file annotations) or stop the process (stop button). The function returns an object similar to the input object 'X' in which the start and end of the sounds have been adjusted.
#' @family test sound alignment
#' @seealso \code{\link{auto_realign}}; \code{\link{find_markers}}; \code{\link{align_test_files}}
#' @examples
#' {
#'   # load example data
#'   data("master_est")
#'
#'   # save example files in working director to recreate a case in which working
#'   # with sound files instead of extended selection tables.
#'   # This doesn't have to be done with your own data as you will
#'   # have them as sound files already.
#'   for (i in unique(test_sounds_est$sound.files)[1:2]) {
#'   writeWave(object = attr(test_sounds_est, "wave.objects")[[i]], file.path(tempdir(), i))
#'   }
#'
#'   # save master file
#'   writeWave(object = attr(master_est, "wave.objects")[[1]], file.path(tempdir(), "master.wav"))
#'
#'   # get marker position
#'   markers <- find_markers(X = master_est, test.files = unique(test_sounds_est$sound.files)[2], 
#'   path = tempdir())
#'
#'   # align all test sounds
#'   alg.tests <- align_test_files(X = master_est, Y = markers)
#'
#'   # add error to alignment
#'   lag <- (as.numeric(as.factor(alg.tests$sound.files)) - 2) / 30
#'   alg.tests$start <- alg.tests$start + lag
#'   alg.tests$end <- alg.tests$end + lag
#'
#'   if(interactive()){
#'   realigned_est <- manual_realign(X = alg.tests, Y = master_est, duration = 2,
#'   ovlp = 50, hop.size = 14, collevels = seq(-140, 0, 5), palette = viridis::mako,
#'   ext.window = FALSE)
#'  }
#' }
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})
#' @references {
#' Araya-Salas, M. (2020). baRulho: baRulho: quantifying degradation of (animal) acoustic signals in R. R package version 1.0.2
#' }

manual_realign <-
  function(X,
           Y,
           hop.size = getOption("hop.size", 11.6),
           wl = getOption("wl", NULL),
           ovlp = getOption("ovlp", 0),
           path = getOption("sound.files.path", "."),
           collevels = seq(-120, 0, 5),
           palette = viridis::viridis,
           duration = 2,
           mar = 0.2,
           step.lengths = c(5, 30),
           flim = NULL,
           label.col = "white",
           ext.window = TRUE,
           width = 10,
           height = 5,
           srt = 0,
           cex = 1,
           fast.spec = TRUE,
           marker = "start_marker",
           grid = 0.2,
           ...) {
    # stop if not run interactively
    if (!interactive()) {
      stop("This function requires an interactive graphic device.")
    }
    
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
    
    # check if markers are found in data
    if (!all(sapply(unique(X$sound.files), function(x)
      sum(X$sound.files == x &
          X$sound.id == marker))))
      .stop("at least 1 sound file does not have 'marker' in 'sound.id' column in 'X'")
    if (!marker %in% Y$sound.id)
      .stop("'marker' not found in 'sound.id' column in 'Y'")
    
    # set external window function
    if (any(Sys.info()[1] == c("Linux", "Windows")))
      extwin <- grDevices::X11 else
      extwin <- grDevices::quartz
    
    # start external graphic device
    if (ext.window)
      extwin(width = width, height = height)
    
    # adjust wl based on hop.size
    wl <- .adjust_wl(wl, X, hop.size, path)
    
    # adjust spectrogram margins
    prev_mar <- par("mar")
    
    # set graph layout
    page_layout <- matrix(
      c(
        0.03,
        0.93,
        0.57,
        1,
        # 1) top spectrogram
        0.03,
        0.93,
        0.14,
        0.57,
        # 2) bottom spectrogram
        0,
        0.05,
        0.14,
        1,
        # 3) frequency label
        0,
        1,
        0,
        0.14,
        # 4) time label
        0.93,
        1,
        0.14,
        0.57,
        # 5) bottom label box
        0.93,
        1,
        0.57,
        1 # 6) top label box
      ),
      ncol = 4,
      byrow = TRUE
    )
    
    # set color for buttoms
    border_col <- "#A0DFB9CC"
    fill_col <- "#FF993399"
    
    # # set panel layout
    invisible(close.screen(all.screens = TRUE))
    suppressWarnings(catch <- split.screen(figs = page_layout))
    
    # reset graphic parameterswhen function is don
    on.exit(par(mar = prev_mar))
    
    # output from master panel
    omp <- .master_panel(Y,
                         path,
                         duration,
                         marker,
                         mar,
                         ovlp,
                         wl,
                         flim,
                         palette,
                         fast.spec,
                         grid,
                         collevels,
                         label.col,
                         srt,
                         cex,
                         ...)
    
    # get name of files in X
    rerec_files <- unique(X$sound.files)
    
    # empty vector to save time change
    step_sum_vector <- vector(length = length(rerec_files))
    names(step_sum_vector) <- rerec_files
    
    # start with first file
    i <- 1
    
    # how to add to each time coordinate (start with nothing)
    step_sum <- 0
    
    # save plot
    prev.plot <- recordPlot()
    
    # negative fake coordinates to get started
    xy <- list(x = -1000, y = -1000)
    
    # loop over re-recorded files
    while (i <= length(rerec_files)) {
      # plot saved plot only if not the first time the while loop runs
      if (i == 1 & step_sum == 0 & xy$x != -1000 & xy$y != -1000) {
        prev.plot
      }
      
      # plot lower panel
      rp <- .responsive_panel(
        X,
        rerec_files,
        i,
        marker,
        step_sum,
        path,
        omp,
        flim,
        collevels,
        ovlp,
        wl,
        palette,
        fast.spec,
        grid,
        label.col,
        border_col,
        fill_col,
        ...
      )
      
      xy <- rp$xy
      
      # While not in next or stop remain
      onn <- .not_next(xy, rp, step_sum, step_sum_vector, step.lengths, i)
      
      xy <- onn$xy
      step_sum <- onn$step_sum
      step_sum_vector <- onn$step_sum_vector
      i <- onn$i
      
      # stop
      if (xy$x > min(rp$xs) &
          xy$x < max(rp$xs) &
          xy$y > min(rp$grYs$stop) & xy$y < max(rp$grYs$stop)) {
        break
      }
    }
    
    # add white transparent color on top of lower spectrogram
    rect(
      par("usr")[1],
      par("usr")[3],
      par("usr")[2],
      par("usr")[4],
      col = adjustcolor("white", 0.5),
      border = NA
    )
    
    # add "stopped by user" and "last sound file" messages
    .stop_by_user(i, xy, rp$xs, rp$grYs, rerec_files)
    
    # adjust time coordinates based on user input
    X <- .adjust_coors(rerec_files, X, step_sum_vector, arg_call = base::match.call())
    
    # close external window
    if (ext.window) {
      Sys.sleep(2)
      dev.off()
    }
    
    return(X)
  }
