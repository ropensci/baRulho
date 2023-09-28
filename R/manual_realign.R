#' Plot spectrograms to check test sound files alignment
#'
#' \code{manual_realign} plots spectrograms to visually inspect alignment precision on test sound files.
#' @usage manual_realign(X, Y, hop.size = getOption("hop.size", 11.6),
#' wl = getOption("wl", NULL), ovlp = getOption("ovlp", 0),
#' path = getOption("sound.files.path", "."), collevels = seq(-120, 0, 5),
#' palette = viridis::viridis, duration = 2, mar = 0.2, step.lengths = c(5, 30),
#' flim = NULL, label.col = "white",  ext.window = TRUE, width = 10, 
#' height = 5, srt = 0, cex = 1, fast.spec = TRUE, 
#' marker = "start_marker", grid = 0.2, ...)
#' @param X Object of class 'data.frame', 'selection_table' or 'extended_selection_table' (the last 2 classes are created by the function \code{\link[warbleR]{selection_table}} from the warbleR package) with the reference to the sounds in the master sound file (typically the output of \code{\link{align_test_files}}). Must contain the following columns: 1) "sound.files": name of the .wav files, 2) "selec": unique selection identifier (within a sound file), 3) "start": start time and 4) "end": end time of selections, 5)  "bottom.freq": low frequency for bandpass, 6) "top.freq": high frequency for bandpass and 7) "sound.id": ID of sounds used to identify counterparts across distances. Each sound must have a unique ID within a distance.
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
#'   # set temporary directory
#'   td <- tempdir()
#'
#'   # load example data
#'   data("master_est")
#'
#'   # save example files in working director to recreate a case in which working
#'   # with sound files instead of extended selection tables.
#'   # This doesn't have to be done with your own data as you will
#'   # have them as sound files already.
#'   for (i in unique(test_sounds_est$sound.files)[1:2]) {
#'     writeWave(object = attr(test_sounds_est, "wave.objects")[[i]], file.path(td, i))
#'   }
#'
#'   # save master file
#'   writeWave(object = attr(master_est, "wave.objects")[[1]], file.path(td, "master.wav"))
#' 
#'   # set path and no progress bar in global options
#'   options(sound.files.path = td, pb = FALSE)
#'
#'   # get marker position
#'   markers <- find_markers(X = master_est, test.files = unique(test_sounds_est$sound.files)[2])
#'
#'   # align all test sounds
#'   alg.tests <- align_test_files(X = master_est, Y = markers, pb = FALSE)
#'   
#'   # add error to alignment
#'   lag <- (as.numeric(as.factor(alg.tests$sound.files)) - 2) / 30
#'   alg.tests$start <- alg.tests$start + lag
#'   alg.tests$end <- alg.tests$end + lag
#'   
#'   \dontrun{
#'   realigned_est <- manual_realign(X = alg.tests, Y = master_df, duration = 2, 
#'   ovlp = 50, hop.size = 14, collevels = seq(-140, 0, 5), palette = mako, 
#'   ext.window = F)
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
      arguments[[i]] <- get(i)
    }

    # check each arguments
    check_results <-
      check_arguments(fun = arguments[[1]], args = arguments)

    # report errors
    report_assertions2(check_results)

    # check that if they are extended selection table must be by song
    if (is_extended_selection_table(X)) {
      if (!attributes(X)$by.song[[1]]) {
        stop2("'X': if is an extended selection table must be created 'by song'")
      }
    }

    if (is_extended_selection_table(Y)) {
      if (!attributes(Y)$by.song[[1]]) {
        stop2("'Y': if is an extended selection table must be created 'by song'")
      }
    }

    # check if markers are found in data
    if (!marker %in% X$sound.id) stop2("'marker' not found in 'sound.id' column in 'X'")
    if (!marker %in% Y$sound.id) stop2("'marker' not found in 'sound.id' column in 'Y'")

    # set external window function
    if (any(Sys.info()[1] == c("Linux", "Windows"))) extwin <- grDevices::X11 else extwin <- grDevices::quartz

    # start external graphic device
    if (ext.window) extwin(width = width, height = height)

    # get sampling rate
    sampling_rate <-
      warbleR::read_sound_file(
        X = X,
        index = 1,
        path = path,
        header = TRUE
      )$sample.rate

    # adjust wl based on hope.size
    if (is.null(wl)) {
      wl <- round(
        sampling_rate * hop.size / 1000,
        0
      )
    }

    # make wl even if odd
    if (!(wl %% 2) == 0) {
      wl <- wl + 1
    }

    # adjust spectrogram margins
    prev_mar <- par("mar")

    # set graph layout
    page_layout <- matrix(c(
      0.03, 0.93, 0.57, 1, # 1) top spectrogram
      0.03, 0.93, 0.14, 0.57, # 2) bottom spectrogram
      0, 0.05, 0.14, 1, # 3) frequency label
      0, 1, 0, 0.14, # 4) time label
      0.93, 1, 0.14, 0.57, # 5) bottom label box
      0.93, 1, 0.57, 1 # 6) top label box
    ), ncol = 4, byrow = TRUE)

    # set color for buttoms
    border_col <- "#A0DFB9CC"
    fill_col <- "#FF993399"

    # # set panel layout
    invisible(close.screen(all.screens = TRUE))
    suppressWarnings(catch <- split.screen(figs = page_layout))

    # reset graphic parameterswhen function is don
    on.exit(par(mar = prev_mar))

    # check duration
    # import sound data
    master_wave_info <- read_wave(
      index = 1,
      Y,
      header = TRUE,
      path = path
    )

    # fix to if higher than sound file duration
    master_dur <- master_wave_info$samples / master_wave_info$sample.rate

    # sto if duration of spectrogram larger than master sound file
    if (master_dur < duration) {
      stop2("'duration' cannot be larger than length of master sound file")
    }

    ## plot master spectrogram ####
    # set start and end of clip to plot
    if (Y$start[Y$sound.id == marker] - mar >= 0) {
      from <- Y$start[Y$sound.id == marker] - mar
    } else {
      from <- 0
    }
    to <- from + duration

    # set margin at the end
    if (to - Y$start[Y$sound.id == marker] < mar) {
      print("asdasd")
      to <- Y$start[Y$sound.id == marker] + mar
      from <- to - duration
    }

    # change to if larger than sound file
    if (to > master_dur) {
      to <- master_dur
      from <- to - duration
    }

    # set margins from the start of marker to be used by spectrograms of test sounds
    left_mar <- Y$start[Y$sound.id == marker] - from
    right_mar <- to - Y$start[Y$sound.id == marker]


    # import sound data
    master_wave <- read_wave(
      Y,
      index = 1,
      from = from,
      to = to,
      path = path
    )

    # frequency label
    screen(3)
    par(mar = c(0, 0, 0, 0), new = TRUE)

    plot(
      1,
      frame.plot = FALSE,
      type = "n",
      yaxt = "n",
      xaxt = "n"
    )

    text(
      x = 0.8,
      y = 1,
      "Frequency (kHz)",
      srt = 90,
      cex = 1.4
    )

    # top box
    screen(6)
    par(
      mar = c(0, 0, 0, 0),
      new = TRUE
    )

    plot(
      1,
      frame.plot = FALSE,
      type = "n",
      yaxt = "n",
      xaxt = "n"
    )

    # add background color
    rect(0, 0, 2, 2, col = "white", border = NA)
    rect(0, 0, 2, 2, col = "#366A9FFF", border = NA)

    top_lab <- vapply(as.character(Y$sound.files[1]), function(x) {if (nchar(x) > 20) paste0(substr(x, 0, 20), "\n", substr(x, 21, nchar(x))) else x}, FUN.VALUE = character(1))

    text(
      x = 1,
      y = 1,
      top_lab,
      cex = 1.2,
      col = "white",
      font = 2,
      srt = 270
    )
    box()

    # top spectrogram
    screen(1)
    par(mar = c(0, 2, 0, 0))

    # set flim (mostly for sound.id labels below)
    if (is.null(flim)) {
      flim <-
        c(0, master_wave@samp.rate / 2000.1)
    } # use 2000.1 to avoid errors at the highest of nyquist frequency

    # plot spectrogram
    warbleR:::spectro_wrblr_int2(
      wave = master_wave,
      collevels = collevels,
      ovlp = ovlp,
      wl = wl,
      flim = flim,
      palette = palette,
      grid = FALSE,
      axisX = FALSE,
      axisY = FALSE,
      flab = "",
      tlab = "",
      fast.spec = fast.spec,
      ...
    )

    # plot grid
    if (grid > 0)
      abline(v = seq(grid, duration(master_wave), by = grid), lty = 4)
    
    # add Freq axis
    at_freq <-
      pretty(seq(0, master_wave@samp.rate / 2000, length.out = 6)[-6], n = 6)
    axis(2,
      at = at_freq,
      labels =
        c("", at_freq[-1])
    )


    # add dotted lines
    abline(
      v = c(Y$start - from, Y$end - from),
      col = label.col,
      lty = 3,
      lwd = 1.5
    )

    # add sound.id labels
    # position of sound.id labels in the freq axis
    y_pos <- flim[2] - 2 * ((flim[2] - flim[1]) / 12)

    # plot sound id labels
    text(
      labels = Y$sound.id,
      x = ((Y$end + Y$start) / 2) - from,
      y = y_pos,
      pos = 3,
      col = label.col,
      srt = srt,
      cex = cex
    )

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

      # get data subset for 1 file
      W <- X[X$sound.files == rerec_files[i], ]

      # bottom box
      screen(5)
      par(
        mar = c(0, 0, 0, 0),
        new = TRUE
      )

      plot(
        1,
        frame.plot = FALSE,
        type = "n",
        yaxt = "n",
        xaxt = "n"
      )

      # add background color
      rect(0, 0, 2, 2, col = "white", border = NA)
      rect(0, 0, 2, 2, col = "#366A9FFF", border = NA)

      bottom_lab <- vapply(as.character(W$sound.files[1]), function(x) {if (nchar(x) > 15) paste0(substr(x, 0, 15), "\n", substr(x, 16, nchar(x))) else x}, FUN.VALUE = character(1))


      # spectrogram title vertical boxes
      text(
        x = 1,
        y = 1,
        bottom_lab,
        cex = 1.2,
        col = "white",
        font = 2,
        srt = 270
      )
      box()

      # plot spectrogram
      screen(2)
      par(mar = c(0, 2, 0, 0))

      from_stp <- W$start[W$sound.id == marker] - left_mar + step_sum
      to_stp <- W$start[W$sound.id == marker] + right_mar + step_sum

      # import sound data
      wave_info <- read_wave(
        index = 1,
        W,
        header = TRUE,
        path = path
      )

      # fix to if higher than sound file duration
      wave_dur <- wave_info$samples / wave_info$sample.rate

      # read wave
      wave <- read_wave(
        W,
        index = 1,
        from = if (from_stp < 0) 0 else from_stp,
        to = if (to_stp > wave_dur) wave_dur else to_stp,
        path = path
      )

      if (from_stp < 0) {
        wave <- seewave::pastew(wave1 = wave, wave2 = tuneR::silence(duration = abs(from_stp), xunit = "time", samp.rate = wave@samp.rate), output = "Wave")
      }

      if (to_stp > wave_dur) {
        wave <- seewave::pastew(wave2 = wave, wave1 = tuneR::silence(duration = to_stp - wave_dur, xunit = "time", samp.rate = wave@samp.rate), output = "Wave")
      }

      # set flim (mostly for sound.id labels below)
      if (is.null(flim)) {
        flim <-
          c(0, wave@samp.rate / 2000.1)
      } # use 2000.1 to avoid errors at the highest of nyquist frequency

      # plot spectrogram
      warbleR:::spectro_wrblr_int2(
        wave = wave,
        collevels = collevels,
        ovlp = ovlp,
        wl = wl,
        flim = flim,
        palette = palette,
        grid = FALSE,
        axisX = FALSE,
        axisY = FALSE,
        flab = "",
        tlab = "",
        fast.spec = fast.spec,
        ...
      )

      # plot grid
      if (grid > 0)
        abline(v = seq(grid, duration(master_wave), by = grid), lty = 4)
      
      # progress bar
      prct <- (i / length(rerec_files))
      y <- grconvertY(y = 0.99, from = "npc", to = "user")
      lines(x = c(0, grconvertX(x = prct - 0.02, from = "npc", to = "user")), y = rep(y, 2), lwd = 7, col = adjustcolor("#E37222", alpha.f = 0.6), xpd = TRUE)
      text(x = grconvertX(x = prct, from = "npc", to = "user"), y = y, xpd = TRUE, labels = paste0(floor(prct * 100), "%"), col = "#E37222", cex = 0.8)

      # add dotted lines on the position of test sounds
      abline(
        v = c(W$start - from_stp + step_sum, W$end - from_stp + step_sum),
        col = label.col,
        lty = 3,
        lwd = 1.5
      )

      # add Freq axis
      at_freq <-
        pretty(seq(0, wave@samp.rate / 2000, length.out = 6)[-6], n = 6)
      axis(2,
        at = at_freq,
        labels = c(at_freq[-length(at_freq)], "")
      )

      # negative fake coordinates to get started
      xy <- list(x = -1000, y = -1000)

      # plot time ticks
      if (i == 1 & step_sum == 0 & xy$x == -1000 & xy$y == -1000) {
        at_time <-
          pretty(seq(0, duration(wave), length.out = 7)[-7], n = 7)
        axis(1,
          at = at_time,
          labels = at_time
        )

        mtext(text = "Time (s)", side = 1, line = 2.2, cex = 1.4)
      }

      # add buttons
      xs <- grconvertX(x = c(0.92, 0.92, 0.99, 0.99), from = "npc", to = "user")
      labels <- c("stop", "long left", "long right", "short left", "short right", "next", "reset")

      # relative position of buttons in y axis
      cpy <- seq(0.93, 0.07, length.out = length(labels))

      mrg <- (cpy[1] - cpy[2]) / 3

      # mid position ofbuttons in c(0, 1) range
      ys <- c(-mrg, mrg, mrg, -mrg)

      grYs <- lapply(seq_len(length(labels)), function(x) {
        grY <- grconvertY(y = cpy[x] - ys, from = "npc", to = "user")
        polygon(x = xs, y = grY, border = border_col, col = fill_col, lwd = 2)
        # "#440154FF" "#482878FF" "#3E4A89FF" "#31688EFF" "#26828EFF" "#1F9E89FF"
        # [7] "#35B779FF" "#6DCD59FF" "#B4DE2CFF" border_col
        # plot symbols
        if (labels[x] == "stop") {
          points(x = mean(xs), y = mean(grY), pch = 15, cex = 1.5, col = border_col)
        }

        if (labels[x] == "long right") {
          text(x = mean(xs), y = mean(grY), labels = ">>", cex = 1.2, font = 2, col = border_col)
        }

        if (labels[x] == "long left") {
          text(x = mean(xs), y = mean(grY), labels = "<<", cex = 1.2, font = 2, col = border_col)
        }

        if (labels[x] == "short right") {
          text(x = mean(xs), y = mean(grY), labels = ">", cex = 1.2, font = 2, col = border_col)
        }

        if (labels[x] == "short left") {
          text(x = mean(xs), y = mean(grY), labels = "<", cex = 1.2, font = 2, col = border_col)
        }

        if (labels[x] == "next") {
          text(x = mean(xs), y = mean(grY), labels = "next", cex = 1.2, font = 2, col = border_col)
        }

        if (labels[x] == "reset") {
          text(x = mean(xs), y = mean(grY), labels = "reset", cex = 1.2, font = 2, col = border_col)
        }

        return(grY)
      })

      names(grYs) <- labels

      # While not in next or stop remain
      while (!(xy$x > min(xs) & xy$x < max(xs) & xy$y > min(grYs$stop) & xy$y < max(grYs$stop)) & !(xy$x > min(xs) & xy$x < max(xs) & xy$y > min(grYs$`next`) & xy$y < max(grYs$`next`))) {
        xy <- locator(n = 1, type = "n")

        # if reset
        if (xy$x > min(xs) & xy$x < max(xs) & xy$y > min(grYs$reset) & xy$y < max(grYs$reset)) {
          step_sum <- 0
          step_sum_vector[i] <- step_sum
          break
        }
        # if next
        if (xy$x > min(xs) & xy$x < max(xs) & xy$y > min(grYs$`next`) & xy$y < max(grYs$`next`)) {
          i <- i + 1
          step_sum <- 0
          step_sum_vector[i] <- step_sum
          break
        }

        # if long left
        if (xy$x > min(xs) & xy$x < max(xs) & xy$y > min(grYs$`long left`) & xy$y < max(grYs$`long left`)) {
          step_sum <- step_sum + max(step.lengths / 1000)
          step_sum_vector[i] <- step_sum
          break
        }

        # if long right
        if (xy$x > min(xs) & xy$x < max(xs) & xy$y > min(grYs$`long right`) & xy$y < max(grYs$`long right`)) {
          step_sum <- step_sum - max(step.lengths / 1000)
          step_sum_vector[i] <- step_sum
          break
        }

        # if short left
        if (xy$x > min(xs) & xy$x < max(xs) & xy$y > min(grYs$`short left`) & xy$y < max(grYs$`short left`)) {
          step_sum <- step_sum + min(step.lengths / 1000)
          step_sum_vector[i] <- step_sum
          break
        }

        # if short right
        if (xy$x > min(xs) & xy$x < max(xs) & xy$y > min(grYs$`short right`) & xy$y < max(grYs$`short right`)) {
          step_sum <- step_sum - min(step.lengths / 1000)
          step_sum_vector[i] <- step_sum
          break
        }

        # if stop
        if (xy$x > min(xs) & xy$x < max(xs) & xy$y > min(grYs$stop) & xy$y < max(grYs$stop)) {
          break
        }
      }

      # stop
      if (xy$x > min(xs) & xy$x < max(xs) & xy$y > min(grYs$stop) & xy$y < max(grYs$stop)) {
        break
      }
    }

    # add white transparent color on top of lower spectrogram
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = adjustcolor("white", 0.5), border = NA)

    # add stopped by user message
    if (xy$x > min(xs) & xy$x < max(xs) & xy$y > min(grYs$stop) & xy$y < max(grYs$stop) & i != length(rerec_files)) {
      text(
        x = mean(par("usr")[1:2]),
        y = mean(par("usr")[3:4]),
        labels = "Stopped by user",
        cex = 1.2,
        col = "black",
        font = 2
      )
    }
    # add last sound file message
    if (xy$x > min(xs) & xy$x < max(xs) & xy$y > min(grYs$`next`) & xy$y < max(grYs$`next`) | xy$x > min(xs) & xy$x < max(xs) & xy$y > min(grYs$stop) & xy$y < max(grYs$stop) & i == length(rerec_files)) {
      text(
        x = mean(par("usr")[1:2]),
        y = mean(par("usr")[3:4]),
        labels = "All sound files were re-aligned",
        cex = 1.2,
        col = "black",
        font = 2
      )
    }
    # adjust time coordinates based on user input
    for (u in rerec_files) {
      X$start[X$sound.files == u] <- X$start[X$sound.files == u] + step_sum_vector[u]
      X$end[X$sound.files == u] <- X$end[X$sound.files == u] + step_sum_vector[u]

      # fix internally for extended selection tables
      if (is_selection_table(X) | is_extended_selection_table(X)) {
        attributes(X)$check.res$start[attributes(X)$check.res$sound.files == u] <- attributes(X)$check.res$start[X$sound.files == u] + step_sum_vector[u]
        attributes(X)$check.res$end[attributes(X)$check.res$sound.files == u] <- attributes(X)$check.res$end[X$sound.files == u] + step_sum_vector[u]

        # fix call attribute
        attributes(X)$call <- base::match.call()
      }
    }

    # close external window
    if (ext.window) {
      Sys.sleep(2)
      dev.off()
    }

    return(X)
  }
