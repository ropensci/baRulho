#' Create a master sound file
#'
#' \code{master_sound_file} creates a master sound file to be used in playback experiments related to sound degradation.
#' @usage master_sound_file(X, file.name, dest.path = getOption("dest.path", "."),
#' overwrite = FALSE, delay = 1, gap.duration = 1, amp.marker = 2, flim = c(0, 4), 
#' cex = 14, path = getOption("sound.files.path", "."))
#' @param X Object of class 'data.frame', 'selection_table' or 'extended_selection_table' (the last 2 classes are created by the function \code{\link[warbleR]{selection_table}} from the warbleR package) with the reference to the sounds in the master sound file. Must contain the following columns: 1) "sound.files": name of the .wav files, 2) "selec": unique selection identifier (within a sound file), 3) "start": start time and 4) "end": end time of selections, 5)  "bottom.freq": low frequency for bandpass and 6) "top.freq": high frequency for bandpass. An optional 'sound.id' column can be included to use a custom label for each sound in the output. This column must contain a unique id for each sound (labels cannot repeated). If not supplied the function will make it by combining the sound file and selection columns.
#' @param file.name Character string indicating the name of the sound file.
#' @param dest.path Character string containing the directory path where the sound file will be saved.
#' If \code{NULL} (default) then the current working directory will be used instead.
#' @param overwrite Logical argument to determine if the function will overwrite any existing sound file with the same file name. Default is the current working directory.
#' @param delay Numeric vector of length 1 to control the duration (in s) of a silence gap at the beginning (and at the end) of the sound file. This can be useful to allow some time at the start of the playback experiment. Default is 1.
#' @param gap.duration Numeric vector of length 1 to control the duration (in s) of silence gaps to be placed in between sounds. Default is 1 s.
#' @param amp.marker Numeric vector of length 1 to use as a constant to amplify markers amplitude. This is useful to increase the amplitude of markers in relation to those of sounds, so it is picked up at further distances. Default is 2.
#' @param flim Numeric vector of length 2 to control the frequency range in which the markers would be found. If \code{NULL} markers would be display across the whole frequency range. Default is c(0, 4).
#' @param cex Numeric vector of length 1 indicating the font size for the start and end markers. Default is 14.
#' @param path Character string containing the directory path where the sound files are found. Only needed when 'X' is not an extended selection table.
#' @return A .wav file in 'path' as well as a data frame in the R environment with the annotations (i.e. time position) of sounds in the master sound file and an additional column 'sound.id' that provides a unique id for each sound in the sound file. This is useful for identifying/labeling sounds in test (re-recorded) sound files for downstream analyses.
#' @export
#' @name master_sound_file
#' @details The function is intended to simplify the creation of master sound files for playback experiments in sound degradation studies. The function clips sounds from sound files (or wave objects from extended selection tables) and concatenates them in a single sound file. The function also adds acoustic markers at the start and end of the playback that can be used to time-sync test (re-recorded) sounds to facilitate the streamlining of degradation quantification.
#' @examples {
#'   # load example data from warbleR
#'   data(list = c(
#'     "Phae.long1", "Phae.long2", "Phae.long3", "Phae.long4",
#'     "lbh_selec_table"
#'   ))
#'
#'   # save sound files to temporary folder
#'   writeWave(Phae.long1, file.path(tempdir(), "Phae.long1.wav"))
#'   writeWave(Phae.long2, file.path(tempdir(), "Phae.long2.wav"))
#'   writeWave(Phae.long3, file.path(tempdir(), "Phae.long3.wav"))
#'   writeWave(Phae.long4, file.path(tempdir(), "Phae.long4.wav"))
#'
#'   # make an extended selection table
#'   est <- selection_table(
#'     X = lbh_selec_table, extended = TRUE, confirm.extended = FALSE,
#'     path = tempdir()
#'   )
#'
#'   # create master sound file
#'   master.sel.tab <- master_sound_file(
#'     X = est, file.name = "example_master",
#'     dest.path = tempdir(), gap.duration = 0.3
#'   )
#'
#' \dontrun{
#'   # the following code exports the selection table to Raven
#'   # using the Rraven package
#'   Rraven::exp_raven(master.sel.tab, path = tempdir(),
#'   file.name = "example_master_selection_table")
#'   }
#' }
#'
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})
#' @seealso \code{\link[Rraven]{exp_raven}}
#' @references {
#' Araya-Salas, M. (2020). baRulho: baRulho: quantifying degradation of (animal) acoustic signals in R. R package version 1.0.2
#' }

master_sound_file <-
  function(X,
           file.name,
           dest.path = getOption("dest.path", "."),
           overwrite = FALSE,
           delay = 1,
           gap.duration = 1,
           amp.marker = 2,
           flim = c(0, 4),
           cex = 14,
           path = getOption("sound.files.path", ".")) {
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

    # get sampling rate
    sampling_rate <-
      warbleR::read_sound_file(
        X = X,
        index = 1,
        path = path,
        header = TRUE
      )$sample.rate


    # check if ghost script is installed
    gsexe <- tools::find_gs_cmd()

    # warning if ghostscript not found
    if (!nzchar(gsexe)) {
      message2(paste(cli::cli_text(colortext(
        "warning: GhostScript was not found. It produces clearer start and end marker. You can download it from {.url https://ghostscript.com}. ",
        as = "red"
      )), colortext("Using 'png()' instead.", as = "red")))
    }

    # set frequency range for markers
    if (is.null(flim)) {
      flim <- c(0, attr(X, "check.results")$sample.rate[1] / 2)
    }

    # save par settings
    oldpar <- par(no.readonly = TRUE)
    on.exit(try(suppressWarnings(par(oldpar)), silent = TRUE))

    # save image of start marker in temporary directory
    if (!nzchar(gsexe)) {
      grDevices::png(
        filename = file.path(tempdir(), "strt_mrkr-img.png"),
        pointsize = 10
      )
    }

    # remove margins in graphic device
    par(mar = rep(0, 4))

    # empty plot
    plot(
      0,
      type = "n",
      axes = FALSE,
      ann = FALSE,
      xlim = c(0, 1),
      ylim = c(0, 1)
    )

    # add text
    text(
      x = 0.5,
      y = 0.5,
      labels = "*start+",
      cex = cex,
      font = 2
    )

    # save image of start marker in temporary directory
    if (nzchar(gsexe)) {
      grDevices::dev2bitmap(file.path(tempdir(), "strt_mrkr-img.png"),
        type = "pngmono",
        res = 30
      )
    }

    # close graph
    dev.off()

    # make image of start marker
    strt_mrkr <-
      warbleR::image_to_wave(
        file = file.path(tempdir(), "strt_mrkr-img.png"),
        plot = FALSE,
        flim = flim,
        samp.rate = sampling_rate / 1000
      )

    # remove image file
    unlink(file.path(tempdir(), "strt_mrkr-img.png"))

    # save image of end marker in temporary directory
    if (!nzchar(gsexe)) {
      grDevices::png(
        filename = file.path(tempdir(), "end_mrkr-img.png"),
        pointsize = 10
      )
    }

    # remove margins in graphic device
    par(mar = rep(0, 4))

    # empty plot
    plot(
      0,
      type = "n",
      axes = FALSE,
      ann = FALSE,
      xlim = c(0, 1),
      ylim = c(0, 1)
    )

    # add text
    text(
      x = 0.5,
      y = 0.5,
      labels = "+end*-",
      cex = cex,
      font = 2
    )

    # save image of end marker in temporary directory
    if (nzchar(gsexe)) {
      dev2bitmap(file.path(tempdir(), "end_mrkr-img.png"),
        type = "pngmono",
        res = 30
      )
    }

    # close graph
    dev.off()

    # conver to wave both
    end_mrkr <-
      warbleR::image_to_wave(
        file = file.path(tempdir(), "end_mrkr-img.png"),
        plot = FALSE,
        flim = flim,
        samp.rate = sampling_rate / 1000
      )

    # remove image file
    unlink(file.path(tempdir(), "end_mrkr-img.png"))

    # remove plots
    nll <- try(dev.off(), silent = TRUE)

    # output wave object
    strt_mrkr <- tuneR::normalize(strt_mrkr)
    end_mrkr <- tuneR::normalize(end_mrkr)

    # frequency range of markers
    strt_mrkr_freq <-
      warbleR::freq_range_detec(strt_mrkr, fsmooth = 0.2, plot = FALSE)
    end_mrkr_freq <-
      warbleR::freq_range_detec(end_mrkr, fsmooth = 0.2, plot = FALSE)

    # amplify markers
    strt_mrkr@left <- strt_mrkr@left * amp.marker
    end_mrkr@left <- end_mrkr@left * amp.marker

    # save duration of markers for creating selection table
    dur_strt_mrkr <- seewave::duration(strt_mrkr)
    dur_end_mrkr <- seewave::duration(end_mrkr)

    # reset margins
    par(mar = c(5, 4, 4, 2) + 0.1)

    # add delay at the beggining
    if (delay > 0) {
      strt_mrkr <-
        seewave::addsilw(
          strt_mrkr,
          d = delay,
          output = "Wave",
          at = "start",
          f = sampling_rate
        )
    }

    # add gap to start marker
    strt_mrkr <-
      seewave::addsilw(
        strt_mrkr,
        d = gap.duration,
        output = "Wave",
        at = "end",
        f = sampling_rate
      )

    # add columns to attach durations
    X$pb.start <- NA

    # add start marker duration
    X$pb.start[1] <- seewave::duration(strt_mrkr)

    # read first selection
    plbck <- warbleR::read_sound_file(X, index = 1, path = path)

    # duration first selection
    dr1 <- seewave::duration(plbck)

    # add gap
    plbck <-
      seewave::addsilw(plbck,
        d = gap.duration,
        output = "Wave",
        at = "end"
      )

    # normalize
    plbck <- tuneR::normalize(plbck)

    # add start marker
    plbck <- seewave::pastew(plbck, strt_mrkr, output = "Wave")

    # add end column
    X$pb.end <- NA

    # add duration of first selection
    X$pb.end[1] <- X$pb.start[1] + dr1

    # concatenate all selection with a loop
    for (i in 2:nrow(X))
    {
      # read waves
      wv <- warbleR::read_sound_file(X, index = i, path = path)

      # save duration in sel tab
      X$pb.start[i] <- seewave::duration(plbck)
      X$pb.end[i] <-
        seewave::duration(plbck) + seewave::duration(wv)

      # add gaps
      wv <-
        seewave::addsilw(
          wv,
          d = gap.duration,
          output = "Wave",
          at = "end",
          f = sampling_rate
        )

      # normalize
      wv <- tuneR::normalize(wv)

      # add to master playback
      plbck <- seewave::pastew(wv, plbck, output = "Wave")
    }

    # add end marker
    plbck <- seewave::pastew(end_mrkr, plbck, output = "Wave")

    # margin range for selections on markers
    mar.f <- (flim[2] - flim[1]) / 3

    # add .wav at the end of file.name if not included
    if (!grepl("\\.wav$", file.name, ignore.case = TRUE)) {
      file.name <- paste0(file.name, ".wav")
    }

    # create selection table
    sel.tab <- data.frame(
      sound.files = file.name,
      selec = 1:(nrow(X) + 2),
      start = c(delay, X$pb.start, X$pb.end[nrow(X)] + gap.duration),
      end = c(
        delay + dur_strt_mrkr,
        X$pb.end,
        length(plbck@left) / sampling_rate
      )
    )

    # add bottom freq info
    if (!is.null(X$bottom.freq)) {
      sel.tab$bottom.freq <-
        c(
          strt_mrkr_freq$bottom.freq,
          X$bottom.freq,
          end_mrkr_freq$bottom.freq
        )
    } else {
      sel.tab$bottom.freq <-
        c(
          strt_mrkr_freq$bottom.freq,
          rep(NA, nrow(X)),
          end_mrkr_freq$bottom.freq
        )
    }


    # add top freq info
    if (!is.null(X$top.freq)) {
      sel.tab$top.freq <-
        c(
          strt_mrkr_freq$top.freq,
          X$top.freq,
          end_mrkr_freq$top.freq
        )
    } else {
      sel.tab$top.freq <-
        c(
          strt_mrkr_freq$top.freq,
          rep(NA, nrow(X)),
          end_mrkr_freq$top.freq
        )
    }

    # add start & end markers
    sel.tab$sound.id <- if (!is.null(X$sound.id)) {
      c("start_marker", X$sound.id, "end_marker")
    } else if (is_extended_selection_table(X)) {
      c("start_marker", X$sound.files, "end_marker")
    } else {
      c(
        "start_marker",
        paste(X$sound.files, X$selec, sep = "-"),
        "end_marker"
      )
    }

    # add delay at the end
    if (delay > 0) {
      plbck <-
        seewave::addsilw(plbck,
          d = delay,
          output = "Wave",
          at = "end"
        )
    }

    # normalize whole master playback
    plbck <- tuneR::normalize(plbck, unit = "16")

    # save master playback
    tuneR::writeWave(plbck, file.path(dest.path, file.name), extensible = FALSE)

    # message to let know users the file has been saved
    ohun:::message2(
      paste0(
        "The file ",
        file.name,
        " has been saved in the directory path '",
        normalizePath(dest.path),
        "'"
      )
    )

    return(sel.tab)
  }
