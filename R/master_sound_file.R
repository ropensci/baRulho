#' Create a master sound file
#'
#' \code{master_sound_file} creates a master sound file to be used in playback experiments related to sound degradation.
#' @inheritParams template_params
#' @param X Object of class 'data.frame', 'selection_table' or 'extended_selection_table' (the last 2 classes are created by the function \code{\link[warbleR]{selection_table}} from the warbleR package) with the reference to the model sounds. Must contain the following columns: 1) "sound.files": name of the .wav files, 2) "selec": unique selection identifier (within a sound file), 3) "start": start time and 4) "end": end time of selections, 5)  "bottom.freq": low frequency for bandpass and 6) "top.freq": high frequency for bandpass. An optional 'sound.id' column can be included to use a custom label for each sound in the output. This column must contain a unique id for each sound (labels cannot repeated). If not supplied the function will make it by combining the sound file and selection columns.
#' @param file.name Character string indicating the name of the sound file.
#' @param dest.path Character string containing the directory path where the sound file will be saved.
#' If \code{NULL} (default) then the current working directory will be used instead.
#' @param overwrite Logical argument to determine if the function will overwrite any existing sound file with the same file name. Default is the current working directory.
#' @param delay Numeric vector of length 1 to control the duration (in s) of a silence gap at the beginning (and at the end) of the sound file. This can be useful to allow some time at the start of the playback experiment. Default is 1.
#' @param gap.duration Numeric vector of length 1 to control the duration (in s) of silence gaps to be placed in between sounds. Default is 1 s.
#' @param amp.marker Numeric vector of length 1 to use as a constant to amplify markers amplitude. This is useful to increase the amplitude of markers in relation to those of sounds, so it is picked up at further distances. Default is 2.
#' @param flim Numeric vector of length 2 to control the (approximate) frequency range in which the markers would be found. If \code{NULL} markers would span across the entire frequency range. Default is c(0, 4). 
#' @param cex Numeric vector of length 1 indicating the font size for the start and end markers. Default is 14.
#' @return A .wav file in 'path' as well as a data frame in the R environment with the annotations (i.e. time position) of sounds in the master sound file and an additional column 'sound.id' that provides a unique id for each sound in the sound file. This is useful for identifying/labeling sounds in test (re-recorded) sound files for downstream analyses.
#' @export
#' @name master_sound_file
#' @details The function is intended to simplify the creation of master sound files for playback experiments in sound degradation studies. The function clips sounds from sound files (or wave objects from extended selection tables) and concatenates them in a single sound file. The function also adds acoustic markers at the start and end of the playback that can be used to time-sync test (re-recorded) sounds to facilitate the streamlining of degradation quantification. There is no predefined limit to the duration of the output master sound file, although long this be constrained by computer memory. As a reference, master sound files of up to 10 min have been created in a 16GB RAM laptop computer.
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
#'     X = lbh_selec_table, extended = TRUE,
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
#' @family prepare acoustic data
#' @references {
#' Araya-Salas, M., Grabarczyk, E. E., Quiroz-Oliva, M., Garcia-Rodriguez, A., & Rico-Guevara, A. (2025). Quantifying degradation in animal acoustic signals with the R package baRulho. Methods in Ecology and Evolution, 00, 1-12. https://doi.org/10.1111/2041-210X.14481
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
      .check_arguments(fun = arguments[[1]], args = arguments)
    
    # report errors
    .report_assertions(check_results)
    
    # get sampling rate
    sampling_rate <-
      warbleR::read_sound_file(
        X = X,
        index = 1,
        path = path,
        header = TRUE
      )$sample.rate
    
    mrkrs <- .make_markers(X, flim, sampling_rate, cex)  
    
    # output wave object
    strt_mrkr <- tuneR::normalize(mrkrs$strt_mrkr)
    end_mrkr <- tuneR::normalize(mrkrs$end_mrkr)
    
    # frequency range of markers
    strt_mrkr_freq <-
      warbleR::freq_range_detec(
        strt_mrkr,
        fsmooth = 0.2,
        plot = FALSE,
        dB.threshold = 20
      )
    
    # fix bottom freq if NA
    if (is.na(strt_mrkr_freq$bottom.freq))
      strt_mrkr_freq$bottom.freq <- flim[1]
    
    # fix top freq if NA
    if (is.na(strt_mrkr_freq$top.freq))
      strt_mrkr_freq$top.freq <- flim[2]
    
    
    end_mrkr_freq <-
      warbleR::freq_range_detec(
        end_mrkr,
        fsmooth = 0.2,
        plot = FALSE,
        dB.threshold = 20
      )
    
    # fix bottom freq if NA
    if (is.na(end_mrkr_freq$bottom.freq))
      end_mrkr_freq$bottom.freq <- flim[1]
    
    # fix top freq if NA
    if (is.na(end_mrkr_freq$top.freq))
      end_mrkr_freq$top.freq <- flim[2]
    
    # amplify markers
    strt_mrkr@left <- strt_mrkr@left * amp.marker
    end_mrkr@left <- end_mrkr@left * amp.marker
    
    # save duration of markers for creating selection table
    dur_strt_mrkr <- seewave::duration(strt_mrkr)
    dur_end_mrkr <- seewave::duration(end_mrkr)
    
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
                       at = "end")
    
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
      selec = seq_along(nrow(X) + 2),
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
        c(strt_mrkr_freq$bottom.freq,
          X$bottom.freq,
          end_mrkr_freq$bottom.freq)
    } else {
      sel.tab$bottom.freq <-
        c(strt_mrkr_freq$bottom.freq,
          rep(NA, nrow(X)),
          end_mrkr_freq$bottom.freq)
    }
    
    
    # add top freq info
    if (!is.null(X$top.freq)) {
      sel.tab$top.freq <-
        c(strt_mrkr_freq$top.freq,
          X$top.freq,
          end_mrkr_freq$top.freq)
    } else {
      sel.tab$top.freq <-
        c(strt_mrkr_freq$top.freq,
          rep(NA, nrow(X)),
          end_mrkr_freq$top.freq)
    }
    
    # add start & end markers
    sel.tab$sound.id <- if (!is.null(X$sound.id)) {
      c("start_marker", X$sound.id, "end_marker")
    } else if (warbleR::is_extended_selection_table(X)) {
      c("start_marker", X$sound.files, "end_marker")
    } else {
      c("start_marker",
        paste(X$sound.files, X$selec, sep = "-"),
        "end_marker")
    }
    
    # add delay at the end
    if (delay > 0) {
      plbck <-
        seewave::addsilw(plbck,
                         d = delay,
                         output = "Wave",
                         at = "end")
    }
    
    # normalize whole master playback
    plbck <- tuneR::normalize(plbck, unit = "16")
    
    # save master playback
    tuneR::writeWave(plbck, file.path(dest.path, file.name), extensible = FALSE)
    
    # message to let know users the file has been saved
    .message(
      paste0(
        "The file ",
        file.name,
        " has been saved in the directory path '",
        normalizePath(dest.path),
        "'"
      )
    )
    
    # add a unique selec id to each annotation
    sel.tab$selec <- seq_len(nrow(sel.tab))
    
    return(sel.tab)
  }
