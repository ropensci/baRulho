#' Find acoustic markers on test sound files
#'
#' \code{find_markers} find acoustic markers on test (re-recorded) sound files using spectrographic cross-correlation.
#' @inheritParams template_params
#' @param X Object of class 'data.frame', 'selection_table' or 'extended_selection_table' (the last 2 classes are created by the function \code{\link[warbleR]{selection_table}} from the warbleR package) with the reference to the sounds in the master sound file. Must contain the following columns: 1) "sound.files": name of the .wav files, 2) "selec": unique selection identifier (within a sound file), 3) "start": start time, 4) "end": end time of selections and 5) "sound.id": unique identifier for each of the annotated sounds in 'X'. Columns for 'top.freq', 'bottom.freq' and 'channel' are optional. The acoustic start and end markers (added by \code{\link{master_sound_file}}) should be labeled as "start_marker" and "end_marker" respectively. Required.
#' @param markers Character vector with the name of the annotations (as in the column 'sound.id') to be used as templates for cross-correlation. Default is \code{c("start_marker", "end_marker")}. Using more than one marker is recommended as the time difference between their position can be used to evaluate the precision of the detection (see 'Value' section).
#' @param test.files Character vector of length 1 with the name(s) of the test (re-recorded) file(s) in which to search for the marker(s). If not supplied all sound files in 'path' are used instead.
#' @param path Character string containing the directory path where test (re-recorded) sound files are found.
#' @param ...	Additional arguments to be passed to \code{\link[ohun]{template_correlator}} for setting cross-correlation parameters (e.g. 'wl', 'ovlp', etc).
#' @return A data frame with test file names, marker id, maximum cross-correlation score for each marker and the start and end where it was detected. If two or more markers are used the function computes an additional column, 'time.mismatch', that compares the time difference between the two markers in the test-files against that in the master sound file. In a perfect detection the value must be 0.
#' @export
#' @name find_markers
#' @details The function takes a master sound file's reference data ('X') and finds the position of acoustics markers ('markers' argument, included as selections in 'X') in the re-recorded sound files. This is used to align signals found in re-recorded sound files according to a master sound file referenced in 'X'. The position of the markers is determined as the highest spectrogram cross-correlation value for each marker using the functions \code{\link[ohun]{template_correlator}} and \code{\link[ohun]{template_detector}}. \strong{Make sure the master sound file (that referred to in 'X') is found in the same folder than the re-recorded sound files}. Take a look at the package vignette for information on how to incorporate this function into a sound degradation analysis workflow. In cases in which markers are not correctly detected editing test sound files to remove audio segments with no target sounds (before the start marker and after the end marker) can improve performance. Using a low 'hop.size' or window length 'wl' (used internally by \code{\link[ohun]{template_correlator}}) can help to improve precision. Other spectrogram types (argument 'type' in \code{\link[ohun]{template_correlator}}) can sometimes show better performance when markers are highly degraded. If frequency range columns are included ('bottom.freq' and 'top.freq', in kHz) cross-correlation will be run on those frequency ranges. All templates must have the same sampling rate and both templates and 'files' (in which to find templates) must also have the same sampling rate.
#' @family test sound alignment
#' @seealso \code{\link{manual_realign}}; \code{\link{auto_realign}}; \code{\link{align_test_files}}; \code{\link{master_sound_file}}
#' @examples {
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
#' }
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})
#' @references
#' Araya-Salas, M., Grabarczyk, E. E., Quiroz-Oliva, M., Garcia-Rodriguez, A., & Rico-Guevara, A. (2025). Quantifying degradation in animal acoustic signals with the R package baRulho. Methods in Ecology and Evolution, 00, 1-12. https://doi.org/10.1111/2041-210X.14481

find_markers <-
  function(X,
           markers = c("start_marker", "end_marker"),
           test.files = NULL,
           path = getOption("sound.files.path", "."),
           pb = getOption("pb", TRUE),
           cores = getOption("mc.cores", 1),
           ...) {
    # check arguments
    arguments <- as.list(base::match.call())
    
    # add objects to argument names
    for (i in names(arguments)[-1]) {
      # use try to avoid errors with argumets from dots (...)
      try(arguments[[i]] <- get(i), silent = TRUE)
    }
    
    # check each arguments
    check_results <-
      .check_arguments(fun = arguments[[1]], args = arguments)
    
    # report errors
    .report_assertions(check_results)
    
    if (!all(markers %in% X$sound.id))
      .stop("at least one value in 'markers' not found in 'sound.id' column in 'X'") else
      X <- X[X$sound.id %in% markers,]
    
    # get sound files in path
    files_in_path <-
      list.files(path = path,
                 pattern = "\\.wav$",
                 ignore.case = TRUE)
    
    # check if there are files
    if (length(files_in_path) == 0) {
      .stop("No .wav files found in 'path'")
    }
    
    # check for master sound file
    if (!any(files_in_path %in% unique(X$sound.files))) {
      .stop(
        "sound file referenced in 'X' not found in 'path' (make sure you put the master sound file in the same folder than the re-recorded files)"
      )
    }
    
    # remove sound files
    if (!is.null(test.files)) {
      if (!all(test.files %in% files_in_path)) {
        .stop("Not all 'test.files' were found in 'path'")
      }
    } else {
      # remove master sound file
      test.files <-
        files_in_path[!files_in_path %in% unique(X$sound.files)]
    }
    
    # get metadata of sound files to get sampling rates
    wi <-
      warbleR::info_sound_files(
        path = path,
        parallel = 1,
        pb = FALSE,
        skip.error = TRUE,
        files = c(test.files, as.character(unique(X$sound.files)))
      )
    
    if (length(unique(wi$sample.rate)) > 1) {
      .stop(
        "Not all sound files share the same sampling rate (check wave properties with warbleR::info_sound_files())"
      )
    }
    
    # run cross correlation
    if (pb)  {
   warbleR:::.update_progress(total = 2)
    }
    
    
    templ_corrs <-
      ohun::template_correlator(
        templates = X,
        files = test.files,
        path = path,
        cores = cores,
        pb = pb,
        ...
      )
    
    # find peaks
    pks <-
      as.data.frame(
        ohun::template_detector(
          template.correlations = templ_corrs,
          cores = cores,
          threshold = 0.0001,
          pb = pb
        )
      )
    
    pks <-
      pks[stats::ave(
        x = -pks$scores,
        as.factor(pks$sound.files),
        as.factor(pks$template),
        FUN = rank
      ) <= 1,]
    
    # rename markers
    pks$marker <-
      vapply(pks$template, function(x)
        as.character(X$sound.id[paste(X$sound.files, X$selec, sep = "-") == x]), FUN.VALUE = character(length = 1))
    
    # fix row labels and selec labels
    pks$selec <- rownames(pks) <- seq_len(nrow(pks))
    
    
    # check that the distance between markers is similar in re-recorded files compare to the master
    # if all sound files have two markers and two markers where used
    if (all(tapply(X$sound.id, X$sound.files, function(y)
      sum(grepl("marker", y))) >= 2) & length(markers) > 1) {
      time_dist_markers <-
        vapply(test.files, function(x)
          pks$end[pks$marker == markers[2] &
                    pks$sound.files == x] - pks$start[pks$marker == markers[1] &
                                                        pks$sound.files == x], FUN.VALUE = numeric(length = 1L))
      
      # subtract duration in master sound file
      time_dist_markers <-
        time_dist_markers - (X$end[X$sound.id == markers[2]] - X$start[X$sound.id == markers[1]])
      
      # add results
      pks$time.mismatch <- NA
      pks$time.mismatch[pks$marker == markers[1]] <-
        time_dist_markers
    }
    
    # remove template column
    pks$template <- NULL
    
    return(pks)
  }
