#' Align test sound files
#'
#' @description
#' `align_test_files()` aligns test (re-recorded) sound files. It uses the
#' position of acoustic markers found by [find_markers()] to infer the
#' position of all other sounds referenced in a master sound file,
#' producing an aligned selection table for the re-recorded files.
#'
#' @inheritParams template_params
#' @param X Object of class `data.frame`, `selection_table`, or
#'   `extended_selection_table` (the last 2 classes are created by
#'   [warbleR::selection_table()] from the **warbleR** package) with the
#'   master sound file annotations. This should be the same data used
#'   for finding the position of markers in [find_markers()]. It should
#'   also contain a `sound.id` column that will be used to label
#'   re-recorded sounds according to their counterpart in the master
#'   sound file.
#' @param Y Object of class `data.frame` with the output of
#'   [find_markers()]. This object contains the position of markers in
#'   the re-recorded sound files. If more than one marker is supplied
#'   for a sound file, only the one with the highest correlation score
#'   (`scores` column in `Y`) is used.
#' @param path Character string containing the directory path where
#'   test (re-recorded) sound files are found.
#' @param by.song Logical argument to indicate if the extended selection
#'   table should be created by song (see the `by.song` argument of
#'   [warbleR::selection_table()]). Default `TRUE`.
#' @param marker Character string to define whether a `"start"` or
#'   `"end"` marker would be used for aligning re-recorded sound files.
#'   Default `NULL`. Deprecated.
#' @param ... Additional arguments to be passed to
#'   [warbleR::selection_table()] for customizing the extended selection
#'   table.
#'
#' @details
#' The function aligns sounds found in re-recorded sound files
#' (referenced in `Y`) according to a master sound file (referenced in
#' `X`). If more than one marker is supplied for a sound file only the
#' one with the highest correlation score (`scores` column in `Y`) is
#' used. The function outputs an `extended selection table` by default.
#'
#' @return
#' An object of the same class as `X` with the aligned sounds from the
#' test (re-recorded) sound files.
#'
#' @seealso [manual_realign()] and [auto_realign()], for fixing small
#'   remaining misalignments; [find_markers()], for locating markers in
#'   the first place; and [plot_aligned_sounds()], for visually checking
#'   the result.
#'
#' @export
#' @name align_test_files
#' @family test sound alignment
#' @examples
#' {
#'   # load example data
#'   data(list = c("master_est", "test_sounds_est"))
#'
#'   # save example files in working director to recreate a case in which working
#'   # with sound files instead of extended selection tables.
#'   # This doesn't have to be done with your own data as you will
#'   # have them as sound files already.
#'   for (i in unique(test_sounds_est$sound.files)[1:2]) {
#'     writeWave(object = attr(test_sounds_est, "wave.objects")[[i]], 
#'               file.path(tempdir(), i))
#'   }
#'
#'   # save master file
#'   writeWave(object = attr(master_est, "wave.objects")[[1]], 
#'         file.path(tempdir(), "master.wav"))
#'
#'   # get marker position for the first test file
#'     markers <- find_markers(X = master_est,
#'     test.files = unique(test_sounds_est$sound.files)[1],
#'     path = tempdir())
#'
#'   # align all test sounds
#'   alg.tests <- align_test_files(X = master_est, Y = markers, 
#'   path = tempdir())
#' }
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})
#' @references 
#' Araya-Salas, M., Grabarczyk, E. E., Quiroz-Oliva, M., Garcia-Rodriguez, A., & Rico-Guevara, A. (2025). Quantifying degradation in animal acoustic signals with the R package baRulho. Methods in Ecology and Evolution, 00, 1-12. https://doi.org/10.1111/2041-210X.14481

align_test_files <-
  function(X,
           Y,
           path = getOption("sound.files.path", "."),
           by.song = TRUE,
           marker = NULL,
           cores = getOption("mc.cores", 1),
           pb = getOption("pb", TRUE),
           ...) {
    # check arguments
    check_results <- .check_arguments(
      fun = "align_test_files",
      args = list(
        X = X,
        Y = Y,
        path = path,
        by.song = by.song,
        marker = marker,
        cores = cores,
        pb = pb
      )
    )
    
    # report errors
    .report_assertions(check_results)
    
    # if more than one marker per test files then keep only the marker with the highest score
    if (any(table(Y$sound.files) > 1)) {
      Y <-
        Y[stats::ave(x = -Y$scores,
                     as.factor(Y$sound.files),
                     FUN = rank) <= 1,]
    }
    
      
    # set clusters for windows OS
    if (Sys.info()[1] == "Windows" & cores > 1) {
      cl <- parallel::makePSOCKcluster(cores)
    } else {
      cl <- cores
    }
    
    # align each file
    out <-
      warbleR:::.pblapply(
        pbar = pb,
        X = seq_len(nrow(Y)),
        cl = cl,
        message = "aligning test sounds", 
        current = 1,
        total = if (warbleR::is_extended_selection_table(X)) 3 else 1,
        FUN = function(y) {
          # compute start and end as the difference in relation to the marker position in the master sound file
          start <-
            X$start + (Y$start[y] - X$start[X$sound.id == Y$marker[y]])
          end <-
            X$end + (Y$start[y] - X$start[X$sound.id == Y$marker[y]])
          
          # make data frame
          W <-
            data.frame(
              sound.files = Y$sound.files[y],
              selec = seq_along(start),
              start,
              end,
              bottom.freq = X$bottom.freq,
              top.freq = X$top.freq,
              sound.id = X$sound.id,
              marker = Y$marker[y]
            )
          
          
          # re set rownames
          rownames(W) <- seq_len(nrow(W))
          
          return(W)
        }
      )
    
    # put data frames togheter
    sync.sls <- do.call(rbind, out)
    
    
    # check if any selection exceeds length of recordings
    wvdr <-
      duration_sound_files(path = path, files = unique(sync.sls$sound.files))
    
    # add duration to data frame
    sync.sls <- merge(sync.sls, wvdr)
    
    # start empty vector to add name of problematic files
    problematic_files <- character()
    
    if (any(sync.sls$end > sync.sls$duration)) {
      .warning(
        x = paste(
          sum(sync.sls$end > sync.sls$duration),
          "selection(s) exceeded sound file length and were removed (run getOption('baRulho')$files_to_check_align_test_files to see which test files were involved)"
        )
      )
      
      problematic_files <-
        append(problematic_files, unique(sync.sls$sound.files[sync.sls$end > sync.sls$duration]))
      
      # remove exceeding selections
      sync.sls <- sync.sls[!sync.sls$end > sync.sls$duration,]
    }
    
    if (any(sync.sls$start < 0)) {
      .warning(
        x = paste(
          sum(sync.sls$start < 0),
          "selection(s) were absent at the start of the files (negative start values) and were removed (run getOption('baRulho')$files_to_check_align_test_files to see which test files were involved)"
        )
      )
      
      problematic_files <-
        append(problematic_files, unique(sync.sls$sound.files[sync.sls$start < 0]))
      
      # remove exceeding selections
      sync.sls <- sync.sls[sync.sls$start >= 0,]
    }
    
    on.exit(options(baRulho = list(
      files_to_check_align_test_files = unique(problematic_files)
    )))
    
    # remove duration column and marker
    sync.sls$duration <- NULL
    
    if (warbleR::is_extended_selection_table(X)) {

      if (by.song)
        # if by song add a numeric column to represent sound files
      {
        sync.sls$song <- as.numeric(as.factor(sync.sls$sound.files))
        by.song <- "song"
      } else {
        by.song <- NULL
      } # rewrite by song as null
      

      warbleR:::.update_progress(current = 1, total = 3)
      
      sync.sls <-
        selection_table(
          sync.sls,
          extended = TRUE,
          path = path,
          by.song = by.song,
          pb = pb,
          verbose = pb,
          ...
        )
      
      # remove song column
      sync.sls$song <- NULL
      
      # fix call attribute
      attributes(sync.sls)$call <- base::match.call()
    }
    
    return(sync.sls)
  }
