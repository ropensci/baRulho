#' Align test sound files
#'
#' \code{align_test_files} aligns test (re-recorded) sound files.
#' @inheritParams template_params
#' @param X object of class 'data.frame', 'selection_table' or 'extended_selection_table' (the last 2 classes are created by the function \code{\link[warbleR]{selection_table}} from the warbleR package). This should be the same data than that was used for finding the position of markers in \code{\link{find_markers}}. It should also contain a 'sound.id' column that will be used to label re-recorded sounds according to their counterpart in the master sound file.
#' @param Y object of class 'data.frame' with the output of \code{\link{find_markers}}. This object contains the position of markers in the re-recorded sound files. If more than one marker is supplied for a sound file only the one with the highest correlation score ('scores' column in 'X') is used.
#' @param path Character string containing the directory path where test (re-recorded) sound files are found.
#' @param by.song Logical argument to indicate if the extended selection table should be created by song (see 'by.song' \code{\link[warbleR]{selection_table}} argument). Default is \code{TRUE}.
#' @param marker Character string to define whether a "start" or "end" marker would be used for aligning re-recorded sound files. Default is \code{NULL}. DEPRECATED.
#' @param ...	Additional arguments to be passed to \code{\link[warbleR]{selection_table}} for customizing extended selection table.
#' @return An object of the same class than 'X' with the aligned sounds from test (re-recorded) sound files.
#' @export
#' @name align_test_files
#' @details The function aligns sounds found in re-recorded sound files (referenced in 'Y') according to a master sound file (referenced in 'X'). If more than one marker is supplied for a sound file only the one with the highest correlation score ('scores' column in 'X') is used. The function outputs an 'extended selection table' by default. 
#' @family test sound alignment
#' @seealso \code{\link{manual_realign}}; \code{\link{find_markers}}; \code{\link{plot_aligned_sounds}}
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
#' @references {
#' Araya-Salas M., E. Grabarczyk, M. Quiroz-Oliva, A. Garcia-Rodriguez, A. Rico-Guevara. (2023), baRulho: an R package to quantify degradation in animal acoustic signals .bioRxiv 2023.11.22.568305.
#' }

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
    
    # if more than one marker per test files then keep only the marker with the highest score
    if (any(table(Y$sound.files) > 1)) {
      Y <-
        Y[stats::ave(x = -Y$scores,
                     as.factor(Y$sound.files),
                     FUN = rank) <= 1,]
    }
    
    if (pb & warbleR::is_extended_selection_table(X)) {
      write(file = "",
            x = paste0("Aligning test sound files (step 1 out of 2):"))
    }
    
    # set clusters for windows OS
    if (Sys.info()[1] == "Windows" & cores > 1) {
      cl <- parallel::makePSOCKcluster(getOption("cl.cores", cores))
    } else {
      cl <- cores
    }
    
    # align each file
    out <-
      warbleR:::pblapply_wrblr_int(
        pbar = pb,
        X = seq_len(nrow(Y)),
        cl = cl,
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
      wavdur(path = path, files = unique(sync.sls$sound.files))
    
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
      if (pb) {
        write(
          file = "",
          x = paste0("Creating extended selection table (step 2 out of 2):")
        )
      }
      if (by.song)
        # if by song add a numeric column to represent sound files
      {
        sync.sls$song <- as.numeric(as.factor(sync.sls$sound.files))
        by.song <- "song"
      } else {
        by.song <- NULL
      } # rewrite by song as null
      
      sync.sls <-
        selection_table(
          sync.sls,
          extended = TRUE,
          confirm.extended = FALSE,
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
