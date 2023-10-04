#' Set reference for test sounds
#'
#' \code{set_reference_sounds} set rows to be used as reference for each test sound.
#' @usage set_reference_sounds(X, method = getOption("method", 1),
#' cores =  getOption("mc.cores", 1), pb =  getOption("pb", TRUE),
#' path = getOption("sound.files.path", "."))
#' @param X Object of class 'data.frame', 'selection_table' or 'extended_selection_table' (the last 2 classes are created by the function \code{\link[warbleR]{selection_table}} from the warbleR package) with the reference to the test sounds . Must contain the following columns: 1) "sound.files": name of the .wav files, 2) "selec": unique selection identifier (within a sound file), 3) "start": start time and 4) "end": end time of selections, 5)  "bottom.freq": low frequency for bandpass, 6) "top.freq": high frequency for bandpass, 7) "sound.id": ID of sounds used to identify counterparts across distances (and transects if more than 1) and 8) "distance": distance (numeric) at which each test sound was re-recorded. A 'transect' column labeling those sounds recorded in the same transect is required if 'method = 2'. 'X' can only have 1 copy for any given sound id in a distance or a transect-distance combination (if column 'transect' is supplied). In addition, 'selec' column values in 'X' cannot be duplicated within a sound file ('sound.files' column) as this combination is used to refer to specific rows in the output 'reference' column.
#' @param cores Numeric vector of length 1. Controls whether parallel computing is applied by specifying the number of cores to be used. Default is 1 (i.e. no parallel computing).
#' @param pb Logical argument to control if progress bar is shown. Default is \code{TRUE}.
#' @param method Numeric vector of length 1 to indicate the 'experimental design' for measuring degradation. Two methods are available:
#' \itemize{
#' \item \code{1}: compare sounds (by 'sound.id') with their counterpart that was recorded at the closest distance to source (e.g. compare a sound recorded at 5m, 10m and 15m with its counterpart recorded at 1m). This is the default method. The function will try to use references from the same transect. However, if there is another test sound from the same 'sound.id' at a shorter distance in other transects, it will be used as reference instead.
#' \item \code{2}: compare all sounds with their counterpart recorded at the distance immediately before within a transect (e.g. a sound recorded at 10m compared with the same sound recorded at 5m, then sound recorded at 15m compared with same sound recorded at 10m and so on). 'transect' column in 'X' is required.
#' }
#' @param path Character string containing the directory path where the sound files are found. Only needed when 'X' is not an extended selection table.
#' @details This function adds a 'reference' column defining which sounds will be used by other functions as reference. Two methods are available (see 'methods' argument description). Note that if users want to define their own reference sound this can be set manually. If so  NAs must be used to indicate rows to be ignored. References must be indicated as a the combination of the 'sound.files' and 'selec' column. For instance, '10m.wav-1' indicates that the row in which the 'selec' column is '1' and the sound file is '10m.wav' should be used as reference. The function also checks that the information 'X' is in the right format so it wont produce errors in downstream analysis (see 'X' argument description for details on format). The function will ignore rows with 'sound.id' 'ambient', 'start_marker' and 'end_marker'
#' @return An object similar to 'X' with one additional column, 'reference', with the ID of the sounds to be used as reference by degradation-quantifying functions in downstream analyses. The ID is created as \code{paste(X$sound.files, X$selec, sep = "-")}.
#' @export
#' @family quantify degradation
#' @seealso \code{\link[warbleR]{check_wavs}}, \code{\link[warbleR]{check_sels}}
#' @export
#' @name set_reference_sounds
#' @export
#' @examples{
#'   # load example data
#'   data("test_sounds_est")
#'
#' # save wav file examples
#' X <- test_sounds_est[test_sounds_est$sound.files != "master.wav", ]
#'
#'
#' # method 1
#' Y <- set_reference_sounds(X = X)
#'
#' # method 2
#' Y <- set_reference_sounds(X = X, method = 2)
#' }
#' @references {Araya-Salas, M., & Smith-Vidaurre, G. (2017). warbleR: An R package to streamline analysis of animal acoustic signals. Methods in Ecology and Evolution, 8(2), 184-191.}
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})
# last modification on jul-5-2016 (MAS)

set_reference_sounds <-
  function(X,
           method = getOption("method", 1),
           cores = getOption("mc.cores", 1),
           pb = getOption("pb", TRUE),
           path = getOption("sound.files.path", ".")) {
    # if (!is.null(X$reference))
    #   stop2("Column 'reference' already found in 'X'. Must be removed first.")
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

    # add sound file selec colums to X (weird column name so it does not overwrite user columns)
    X$.sgnl.temp <- paste(X$sound.files, X$selec, sep = "-")

    ## logical to control whether references are added
    add_refs <- TRUE

    if (!is_extended_selection_table(X) & !is_selection_table(X)) {
      # print message
      if (pb) {
        write(file = "", x = "Checking annotations (step 1 out of 2):")
      }

      X_check <-
        warbleR::check_sels(
          X = X,
          parallel = cores,
          verbose = FALSE,
          pb = pb,
          fix.selec = FALSE,
          path = path
        )

      if (any(X_check$check.res != "OK")) {
        options("warn" = 1)
        warning2(x = "Reference cannot be added as 'X' doesn't have the right format. Check the column 'check.res' in the output data frame to spot problematic annotations and their issues. See warbleR's check_sels() function documentation for checking details.")
        X <- X_check

        add_refs <- FALSE
      }
    }

    if (add_refs) {
      # make it a regular data frame (no est)
      X2 <- as.data.frame(X)

      # set clusters for windows OS
      if (Sys.info()[1] == "Windows" & cores > 1) {
        cl <- parallel::makePSOCKcluster(getOption("cl.cores", cores))
      } else {
        cl <- cores
      }

      # add column with names of the reference sounds to be compared against
      if (is.null(X$reference)) {
        if (pb & !is_extended_selection_table(X) & !is_selection_table(X)) {
          write(file = "", x = "Computing references (step 2 out of 2):")
        }

        X$reference <-
          unlist(warbleR:::pblapply_wrblr_int(X = X2$.sgnl.temp, pbar = pb, cl = cl, function(x, meth = method) {
            # extract for single sound id and order by distance
            Y <-
              X2[X2$sound.id == X$sound.id[X2$.sgnl.temp == x], , drop = FALSE]

            if (!Y$sound.id[1] %in% c("ambient", "start_marker", "end_marker")) {
              # Order by distance
              Y <- Y[order(Y$distance), ]

              # method 1 compare to closest distance to source
              if (meth == 1) {
                # if column transect is found select the lowest distance in that trasnect
                if (!is.null(X2$transect)) {
                  W <-
                    X2[X2$transect == X2$transect[X2$.sgnl.temp == x] &
                      X2$sound.id == X$sound.id[X2$.sgnl.temp == x], , drop = FALSE]

                  # if there is another distance that is shorter in other transects for that signal, use that distance
                  z <- if (min(W$distance) <= min(Y$distance)) {
                    W$.sgnl.temp[which.min(W$distance)]
                  } else {
                    Y$.sgnl.temp[which.min(Y$distance)]
                  }
                } else {
                  z <- Y$.sgnl.temp[which.min(Y$distance)]
                }
              } else {
                # if method 2
                # get those from the same transect and same sound id
                W <-
                  X2[X2$transect == X2$transect[X2$.sgnl.temp == x] &
                    X2$sound.id == X$sound.id[X2$.sgnl.temp == x], , drop = FALSE]
                W <- W[order(W$distance), ]

                # if not the first row then the previous row
                if (W$.sgnl.temp[1] != x) {
                  z <- W$.sgnl.temp[which(W$.sgnl.temp == x) - 1]
                } else {
                  # else the first row
                  z <- x
                }
              }

              # set reference to NA if is the same than the current row
              if (z == x) {
                z <- NA
              }
            } else {
              z <- NA
            }

            return(z)
          }))
      }

      # remove temp column
      X$.sgnl.temp <- NULL

      if (all(is.na(X$reference))) {
        stop2("No single reference could be determined")
      }
    }

    return(X)
  }
