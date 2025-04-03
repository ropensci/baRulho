#' Find a segment of ambient noise to be used as reference
#'
#' \code{spot_ambient_noise} finds a segment of ambient noise to be used as reference by other functions.
#' @inheritParams template_params
#' @param X Object of class 'data.frame', or 'selection_table' (a class are created by the function \code{\link[warbleR]{selection_table}} from the warbleR package) with the test sound files' annotations ('extended_selection_table' are not supported). Must contain the following columns: 1) "sound.files": name of the .wav files, 2) "selec": unique selection identifier (within a sound file), 3) "start": start time and 4) "end": end time of selections, 5)  "bottom.freq": low frequency for bandpass, 6) "top.freq": high frequency for bandpass, 7) "sound.id": ID of sounds used to identify counterparts across distances/transects. 'selec' column values in 'X' cannot be duplicated within a sound file ('sound.files' column) as this combination is used to refer to specific rows.
#' @param length Numeric. Length (in s) of the segments to be used as ambient noise. Must be supplied. Default is \code{NULL}.
#' @param ovlp Numeric vector of length 1 specifying the percentage of overlap between two
#'   consecutive segments. Default is 0. Can be set globally for the current R session via the "ovlp" option (see \code{\link[base]{options}}).
#' @param fun Function to be applied to select the segment to be used as ambient noise. It must be a function that takes a numeric vector (peak sound pressure level values for each candidate segment) and a single value with the index of the value to keep. Default is \code{function(x) which.min(abs(x - mean(x)))}.
#' @details This function finds a segment of ambient noise to be used as reference by other functions. The function first finds candidate segments that do not overlap with annotated sounds in 'X'. Then, it calculates the peak sound pressure level (SPL) of each candidate segment and applies the function supplied by the argument 'fun' to select a single segment. By default 'fun' searches for the segment with the closest value to the mean peak SPL across all candidate segments. Ambient noise annotations are added as a new row in 'X'. Ambient noise annotations are used by the functions \code{\link{signal_to_noise_ratio}} and \code{\link{noise_profile}} to determine background noise levels. Note that this function does not work with annotations in 'extended_selection_table' format.
#' @return An object similar to 'X' with one additional row for each sound file, containing the selected 'ambient' reference. 
#' @family prepare acoustic data
#' @seealso \code{\link{signal_to_noise_ratio}}, \code{\link{noise_profile}} 
#' @export
#' @name spot_ambient_noise
#' @examples {
#' # set temporary directory
#' td <- tempdir()  
#' # load example data
#' data("test_sounds_est")

#' ########## save acoustic data (This doesn't have to be done 
#' # with your own data as you will have them as sound files already.)
#'  # save example files in working director 
#' for (i in unique(test_sounds_est$sound.files)[1:2]) {
#'  writeWave(object = attr(test_sounds_est, "wave.objects")[[i]],
#'            file.path(tempdir(), i))
#' }
#' test_sounds_df <- as.data.frame(test_sounds_est)
#' test_sounds_df <- test_sounds_df[test_sounds_df$sound.id != "ambient", ]
#' test_sounds_df <- 
#'  test_sounds_df[test_sounds_df$sound.files %in% 
#'   unique(test_sounds_est$sound.files)[1:2], ]
#' ####
#' # closest to mean (default)
#' spot_ambient_noise(X = test_sounds_df, path = td, length = 0.12, ovlp = 20)
#' 
#' # min peak
#' spot_ambient_noise(X = test_sounds_df, path = td, length = 0.12, ovlp = 20, fun = which.min)
#' }
#' @references {
#' #' Araya-Salas, M., Grabarczyk, E. E., Quiroz-Oliva, M., Garcia-Rodriguez, A., & Rico-Guevara, A. (2025). Quantifying degradation in animal acoustic signals with the R package baRulho. Methods in Ecology and Evolution, 00, 1-12. https://doi.org/10.1111/2041-210X.14481
#' Araya-Salas, M., & Smith-Vidaurre, G. (2017). warbleR: An R package to streamline analysis of animal acoustic signals. Methods in Ecology and Evolution, 8(2), 184-191.
#' }
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})

spot_ambient_noise <-
  function(X,
           cores = getOption("mc.cores", 1),
           pb = getOption("pb", TRUE),
           path = getOption("sound.files.path", "."),
           length = NULL,
           ovlp = 0, 
           fun = function(x) which.min(abs(x - mean(x)))) {
    
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
    
      
      # set clusters for windows OS
      if (Sys.info()[1] == "Windows" & cores > 1) {
        cl <- parallel::makePSOCKcluster(cores)
      } else {
        cl <- cores
      }
      
    # split by sound file
    X_split <- split(X, X$sound.files)
    
    output_anns_list <- warbleR:::.pblapply(
              X = X_split,
              pbar = pb,
              cl = cl,
              message = "Spotting candidate segments",
              current = 1, 
              total = 1,
              .spot_ambient_noise,
              length = length,
              ovlp = ovlp,
              path = path, 
              fun = fun,
              cores = cores
              )
    
    # combined results
    output_anns <- do.call(rbind, output_anns_list)
    
    # fix row names
    rownames(output_anns) <- seq_len(nrow(output_anns))
    
    # inform that some didn't find a non-overlaping segment
    if (sum(output_anns$sound.id == "ambient") < length(unique(output_anns$sound.files))) {
      .message("At least 1 sound file did not have a segment not overlapping with annotated sounds (try lower 'length' values)")
    }
    return(output_anns)
  }
