#' Measure full spectrum sound noise profiles
#'
#' \code{noise_profile} Measure full spectrum sound pressure levels (i.e. noise profiles) in sound files or extended selection tables.
#' @inheritParams template_params
#' @param X Object of class 'data.frame', 'selection_table' or 'extended_selection_table' (the last 2 classes are created by the function \code{\link[warbleR]{selection_table}} from the warbleR package) with the reference to the test sounds . Must contain the following columns: 1) "sound.files": name of the .wav files, 2) "selec": unique selection identifier (within a sound file), 3) "start": start time and 4) "end": end time of selections, 5)  "bottom.freq": low frequency for bandpass, 6) "top.freq": high frequency for bandpass and 7) "sound.id": ID of sounds used to identify counterparts across distances (needed for "custom" noise reference, see "noise.ref" argument). Default is \code{NULL}.
#' @param files Character vector with names of wave files to be analyzed. Files must be found in 'path' supplied (or in the working directory if 'path' is not supplied). Default is \code{NULL}.
#' @param mar numeric vector of length 1. Specifies the margins adjacent to
#'   the start and end points of selection over which to measure ambient noise. Required if 'X' is supplied and ignored if not supplied. Default is \code{NULL}.
#' @param noise.ref Character vector of length 1 to determined which noise segment must be used for measuring ambient noise. Ignored if 'X' is not supplied. Two options are available:
#' \itemize{
#' \item \code{adjacent}: measure ambient noise right before the sound (using argument 'mar' to define duration of ambient noise segments).
#' \item \code{custom}: measure ambient noise segments referenced in the selection table (labeled as 'ambient' in the 'sound.id' column).
#' }
#' @param parallel DEPRECATED. Use 'cores' instead.
#' @param bp Numeric vector of length 2 giving the lower and upper limits of a frequency bandpass filter (in kHz). Default is \code{NULL}.
#' @param hop.size A numeric vector of length 1 specifying the time window duration (in ms). Default is 1 ms, which is equivalent to ~45 wl for a 44.1 kHz sampling rate. Ignored if 'wl' is supplied.
#' @param wl A numeric vector of length 1 specifying the window length of the spectrogram, default
#' is NULL. Ignored if \code{bp = NULL}. If supplied, 'hop.size' is ignored.
#' Note that lower values will increase time resolution, which is more important for amplitude ratio calculations.
#' @param PSD Logical to control whether the Probability Mass Function (the probability distribution of frequencies). See \code{\link[seewave]{meanspec}}. Default is \code{FALSE}.
#' @param norm Logical to control whether amplitude values are normalized (divided by the maximum) so the highest value is 1. See \code{\link[seewave]{meanspec}}. Default is \code{TRUE}.
#' @param dB A character string of length 1 specifying the type dB to return: "max0" for a maximum dB value at 0, "A", "B", "C", "D", and "ITU" for common dB weights. See \code{\link[seewave]{meanspec}}. Default is \code{"A"}.
#' @param averaged Logical to control if frequency spectra are averaged within a sound file. Default is \code{TRUE}.
#' @return A data frame containing the frequency spectra for each sound file or wave object (if 'X' is supplied and is of class 'extended.selection.table').
#' @export
#' @name noise_profile
#' @details The function estimates full spectrum sound pressure levels (i.e. noise profiles) of ambient noise. This can be done on data frames/(extended) selection tables (using the segments containing no target sound or the 'ambient' sound id) or over complete sound files in the working directory (or path supplied). The function uses \code{\link[seewave]{meanspec}} internally to calculate frequency spectra.
#' @examples {
#'   # load example data
#'   data("test_sounds_est")
#'
#'   # measure on custom noise reference
#'   noise_profile(X = test_sounds_est, mar = 0.01, pb = FALSE, noise.ref = "custom")
#'
#'   # remove noise selections so noise is measured right before the signals
#'   pe <- test_sounds_est[test_sounds_est$sound.id != "ambient", ]
#'
#'   noise_profile(X = pe, mar = 0.01, pb = FALSE, noise.ref = "adjacent")
#' }
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})
#' @family miscellaneous
#' @seealso \code{\link{excess_attenuation}}
#' @references {
#' Araya-Salas, M. (2020). baRulho: baRulho: quantifying degradation of (animal) acoustic signals in R. R package version 1.0.2.
#'
#' }

noise_profile <-
  function(X = NULL,
           files = NULL,
           mar = NULL,
           noise.ref = "adjacent",
           parallel = NULL,
           cores = getOption("mc.cores", 1),
           pb = getOption("pb", TRUE),
           path = getOption("sound.files.path", "."),
           bp = NULL,
           hop.size = getOption("hop.size", 1),
           wl = getOption("wl", NULL),
           PSD = FALSE,
           norm = TRUE,
           dB = "A",
           averaged = TRUE) {
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
    
    # if X is not supplied modify noise reference
    if (!is.null(X)) {
      
      # keep only 'ambient' selections
      if (noise.ref == "custom") {
        X <- X[X$sound.id == "ambient", ]
      }
      if (noise.ref == "adjacent" &
          is.null(mar)) {
        .stop("'mar' must be supplied when 'noise.ref == 'adjacent''")
      }
    } else
      # if  no  files and no X get files in path
      if (is.null(files)) {
        # check path to working directory
        if (is.null(path)) {
          path <- getwd()
        } else if (!dir.exists(path)) {
          .stop("'path' provided does not exist")
        } else {
          path <- normalizePath(path)
        }
        
        files <-
          list.files(path = path,
                     pattern = "\\.wav$",
                     ignore.case = TRUE)
        
        if (length(files) == 0) {
          .stop("No files found in working directory")
        }
      }
    
    # check files
    if (!is.null(files)) {
      
      # created selection table from sound files
      X <-
        warbleR::selection_table(whole.recs = TRUE,
                                 pb = FALSE,
                                 path = path)
      
      # filter sound files in files
      X <- X[X$sound.files %in% files,]
      
      # add sound column
      X$sound.id <- "ambient"
      
      # set noise.ref to ambient so the whole sound file is measured
      noise.ref <- "custom"
    }
    
    # adjust wl based on hop.size
    wl <- .adjust_wl(wl, X, hop.size, path)
      
    # set clusters for windows OS
    if (Sys.info()[1] == "Windows" & cores > 1) {
      cl <-
        parallel::makePSOCKcluster(getOption("cl.cores", cores))
    } else {
      cl <- cores
    }
    
    # calculate STR
    noise.profiles <-
      warbleR:::pblapply_wrblr_int(
        pbar = pb,
        seq_len(nrow(X)),
        cl = cl,
        .noise_profile,
        Y = X,
        noise.ref = noise.ref,
        mar = mar,
        path = path,
        wl = wl,
        PSD = PSD,
        bp = bp,
        norm = norm,
        dB = dB
      )
    
    # get numbers of rows
    rws <- vapply(noise.profiles, nrow, FUN.VALUE = numeric(1))
    
    # make all the same length if noise.ref is adjacent
    if (length(unique(rws)) > 1 & noise.ref == "adjacent") {
      noise.profiles <- .same_length_noise(noise.profiles, rws)  
    }
    
    # put together in 1 data frame
    noise.profile <- do.call(rbind, noise.profiles)
    
    # get mean by sound file
    if (averaged) {
      noise.profile <-
        aggregate(amp ~ sound.files + freq, data = noise.profile, FUN = mean)
    }
    
    # sort by sound file and freq
    noise.profile <-
      noise.profile[order(noise.profile$sound.files, noise.profile$freq),]
    
    rownames(noise.profile) <- seq_len(nrow(noise.profile))
    
    return(noise.profile)
  }
