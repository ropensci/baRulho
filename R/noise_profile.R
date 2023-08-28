#' Measure full spectrum sound noise profiles
#'
#' \code{noise_profile} Measure full spectrum sound pressure levels (i.e. noise profiles) in sound files or extended selection tables.
#' @usage noise_profile(X = NULL, files = NULL, mar = NULL,
#' noise.ref = "adjacent", parallel = NULL, cores = getOption("mc.cores", 1),
#' pb = getOption("pb", TRUE), path = getOption("sound.files.path", "."),
#' bp = NULL, hop.size = getOption("hop.size", 1), wl = getOption("wl", NULL),
#' PSD = FALSE, norm = TRUE, dB = "A", averaged = TRUE)
#' @param X Object of class 'data.frame', 'selection_table' or 'extended_selection_table' (the last 2 classes are created by the function \code{\link[warbleR]{selection_table}} from the warbleR package) with the reference to the sounds in the master sound file. Must contain the following columns: 1) "sound.files": name of the .wav files, 2) "selec": unique selection identifier (within a sound file), 3) "start": start time and 4) "end": end time of selections, 5)  "bottom.freq": low frequency for bandpass, 6) "top.freq": high frequency for bandpass and 7) "sound.id": ID of sounds used to identify counterparts across distances (needed for "custom" noise reference, see "noise.ref" argument). Default is \code{NULL}.
#' @param files Character vector with names of wave files to be analyzed. Files must be found in 'path' supplied (or in the working directory if 'path' is not supplied). Default is \code{NULL}.
#' @param mar numeric vector of length 1. Specifies the margins adjacent to
#'   the start and end points of selection over which to measure ambient noise. Required if 'X' is supplied and ignored if not supplied. Default is \code{NULL}.
#' @param noise.ref Character vector of length 1 to determined which noise segment must be used for measuring ambient noise. Ignored if 'X' is not supplied. Two options are available:
#' \itemize{
#' \item \code{adjacent}: measure ambient noise right before the sound (using argument 'mar' to define duration of ambient noise segments).
#' \item \code{custom}: measure ambient noise segments referenced in the selection table (labeled as 'ambient' in the 'sound.id' column).
#' }
#' @param parallel DEPRECATED. Use 'cores' instead.
#' @param cores Numeric vector of length 1. Controls whether parallel computing is applied by specifying the number of cores to be used. Default is 1 (i.e. no parallel computing).
#' @param pb Logical argument to control if progress bar is shown. Default is \code{TRUE}.
#' @param path Character string containing the directory path where the sound files are located.
#' If \code{NULL} (default) then the current working directory is used.
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
#' @examples {{ # load example data
#'   data("test_sounds_est")
#'
#'   # measure on custom noise reference
#'   noise_profile(X = test_sounds_est, mar = 0.01, pb = FALSE, noise.ref = "custom")
#'
#'   # remove noise selections so noise is measured right before the signals
#'   pe <- test_sounds_est[test_sounds_est$sound.id != "ambient", ]
#'
#'   noise_profile(X = pe, mar = 0.01, pb = FALSE, noise.ref = "adjacent") }}
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})
#' @seealso \code{\link{excess_attenuation}}
#' @references {
#' Araya-Salas, M. (2020). baRulho: baRulho: quantifying habitat-induced degradation of (animal) acoustic signals in R. R package version 1.0.2.
#'
#' }

noise_profile <-
  function(X = NULL,
           files = NULL,
           mar = NULL,
           noise.ref = "adjacent",
           parallel = NULL, cores = getOption("mc.cores", 1),
           pb = getOption("pb", TRUE),
           path = getOption("sound.files.path", "."),
           bp = NULL,
           hop.size = getOption("hop.size", 1),
           wl = getOption("wl", NULL),
           PSD = FALSE,
           norm = TRUE,
           dB = "A",
           averaged = TRUE) {
    if (!is.null(X)) {
      # invert selections so gaps become selections instead if noise.ref != ambient
      if (noise.ref == "custom" &
        !any(X$sound.id == "ambient")) {
        stop2("'noise.ref = custom' but no 'ambient' label found in 'sound.id' column ")
      }

      # keep only 'ambient' selections
      if (noise.ref == "custom") {
        X <- X[X$sound.id == "ambient", ]
      }

      if (noise.ref == "adjacent" &
        is.null(mar)) {
        stop2("'mar' must be supplied when 'noise.ref == 'adjacent''")
      }
    } else
    # if  no  files and no X get files in path
    if (is.null(files)) {
      # check path to working directory
      if (is.null(path)) {
        path <- getwd()
      } else if (!dir.exists(path)) {
        stop2("'path' provided does not exist")
      } else {
        path <- normalizePath(path)
      }

      files <-
        list.files(
          path = path,
          pattern = "\\.wav$",
          ignore.case = TRUE
        )

      if (length(files) == 0) {
        stop2("No files found in working directory (alternatively supply 'X')")
      }
    }

    # check files
    if (!is.null(files)) {
      if (any(!file.exists(file.path(path, files)))) {
        stop2(paste(paste(files[!file.exists(files)], collapse = "/"), "was (were) not found"))
      }

      # created selection table from sound files
      X <-
        warbleR::selection_table(
          whole.recs = TRUE,
          pb = FALSE,
          path = path
        )

      # filter sound files in files
      X <- X[X$sound.files %in% files, ]

      # add sound column
      X$sound.id <- "ambient"

      # set noise.ref to ambient so the whole sound file is measured
      noise.ref <- "custom"
    }

    # If cores is not numeric
    if (!is.numeric(cores)) {
      stop2("'cores' must be a numeric vector of length 1")
    }

    if (any(!(cores %% 1 == 0), cores < 1)) {
      stop2("'cores' should be a positive integer")
    }

    # hopsize
    if (!is.numeric(hop.size) |
      hop.size < 0) {
      stop2("'hop.size' must be a positive number")
    }

    # adjust wl based on hope.size
    if (is.null(wl)) {
      wl <-
        round(
          read_sound_file(
            X,
            index = 1,
            header = TRUE,
            path = path
          )$sample.rate * hop.size / 1000,
          0
        )
    }

    # make wl even if odd
    if (!(wl %% 2) == 0) {
      wl <- wl + 1
    }

    # set clusters for windows OS
    if (Sys.info()[1] == "Windows" & cores > 1) {
      cl <-
        parallel::makePSOCKcluster(getOption("cl.cores", cores))
    } else {
      cl <- cores
    }

    # calculate STR
    noise.profiles <-
      warbleR:::pblapply_wrblr_int(pbar = pb, seq_len(nrow(X)), cl = cl, function(y) {
        # extract  complete sound file for custom or files in folder
        if (noise.ref == "custom") {
          noise.wv <-
            warbleR::read_sound_file(
              X = X,
              index = y,
              from = 0,
              to = Inf,
              path = path
            )
        }

        if (noise.ref == "adjacent") {
          # reset time coordinates of sounds if lower than 0 o higher than duration
          stn <- X$start[y] - mar

          if (stn < 0) {
            stn <- 0
          }

          # read ambient noise
          noise.wv <-
            warbleR::read_sound_file(
              X = X,
              index = y,
              from = stn,
              to = X$start[y],
              path = path
            )
        }

        # mean spec
        mspc <-
          meanspec(
            wave = noise.wv,
            f = noise.wv@samp.rate,
            plot = FALSE,
            wl = wl,
            ovlp = 0,
            PSD = PSD,
            PMF = FALSE,
            norm = norm,
            dB = dB
          )

        # name columns
        colnames(mspc) <- c("freq", "amp")

        # add sound file name
        mspc <-
          data.frame(
            sound.files = X$sound.files[y],
            selec = X$selec[y],
            mspc
          )

        # add band-pass frequency filter
        if (!is.null(bp)) {
          mspc <- mspc[mspc$freq >= bp[1] & mspc$freq <= bp[2], ]
        }

        return(mspc)
      })

    # get numbers of rows
    rws <- sapply(noise.profiles, nrow)

    # make all the same length if noise.ref is adjacent
    if (length(unique(rws)) > 1 & noise.ref == "adjacent") {
      # gt freq range of minimum
      fr.range <- range(noise.profiles[[which.min(rws)]]$frequency)

      # interpolate so all have the same number of frequency bins
      noise.profiles <- lapply(noise.profiles, function(Y) {
        # interpolate
        Yappr <- approx(
          x = Y$freq,
          y = Y$amp,
          xout = seq(
            from = fr.range[1],
            to = fr.range[2],
            length.out = min(rws)
          ),
          method = "linear"
        )

        Ydf <-
          data.frame(
            sound.files = Y$sound.files[1],
            selec = Y$selec[1],
            freq = Yappr$x,
            amp = Yappr$y
          )

        return(Ydf)
      })
    }

    # put together in 1 data frame
    noise.profile <- do.call(rbind, noise.profiles)

    # get mean by sound file
    if (averaged) {
      noise.profile <-
        aggregate(amp ~ sound.files + freq, data = noise.profile, FUN = mean)
    }

    # sort by sound file and freq
    noise.profile <- noise.profile[order(noise.profile$sound.files, noise.profile$freq), ]

    rownames(noise.profile) <- seq_len(nrow(noise.profile))

    return(noise.profile)
  }
