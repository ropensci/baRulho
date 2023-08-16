#' Measure attenuation as signal-to-noise ratio
#'
#' \code{signal_to_noise_ratio} measures attenuation as signal-to-noise ratio of sounds referenced in an extended selection table.
#' @usage signal_to_noise_ratio(X, mar, parallel = NULL, cores = getOption("mc.cores", 1),
#' pb = getOption("pb", TRUE), eq.dur = FALSE, noise.ref = "adjacent", type = 1,
#' bp = 'freq.range', output = NULL, hop.size = getOption("hop.size", 1),
#' wl = getOption("wl", NULL), ovlp = getOption("ovlp", 0),
#' path = getOption("sound.files.path", "."))
#' @param X Object of class 'data.frame', 'selection_table' or 'extended_selection_table' (the last 2 classes are created by the function \code{\link[warbleR]{selection_table}} from the warbleR package) with the reference to the sounds in the master sound file. Must contain the following columns: 1) "sound.files": name of the .wav files, 2) "selec": unique selection identifier (within a sound file), 3) "start": start time and 4) "end": end time of selections, 5)  "bottom.freq": low frequency for bandpass, 6) "top.freq": high frequency for bandpass and 7) "sound.id": ID of sounds used to identify counterparts across distances (only needed for "custom" noise reference, see "noise.ref" argument).
#' @param mar numeric vector of length 1. Specifies the margins adjacent to
#'   the start and end points of selection over which to measure ambient noise.
#' @param parallel DEPRECATED. Use 'cores' instead.
#' @param cores Numeric vector of length 1. Controls whether parallel computing is applied by specifying the number of cores to be used. Default is 1 (i.e. no parallel computing).
#' @param pb Logical argument to control if progress bar is shown. Default is \code{TRUE}.
#' @param eq.dur Logical. Controls whether the ambient noise segment that is measured has the same duration
#' to that of the sound (if \code{TRUE}. Default is \code{FALSE}). If \code{TRUE} then 'mar' and 'noise.ref' arguments are ignored.
#' @param noise.ref Character vector of length 1 to determined which noise segment must be used for measuring ambient noise. Two options are available:
#' \itemize{
#' \item \code{adjacent}: measure ambient noise right before the sound (using argument 'mar' to define duration of ambient noise segments).
#' \item \code{custom}: measure ambient noise segments referenced in the selection table (labeled as 'ambient' in the 'sound.id' column). Those segments will be used to apply the same ambient noise reference to all sounds in a sound file. Therefore, at least one 'ambient' selection for each sound file must be provided. If several 'ambient' selections by sound file are supplied, then the root mean square of the amplitude envelope will be averaged across those selections.
#' }
#' @param type Numeric vector of length 1. Selects the formula to be used to calculate the signal-to-noise ratio (S = signal
#' , N = background noise):
#' \itemize{
#' \item \code{1}: ratio of S amplitude envelope root mean square to N amplitude envelope root mean square
#'  (\code{20 * log10(rms(env(S))/rms(env(N)))}) as described by Darden (2008).
#' \item \code{2}: ratio of the difference between S amplitude envelope root mean square and N amplitude envelope root mean square to N amplitude envelope root mean square (\code{20 * log10((rms(env(S)) - rms(env(N)))/rms(env(N)))}, as described by Dabelsteen et al (1993).
#' }
#' @param bp Numeric vector of length 2 giving the lower and upper limits of a frequency bandpass filter (in kHz). Alternatively, when set to 'freq.range' (default) the function uses the 'bottom.freq' and 'top.freq' as the bandpass range.
#' @param output DEPRECATED. Now the output format mirrors the class of the input 'X'.
#' @param hop.size A numeric vector of length 1 specifying the time window duration (in ms). Default is 1 ms, which is equivalent to ~45 wl for a 44.1 kHz sampling rate. Ignored if 'wl' is supplied.
#' @param wl A numeric vector of length 1 specifying the window length of the spectrogram, default
#' is NULL. Ignored if \code{bp = NULL}. If supplied, 'hop.size' is ignored. Note that lower values will increase time resolution, which is more important for amplitude ratios calculations.
#' @param ovlp Numeric vector of length 1 specifying the percent overlap between two
#'   consecutive windows, as in \code{\link[seewave]{spectro}}. Default is 0. Only used for bandpass filtering.
#' @param path Character string containing the directory path where the sound files are found. Only needed when 'X' is not an extended selection table.
#' @return Object 'X' with an additional column, 'signal.to.noise.ratio',
#' with the signal-to-noise ratio values.
#' @export
#' @name signal_to_noise_ratio
#' @details Signal-to-noise ratio (SNR) measures sound amplitude level in relation to ambient noise. Noise is measured on the background noise immediately before the test sound. A general margin in which ambient noise will be measured must be specified. Alternatively, a selection of ambient noise can be used as reference (see 'noise.ref' argument). When margins overlap with another sound nearby, SNR will be inaccurate, so margin length should be carefully considered. Any SNR less than or equal to one suggests background noise is equal to or overpowering the sound. The function will measure signal-to-noise ratio within the supplied frequency range (e.g. bandpass) of the reference signal ('bottom.freq' and 'top.freq' columns in 'X') by default (that is, when \code{bp = 'freq.range'}.
#' @examples
#' {
#'   # load example data
#'   data("degradation_est")
#'
#'   # create subset of data with only re-recorded files
#'   rerecorded_est <- degradation_est[degradation_est$sound.files != "master.wav", ]
#'
#'   # using measure ambient noise reference selections
#'   signal_to_noise_ratio(X = rerecorded_est, mar = 0.05, noise.ref = "custom")
#'
#'   # remove ambient selections
#'   rerecorded_est <- rerecorded_est[rerecorded_est$sound.id != "ambient", ]
#'
#'   # using margin for ambient noise of 0.05 and adjacent measure ambient noise reference
#'   signal_to_noise_ratio(X = rerecorded_est, mar = 0.05, noise.ref = "adjacent")
#' }
#'
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})
#' @seealso \code{\link{excess_attenuation}}
#' @references {
#' Araya-Salas, M. (2020). baRulho: baRulho: quantifying habitat-induced degradation of (animal) acoustic signals in R. R package version 1.0.2
#'
#' Dabelsteen, T., Larsen, O. N., & Pedersen, S. B. (1993). Habitat-induced degradation of sound signals: Quantifying the effects of communication sounds and bird location on blur ratio, excess attenuation, and signal-to-noise ratio in blackbird song. The Journal of the Acoustical Society of America, 93(4), 2206.
#'
#' Darden, SK, Pedersen SB, Larsen ON, & Dabelsteen T. (2008). Sound transmission at ground level in a short-grass prairie habitat and its implications for long-range communication in the swift fox *Vulpes velox*. The Journal of the Acoustical Society of America, 124(2), 758-766.
#' }
# last modification on nov-01-2019 (MAS)

signal_to_noise_ratio <-
  function(X,
           mar = NULL,
           parallel = NULL,
           cores = getOption("mc.cores", 1),
           pb = getOption("pb", TRUE),
           eq.dur = FALSE,
           noise.ref = "adjacent",
           type = 1,
           bp = "freq.range",
           output = NULL,
           hop.size = getOption("hop.size", 1),
           wl = getOption("wl", NULL),
           ovlp = getOption("ovlp", 0),
           path = getOption("sound.files.path", ".")) {
    # check arguments
    arguments <- as.list(base::match.call())

    # add objects to argument names
    for (i in names(arguments)[-1]) {
      arguments[[i]] <- get(i)
    }

    # check each arguments
    check_results <- check_arguments(fun = arguments[[1]], args = arguments)

    # report errors
    checkmate::reportAssertions(check_results)

    # get sampling rate
    sampling_rate <- warbleR::read_sound_file(X = X, index = 1, path = path, header = TRUE)$sample.rate

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

    # check sound.id column
    if (is.null(X$sound.id)) {
      if (noise.ref == "custom") {
        stop2("'sound.id' required when 'noise.ref == 'custom''")
      }

      X$sound.id <- "no.sound.id.column"
    }

    # check if 'ambient' is found in  sound.id column
    if (!any(X$sound.id %in% "ambient") &
      noise.ref == "custom") {
      stop2(
        "'ambient' selections must be contained in 'X' and labeled in a 'sound.id' column as 'ambient' when 'noise.ref = 'custom'"
      )
    }

    if (noise.ref == "custom" &
      any(sapply(unique(X$sound.files), function(x) {
        sum(X$sound.files == x &
          X$sound.id == "ambient")
      }) == 0)) {
      stop2(
        "Each sound file referenced in 'X' must have at least 1 'ambient' selection when 'noise.ref = 'custom'"
      )
    }

    # set clusters for windows OS
    if (Sys.info()[1] == "Windows" & cores > 1) {
      cl <-
        parallel::makePSOCKcluster(getOption("cl.cores", cores))
    } else {
      cl <- cores
    }

    # calculate all envelops with a apply function
    envs <-
      warbleR:::pblapply_wrblr_int(
        X = seq_len(nrow(X)),
        pbar = pb,
        cl = cl,
        FUN = function(y) {
          if (noise.ref == "custom") {
            # read sound clip
            signal <- warbleR::read_sound_file(X = X, index = y, path = path)

            # add band-pass frequency filter
            if (!is.null(bp)) {
              # filter to bottom and top freq range
              if (bp[1] == "freq.range") {
                bp <- c(X$bottom.freq[y], X$top.freq[y])
              }

              signal <-
                seewave::ffilter(
                  signal,
                  f = sampling_rate,
                  from = bp[1] * 1000,
                  ovlp = 0,
                  to = bp[2] * 1000,
                  bandpass = TRUE,
                  wl = wl,
                  output = "Wave"
                )
            }

            # get RMS for signal
            sig.env <-
              seewave::env(signal,
                f = sampling_rate,
                envt = "abs",
                plot = FALSE
              )

            bg.env <- NA
          }

          if (noise.ref == "adjacent") {
            # set margin to half of signal duration
            if (eq.dur) {
              mar <-
                (X$end[y] - X$start[y])
            } else if (is.null(mar)) {
              stop2("'mar' must be provided when 'eq.dur = FALSE'")
            }

            # Read sound files to get sample rate and length
            r <-
              warbleR::read_sound_file(
                X = X,
                index = y,
                from = X$start[y] - mar,
                to = X$end[y],
                header = TRUE,
                path = path
              )

            # reset time coordinates of sounds if lower than 0 o higher than duration
            stn <- X$start[y] - mar
            enn <- X$end[y] + mar
            mar1 <- mar

            if (stn < 0) {
              mar1 <- mar1 + stn
              stn <- 0
            }

            mar2 <- mar1 + X$end[y] - X$start[y]

            if (enn > r$samples / sampling_rate) {
              enn <- r$samples / sampling_rate
            }

            # read sound and margin
            noise_sig <-
              warbleR::read_sound_file(
                X = X,
                index = y,
                from = stn,
                to = enn,
                path = path
              )

            # add band-pass frequency filter
            if (!is.null(bp)) {
              # filter to bottom and top freq range
              if (bp[1] == "freq.range") {
                bp <- c(X$bottom.freq[y], X$top.freq[y])
              }

              noise_sig <-
                seewave::ffilter(
                  wave = noise_sig,
                  f = sampling_rate,
                  from = bp[1] * 1000,
                  ovlp = ovlp,
                  to = bp[2] * 1000,
                  bandpass = TRUE,
                  wl = wl,
                  output = "Wave"
                )
            }


            # read clip with sound
            signal <-
              seewave::cutw(noise_sig,
                from = mar1,
                to = mar2,
                f = sampling_rate
              )

            # get envelop for sound
            sig.env <-
              seewave::env(signal,
                f = sampling_rate,
                envt = "hil",
                plot = FALSE
              )

            # cut ambient noise before sound
            noise1 <-
              seewave::cutw(noise_sig,
                from = 0,
                to = mar1,
                f = sampling_rate
              )

            # get envelop for background noise
            bg.env <-
              seewave::env(noise1,
                f = sampling_rate,
                envt = "hil",
                plot = FALSE
              )
          }
          return(list(sig.env = sig.env, bg.env = bg.env))
        }
      )

    # add sound file selec column and names to envelopes (weird column name so it does not overwrite user columns)
    X$TEMP....y <-
      names(envs) <- paste(X$sound.files, X$selec, sep = "-")

    # calculate SNR
    X$signal.to.noise.ratio <- sapply(seq_len(nrow(X)), function(y) {
      if (X$sound.id[y] != "ambient") {
        suppressWarnings({
          # sound RMS
          sig_RMS <- seewave::rms(envs[[X$TEMP....y[y]]]$sig.env)

          # get reference ambient noise RMS
          if (noise.ref == "adjacent") {
            bg_RMS <- seewave::rms(envs[[X$TEMP....y[y]]]$bg.env)
          } else {
            # get envelopes from ambient selections
            bg_envs <-
              sapply(envs[X$TEMP....y[X$sound.files == X$sound.files[y] &
                X$sound.id == "ambient"]], "[", "sig.env")

            # get mean RMS from combined envelopes
            bg_RMS <- seewave::rms(unlist(sapply(bg_envs, as.vector)))
          }

          # Calculate signal-to-noise ratio
          if (type == 1) {
            snr <- 20 * log10(sig_RMS / bg_RMS)
          }

          if (type == 2) {
            snr <- 20 * log10((sig_RMS - bg_RMS) / bg_RMS)
          }
        })
      } else {
        snr <- NA
      } # return NA if current row is noise

      return(snr)
    })

    # remove temporary column
    X$TEMP....y <- NULL

    # remove sound.id column
    if (X$sound.id[1] == "no.sound.id.column") {
      X$sound.id <- NULL
    }

    # fix call if not a data frame
    if (!is.data.frame(X)) {
      attributes(X)$call <-
        base::match.call()
    } # fix call attribute


    return(X)
  }
