#' Measure reverberations as tail-to-signal ratio
#'
#' \code{tail_to_signal_ratio} measures reverberations as tail-to-signal ratio of sounds referenced in an extended selection table.
#' @inheritParams template_params
#' @param X Object of class 'data.frame', 'selection_table' or 'extended_selection_table' (the last 2 classes are created by the function \code{\link[warbleR]{selection_table}} from the warbleR package) with the reference to the test sounds . Must contain the following columns: 1) "sound.files": name of the .wav files, 2) "selec": unique selection identifier (within a sound file), 3) "start": start time and 4) "end": end time of selections, 5)  "bottom.freq": low frequency for bandpass and 6) "top.freq": high frequency for bandpass.
#' @param mar numeric vector of length 1. Specifies the margins adjacent to
#'   end of the sound over which to measure tail power.
#' @param hop.size A numeric vector of length 1 specifying the time window duration (in ms). Default is 1 ms, which is equivalent to ~45 wl for a 44.1 kHz sampling rate. Ignored if 'wl' is supplied.
#' @param tsr.formula Integer vector of length 1. Determine the formula to be used to calculate the tail-to-signal ratio (S = signal, T = tail, N = background noise):
#' \itemize{
#' \item \code{1}: ratio of T amplitude envelope quadratic mean to S amplitude envelope quadratic mean
#'  (\code{rms(env(T))/rms(env(S))}) as described by Dabelsteen et al. (1993).
#' \item \code{2}: ratio of T amplitude envelope quadratic mean to N amplitude envelope quadratic mean (\code{rms(env(T))/rms(env(N))}). N is measure in the margin right before the sound. So tsr.formula 2 actually measures tail-to-noise ratio.
#' }
#' @param wl A numeric vector of length 1 specifying the window length of the spectrogram, default
#' is NULL. Ignored if \code{bp = NULL}. If supplied, 'hop.size' is ignored.
#' Note that lower values will increase time resolution, which is more important for amplitude calculations.
#' @param ovlp Numeric vector of length 1 specifying the percentage of overlap between two
#'   consecutive windows, as in \code{\link[seewave]{spectro}}. Default is 0. Only used for bandpass filtering.
#' @return Object 'X' with an additional column, 'tail.to.signal.ratio',
#' with the tail-to-signal ratio values (in dB).
#' @export
#' @name tail_to_signal_ratio
#' @details Tail-to-signal ratio (TSR) measures the ratio of power in the tail of reverberations to that in the test sound. A general margin in which reverberation tail will be measured must be specified. The function will measure TSR within the supplied frequency range (e.g. bandpass) of the reference sound ('bottom.freq' and 'top.freq' columns in 'X'). Two methods for computing reverberations are provided (see 'tsr.formula' argument). Note that 'tsr.formula' 2 is not equivalent to the original description of TSR in Dabelsteen et al. (1993) and  is better referred to as tail-to-noise ratio. Tail-to-signal ratio values are typically negative as signals tend to have higher power than that in the reverberating tail. TSR can be ~0 when both tail and signal have very low amplitude.
#' @examples {
#'   # load example data
#'
#'   data("test_sounds_est")
#'
#'   # set global options
#'   options(pb = FALSE)
#'
#'   # using margin for noise of 0.01
#'   tsr <- tail_to_signal_ratio(X = test_sounds_est, mar = 0.01)
#'
#'   # use tsr.formula 2 which is equivalent to tail-to-noise ratio
#'   tsr <- tail_to_signal_ratio(X = test_sounds_est, mar = 0.01, tsr.formula = 2)
#' }
#'
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})
#' @family quantify degradation
#' @seealso \code{\link{excess_attenuation}}
#' @references {
#' Araya-Salas M., E. Grabarczyk, M. Quiroz-Oliva, A. Garcia-Rodriguez, A. Rico-Guevara. (2023), baRulho: an R package to quantify degradation in animal acoustic signals .bioRxiv 2023.11.22.568305.
#'
#' Darden, SK, Pedersen SB, Larsen ON, & Dabelsteen T. (2008). Sound transmission at ground level in a short-grass prairie habitat and its implications for long-range communication in the swift fox *Vulpes velox*. The Journal of the Acoustical Society of America, 124(2), 758-766.
#'
#' Mathevon, N., Dabelsteen, T., & Blumenrath, S. H. (2005). Are high perches in the blackcap Sylvia atricapilla song or listening posts? A sound transmission study. The Journal of the Acoustical Society of America, 117(1), 442-449.
#' }

tail_to_signal_ratio <- function(X,
                                 mar,
                                 cores = getOption("mc.cores", 1),
                                 pb = getOption("pb", TRUE),
                                 tsr.formula = 1,
                                 bp = "freq.range",
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
  X$tail.to.signal.ratio <-
    unlist(warbleR:::pblapply_wrblr_int(X = seq_len(nrow(X)), pbar = pb, cl = cl, function(y) {
      if (X$sound.id[y] != "ambient") {
        # Read sound files to get sample rate and length
        r <-
          warbleR::read_sound_file(
            X = X,
            index = y,
            from = 0,
            to = Inf,
            header = TRUE,
            path = path
          )
        
        
        # reset time coordinates of sounds if higher than duration
        enn <- X$end[y] + mar
        if (enn > r$samples / sampling_rate) {
          enn <- r$samples / sampling_rate
        }
        
        # read sound and margin
        tail.wv <-
          warbleR::read_sound_file(
            X = X,
            index = y,
            from = X$end[y],
            to = enn,
            path = path
          )
        
        # read sound
        if (tsr.formula == 1) {
          signal <-
            warbleR::read_sound_file(X = X,
                                     index = y,
                                     path = path)
        }
        
        # read background noise right before the sound
        if (tsr.formula == 2) {
          signal <-
            warbleR::read_sound_file(
              X = X,
              index = y,
              from = X$start[y] - mar,
              to = X$start[y],
              path = path
            )
        }
        
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
              ovlp = ovlp,
              to = bp[2] * 1000,
              bandpass = TRUE,
              wl = wl,
              output = "Wave"
            )
          
          tail.wv <-
            seewave::ffilter(
              tail.wv,
              f = sampling_rate,
              from = bp[1] * 1000,
              ovlp = ovlp,
              to = bp[2] * 1000,
              bandpass = TRUE,
              wl = wl,
              output = "Wave"
            )
        }
        
        
        # get RMS for sound (or noise if tsr.formula 2)
        sig.env <-
          seewave::env(signal,
                       f = sampling_rate,
                       envt = "hil",
                       plot = FALSE)
        
        # get RMS for background noise
        tail.env <-
          seewave::env(tail.wv,
                       f = sampling_rate,
                       envt = "hil",
                       plot = FALSE)
        
        # sound (or noise) RMS
        sig_RMS <- seewave::rms(sig.env)
        
        # convert to 0.0001 if sig_RMS is 0 to avoid errors in STR measurements
        if (sig_RMS == 0)
          sig_RMS <- 0.0001
        
        # get tail RMS
        tail_RMS <- seewave::rms(tail.env)
        
        # convert to 0.0001 if tail_RMS is 0 to avoid errors in STR measurements
        if (tail_RMS == 0)
          tail_RMS <- 0.0001
        
        # Calculate tail.to.signal ratio
        str <- tail_RMS / sig_RMS
        
        str <- 10 * log(str)
      } else {
        str <- NA
      }
      
      return(str)
    }))
  
  # fix call if not a data frame
  if (!is.data.frame(X)) {
    attributes(X)$call <-
      base::match.call()
  } # fix call attribute
  
  # remove sound.id column if was not  found in X
  if (X$sound.id[1] == "no.sound.id.column") {
    X$sound.id <- NULL
  }
  
  return(X)
}
