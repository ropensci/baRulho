#' Measure attenuation as signal-to-noise ratio
#'
#' \code{signal_to_noise_ratio} measures attenuation as signal-to-noise ratio of sounds referenced in an extended selection table.
#' @inheritParams template_params
#' @param X Object of class 'data.frame', 'selection_table' or 'extended_selection_table' (the last 2 classes are created by the function \code{\link[warbleR]{selection_table}} from the warbleR package) with the test sound files' annotations (typically the output of \code{\link{align_test_files}}). Must contain the following columns: 1) "sound.files": name of the .wav files, 2) "selec": unique selection identifier (within a sound file), 3) "start": start time and 4) "end": end time of selections, 5)  "bottom.freq": low frequency for bandpass, 6) "top.freq": high frequency for bandpass and 7) "sound.id": ID of sounds used to identify counterparts across distances (only needed for "custom" noise reference, see "noise.ref" argument). If the "sound.id" column is supplied then SNR is only computed for those rows with a sound.id different from "ambient", "start_marker" or "end_marker". 
#' @param mar numeric vector of length 1. Specifies the margins adjacent to
#'   the start point of the annotation over which to measure ambient noise.
#' @param eq.dur Logical. Controls whether the ambient noise segment that is measured has the same duration
#' to that of the sound (if \code{TRUE}. Default is \code{FALSE}). If \code{TRUE} then 'mar' and 'noise.ref' arguments are ignored.
#' @param snr.formula Integer vector of length 1. Selects the formula to be used to calculate the signal-to-noise ratio (S = signal
#' , N = background noise):
#' \itemize{
#' \item \code{1}: ratio of S amplitude envelope root mean square to N amplitude envelope root mean square
#'  (\code{20 * log10(rms(env(S))/rms(env(N)))}) as described by Darden (2008).
#' \item \code{2}: ratio of the difference between S amplitude envelope root mean square and N amplitude envelope root mean square to N amplitude envelope root mean square (\code{20 * log10((rms(env(S)) - rms(env(N)))/rms(env(N)))}, as described by Dabelsteen et al (1993).
#' }
#' @param hop.size A numeric vector of length 1 specifying the time window duration (in ms). Default is 1 ms, which is equivalent to ~45 wl for a 44.1 kHz sampling rate. Ignored if 'wl' is supplied. Can be set globally for the current R session via the "hop.size" option (see \code{\link[base]{options}}).
#' @param wl A numeric vector of length 1 specifying the window length of the spectrogram, default
#' is NULL. Ignored if \code{bp = NULL}. If supplied, 'hop.size' is ignored. Note that lower values will increase time resolution, which is more important for amplitude ratios calculations.
#' @param ovlp Numeric vector of length 1 specifying the percentage of overlap between two
#'   consecutive windows, as in \code{\link[seewave]{spectro}}. Default is 0. Only used for bandpass filtering. Can be set globally for the current R session via the "ovlp" option (see \code{\link[base]{options}}).
#' @return Object 'X' with an additional column, 'signal.to.noise.ratio',
#' with the signal-to-noise ratio values (in dB).
#' @export
#' @name signal_to_noise_ratio
#' @details Signal-to-noise ratio (SNR) measures sound amplitude level in relation to ambient noise. Noise is measured on the background noise immediately before the test sound. A general margin in which ambient noise will be measured must be specified. Alternatively, a selection of ambient noise can be used as reference (see 'noise.ref' argument). When margins overlap with another sound nearby, SNR will be inaccurate, so margin length should be carefully considered. Any SNR less than or equal to one suggests background noise is equal to or overpowering the sound. The function will measure signal-to-noise ratio within the supplied frequency range (e.g. bandpass) of the reference signal ('bottom.freq' and 'top.freq' columns in 'X') by default (that is, when \code{bp = 'freq.range'}. SNR can be ~0 when both tail and signal have very low amplitude.
#' @examples {
#'   # load example data
#'   data("test_sounds_est")
#'
#'   # using measure ambient noise reference selections
#'   signal_to_noise_ratio(X = test_sounds_est, mar = 0.05, noise.ref = "custom")
#'
#'   # using margin for ambient noise of 0.05 and adjacent measure ambient noise reference
#'   signal_to_noise_ratio(X = test_sounds_est, mar = 0.05, noise.ref = "adjacent")
#' }
#'
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})
#' @family quantify degradation
#' @seealso \code{\link{excess_attenuation}}
#' @references 
#' Araya-Salas, M., Grabarczyk, E. E., Quiroz-Oliva, M., Garcia-Rodriguez, A., & Rico-Guevara, A. (2025). Quantifying degradation in animal acoustic signals with the R package baRulho. Methods in Ecology and Evolution, 00, 1-12. https://doi.org/10.1111/2041-210X.14481
#' Holland J, Dabelsteen T, Pedersen SB, Paris AL (2001) Potential ranging cues contained within the energetic pauses of transmitted wren song. Bioacoustics 12(1):3-20.
#' Darden, SK, Pedersen SB, Larsen ON, & Dabelsteen T. (2008). Sound transmission at ground level in a short-grass prairie habitat and its implications for long-range communication in the swift fox *Vulpes velox*. The Journal of the Acoustical Society of America, 124(2), 758-766.

signal_to_noise_ratio <-
  function(X,
           mar = NULL,
           cores = getOption("mc.cores", 1),
           pb = getOption("pb", TRUE),
           eq.dur = FALSE,
           noise.ref = c("adjacent", "custom"),
           snr.formula = 1,
           bp = "freq.range",
           hop.size = getOption("hop.size", 1),
           wl = getOption("wl", NULL),
           ovlp = getOption("ovlp", 0),
           path = getOption("sound.files.path", ".")) {
    
    # assign a value to noise.ref
    noise.ref <- rlang::arg_match(noise.ref)
    
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
    
    # check sound.id column
    if (is.null(X$sound.id)) {
      if (noise.ref == "custom") {
        .stop("'sound.id' required when 'noise.ref == 'custom''")
      }
      
      X$sound.id <- "no.sound.id.column"
    }
    
    # check if 'ambient' is found in  sound.id column
    if (!any(X$sound.id %in% "ambient") &
        noise.ref == "custom") {
      .stop(
        "'ambient' selections must be contained in 'X' and labeled in a 'sound.id' column as 'ambient' when 'noise.ref = 'custom'"
      )
    }
    
    if (noise.ref == "custom" & warbleR::is_extended_selection_table(X)){
      if (!attributes(X)$by.song[[1]]) {
        .stop(
          "if 'X' is an extended selection table, it should be created 'by.song' when 'noise.ref = 'custom'' (see 'https://marce10.github.io/warbleR/articles/b_annotation_data_format.html#by-element-vs-by-song-extended-selection-tables')"
        )
      }
    }
    
    if (noise.ref == "custom" &
        any(vapply(unique(X$sound.files), function(x) {
          sum(X$sound.files == x &
              X$sound.id == "ambient")
        }, FUN.VALUE = numeric(1)) == 0)) {
      .stop(
        "Each sound file referenced in 'X' must have at least 1 'ambient' selection when 'noise.ref = 'custom'"
      )
    }
    
    # 'mar' is needed when not using equal duration
    if (!eq.dur & is.null(mar) & noise.ref == "adjacent") {
      .stop("'mar' must be supplied when 'eq.dur = FALSE'")
    }
    
    # set clusters for windows OS
    if (Sys.info()[1] == "Windows" & cores > 1) {
      cl <-
        parallel::makePSOCKcluster(cores)
    } else {
      cl <- cores
    }
    
    # get index number
    # exclude "ambient" sounds if "sound.id" column was supplied
    target_rows <- if (!is.null(X$sound.id)){
      which(!X$sound.id %in% c("ambient", "start_marker", "end_marker")) } else {
        seq_len(nrow(X))
      }
    
    # target rows for rms differ when reference is "ambient"
    target_rows_rms <- if (noise.ref == "custom") {
      which(!X$sound.id %in% c("start_marker", "end_marker"))
    } else {
      target_rows
      }

    
    # calculate all RMS of envelopes with a apply function
    rms_list <-
      warbleR:::.pblapply(
        X = target_rows_rms,
        pbar = pb,
        cl = cl,
        message = "computing signal-to-noise ratio",
        total = 1,
        FUN = .rms,
        Y = X,
        mar = mar,
        noise.ref = noise.ref,
        sampling_rate = sampling_rate,
        wl = wl,
        path = path,
        eq.dur = eq.dur,
        ovlp = ovlp,
        bp = bp
      )
    
    # add sound file selec column and names to envelopes (weird column name so it does not overwrite user columns)
    X$.y <- paste(X$sound.files, X$selec, sep = "-")
    names(rms_list) <- paste(X$sound.files[target_rows_rms], X$selec[target_rows_rms], sep = "-")
    
    # calculate SNR
    X$signal.to.noise.ratio[target_rows] <-
      vapply(
        X = target_rows,
        FUN = .snr,
        W = X,
        noise.ref = noise.ref,
        type = snr.formula,
        rms_list = rms_list,
        FUN.VALUE = numeric(1)
      )
    
    # remove temporary column
    X$.y <- NULL
    
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
