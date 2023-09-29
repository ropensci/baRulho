#' Measure blur ratio in the frequency domain
#'
#' \code{spectrum_blur_ratio} measures blur ratio of frequency spectra from sounds referenced in an extended selection table.
#' @usage spectrum_blur_ratio(X, parallel = NULL, cores = getOption("mc.cores", 1),
#' pb = getOption("pb", TRUE),
#' spec.smooth = getOption("spec.smooth", 5),
#' output = NULL, spectra = FALSE, res = 150,
#' hop.size = getOption("hop.size", 11.6), wl = getOption("wl", NULL),
#'  ovlp = getOption("ovlp", 70), path = getOption("sound.files.path", "."))
#' @param X The output of \code{\link{set_reference_sounds}} which is an object of class 'data.frame', 'selection_table' or 'extended_selection_table' (the last 2 classes are created by the function \code{\link[warbleR]{selection_table}} from the warbleR package) with the reference to the sounds in the master sound file. Must contain the following columns: 1) "sound.files": name of the .wav files, 2) "selec": unique selection identifier (within a sound file), 3) "start": start time and 4) "end": end time of selections, 5)  "bottom.freq": low frequency for bandpass, 6) "top.freq": high frequency for bandpass, 7) "sound.id": ID of sounds used to identify counterparts across distances and 8) "reference": identity of sounds to be used as reference for each test sound (row). See \code{\link{set_reference_sounds}} for more details on the structure of 'X'.
#' @param parallel DEPRECATED. Use 'cores' instead.
#' @param cores Numeric vector of length 1. Controls whether parallel computing is applied by specifying the number of cores to be used. Default is 1 (i.e. no parallel computing).
#' @param pb Logical argument to control if progress bar is shown. Default is \code{TRUE}.
#' @param spec.smooth Numeric vector of length 1 determining the length of the sliding window used for a sum smooth for power spectrum calculation (in kHz). Default is 5.
#' @param output DEPRECATED. Now the output format mirrors the class of the input 'X'. To obtain the power spectra use 'spectra = TRUE'.
#' @param spectra Logical to control if power spectra are returned (as attributes). Default is \code{FALSE}.
#' @param res Numeric argument of length 1. Controls image resolution. Default is 150 (faster) although 300 - 400 is recommended for publication/presentation quality.
#' @param hop.size A numeric vector of length 1 specifying the time window duration (in ms). Default is 11.6 ms, which is equivalent to 512 wl for a 44.1 kHz sampling rate. Ignored if 'wl' is supplied.
#' @param wl A numeric vector of length 1 specifying the window length of the spectrogram, default
#' is NULL. If supplied, 'hop.size' is ignored. Applied to both spectra and spectrograms on image files.
#' @param ovlp Numeric vector of length 1 specifying the percent overlap between two
#'   consecutive windows, as in \code{\link[seewave]{spectro}}. Default is 70. Applied to both spectra and spectrograms on image files.
#' @param path Character string containing the directory path where the sound files are found. Only needed when 'X' is not an extended selection table.
#' @return Object 'X' with an additional column, 'spectrum.blur.ratio', containing the computed spectrum blur ratio values. If \code{spectra = TRUE} the output would include power spectra for all sounds as attributes ('attributes(X)$spectra').
#' @export
#' @name spectrum_blur_ratio
#' @details Spectral blur ratio measures the degradation of sound as a function of the change in sound power in the frequency domain, analogous to the blur ratio proposed by Dabelsteen et al (1993) for the time domain (and implemented in \code{\link{blur_ratio}}). Low values indicate low degradation of sounds. The function measures the blur ratio of spectra from sounds in which a reference playback has been re-recorded at different distances. Spectral blur ratio is measured as the mismatch between power spectra (expressed as probability density functions) of the reference sound and the re-recorded sound. The function compares each sound type to the corresponding reference sound. The 'sound.id' column must be used to tell the function to only compare sounds belonging to the same category (e.g. song-types). Two methods for setting the experimental design are provided. All wave objects in the extended selection table must have the same sampling rate so the length of spectra is comparable. The function uses \code{\link[seewave]{spec}} internally to compute power spectra.
#' @examples {
#'   # load example data
#'   data("test_sounds_est")
#'   
#' # add reference to X
#' X <- set_reference_sounds(X = test_sounds_est)
#'
#'   # get spetrum blur ratio
#'   spectrum_blur_ratio(X = X)
#'
#'   # using method 2
#'    X <- set_reference_sounds(X = test_sounds_est, method = 2)
#'   spectrum_blur_ratio(X = X)
#'
#'   # get power spectra
#'   sbr <- spectrum_blur_ratio(X = X, spectra = TRUE)
#'   spctr <- attributes(sbr)$spectra
#'
#'   # make distance a factor for plotting
#'   spctr$distance <- as.factor(spctr$distance)
#'
#'   \dontrun{
#'   # plot
#'   ggplot(spctr[spctr$freq > 0.3, ], aes(y = amp, x = freq,
#'   col = distance)) +
#'   geom_line() +
#'   facet_wrap(~sound.id) +
#'   scale_color_viridis_d(alpha = 0.7) +
#'   labs(x = "Frequency (kHz)", y = "Amplitude (PMF)") +
#'   coord_flip() +
#'   theme_classic()
#'   }
#' }
#'
#' @seealso \code{\link{blur_ratio}}
#' @family quantify degradation
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})
#' @references {
#' Dabelsteen, T., Larsen, O. N., & Pedersen, S. B. (1993). Habitat-induced degradation of sound signals: Quantifying the effects of communication sounds and bird location on blur ratio, excess attenuation, and signal-to-noise ratio in blackbird song. The Journal of the Acoustical Society of America, 93(4), 2206.
#'
#' Araya-Salas, M. (2020). baRulho: baRulho: quantifying degradation of (animal) acoustic signals in R. R package version 1.0.2
#' }

spectrum_blur_ratio <-
  function(X,
           parallel = NULL,
           cores = getOption("mc.cores", 1),
           pb = getOption("pb", TRUE),
           spec.smooth = getOption("spec.smooth", 5),
           output = NULL,
           spectra = FALSE,
           res = 150,
           hop.size = getOption("hop.size", 11.6),
           wl = getOption("wl", NULL),
           ovlp = getOption("ovlp", 70),
           path = getOption("sound.files.path", ".")) {
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
    
    # total number of steps depending on whether envelopes are returned
    steps <- if (spectra) 3 else 2
    
    # get sampling rate
    sampling_rate <-
      warbleR::read_sound_file(
        X = X,
        index = 1,
        path = path,
        header = TRUE
      )$sample.rate
    
    
    # adjust wl based on hope.size
    if (is.null(wl)) {
      wl <-
        round(sampling_rate * hop.size / 1000,
              0)
    }
    
    # make wl even if odd
    if (!(wl %% 2) == 0) {
      wl <- wl + 1
    }
    
    # add sound file selec colums to X (weird column name so it does not overwrite user columns)
    X$.sgnl.temp <- paste(X$sound.files, X$selec, sep = "-")
    
    # get names of envelopes involved (those as test with reference or as reference)
    target_sgnl_temp <-
      unique(c(X$.sgnl.temp[!is.na(X$reference)], X$reference[!is.na(X$reference)]))
    
    # set clusters for windows OS
    if (Sys.info()[1] == "Windows" & cores > 1) {
      cl <-
        parallel::makePSOCKcluster(getOption("cl.cores", cores))
    } else {
      cl <- cores
    }
    
    # print message
    if (pb)
      write(
        file = "",
        x = paste0("Computing power spectra (step 1 out of ", steps, "):")
      )
    
    # calculate all spectra apply function
    specs <-
      warbleR:::pblapply_wrblr_int(
        pbar = pb,
        X = target_sgnl_temp,
        cl = cl,
        FUN = function(y,
                       ssmth = spec.smooth,
                       wln = wl,
                       Q = X,
                       pth = path) {
          spctr_FUN(
            y,
            spec.smooth = ssmth,
            wl = wln,
            X = Q,
            path = pth
          )
        }
      )
    
    # add sound file selec names to spectra
    names(specs) <- target_sgnl_temp
    
    # print second message
    if (pb)
      write(
        file = "",
        x = paste0("Computing spectrum blur ratio (step 2 out of ", steps, "):")
      )
    
    # get blur ratio
    # calculate all spectra apply function
    X$spectrum.blur.ratio <-
      unlist(warbleR:::pblapply_wrblr_int(
        X = seq_len(nrow(X)),
        pbar = pb,
        cl = cl,
        FUN = function(x, 
                       Q = X,
                       wle = wl,
                       ovp = ovlp,
                       spcs = specs,
                       sr = sampling_rate) {
          blur_sp_FUN(x,
                      X = Q,
                      ovlp = ovp,
                      wl = wle,
                      specs = spcs,
                      sampling_rate = sr)
        }
      ))
    
    # remove temporal columns
    X$.sgnl.temp <- NULL
    
    if (pb)
      write(file = "", x = "Saving envelopes (step 3 out of 3):")
    
    
    # convert to list instead of extended selection table, add envelopes
    if (spectra) {
      spec.dfs <-
        warbleR:::pblapply_wrblr_int(
          X = seq_along(specs),
          cl = cl,
          pbar = pb,
          FUN = function(y) {
            # extract 1 envelope
            x <- specs[[y]]
            
            # convert envelopes to PMF (probability mass function)
            x[, 2] <- x[, 2] / sum(x[, 2])
            
            # put in data framme
            out <-
              data.frame(
                sound = names(specs)[y],
                sound.id = X$sound.id[paste(X$sound.files, X$selec, sep = "-") == names(specs)[y]],
                distance = X$distance[paste(X$sound.files, X$selec, sep = "-") == names(specs)[y]],
                freq = x[, 1],
                amp = x[, 2]
              )
            
            return(out)
          }
        )
      
      # put together in a single data frame
      spec.df <- do.call(rbind, spec.dfs)
      
      # add envelopes as attributes
      attributes(X)$spectra <- spec.df
    }
    
    # return data frame
    if (is_extended_selection_table(X) | is_selection_table(X)) {
      attributes(X)$call <- base::match.call()
    }
    return(X)
  }
