#' Measure blur ratio in the time domain
#'
#' \code{blur_ratio} measures blur ratio in sounds referenced in an extended selection table.
#' @usage blur_ratio(X, parallel = NULL, cores = getOption("mc.cores", 1),
#' pb = getOption("pb", TRUE),  
#' env.smooth = getOption("env.smooth", 200),
#' output = NULL, envelopes = FALSE, 
#' hop.size = getOption("hop.size", 11.6), wl = getOption("wl", NULL),
#' ovlp = getOption("ovlp", 70), path = getOption("sound.files.path", "."))
#' @param X The output of \code{\link{set_reference_sounds}} which is an object of class 'data.frame', 'selection_table' or 'extended_selection_table' (the last 2 classes are created by the function \code{\link[warbleR]{selection_table}} from the warbleR package) with the reference to the sounds in the master sound file. Must contain the following columns: 1) "sound.files": name of the .wav files, 2) "selec": unique selection identifier (within a sound file), 3) "start": start time and 4) "end": end time of selections, 5)  "bottom.freq": low frequency for bandpass, 6) "top.freq": high frequency for bandpass, 7) "sound.id": ID of sounds used to identify counterparts across distances and 8) "reference": identity of sounds to be used as reference for each test sound (row). See \code{\link{set_reference_sounds}} for more details on the structure of 'X'.
#' @param parallel DEPRECATED. Use 'cores' instead.
#' @param cores Numeric vector of length 1. Controls whether parallel computing is applied by specifying the number of cores to be used. Default is 1 (i.e. no parallel computing).
#' @param pb Logical argument to control if progress bar is shown. Default is \code{TRUE}.
#' @param env.smooth Numeric vector of length 1 determining the length of the sliding window (in amplitude samples) used for a sum smooth for amplitude envelope calculation (used internally by \code{\link[seewave]{env}}). Default is 200.
#' @param output DEPRECATED. Now the output format mirrors the class of the input 'X'. To obtain the amplitude envelopes use 'envelopes = TRUE'.
#' @param envelopes Logical to control if envelopes are returned (as attributes, 'attributes(X)$envelopes'). Default is \code{FALSE}.
#' @param hop.size A numeric vector of length 1 specifying the time window duration (in ms). Default is 11.6 ms, which is equivalent to 512 wl for a 44.1 kHz sampling rate. Ignored if 'wl' is supplied.
#' @param wl A numeric vector of length 1 specifying the window length of the spectrogram, default
#' is NULL. If supplied, 'hop.size' is ignored.
#' @param ovlp Numeric vector of length 1 specifying the percent overlap between two
#'   consecutive windows, as in \code{\link[seewave]{spectro}}. Default is 70. Used for applying bandpass filtering.
#' @param path Character string containing the directory path where the sound files are found. Only needed when 'X' is not an extended selection table.
#' @return Object 'X' with two additional columns, 'reference' and 'blur.ratio', containing containing the id of the sound used as reference and the computed blur ratio values, respectively. If \code{envelopes = TRUE} the output would include amplitude envelopes for all sounds as attributes ('attributes(X)$envelopes').
#' @export
#' @name blur_ratio
#' @details Blur ratio measures the degradation of sound as a change in sound power in the time domain as described by Dabelsteen et al (1993). Low values indicate low degradation of sounds. The function measures the blur ratio on sounds in which a reference playback has been re-recorded at different distances. Blur ratio is measured as the mismatch between amplitude envelopes (expressed as probability mass functions) of the reference sound and the re-recorded sound. By converting envelopes to probability mass functions the effect of energy attenuation is removed, focusing the analysis on the modification of the envelope shape. The function compares each sound to the corresponding reference sound within the supplied frequency range (e.g. bandpass) of the reference sound ('bottom.freq' and 'top.freq' columns in 'X'). The 'sound.id' column must be used to tell the function to only compare sounds belonging to the same category (e.g. song-types). Two methods for setting the experimental design are provided. All wave objects in the extended selection table must have the same sampling rate so the length of envelopes is comparable.
#' @seealso \code{\link{envelope_correlation}}, \code{\link{spectrum_blur_ratio}}
#' @examples {
#'   # load example data
#'   data("test_sounds_est")
#'   
#'  # add reference to X
#'  X <- set_reference_sounds(X = test_sounds_est)
#'   blur_ratio(X = X)
#'
#'   # using method 2
#' X <- set_reference_sounds(X = test_sounds_est, method = 2)
#'   # blur_ratio(X = X)
#'
#'   # get envelopes
#'   br <- blur_ratio(X = X, envelopes = TRUE)
#'   envs <- attributes(br)$envelopes
#'
#'   # make distance a factor for plotting
#'   envs$distance <- as.factor(envs$distance)
#'
#'   \dontrun{
#'   # plot
#'   ggplot(envs, aes(x= time, y = amp, col = distance)) +
#'   geom_line() + facet_wrap(~ sound.id) +
#'   scale_color_manual(values = viridis::viridis(4)) +
#'   labs(x = "Time (s)", y = "Amplitude (PMF)") +
#'   theme_classic()
#'   }
#' }
#'
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})
#' @references {
#' Dabelsteen, T., Larsen, O. N., & Pedersen, S. B. (1993). Habitat-induced degradation of sound signals: Quantifying the effects of communication sounds and bird location on blur ratio, excess attenuation, and signal-to-noise ratio in blackbird song. The Journal of the Acoustical Society of America, 93(4), 2206.
#'
#' Araya-Salas, M. (2020). baRulho: baRulho: quantifying habitat-induced degradation of (animal) acoustic signals in R. R package version 1.0.2
#' }

blur_ratio <-
  function(X,
           parallel = NULL,
           cores = getOption("mc.cores", 1),
           pb = getOption("pb", TRUE),
           env.smooth = getOption("env.smooth", 200),
           output = NULL,
           envelopes = FALSE,
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
    steps <- if (envelopes) 3 else 2
    
    # get sampling rate assuming is the same for all sound files
    sampling.rate <- read_sound_file(
      X,
      index = 1,
      header = TRUE,
      path = path
    )$sample.rate
    
    # adjust wl based on hope.size
    if (is.null(wl)) {
      wl <-
        round(
          sampling.rate * hop.size / 1000,
          0
        )
    }
    
    # make wl even if odd
    if (!(wl %% 2) == 0)
      wl <- wl + 1
    
    # set clusters for windows OS
    if (Sys.info()[1] == "Windows" & cores > 1) {
      cl <- parallel::makePSOCKcluster(getOption("cl.cores", cores))
    } else {
      cl <- cores
    }
    
    # add sound file selec colums to X (weird column name so it does not overwrite user columns)
    X$.sgnl.temp <- paste(X$sound.files, X$selec, sep = "-")
    
    # get names of envelopes involved (those as test with reference or as reference)
    target_sgnl_temp <- unique(c(X$.sgnl.temp[!is.na(X$reference)], X$reference[!is.na(X$reference)]))
    
    # print message
    if (pb)
      write(file = "", x = paste0("Computing amplitude envelopes (step 1 out of ", steps, "):"))
    
    # calculate all envelops apply function
    envs <-
      warbleR:::pblapply_wrblr_int(
        pbar = pb,
        X = target_sgnl_temp,
        cl = cl,
        FUN = function(x, ssmth = env.smooth, ovl = ovlp, Q = X, wln = wl, pth = path){
                       env_FUN(X = Q, y = x, env.smooth = ssmth, ovlp = ovl, wl = wln, path = pth) 
}
)
    
    # add sound file selec column as names to envelopes
    names(envs) <- target_sgnl_temp
    
    # set options to run loops
    if (pb)
      write(file = "", x = paste0("Computing blur ratio (step 2 out of ", steps, "):"))
    
    # get blur ratio
    # calculate all envelops apply function
    X$blur.ratio <-
      unlist(warbleR:::pblapply_wrblr_int(
        pbar = pb,
        X = seq_len(nrow(X)),
        cl = cl,
        FUN = function(x,
                       Q = X,
                       nvs = envs,
                       wle = wl,
                       ovp = ovlp,
                       sr = sampling.rate) {
          blur_FUN(
            x,
            X = Q,
            envs = nvs,
            ovlp = ovp,
            wl = wle,
            sampling.rate = sr
          )
        }
      )
      )
    

    # remove temporal column
    X$.sgnl.temp <- NULL
    
    # convert to list instead of extended selection table, add envelopes
    if (envelopes) {
      
      if (pb)
        write(file = "", x = "Saving envelopes (step 3 out of 3):")
      
      # get envelopes in a data frame
      env.dfs <-
        warbleR:::pblapply_wrblr_int(pbar = pb, 1:length(envs), cl = cl, function(y) {
          # extract 1 envelope
          x <- envs[[y]]
          
          # convert envelopes to PMF (probability mass function)
          x <- x / sum(x)
          
          # put in data framme
          out <-
            data.frame(
              sound = names(envs)[y],
              sound.id = X$sound.id[paste(X$sound.files, X$selec, sep = "-") == names(envs)[y]],
              distance = X$distance[paste(X$sound.files, X$selec, sep = "-") == names(envs)[y]],
              time = seq(
                from = 0,
                to = length(x) / (attr(X, "check.results")$sample.rate[1] * 1000),
                along.with = x
              ),
              amp = x
            )
          
          return(out)
        })
      
      # put together in a single data frame
      env.df <- do.call(rbind, env.dfs)
      
      # add envelopes as attributes
      attributes(X)$envelopes <- env.df
    }
    
    # return data frame
    if (is_extended_selection_table(X) | is_selection_table(X)) {
      attributes(X)$call <- base::match.call()
    }
    return(X)
  }
