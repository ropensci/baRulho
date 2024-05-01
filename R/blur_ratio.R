#' Measure blur ratio in the time domain
#'
#' \code{blur_ratio} measures blur ratio in sounds referenced in an extended selection table.
#' @inheritParams template_params
#' @param env.smooth Numeric vector of length 1 determining the length of the sliding window (in amplitude samples) used for a sum smooth for amplitude envelope calculation (used internally by \code{\link[seewave]{env}}). Default is 200.
#' @param envelopes Logical to control if envelopes are returned (as attributes, 'attributes(X)$envelopes'). Default is \code{FALSE}.
#' @param ovlp Numeric vector of length 1 specifying the percentage of overlap between two
#'   consecutive windows, as in \code{\link[seewave]{spectro}}. Default is 70. Used for applying bandpass filtering.
#' @param n.samples Numeric vector of length 1 specifying the number of amplitude samples to use for representing amplitude envelopes. Default is 100. If null the raw amplitude envelope is used (note that this can result in high RAM memory usage for large data sets). Amplitude envelope values are interpolated using \code{\link[stats]{approx}}.
#' @return Object 'X' with an additional column,  'blur.ratio', containing the computed blur ratio values. If \code{envelopes = TRUE} the output would include amplitude envelopes for all sounds as attributes ('attributes(X)$envelopes').
#' @export
#' @name blur_ratio
#' @details Blur ratio measures the degradation of sound as a change in sound power in the time domain as described by Dabelsteen et al (1993). Low values indicate low degradation of sounds. The function measures the blur ratio on sounds in which a reference playback has been re-recorded at different distances. Blur ratio is measured as the mismatch between amplitude envelopes (expressed as probability mass functions) of the reference sound and the re-recorded sound. By converting envelopes to probability mass functions the effect of energy attenuation is removed, focusing the analysis on the modification of the envelope shape. The function compares each sound to the corresponding reference sound within the supplied frequency range (e.g. bandpass) of the reference sound ('bottom.freq' and 'top.freq' columns in 'X'). The 'sound.id' column must be used to tell the function to only compare sounds belonging to the same category (e.g. song-types). Two methods for setting the experimental design are provided. All wave objects in the extended selection table must have the same sampling rate so the length of envelopes is comparable.
#' @family quantify degradation
#' @seealso \code{\link{envelope_correlation}}, \code{\link{spectrum_blur_ratio}}
#' @examples {
#'   # load example data
#'   data("test_sounds_est")
#'
#'  # add reference to X
#'   X <- set_reference_sounds(X = test_sounds_est)
#'   blur_ratio(X = X)
#'
#'   # using method 2
#' X <- set_reference_sounds(X = test_sounds_est, method = 2)
#'   blur_ratio(X = X)
#'
#'   # get envelopes
#'   br <- blur_ratio(X = X, envelopes = TRUE)
#'   envs <- attributes(br)$envelopes
#'
#'   # make distance a factor for plotting
#'   envs$distance <- as.factor(envs$distance)
#'
#'   
#'   # plot
#'   rlang::check_installed("ggplot2")
#'   library(ggplot2)
#'   
#'   ggplot(envs, aes(x= time, y = amp, col = distance)) +
#'   geom_line() + facet_wrap(~ sound.id) +
#'   scale_color_viridis_d() +
#'   labs(x = "Time (s)", y = "Amplitude (PMF)") +
#'   theme_classic()
#'   
#' }
#'
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})
#' @references {
#' Dabelsteen, T., Larsen, O. N., & Pedersen, S. B. (1993). Habitat-induced degradation of sound signals: Quantifying the effects of communication sounds and bird location on blur ratio, excess attenuation, and signal-to-noise ratio in blackbird song. The Journal of the Acoustical Society of America, 93(4), 2206.
#'
#' Araya-Salas M., E. Grabarczyk, M. Quiroz-Oliva, A. Garcia-Rodriguez, A. Rico-Guevara. (2023), baRulho: an R package to quantify degradation in animal acoustic signals .bioRxiv 2023.11.22.568305.
#' }

blur_ratio <-
  function(X,
           cores = getOption("mc.cores", 1),
           pb = getOption("pb", TRUE),
           env.smooth = getOption("env.smooth", 200),
           envelopes = FALSE,
           hop.size = getOption("hop.size", 11.6),
           wl = getOption("wl", NULL),
           ovlp = getOption("ovlp", 70),
           n.samples = 100,
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
    
    # total number of steps depending on whether envelopes are returned
    steps <- if (envelopes)
      3
    else
      2
    
    # get sampling rate assuming is the same for all sound files
    sampling.rate <- read_sound_file(X,
                                     index = 1,
                                     header = TRUE,
                                     path = path)$sample.rate
    
    # adjust wl based on hop.size
    wl <- .adjust_wl(wl, X, hop.size, path)
    
    
    # set clusters for windows OS
    if (Sys.info()[1] == "Windows" & cores > 1) {
      cl <- parallel::makePSOCKcluster(getOption("cl.cores", cores))
    } else {
      cl <- cores
    }
    
    # add sound file selec colums to X (weird column name so it does not overwrite user columns)
    X$.sgnl.temp <- paste(X$sound.files, X$selec, sep = "-")
    
    # get names of envelopes involved (those as test with reference or as reference)
    target_sgnl_temp <-
      unique(c(X$.sgnl.temp[!is.na(X$reference)], X$reference[!is.na(X$reference)]))
    
    # print message
    if (pb)
      write(
        file = "",
        x = paste0("Computing amplitude envelopes (step 1 out of ", steps, "):")
      )
    
    # calculate all envelops apply function
    envs <-
      warbleR:::pblapply_wrblr_int(
        pbar = pb,
        X = target_sgnl_temp,
        cl = cl,
        FUN = function(x,
                       ssmth = env.smooth,
                       ovl = ovlp,
                       Q = X,
                       wln = wl,
                       pth = path,
                       n.samp = n.samples) {
          .env(
            X = Q,
            y = x,
            env.smooth = ssmth,
            ovlp = ovl,
            wl = wln,
            path = pth,
            n.samples = n.samp
          )
        }
      )
    
    # add sound file selec column as names to envelopes
    names(envs) <- target_sgnl_temp
    
    # set options to run loops
    if (pb)
      write(file = "",
            x = paste0("Computing blur ratio (step 2 out of ", steps, "):"))
    
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
          .blur(
            x,
            X = Q,
            envs = nvs,
            ovlp = ovp,
            wl = wle,
            sampling.rate = sr
          )
        }
      ))
    
    
    # remove temporal column
    X$.sgnl.temp <- NULL
    
    # convert to list instead of extended selection table, add envelopes
    if (envelopes) {
      if (pb)
        write(file = "", x = "Saving envelopes (step 3 out of 3):")
      
      # get envelopes in a data frame
      env.dfs <-
        warbleR:::pblapply_wrblr_int(pbar = pb, seq_along(envs), cl = cl, function(y) {
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
    if (warbleR::is_extended_selection_table(X) | is_selection_table(X)) {
      attributes(X)$call <- base::match.call()
    }
    return(X)
  }
