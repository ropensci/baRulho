#' Create synthetic sounds
#'
#' \code{synth_sounds} create synthetic sounds
#' @inheritParams template_params
#' @param replicates Numeric vector of length 1 indicating the number of replicates for each treatment combination. Default is 1. Useful for measuring variation in transmission parameters.
#' @param frequencies Numeric vector with the different frequencies (in seconds) to synthesize. A Brownian bridge motion stochastic process (\code{diff.fun == "BB"}) is used to simulate frequency modulation (see \code{\link[warbleR]{simulate_songs}}).
#' @param durations Numeric vector with the different durations (in seconds) to synthesize.
#' @param nharmonics Numeric vector of length 1 specifying the number of harmonics to simulate. 1 indicates that only the fundamental
#' frequency harmonic will be simulated.
#' @param fm Logical to control if both frequency modulated sounds and pure tones (i.e. non-modulated sounds) are synthesize. If \code{FALSE} (default) only pure tones are synthesized.
#' @param am Logical to control if both amplitude modulated sounds and non-modulated sounds are synthesize. If \code{FALSE} (default) only non-modulated sounds are synthesized.
#' @param am.amps Numeric vector with the relative amplitude for each time step to simulate amplitude modulation (only applied to the fundamental frequency). The default value (\code{rep(c(1:4, 3:2), length.out = 11)}) has 2 amplitude peaks (although only applied if 'am = TRUE')
#' @param mar Numeric vector with the duration of margins of silence around sounds in seconds. Default is \code{0.05}.
#' @param seed Numeric vector of length 1. This allows users to get the same results in different runs (using \code{\link[base:Random]{se.seed}} internally). Default is \code{NULL}.
#' @param sampling.rate Numeric vector of length 1. Sets the sampling frequency of the wave object (in kHz). Default is 44.1.
#' @param sig2 Numeric vector of length 1 defining the sigma value of the brownian motion model (used for simulating frequency modulation). Default is 0.3.
#' @param shuffle Logical to control if the position of sounds is randomized. Having all sounds from the same treatment in a sequence can be problematic if an environmental noise masks them. Hence 'shuffle' is useful to avoid having sounds from the same treatment next to each other. Default is \code{FALSE}.
#' @param hrm.freqs Numeric vector with the frequencies of the harmonics relative to the fundamental frequency. The default values are c(1/2, 1/3, 2/3, 1/4, 3/4, 1/5, 1/6, 1/7, 1/8, 1/9, 1/10).
#' @param sampling.rate Numeric vector of length 1. Sets the sampling frequency of the wave object (in kHz). Default is 44.1.
#' @return An extended selection table, which can be input into \code{\link{master_sound_file}} to create the .wav file. The table contains columns for each of the varying features a 'treatment' column (useful to tell the acoustic features of each sound) and a 'replicate' column indicating the replicates for each 'treatment'.
#' @family prepare acoustic data
#' @seealso \code{\link[warbleR]{simulate_songs}} from the package warbleR.
#' @export
#' @name synth_sounds
#' @details This function creates synthetic sounds that can be used for playback experiments to understand the link between signal structure and its transmission properties. The function can add variation in signal structure in 5 features:
#' \itemize{
#'    \item \code{frequency}: continuous, argument 'frequencies'.
#'    \item \code{duration}: continuous, argument 'durations'.
#'    \item \code{harmonic structure}: binary (harmonics vs no-harmonics), arguments 'nharmonics' and 'hrm.freqs'.
#'    \item \code{frequency modulation}: variation in fundamental frequency across time. Binary (modulated vs non-modulated), arguments 'fm' and 'sig2'.
#'    \item \code{amplitude modulation}: variation in amplitude across time. Binary (modulated vs non-modulated), arguments 'am' and 'am.amps'.
#' }
#' Sound for all possible combinations of the selected structure dimensions will be synthesized. The output is an extended selection table, which can be input into \code{\link{master_sound_file}} to create the .wav file. The functions uses \code{\link[warbleR]{simulate_songs}} internally for synthesizing individual sounds. A Brownian bridge motion stochastic process (\code{diff.fun == "BB"}) is used to simulate frequency modulation. The output table contains columns for each of the varying features and a 'treatment' column (useful to tell sound from the same combination of features when using replicates).

#' @examples
#' \dontrun{
#'
#' synthetic_est <- synth_sounds(
#'   mar = 0.01,
#'   frequencies = c(1, 2, 3, 5),
#'   durations = 0.1,
#'   fm = TRUE,
#'   am = TRUE,
#'   nharmonics = 4,
#'   shuffle = TRUE,
#'   replicates = 3
#' )
#' }
#'
#' @references {
#' Araya-Salas, M., & Smith-Vidaurre, G. (2017). warbleR: An R package to streamline analysis of animal acoustic signals. Methods in Ecology and Evolution, 8(2), 184-191.
#' }
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})


synth_sounds <-
  function(replicates = 1,
           frequencies,
           durations,
           nharmonics = 1,
           fm = FALSE,
           am = FALSE,
           am.amps = rep(c(1:4, 3:2), length.out = 11),
           mar = 0.05,
           seed = NULL,
           sig2 = 0.3,
           shuffle = FALSE,
           hrm.freqs = c(1 / 2, 1 / 3, 2 / 3, 1 / 4, 3 / 4, 1 / 5, 1 / 6, 1 / 7, 1 /
                           8, 1 / 9, 1 / 10),
           sampling.rate = 44.1,
           pb = TRUE) {
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
    
    
    if (am & length(am.amps) < 2) {
      .stop("if 'am = TRUE' 'am.amps' must have more than 1 value 'length(am.amps) > 1'")
    }
    
    if (any(frequencies >= sampling.rate / 2)) {
      .stop("some frequencies are higher than the nyquist frequency (sampling.rate / 2)")
    }
    
    # set number of steps (default is 10)
    steps <- if (length(am.amps) > 2) {
      length(am.amps)
    } else {
      10
    }
    
    # make all possible combinations
    eg <- expand.grid(
      duration = durations,
      fm = c("no.fm", "fm"),
      am = c("no.am", "am"),
      harmonics = c("no.harmonics", "harmonics"),
      stringsAsFactors = FALSE
    )
    
    # remove levels when treatments not included
    if (!fm) {
      eg <- eg[eg$fm == "no.fm", ]
    }
    
    if (!am) {
      eg <- eg[eg$am == "no.am", ]
    }
    
    if (nharmonics == 1) {
      eg <- eg[eg$harmonics == "no.harmonics", ]
    }
    
    # create temporary directory
    temp_dir <- tempfile()
    dir.create(temp_dir)
    
    # simulate songs
    sim.songs <-
      warbleR:::pblapply_wrblr_int(
        seq_len(nrow(eg)),
        pbar = pb,
        FUN = .sim_song,
        temp_dir = temp_dir,
        eg = eg,
        frequencies = frequencies,
        steps = steps,
        am.amps = am.amps,
        nharmonics = nharmonics,
        mar = mar,
        sig2 = sig2,
        seed = seed,
        hrm.freqs = hrm.freqs,
        sampling.rate = sampling.rate
      )
    
    # name with parameters
    names(sim.songs) <-
      vapply(seq_len(nrow(eg)), function(x) {
        paste(eg[x, ], collapse = "_")
      }, FUN.VALUE = character(1))
    
    # extract select tables
    sim.song.sts <- lapply(sim.songs, function(X) {
      X$selec.table
    })
    
    sim.song.st <- do.call(.bind, sim.song.sts)
    
    if (shuffle) {
      sim.song.st <- sim.song.st[sample(seq_len(nrow(sim.song.st))), ]
    }
    
    if (any(sim.song.st$top.freq > sampling.rate / 2)){
    sim.song.st$top.freq[sim.song.st$top.freq >= sampling.rate / 2] <- (sampling.rate / 2)
    
    .warning("Some sounds had a top frequency higher than the nyquist frequency (sampling.rate / 2) and their frequency range was truncated")
    }
    
    
    # make a single extended selection table for simulation
    sim_sounds_est <-
      selection_table(
        mar = mar,
        X = sim.song.st,
        extended = TRUE,
        pb = FALSE,
        confirm.extended = FALSE,
        path = temp_dir,
        verbose = FALSE
      )
    
    # clean column names
    sim_sounds_est$frequency <- sim_sounds_est$sim.freq
    sim_sounds_est$sim.freq <- NULL
    
    # add treatments
    sim_sounds_est$duration <-
      sim_sounds_est$end - sim_sounds_est$start
    
    if (fm) {
      sim_sounds_est$frequency.modulation <-
        ifelse(grepl("no.fm", sim_sounds_est$sound.files), "tonal", "fm")
    }
    
    if (am) {
      sim_sounds_est$amplitude.modulation <-
        ifelse(grepl("no.am", sim_sounds_est$sound.files), "flat", "am")
    }
    
    if (nharmonics > 1) {
      sim_sounds_est$harmonics <-
        ifelse(grepl("no.harm", sim_sounds_est$sound.files),
               "pure.tone",
               "harmonics")
    }
    
    # replicate treatments
    if (replicates > 1) {
      sim_sounds_est <-
        .rep_synth_sound(
          y = seq_len(replicates),
          sim_sounds_est = sim_sounds_est,
          replicates = replicates
        )
    }
    
  # label columns output table
    sim_sounds_est <- .label_synth_est(sim_sounds_est, durations, frequencies)
    
    # fix call attribute
    attributes(sim_sounds_est)$call <- base::match.call()
    
    return(sim_sounds_est)
}