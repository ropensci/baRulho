#' Create synthetic sounds
#'
#' \code{synth_sounds} create synthetic sounds
#' @usage synth_sounds(replicates = 1, frequencies, durations, nharmonics = 1,
#' fm = FALSE, am = FALSE, am.amps = rep(c(1:4, 3:2), length.out = 11),
#' mar = 0.05, seed = NULL, sig2 = 0.3, shuffle = FALSE,
#' hrm.freqs = c(1/2, 1/3, 2/3, 1/4, 3/4, 1/5, 1/6, 1/7, 1/8, 1/9, 1/10),
#' sampling.rate = 44.1)
#' @param replicates Numeric vector of length 1 indicating the number of replicates for each treatment combination. Default is 1. Useful for measuring variation in transmission parameters.
#' @param frequencies Numeric vector with the different frequencies (in seconds) to synthesize. A Brownian bridge motion stochastic process (\code{diff.fun == "BB"}) is used to simulate frequency modulation (see \code{\link[warbleR]{simulate_songs}}).
#' @param durations Numeric vector with the different durations (in seconds) to synthesize.
#' @param nharmonics Numeric vector of length 1 specifying the number of harmonics to simulate. 1 indicates that only the fundamental
#' frequency harmonic will be simulated.
#' @param fm Logical to control if both frequency modulated sounds and pure tones (i.e. non-modulated sounds) are synthesize. If \code{FALSE} (default) only pure tones are synthesized.
#' @param am Logical to control if both amplitude modulated sounds and non-modulated sounds are synthesize. If \code{FALSE} (default) only non-modulated sounds are synthesized.
#' @param am.amps Numeric vector with the relative amplitude for each time step to simulate amplitude modulation (only applied to the fundamental frequency). The default value (\code{rep(c(1:4, 3:2), length.out = 11)}) has 2 amplitude peaks (although only applied if 'am = TRUE')
#' @param mar Numeric vector with the duration of margins of silence around sounds in seconds. Default is \code{0.05}.
#' @param seed Numeric vector of length 1. This allows users to get the same results in different runs (using \code{\link[base:Random]{set.seed}} internally). Default is \code{NULL}.
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
           sampling.rate = 44.1) {
    if (am & length(am.amps) < 2) {
      stop2("if 'am = TRUE' 'am.amps' must have more than 1 value 'length(am.amps) > 1'")
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
      eg <- eg[eg$fm == "no.fm",]
    }
    
    if (!am) {
      eg <- eg[eg$am == "no.am",]
    }
    
    if (nharmonics == 1) {
      eg <- eg[eg$harmonics == "no.harmonics",]
    }
    
    # create temporary directory
    temp_dir <- tempfile()
    dir.create(temp_dir)
    
    # simulate songs
    sim.songs <-
      warbleR:::pblapply_wrblr_int(seq_len(nrow(eg)), function(x) {
        sm.sng <- warbleR::simulate_songs(
          n = length(frequencies),
          durs = eg$dur[x],
          freqs = frequencies,
          samp.rate = sampling.rate,
          gaps = mar * 3 / 2,
          am.amps = if (eg$am[x] == "no.am") {
            1
          } else {
            am.amps
          },
          harms = if (eg$harm[x] == "no.harm") {
            1
          } else {
            nharmonics
          },
          harm.amps = if (eg$harm[x] == "no.harm") {
            1
          } else {
            nharmonics:1
          },
          diff.fun = if (eg$fm[x] == "fm") {
            "BB"
          } else {
            "pure.tone"
          },
          selec.table = TRUE,
          sig2 = sig2,
          steps = steps,
          file.name = paste(eg[x,], collapse = "_"),
          bgn = 0,
          seed = seed,
          path = temp_dir,
          hrm.freqs = hrm.freqs
        )
        
        # add freq room if pure tone
        if (eg$fm[x] == "no.fm") {
          sm.sng$selec.table$bottom.freq <-
            sm.sng$selec.table$bottom.freq - 0.2
          sm.sng$selec.table$top.freq <-
            sm.sng$selec.table$top.freq + 0.2
        }
        
        sm.sng$selec.table$bottom.freq[sm.sng$selec.table$bottom.freq < 0] <-
          0.1
        
        sm.sng$selec.table$sim.freq <- as.character(frequencies)
        
        return(sm.sng)
      })
    
    
    # name with parameters
    names(sim.songs) <-
      vapply(seq_len(nrow(eg)), function(x) {
        paste(eg[x,], collapse = "_")
      }, FUN.VALUE = character(1))
    
    # extract select tables
    sim.song.sts <- lapply(sim.songs, function(X) {
      X$selec.table
    })
    
    sim.song.st <- do.call(rbind2, sim.song.sts)
    
    if (shuffle) {
      sim.song.st <- sim.song.st[sample(seq_len(nrow(sim.song.st))),]
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
      rep_est_list <- lapply(seq_len(replicates), function(x) {
        Y <- sim_sounds_est
        Y$selec <- x
        attr(Y, "check.results")$selec <- x
        return(Y)
      })
      
      sim_sounds_est <- rep_est_list[[1]]
      
      for (i in 2:replicates) {
        suppressWarnings(sim_sounds_est <-
                           rbind(sim_sounds_est, rep_est_list[[i]]))
      }
    }
    
    # rename sound files
    sim_sounds_est <-
      warbleR::rename_est_waves(X = sim_sounds_est,
                                new.sound.files = paste0("synthetic_sound_", seq_along(unique(
                                  sim_sounds_est$sound.files
                                ))))
    
    sim_sounds_est$old.sound.file.name <- NULL
    
    # add single treatment column
    dur_label <- if (length(durations) > 1) {
      paste0("dur:", sim_sounds_est$duration)
    } else {
      NULL
    }
    freq_label <- if (length(frequencies) > 1) {
      paste0("freq:", sim_sounds_est$frequency)
    } else {
      NULL
    }
    freq_dur_label <- paste(dur_label, freq_label, sep = ";")
    
    sim_sounds_est$treatment <- if (ncol(sim_sounds_est) > 8) {
      sim_sounds_est$treatment <-
        paste(freq_dur_label,
              apply(sim_sounds_est[, 9:ncol(sim_sounds_est)], 1, paste, collapse = ";"),
              sep = ";")
    } else {
      freq_dur_label
    }
    
    # add treatment column
    sim_sounds_est$treatment <-
      gsub("^;", "", sim_sounds_est$treatment)
    
    # add sound id column (a unique identifier for each sound)
    sim_sounds_est$replicate <- 1
    
    for (i in 2:nrow(sim_sounds_est)) {
      sim_sounds_est$replicate[i] <-
        sum(sim_sounds_est$treatment[1:i] == sim_sounds_est$treatment[i])
    }
    
    sim_sounds_est$sound.id <-
      paste(sim_sounds_est$treatment, sim_sounds_est$replicate, sep = "_")
    
    # reset row names
    rownames(sim_sounds_est) <- seq_len(nrow(sim_sounds_est))
    
    # fix call attribute
    attributes(sim_sounds_est)$call <- base::match.call()
    
    return(sim_sounds_est)
  }
