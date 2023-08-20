#' Save multipanel plots with reference and test sounds
#'
#' \code{plot_degradation} creates multipanel plots (as image files) with reference and test sounds by distance and transect.
#' @usage plot_degradation(X, nrow = 4, ssmooth = getOption("ssmooth", 200),
#' hop.size = getOption("hop.size", 11.6), wl = getOption("wl", NULL),
#' ovlp = getOption("ovlp", 70),  path = getOption("sound.files.path", "."),
#' dest.path = getOption("dest.path", "."), cores = getOption("mc.cores", 1),
#' pb = getOption("pb", TRUE), collevels = seq(-120, 0, 5),
#' palette = viridis::viridis, flim = c("-1", "+1"), envelope = TRUE, spectrum = TRUE,
#' heights = c(4, 1), widths = c(5, 1), margins = c(2, 1), row.height = 2, col.width = 2,
#' colors = viridis::mako(4, alpha = 0.3), res = 120, ...)
#' @param X Object of class 'data.frame', 'selection_table' or 'extended_selection_table' (the last 2 classes are created by the function \code{\link[warbleR]{selection_table}} from the warbleR package) with the reference to the sounds in the master sound file. Must contain the following columns: 1) "sound.files": name of the .wav files, 2) "selec": unique selection identifier (within a sound file), 3) "start": start time and 4) "end": end time of selections, 5)  "bottom.freq": low frequency for bandpass, 6) "top.freq": high frequency for bandpass, 7) "sound.id": ID of sounds used to identify counterparts across distances and 8) "distance": distance at which each test sound was re-recorded. Each sound must have a unique ID within a distance. An additional 'transect' column labeling those sounds recorded in the same transect is required if including recordings from more than 1 transect.
#' @param nrow Numeric vector of length 1 with the number of rows per image file. Default is 4. This would be dynamically adjusted if more rows than needed are set.
#' @param ssmooth Numeric vector of length 1 determining the length of the sliding window (in amplitude samples) used for a sum smooth for amplitude envelope and power spectrum calculations (used internally by \code{\link[seewave]{env}}). Default is 200.
#' @param hop.size A numeric vector of length 1 specifying the time window duration (in ms). Default is 11.6 ms, which is equivalent to 512 wl for a 44.1 kHz sampling rate. Ignored if 'wl' is supplied.
#' @param wl A numeric vector of length 1 specifying the window length of the spectrogram, default
#' is NULL. If supplied, 'hop.size' is ignored.
#' @param ovlp Numeric vector of length 1 specifying the percent overlap between two
#'   consecutive windows, as in \code{\link[seewave]{spectro}}. Only used when plotting. Default is 70. Applied to both spectra and spectrograms on image files.
#' @param path Character string containing the directory path where the sound files are found. Only needed when 'X' is not an extended selection table.
#' @param dest.path Character string containing the directory path where the image files will be saved. If NULL (default) then the folder containing the sound files will be used instead.
#' @param cores Numeric vector of length 1. Controls whether parallel computing is applied by specifying the number of cores to be used. Default is 1 (i.e. no parallel computing).
#' @param pb Logical argument to control if progress bar is shown. Default is \code{TRUE}.
#' @param collevels	Numeric vector indicating a set of levels which are used to partition the amplitude range of the spectrogram (in dB) as in \code{\link[seewave]{spectro}}. Default is \code{seq(-120, 0, 5)}.
#' @param palette A color palette function to be used to assign colors in the
#'   plot, as in \code{\link[seewave]{spectro}}. Default is \code{\link[viridis]{viridis}}.
#' @param flim A numeric vector of length 2 indicating the highest and lowest frequency limits (kHz) of the spectrogram, as in \code{\link[seewave]{spectro}}. Default is \code{NULL}. Alternatively, a character vector similar to \code{c("-1", "1")} in which the first number is the value to be added to the minimum bottom frequency in 'X' and the second the value to be added to the maximum top frequency in 'X'. This is computed independently for each sound id so the frequency limit better fits the frequency range of the annotated signals. This is useful when test sounds show marked differences in their frequency ranges.
#' @param envelope Logical to control if envelopes are plotted. Default is \code{TRUE}.
#' @param spectrum Logical to control if power spectra are plotted. Default is \code{TRUE}.
#' @param heights Numeric vector of length 2 to control the relative heights of spectrogram (first number) and amplitude envelope (second number) when \code{envelope = TRUE}. Default is c(4, 1).
#' @param widths Numeric vector of length 2 to control the relative widths of spectrogram (first number) and power spectrum (second number) when \code{spectrum = TRUE}. Default is c(5, 1).
#' @param margins Numeric vector of length 2 to control the relative time of the test sound (first number) and adjacent margins (i.e. adjacent background noise, second number) to be included in the spectrogram \code{spectrum = TRUE}. Default is c(2, 1) which means that each margin next to the sound is half the duration of the sound. Note that all spectrograms will have the same time length so margins will be calculated to ensure all spectrograms match the duration of the spectrogram in the longest sound. As such, this argument controls the margin on the longest sound.
#' @param row.height Numeric vector of length 1 controlling the height (in inches) of sound panels in the output image file. Default is 2.
#' @param col.width Numeric vector of length 1 controlling the width (in inches) of sound panels in the output image file. Default is 2.
#' @param colors Character vector of length 4 containing the colors to be used for the background of column and row title panels (element 1), the color of amplitude envelopes (element 2), the color of power spectra (element 3), and the background color of envelopes and spectra (element 4).
#' @param res Numeric argument of length 1. Controls image resolution. Default is 120 (faster) although 300 - 400 is recommended for publication/presentation quality.
#' @param ... Additional arguments to be passed to the internal spectrogram
#' creating function for customizing graphical output. The function is a modified
#' version of \code{\link[seewave]{spectro}}, so it takes the same arguments.
#' @return One ore more image files with a multipanel figure of spectrograms of test sound by distance, sound id and transect.
#' @export
#' @name plot_degradation
#' @details The function aim to simplify the visual inspection of sound degradation by producing multipanel figures (saved in 'dest.path') containing visualizations of each test sound and its reference. Sounds are sorted by distance and transect (if more than 1). Visualizations include spectrograms, amplitude envelopes and power spectra (the last 2 are optional). Each row includes all the copies of a sound id for a given transect, also including its reference if it comes from another transect.
#' @examples
#' {
#'   # load example data
#'   data("degradation_est")
#'
#'   # create subset of data with only re-recorded files
#'   rerecorded_est <- degradation_est[degradation_est$sound.files != "master.wav", ]
#'
#'   # order so spectrograms from same sound id as close in the graph
#'   rerecorded_est <- rerecorded_est[order(rerecorded_est$sound.id), ]
#'
#'   # set directory to save image files
#'   options(dest.path = tempdir())
#'
#'   # plot degradation spectrograms
#'   plot_degradation(
#'     X = rerecorded_est, nrow = 3, ovlp = 95,
#'     colors = viridis::magma(4, alpha = 0.3),
#'     palette = viridis::magma
#'   )
#'
#'   # using other color palettes
#'   plot_degradation(
#'     X = rerecorded_est, nrow = 3, ovlp = 95,
#'     colors = viridis::magma(4, alpha = 0.3),
#'     palette = viridis::magma
#'   )
#'
#'   # missing some data
#'   plot_degradation(
#'     X = rerecorded_est[-3, ], nrow = 2, ovlp = 95,
#'     colors = viridis::mako(4, alpha = 0.4), palette = viridis::mako, wl = 200
#'   )
#'
#'   # changing marging
#'   plot_degradation(X = rerecorded_est, margins = c(5, 1), nrow = 6, ovlp = 95)
#'
#'   # more rows than needed
#'   plot_degradation(X = rerecorded_est, nrow = 10, ovlp = 90)
#' }
#'
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})
#' @seealso \code{\link{blur_ratio}}, \code{\link{realign_test_sounds}}, \code{\link{plot_align_sounds}}
#' @references {
#' Araya-Salas, M. (2020). baRulho: baRulho: quantifying habitat-induced degradation of (animal) acoustic signals in R. R package version 1.0.2
#' }
#'
plot_degradation <-
  function(X,
           nrow = 4,
           ssmooth = getOption("ssmooth", 200),
           hop.size = getOption("hop.size", 11.6),
           wl = getOption("wl", NULL),
           ovlp = getOption("ovlp", 70),
           path = getOption("sound.files.path", "."),
           dest.path = getOption("dest.path", "."),
           cores = getOption("mc.cores", 1),
           pb = getOption("pb", TRUE),
           collevels = seq(-120, 0, 5),
           palette = viridis::viridis,
           flim = c("-1", "+1"),
           envelope = TRUE,
           spectrum = TRUE,
           heights = c(4, 1),
           widths = c(5, 1),
           margins = c(2, 1),
           row.height = 2,
           col.width = 2,
           colors = viridis::mako(4, alpha = 0.3),
           res = 120,
           ...) {
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

    # set colors
    bg_titles <- colors[1]
    env_fill <- colors[2]
    spc_fill <- colors[3]
    bg_sp_env <- colors[4]

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
      cl <- parallel::makePSOCKcluster(getOption("cl.cores", cores))
    } else {
      cl <- cores
    }

    # remove ambient sounds
    X <- X[X$sound.id != "ambient", ]

    # stop if more than 1 sample per distance but no transect info
    if (any(table(X$sound.id, X$distance)) > 1 &
      is.null(X$transect)) {
      stop2(
        "There are more than 1 test sound per sound.id/distance combination but no 'transect' column to group by transect"
      )
    }

    # stop if there are more than 1 sample for distance sound id and transect combination (only 1 samp)
    if (!is.null(X$transect)) {
      if (any(table(X$sound.id, X$distance, X$transect)) > 1) {
        stop2("There is more than 1 test sound per sound.id/distance/transect combination")
      }
    }

    if (is.null(X$transect)) {
      transects <- X$transect <- 1
    } else {
      transects <- unique(X$transect)
    }

    # add sound file selec column and names to X (weird column name so it does not overwrite user columns)
    if (pb) {
      write(
        file = "",
        x = paste0("Preparing data for analysis (step 1 out of 2):")
      )
    }

    X <- prep_X_bRlo_int(X,
      method = 1,
      cores = cores,
      pb = pb
    )

    X_df <- as.data.frame(X)

    # calculate time and freq ranges based on all recs
    rangs <- lapply(1:nrow(X), function(i) {
      r <- read_sound_file(
        X = X,
        path = path,
        index = i,
        header = TRUE
      )
      f <- r$sample.rate

      # change mar to mar (if provided)
      adj.mar <- (X$end[i] - X$start[i]) * (margins[2] / margins[1])

      t <- c(X$start[i] - adj.mar, X$end[i] + adj.mar)

      if (t[1] < 0) {
        t[1] <- 0
      }

      if (t[2] > r$samples / f) {
        t[2] <- r$samples / f
      }

      return(data.frame(mardur = t[2] - t[1]))
    })

    rangs <- do.call(rbind, rangs)

    # for each sound fix start and end based on margin
    X$mar.end <- X$mar.start <- NA
    for (u in seq_len(nrow(X))) {
      dur <- X$end[u] - X$start[u]
      if (dur < max(rangs$mardur)) {
        X$mar.end[u] <- (max(rangs$mardur) - dur) / 2
        X$mar.start[u] <- (max(rangs$mardur) - dur) / 2
        if (X$start[u] - X$mar.start[u] < 0) {
          X$mar.end[u] <-
            X$mar.end[u] + abs(X$start[u] - X$mar.start[u])
          X$mar.start[u] <- X$start[u]
        }
      }
    }

    # create data subsets (element in list) with all the copies of a sound id in transect, also including its reference if it comes from another transect
    X_df$sound.id.transect <-
      paste(X_df$sound.id, X$transect, sep = "\n")

    soundid_X_list <-
      lapply(unique(X_df$sound.id.transect), function(x) {
        Y <-
          X_df[X_df$.sgnl.temp %in% unique(X_df$reference[X_df$sound.id.transect == x]) |
            X_df$sound.id.transect == x, ]
        Y <- Y[order(Y$distance), ]
        return(Y)
      })

    # add numeric id to each subset
    soundid_X_list <-
      lapply(seq_along(soundid_X_list), function(x) {
        soundid_X_list[[x]]$seq_number <- x
        return(soundid_X_list[[x]])
      })

    # get unique distances and maximum number of columns (distances) for any sound id
    distances <- sort(unique(X$distance))
    ncol <- length(distances)

    # put all into a single data frame
    soundid_X <- do.call(rbind, soundid_X_list)

    # sort by sound id as in X
    # Create a factor with desired order
    factor_order <- factor(soundid_X$sound.id, levels = unique(X$sound.id))

    # Sort
    soundid_X <- soundid_X[order(factor_order), ]

    # add page in which will be printed each subset
    soundid_X$page <-
      as.numeric(cut(soundid_X$seq_number, breaks = nrow * 0:max(soundid_X$seq_number)))

    # add seq number to know which will be plotted together
    soundid_X$sound.id.seq <-
      paste(soundid_X$sound.id, soundid_X$seq_number, sep = "-")

    # adjust nrow if less sound.id.seq than nrow
    if (nrow > length(unique(soundid_X$sound.id.seq))) {
      nrow <- length(unique(soundid_X$sound.id.seq))
    }

    # set image size
    img_width <- ncol * col.width
    img_heigth <- nrow * row.height

    # set basic layout for all pages
    # spectrogram panels
    # vertical lines
    # 0.9 is total space for side panels (freq and sound id labels)
    # 0.5 is space for freq panel and 0.4 for sound id pannel
    verticals <- seq(
      from = (1 / (ncol + 0.9)) * 0.5,
      to = 1 - ((1 / (ncol + 0.9)) * 0.4),
      length.out = ncol + 1
    )

    # left margin
    spec_lf <- rep(verticals[-length(verticals)], times = nrow)

    # right margin
    spec_rg <- rep(verticals[-1], times = nrow)

    # horizontal lines
    # 0.7 is total space for top and bottom panels (time and distance labels)
    # 0.4 is space for time panel and 0.3 for distance pannel

    horizontals <- seq(
      from = 1 - ((1 / (nrow + 0.7)) * 0.3),
      to = (1 / (nrow + 0.7)) * 0.4,
      length.out = nrow + 1
    )

    # bottom margin
    spec_btm <- rep(horizontals[-1], each = ncol)

    # top margin
    spec_tp <- rep(horizontals[-length(horizontals)], each = ncol)

    # put together in a matrix
    spec_m <- cbind(spec_lf, spec_rg, spec_btm, spec_tp)

    # add additional panels for spectra and/or envelopes
    if (spectrum | envelope) {
      spec_m_list <- lapply(seq_len(nrow(spec_m)), function(y) {
        psp <- envel <- spectr <- spec_m[y, ]

        if (envelope) {
          psp[3] <-
            envel[4] <-
            spectr[3] <-
            spectr[3] + ((spectr[4] - spectr[3]) / (heights[1] / heights[2]))
        }

        if (spectrum) {
          psp[2] <-
            envel[1] <-
            spectr[1] <-
            spectr[1] + ((spectr[2] - spectr[1]) / (widths[1] / widths[2]))
        }

        m <-
          matrix(
            c(spectr, if (spectrum) {
              psp
            }, if (envelope) {
              envel
            }),
            ncol = 4,
            byrow = TRUE
          )

        return(m)
      })

      spec_m <- do.call(rbind, spec_m_list)
    }

    ## label panels
    # left margin
    lab_lf <-
      c(rep(max(verticals), each = nrow), verticals[-length(verticals)], 0, 0)

    # right margin
    lab_rg <-
      c(rep(1, each = nrow), verticals[-1], min(verticals), 1)

    # bottom margin
    lab_btm <-
      c(horizontals[-1], rep(max(horizontals), ncol), min(horizontals), 0)

    # top margin
    lab_tp <-
      c(
        horizontals[-length(horizontals)],
        rep(1, ncol),
        max(horizontals),
        min(horizontals)
      )

    # put together in a matrix
    lab_m <- cbind(lab_lf, lab_rg, lab_btm, lab_tp)

    # single data frame with all panels
    page_layout <- rbind(spec_m, lab_m)

    # testing layout screens
    # ss <- split.screen(figs = page_layout)
    # for(i in 1:nrow(page_layout))
    # {screen(i)
    #   par( mar = rep(0, 4))
    #   plot(0.5, xlim = c(0,1), ylim = c(0,1), type = "n", axes = FALSE, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
    #   box()
    #   text(x = 0.5, y = 0.5, labels = i)
    # }
    #
    # close.screen(all.screens = T)

    # reset graphic device on exit
    on.exit(grDevices::dev.off())

    if (pb) {
      write(file = "", x = "Saving plots (step 2 out of 2):")
    }

    out <-
      warbleR:::pblapply_wrblr_int(unique(soundid_X$page), function(x, flm = flim) {
        # extract data subset for a page
        Y <- soundid_X[soundid_X$page == x, ]

        # start graphic device
        warbleR:::img_wrlbr_int(
          filename = paste0("plot_degradation_p", x, ".jpeg"),
          path = dest.path,
          width = img_width,
          height = img_heigth,
          units = "in",
          res = res,
        )

        # set panel layout
        suppressWarnings(catch <- split.screen(figs = page_layout))

        # frequency label
        par(mar = c(0, 0, 0, 0), new = TRUE)

        # activate screen for frequency axis label
        screen((nrow * ncol * (sum(
          c(envelope, spectrum)
        ) + 1)) + nrow + ncol + 1)

        par(mar = c(0, 0, 0, 0), new = TRUE)

        plot(
          1,
          frame.plot = FALSE,
          type = "n",
          yaxt = "n",
          xaxt = "n"
        )

        text(
          x = 0.8,
          y = 1,
          "Frequency (kHz)",
          srt = 90,
          cex = 1.2
        )

        # activate screen for time label
        par(mar = c(0, 0, 0, 0), new = TRUE)

        screen((nrow * ncol * (sum(
          c(envelope, spectrum)
        ) + 1)) + nrow + ncol + 2)

        par(mar = c(0, 0, 0, 0), new = TRUE)

        plot(
          1,
          frame.plot = FALSE,
          type = "n",
          yaxt = "n",
          xaxt = "n"
        )

        text(
          x = 1,
          y = 0.75,
          "Time (s)",
          cex = 1.2
        )


        # get combination of distances and sound id to loop over it
        grd <-
          expand.grid(
            distance = distances,
            sound.id.seq = unique(Y$sound.id.seq)
          )


        # add screen number for spectrogram
        grd$screen <-
          seq(
            from = 1,
            to = nrow * ncol * (sum(c(
              envelope, spectrum
            )) + 1),
            by = sum(spectrum, envelope) + 1
          )[seq_len(nrow(grd))]

        grd$.sgnl.temp <- sapply(seq_len(nrow(grd)), function(o) {
          sgnl <- Y$.sgnl.temp[Y$distance == grd$distance[o] &
            Y$sound.id.seq == grd$sound.id.seq[o]]
          if (length(sgnl) == 0) sgnl <- NA
          return(sgnl)
        })

        # start with empty signal id vector so it can be recorded with prev_sgnl (this is needed for adding frequency labels when the spectrogram at the left side is missing)
        sgnl <- NA

        # loop to create spectrograms
        for (i in grd$screen) {
          prev_sgnl <- sgnl

          # get id of signal
          sgnl <- grd$.sgnl.temp[grd$screen == i]

          # if signal exists
          if (!is.na(sgnl)) {
            # get row index for signals in X
            indx <- which(X$.sgnl.temp == sgnl)

            # read wave
            wave <-
              read_wave(
                X = X,
                index = indx,
                from = X$start[indx] - X$mar.start[indx],
                to = X$end[indx] + X$mar.end[indx],
                path = path
              )

            # set frequency limits
            if (is.character(flm)) {
              fl <-
                c(
                  min(Y$bottom.freq[Y$.sgnl.temp == sgnl]) + as.numeric(flm[1]),
                  max(Y$top.freq[Y$.sgnl.temp == sgnl]) + as.numeric(flm[2])
                )
            } else {
              fl <- flm
            }
            # fix if lower than 0
            if (fl[1] < 0) {
              fl[1] <- 0
            }

            par(mar = c(0, 0, 0, 0), new = TRUE)

            # start screen
            screen(i)
            par(mar = c(0, 0, 0, 0), new = TRUE)
            curr_dist <- grd$distance[grd$screen == i]

            # plot spectrogram
            warbleR:::spectro_wrblr_int2(
              wave = wave,
              palette = palette,
              axisX = FALSE,
              axisY = FALSE,
              grid = FALSE,
              collevels = collevels,
              flim = fl,
              wl = wl,
              ...
            )

            # add vertical lines
            # add dotted lines
            abline(
              v = c(Y$mar.start[Y$.sgnl.temp == sgnl], X$end[indx] - X$start[indx] + Y$mar.start[Y$.sgnl.temp == sgnl]),
              col = "white",
              lty = 3,
              lwd = 1.5
            )

            if (spectrum) {
              # set screen
              i <- i + 1
              screen(i)
              par(
                mar = c(0, 0, 0, 0),
                new = TRUE
              )

              # get power spectrum
              spc <- spec(wave = wave, plot = FALSE)

              # smooth
              spc[, 2] <-
                warbleR::envelope(x = spc[, 2], ssmooth = ssmooth)

              # filter to flim
              spc <- spc[spc[, 1] > fl[1] & spc[, 1] < fl[2], ]

              # set white plot
              plot(
                x = spc[, 2],
                y = spc[, 1],
                type = "l",
                frame.plot = FALSE,
                yaxt = "n",
                xaxt = "n",
                col = spc_fill,
                xaxs = "i",
                yaxs = "i"
              )

              # add background color
              rect(min(spc[, 2]), min(spc[, 1]), max(spc[, 2]), max(spc[, 1]), col = "white", border = NA)
              rect(min(spc[, 2]), min(spc[, 1]), max(spc[, 2]), max(spc[, 1]), col = bg_sp_env, border = NA)

              # add polygon with spectrum shape
              polygon(rbind(c(0, 0), spc[, 2:1]), col = spc_fill, border = NA)

              box()

              par(
                mar = c(0, 0, 0, 0),
                bg =  "#FFFFFF00",
                new = TRUE
              )
            }

            # plot frequency ticks
            if (page_layout[i, 1] == min(page_layout[1:(nrow * ncol * (sum(c(
              envelope, spectrum
            )) + 1)), 1]) |
              is.na(prev_sgnl) & curr_dist > distances[1]) {
              at_freq <-
                pretty(seq(0, wave@samp.rate / 2000, length.out = 10)[-10], n = 10)
              axis(2,
                at = at_freq,
                labels = if (envelope) {
                  at_freq
                } else {
                  c(at_freq[-length(at_freq)], "")
                }
              )
            }

            # plot amplitude envelope
            if (envelope) {
              # set screen
              i <- i + 1
              screen(i)

              # get power spectrum
              envlp <-
                warbleR::envelope(
                  x = seewave::ffilter(
                    wave = wave,
                    from = fl[1] * 1000,
                    to = fl[2] * 1000
                  ),
                  ssmooth = ssmooth
                )

              par(
                mar = c(0, 0, 0, 0),
                new = TRUE
              )


              # set white plot
              plot(
                x = seq(0, duration(wave), along.with = envlp),
                y = envlp,
                type = "l",
                frame.plot = FALSE,
                yaxt = "n",
                xaxt = "n",
                col = spc_fill,
                xaxs = "i",
                yaxs = "i"
              )

              # add background color
              rect(0, min(envlp), duration(wave), max(envlp), col = "white", border = NA)
              rect(0, min(envlp), duration(wave), max(envlp), col = bg_sp_env, border = NA)

              # add polygon with envelope shape
              polygon(
                rbind(c(0, 0), cbind(
                  seq(0, duration(wave), along.with = envlp), envlp
                )),
                col = spc_fill,
                border = NA
              )
              box()
            } else if (spectrum) {
              par(mar = c(0, 0, 0, 0), bg = "#FFFFFF00")
              screen(i - 1)
            }

            # plot time ticks
            at_time <-
              pretty(seq(0, duration(wave), length.out = 10)[-10], n = 4)
            axis(1,
              at = at_time,
              labels = if (spectrum & !anyNA(grd$.sgnl.temp)) {
                at_time
              } else {
                c("", at_time[c(-1, -length(at_time))], "")
              }
            )
          }
        }

        # plot distance labels
        dist_labs <- paste(distances, "m")

        # add labels
        for (e in seq_len(ncol) + ((nrow * (ncol) * (sum(
          c(envelope, spectrum)
        ) + 1)) + nrow)) {
          # activate screen
          screen(e)
          par(
            mar = c(0, 0, 0, 0),
            new = TRUE
          )

          plot(
            1,
            frame.plot = FALSE,
            type = "n",
            yaxt = "n",
            xaxt = "n"
          )

          # add background color
          rect(0, 0, 2, 2, col = "white", border = NA)
          rect(0, 0, 2, 2, col = bg_titles, border = NA)

          text(
            x = 1,
            y = 1,
            dist_labs[e - ((nrow * (ncol) * (sum(
              c(envelope, spectrum)
            ) + 1)) + nrow)],
            cex = 1.2,
            col = "white",
            font = 2
          )
          box()
        }
        # get sound id transect labels
        soundid_transect_labs <-
          unique(Y$sound.id.transect[!Y$.sgnl.temp %in% Y$reference])

        # plot sound id transect labels
        for (e in (seq_len(nrow) + (nrow * ncol * (sum(
          c(envelope, spectrum)
        ) + 1)))[seq_along(soundid_transect_labs)]) {
          # activate screen
          screen(e)
          par(
            mar = c(0, 0, 0, 0),
            new = TRUE
          )

          plot(
            1,
            frame.plot = FALSE,
            type = "n",
            yaxt = "n",
            xaxt = "n"
          )

          # add background color
          rect(0, 0, 2, 2, col = "white", border = NA)
          rect(0, 0, 2, 2, col = bg_titles, border = NA)

          text(
            x = 1,
            y = 1,
            soundid_transect_labs[e - (nrow * ncol * (sum(c(
              envelope, spectrum
            )) + 1))],
            font = 2,
            srt = 270,
            cex = 1.2,
            col = "white"
          )
          box()
        }



        close.screen(all.screens = T)
      })
  }
