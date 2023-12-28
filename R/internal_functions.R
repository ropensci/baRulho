### internal functions not to be called by users ###

# stop function that doesn't print call
.stop <- function(...) {
  stop(..., call. = FALSE)
}

# calculate time and freq ranges based on all recs
.time_freq_range_files <- function(X, i, path, margins) {
  rang_list <- lapply(seq_len(nrow(X)), function(i) {
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
  
  rangs <- do.call(rbind, rang_list)
  
  return(rangs)
}

# It is a modified version of warbleR::find_peaks
# that allows to define internally if progress bar would be used (pbapply::pblapply uses pboptions to do this)
# Find cross-correlation peaks

.find_peaks <-
  function(xc.output,
           cores = getOption("mc.cores", 1),
           cutoff = 0.4,
           pb = getOption("pb", TRUE),
           max.peak = FALSE,
           output = "data.frame") {
    # set clusters for windows OS and no soz
    if (Sys.info()[1] == "Windows" & cores > 1) {
      cl <- parallel::makePSOCKcluster(getOption("cl.cores", cores))
    } else {
      cl <- cores
    }
    
    # loop over scores of each dyad
    pks <-
      warbleR:::pblapply_wrblr_int(
        pbar = pb,
        X = unique(xc.output$scores$dyad),
        cl = cl,
        FUN = function(i) {
          # extract data for a dyad
          dat <- xc.output$scores[xc.output$scores$dyad == i, ]
          
          # check xc.output being a autodetec.output object
          if (!(is(xc.output, "xcorr.output") |
                is(xc.output, "xc.output"))) {
            .stop("'xc.output' must be and object of class 'xcorr.output'")
          }
          
          ## get peaks as the ones higher than previous and following scores
          pks <-
            dat[c(FALSE, diff(dat$score) > 0) &
                  c(rev(diff(rev(dat$score)) > 0), FALSE) &
                  dat$score > cutoff, , drop = FALSE]
          
          # get the single highest peak
          if (max.peak) {
            pks <- dat[which.max(dat$score), , drop = FALSE]
          }
          
          return(pks)
        }
      )
    
    # put results in a data frame
    peaks <- do.call(rbind, pks)
    
    # relabel rows
    if (nrow(peaks) > 0) {
      rownames(peaks) <- seq_len(nrow(peaks))
      
      # remove dyad column
      peaks$dyad <- NULL
      
      #### name as in a warbleR selection table
      # remove selec info at the end
      peaks$sound.files <-
        substr(peaks$sound.files,
               start = 0,
               regexpr("\\-[^\\-]*$", peaks$sound.files) - 1)
      
      #### add start and end
      # add template column to selection table in xc.output
      Y <- xc.output$org.selection.table
      Y$template <- paste(Y$sound.files, Y$selec, sep = "-")
      
      # Y <- Y[Y$template %in% comp_mat[, 1], ]
      
      # add start as time - half duration of template
      peaks$start <- vapply(seq_len(nrow(peaks)), function(i) {
        peaks$time[i] -
          ((Y$end[Y$template == peaks$template[i]] -
              Y$start[Y$template == peaks$template[i]]) / 2)
      }, FUN.VALUE = numeric(1))
      
      # add end as time + half duration of template
      peaks$end <- vapply(seq_len(nrow(peaks)), function(i) {
        peaks$time[i] +
          ((Y$end[Y$template == peaks$template[i]] -
              Y$start[Y$template == peaks$template[i]]) / 2)
      }, FUN.VALUE = numeric(1))
      
      # add selec labels
      peaks$selec <- 1
      
      if (nrow(peaks) > 1) {
        for (i in 2:nrow(peaks)) {
          if (peaks$sound.files[i] == peaks$sound.files[i - 1]) {
            peaks$selec[i] <- peaks$selec[i - 1] + 1
          }
        }
      }
      
      # sort columns in a intuitive order
      peaks <- warbleR::sort_colms(peaks)
      
      # output results
      if (output == "data.frame") {
        return(peaks)
      } else {
        output_list <- list(
          selection.table = peaks,
          scores = xc.output$scores,
          cutoff = cutoff,
          call = base::match.call(),
          spectrogram = xc.output$spectrogram
          # warbleR.version = packageVersion("warbleR")
        )
        
        class(output_list) <- c("list", "find_peaks.output")
        
        return(output_list)
      }
    } else {
      # no detections
      write(file = "", x = "no peaks above cutoff were detected")
      
      return(NULL)
    }
  }

.onAttach <-
  function(libname, pkgname) {
    packageStartupMessage("\nPlease cite 'baRulho' as: \n")
    packageStartupMessage(
      "Araya-Salas, M. (2020), baRulho: quantifying degradation of (animal) acoustic signals in R. R package version 1.0.0"
    )
    
    invisible(TRUE)
  }

# set options when loading package
.onLoad <- function(libname, pkgname) {
  # set options
  options("baRulho_check_args" = TRUE)
  invisible(NULL)
}

# remove options when unloading
.onUnload <- function(libpath) {
  options(baRulho_check_args = NULL)
  invisible(NULL)
}

# warning function that doesn't print call
.warning <- function(x, color = "magenta") {
  warning(.colortext(x, as = color), call. = FALSE)
}

# message function that changes colors
.message <- function(x, color = "black") {
  message(.colortext(x, as = color))
}

# coloring text
.colortext <-
  function(text,
           as = c("red",
                  "blue",
                  "green",
                  "magenta",
                  "cyan",
                  "orange",
                  "black",
                  "silver")) {
    if (.has_color()) {
      unclass(cli::make_ansi_style(.baRulho_style(as))(text))
    } else {
      text
    }
  }

.has_color <- function() {
  cli::num_ansi_colors() > 1
}

.baRulho_style <-
  function(color = c("red",
                     "blue",
                     "green",
                     "magenta",
                     "cyan",
                     "orange",
                     "black",
                     "silver")) {
    type <- match.arg(color)
    
    c(
      red = "red",
      blue = "blue",
      green = "green",
      magenta = "magenta",
      cyan = "cyan",
      orange = "orange",
      black = "black",
      silver = "silver"
    )[[color]]
  }

# internal function to get detection distance from spl and cutoff used in detection_distance()
.detection_distance <-
  function(spl.cutoff,
           spl,
           frequency,
           distance,
           temp = 20,
           rh = 60,
           pa = 101325,
           hab.att.coef = 0.02,
           max.distance = 1000,
           resolution = 0.1) {
    # initial SPL and distance
    L <- spl
    iter_dist <- 0
    
    # loop until SPL is equal or lower than background noise
    while (L - spl.cutoff > 0) {
      iter_dist <- iter_dist + resolution
      att <-
        attenuation(
          frequency = frequency,
          temp = temp,
          dist = iter_dist,
          dist0 = distance,
          rh = rh,
          pa = pa,
          hab.att.coef = hab.att.coef
        )
      L <- spl - att$combined.attenuation
      
      if (iter_dist >= max.distance) {
        iter_dist <- NA
        break
      }
    }
    
    return(iter_dist)
  }
## function to measure blur ratio

.spctr <-
  function(y,
           spec.smooth,
           wl,
           X,
           path,
           meanspc = FALSE,
           ovlp,
           n.bins) {
    # load clip
    clp <- warbleR::read_sound_file(X = X,
                                    index = which(X$.sgnl.temp == y),
                                    path = path)
    
    # calculate spectrum
    clp.spc <- if (meanspc) {
      # mean spec
      meanspec(
        wave = clp,
        f = clp@samp.rate,
        plot = FALSE,
        wl = wl,
        ovlp = ovlp
      )
    } else {
      seewave::spec(
        wave = clp,
        f = clp@samp.rate,
        plot = FALSE,
        wl = wl
      )
    }
    
    # smoothing
    clp.spc[, 2] <-
      warbleR::envelope(x = clp.spc[, 2],
                        ssmooth = spec.smooth)
    
    # thin
    if (!is.null(n.bins)) {
      # reduce size of envelope
      
      # return NA if doesn't have at least 2 non-NA values (need at least two non-NA values to interpolate)
      if (sum(!is.na(clp.spc[, 2])) >= 2) {
        clp.spc_list <-
          stats::approx(
            x = clp.spc[, 1],
            y = clp.spc[, 2],
            n = n.bins,
            method = "linear"
          )
        
        # make it a matrix
        clp.spc <- cbind(clp.spc_list[[1]], clp.spc_list[[2]])
      } else {
        clp.spc <-  cbind(NA, NA)
      }
    }
    
    return(clp.spc)
  }

## function to measure blur ratio
.blur <-
  function(X,
           envs,
           x,
           ovlp,
           wl,
           sampling.rate) {
    # get names of sound and reference
    sgnl <- X$.sgnl.temp[x]
    rfrnc <- X$reference[x]
    
    # if reference is NA return NA
    if (is.na(rfrnc)) {
      bl.rt <- NA
    } else {
      # extract envelope for sound and model
      sgnl.env <- envs[[which(names(envs) == sgnl)]]
      rfrnc.env <- envs[[which(names(envs) == rfrnc)]]
      
      # make them the same length as the shortest one
      if (length(sgnl.env) > length(rfrnc.env)) {
        sgnl.env <- sgnl.env[seq_along(rfrnc.env)]
      }
      if (length(rfrnc.env) > length(sgnl.env)) {
        rfrnc.env <- rfrnc.env[seq_along(sgnl.env)]
      }
      
      # duration (any sampling rate works as they all must have the same sampling rate)
      dur <-
        length(sgnl.env) / sampling.rate
      
      # convert envelopes to PMF (probability mass function)
      rfrnc.pmf <- rfrnc.env / sum(rfrnc.env)
      sgn.pmf <- sgnl.env / sum(sgnl.env)
      
      # get blur ratio as half the sum of absolute differences between envelope PMFs
      bl.rt <- sum(abs(rfrnc.pmf - sgn.pmf)) / 2
    }
    return(bl.rt)
  }


## function to measure spectrum blur ratio
.blur_sp <-
  function(x, X, ovlp, wl, specs, sampling_rate) {
    # get names of sound and reference
    sgnl <- X$.sgnl.temp[x]
    rfrnc <- X$reference[x]
    
    # if reference is NA return NA
    if (is.na(rfrnc)) {
      sp.bl.rt <- NA
    } else {
      # extract spectrum for sound and model
      sgnl.spc <- specs[[which(names(specs) == sgnl)]]
      rfrnc.spc <- specs[[which(names(specs) == rfrnc)]]
      
      # compute blur ratio only if power distribution have data
      if (!is.na(sgnl.spc[1, 1]) & !is.na(rfrnc.spc[1, 1])) {
        # make them the same number of rows
        sgnl.spc <-
          sgnl.spc[1:(min(c(nrow(sgnl.spc), nrow(rfrnc.spc)))), ]
        rfrnc.spc <-
          rfrnc.spc[1:(min(c(nrow(sgnl.spc), nrow(rfrnc.spc)))), ]
        
        # make test the same frequency range as reference
        bp <-
          c(X$bottom.freq[X$.sgnl.temp == rfrnc], X$top.freq[X$.sgnl.temp == rfrnc])
        
        bp <- bp + c(-0.2, 0.2) # add 0.2 kHz buffer
        if (bp[1] < 0) {
          # force 0 if negative
          bp[1] <- 0
        }
        if (bp[2] > ceiling(sampling_rate / 2000) - 1) {
          bp[2] <-
            ceiling(sampling_rate / 2000) - 1
        } # force lower than nyquist freq if higher
        
        # apply bandpass by shrinking freq range and remove freq column based on reference freq bins
        sgnl.spc <-
          sgnl.spc[rfrnc.spc[, 1] > bp[1] &
                     rfrnc.spc[, 1] < bp[2], 2]
        rfrnc.spc <-
          rfrnc.spc[rfrnc.spc[, 1] > bp[1] &
                      rfrnc.spc[, 1] < bp[2], 2]
        
        
        # convert envelopes to PMF (probability mass function)
        rfrnc.pmf <- rfrnc.spc / sum(rfrnc.spc)
        sgnl.pmf <- sgnl.spc / sum(sgnl.spc)
        
        # get blur ratio as half the sum of absolute differences between spectra PMFs
        sp.bl.rt <- sum(abs(rfrnc.pmf - sgnl.pmf)) / 2
      } else {
        sp.bl.rt <- NA
      }
      
    }
    return(sp.bl.rt)
  }

## function to plot  blur ratios
.plot_blur <-
  function(X,
           energy_vectors,
           spectr,
           path,
           dest.path,
           x,
           res,
           ovlp,
           wl,
           collevels,
           palette,
           bp,
           flim,
           colors) {
    # set colors
    ref_col <- colors[1]
    test_col <- colors[2]
    blur_col <- colors[3]
    
    
    # get names of sound and reference
    sgnl <- X$.sgnl.temp[x]
    rfrnc <- X$reference[x]
    
    # if reference is NA return NA
    if (!is.na(rfrnc)) {
      # extract envelope for sound and model
      sgnl.energy <-
        energy_vectors[[which(names(energy_vectors) == sgnl)]]
      rfrnc.energy <-
        energy_vectors[[which(names(energy_vectors) == rfrnc)]]
      
      # make them the same length as the shortest one
      if (length(sgnl.energy) > length(rfrnc.energy)) {
        sgnl.energy <- sgnl.energy[seq_along(rfrnc.energy)]
      }
      if (length(rfrnc.energy) > length(sgnl.energy)) {
        rfrnc.energy <- rfrnc.energy[seq_along(sgnl.energy)]
      }
      
      sampling_rate <- warbleR::read_sound_file(
        X = X,
        index = x,
        path = path,
        header = TRUE
      )$sample.rate
      
      dur <-
        length(sgnl.energy) / sampling_rate
      
      # run band pass
      if (spectr) {
        # make them the same frequency range as reference
        bp <-
          c(X$bottom.freq[X$.sgnl.temp == rfrnc], X$top.freq[X$.sgnl.temp == rfrnc])
        
        bp <- bp + c(-0.2, 0.2) # add 0.2 kHz buffer
        if (bp[1] < 0) {
          # force 0 if negative
          bp[1] <- 0
        }
        if (bp[2] > ceiling(sampling_rate / 2000) - 1) {
          bp[2] <-
            ceiling(sampling_rate / 2000) - 1
        } # force lower than nyquist freq if higher
        
        # apply bandpass by shrinking freq range and remove freq column based on reference freq bins
        sgnl.energy <-
          sgnl.energy[rfrnc.energy[, 1] > bp[1] &
                        rfrnc.energy[, 1] < bp[2], 2]
        rfrnc.energy <-
          rfrnc.energy[rfrnc.energy[, 1] > bp[1] &
                         rfrnc.energy[, 1] < bp[2], 2]
      }
      
      # convert envelopes to PMF (probability mass function)
      rfrnc.pmf <- rfrnc.energy / sum(rfrnc.energy)
      sgn.pmf <- sgnl.energy / sum(sgnl.energy)
      
      # get blur ratio as half the sum of absolute differences between envelope PMFs
      bl.rt <- sum(abs(rfrnc.pmf - sgn.pmf)) / 2
      
      img_name <- paste0(if (spectr) {
        "spectrum_blur_ratio_"
      } else {
        "blur_ratio_"
      },
      X$sound.id[x],
      "-",
      rfrnc,
      "-",
      sgnl,
      ".jpeg")
      
      # plot
      warbleR:::img_wrlbr_int(
        filename = img_name,
        path = dest.path,
        width = 10.16 * 1.5,
        height = 10.16,
        units = "cm",
        res = res
      )
      
      
      # matrix for layout
      page_layout <- matrix(
        c(
          0.06,
          0.4,
          0,
          0.562,
          # bottom left spectrogram
          0.06,
          0.4,
          0.562,
          1,
          # top left spectrogram
          0.4,
          1,
          0,
          1,
          # right pannel with blur ratio
          0,
          0.06,
          0.1,
          1
        ),
        nrow = 4,
        byrow = TRUE
      )
      
      # testing layout screens
      # ss <- split.screen(figs = page_layout)
      # for(i in seq_len(nrow(page_layout)))
      # {screen(i)
      #   par( mar = rep(0, 4))
      #   plot(0.5, xlim = c(0,1), ylim = c(0,1), type = "n", axes = FALSE, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
      #   box()
      #   text(x = 0.5, y = 0.5, labels = i)
      # }
      #
      # save par settings
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))
      
      # close if open any screen
      invisible(close.screen(all.screens = TRUE))
      
      # split screen
      split.screen(page_layout)
      
      ## plot spectros
      
      # index of reference
      rf.indx <-
        which(paste(X$sound.files, X$selec, sep = "-") == rfrnc)
      
      # freq limit of reference
      flm <- c(X$bottom.freq[rf.indx], X$top.freq[rf.indx])
      
      # set frequency limits
      if (is.character(flim)[1]) {
        flm <-
          c(flm[1] + as.numeric(flim[1]),
            flm[2] + as.numeric(flim[2]))
      } else {
        flm <- flim
      }
      # fix if lower than 0
      if (flm[1] < 0) {
        flm[1] <- 0
      }
      
      #####
      # end for sound and reference
      rf.info <-
        warbleR::read_sound_file(
          X = X,
          index = rf.indx,
          header = TRUE,
          path = path
        )
      rf.dur <- rf.info$samples / rf.info$sample.rate
      
      # fix upper frequency in flim
      if (flm[2] > rf.info$sample.rate / 2000) {
        flm[2] <- rf.info$sample.rate / 2000
      }
      
      sgnl.info <-
        warbleR::read_sound_file(
          X = X,
          index = x,
          header = TRUE,
          path = path
        )
      sgnl.dur <- sgnl.info$samples / sgnl.info$sample.rate
      
      # calculate margin for spectrogram, before and after
      mar.rf.af <-
        mar.rf.bf <- (X$end[rf.indx] - X$start[rf.indx]) / 4
      
      # start for sound and reference
      strt.sgnl <- X$start[x] - mar.rf.bf
      if (strt.sgnl < 0) {
        strt.sgnl <- 0
      }
      strt.rf <- X$start[rf.indx] - mar.rf.bf
      if (strt.rf < 0) {
        strt.rf <- 0
      }
      
      end.sgnl <- X$end[x] + mar.rf.af
      if (end.sgnl > sgnl.dur) {
        end.sgnl <- sgnl.dur
      }
      end.rf <- X$end[rf.indx] + mar.rf.af
      if (end.rf > rf.dur) {
        end.rf <- rf.dur
      }
      
      # extract clip reference and sound
      clp.sgnl <-
        warbleR::read_sound_file(
          X = X,
          index = x,
          from = strt.sgnl,
          to = end.sgnl,
          path = path
        )
      clp.rfnc <-
        warbleR::read_sound_file(
          X = X,
          index = rf.indx,
          from = strt.rf,
          to = end.rf,
          path = path
        )
      
      
      # frequency axis for spectrograms
      screen(4)
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
        y = 1,
        "Frequency (kHz)",
        srt = 90,
        cex = 1.2
      )
      
      
      # sound at bottom left
      screen(1)
      par(mar = c(3, 2, 0.15, 0.3))
      
      warbleR:::spectro_wrblr_int2(
        wave = clp.sgnl,
        f = clp.sgnl@samp.rate,
        flim = flm,
        axisX = FALSE,
        axisY = FALSE,
        tlab = NULL,
        flab = NULL,
        main = NULL,
        grid = FALSE,
        rm.zero = TRUE,
        cexaxis = 1,
        add = TRUE,
        ovlp = ovlp,
        wl = wl,
        collevels = collevels,
        palette = palette
      )
      
      at_freq <-
        pretty(seq(0, clp.sgnl@samp.rate / 2000, length.out = 10)[-10], n = 10)
      axis(2,
           at = at_freq,
           labels =
             c(at_freq[-length(at_freq)], ""))
      
      # plot time ticks
      at_time <-
        pretty(seq(0, duration(clp.sgnl), length.out = 10)[-10], n = 4)
      axis(1,
           at = at_time,
           labels = c(at_time[c(-length(at_time))], ""))
      
      # add x axis label
      mtext(text = "Time (s)",
            side = 1,
            line = 2)
      
      
      
      # lines showing position of sound
      abline(
        v = if (!spectr) {
          c(mar.rf.bf, X$end[x] - X$start[x] + mar.rf.bf)
        } else {
          NULL
        },
        h = if (spectr) {
          bp
        } else {
          NULL
        },
        col =
          "white",
        lty = 2
      )
      
      # add box with sound color
      box(col = adjustcolor(test_col, 0.6), lwd = 3)
      
      # reference at top left
      screen(2)
      par(mar = c(0, 2, 0.3, 0.3))
      
      warbleR:::spectro_wrblr_int2(
        wave = clp.rfnc,
        f = clp.rfnc@samp.rate,
        flim = flm,
        axisX = FALSE,
        axisY = FALSE,
        tlab = NULL,
        flab = NULL,
        main = NULL,
        grid = FALSE,
        rm.zero = TRUE,
        cexaxis = 1,
        add = TRUE,
        ovlp = ovlp,
        wl = wl,
        collevels = collevels,
        palette = palette
      )
      
      
      # add box with reference color
      box(col = adjustcolor(ref_col, 0.6), lwd = 3)
      
      # lines showing position of sound
      abline(
        v = if (!spectr) {
          c(mar.rf.bf, X$end[x] - X$start[x] + mar.rf.bf)
        } else {
          NULL
        },
        h = if (spectr) {
          bp
        } else {
          NULL
        },
        col =
          "white",
        lty = 2
      )
      
      at_freq <-
        pretty(seq(0, clp.rfnc@samp.rate / 2000, length.out = 10)[-10], n = 10)
      axis(2,
           at = at_freq,
           labels =
             c(at_freq[-length(at_freq)], ""))
      
      # plot envelopes
      screen(3)
      
      # set image margins
      par(mar = c(4, 1, 4, if (spectr) {
        4.2
      } else {
        3.2
      }))
      
      # plot envelope
      if (!spectr) {
        # time values for plots
        time.vals <- seq(0, dur, length.out = length(sgnl.energy))
        
        # reference envelope first
        plot(
          time.vals,
          rfrnc.pmf,
          type = "l",
          xlab = "",
          ylab = "",
          col = adjustcolor(ref_col, 0.7),
          ylim = c(min(rfrnc.pmf, sgn.pmf), max(rfrnc.pmf, sgn.pmf) * 1.1),
          cex.main = 0.8,
          lwd = 1.4,
          yaxt = "n",
          xaxs = "i",
          yaxs = "i"
        )
        
        # add background color
        rect(
          par("usr")[1],
          par("usr")[3],
          par("usr")[2],
          par("usr")[4],
          col = adjustcolor("#DEF5E5FF", 0.4),
          border = NA
        )
        
        # white envelope polygon
        # add 0s at star and end so polygon doesnt twist
        sgn.pmf[c(1, nrow(sgn.pmf))] <- 0
        
        # add polygon with envelope shape
        polygon(
          x = time.vals,
          y = sgn.pmf,
          col = "white",
          border = NA
        )
        
        # add x axis label
        mtext(text = "Time (s)",
              side = 1,
              line = 2.5)
        
        # add title
        mtext(
          text = paste("Sound ID:", X$sound.id[x]),
          side = 3,
          line = 3,
          cex = 1
        )
        mtext(
          text = paste("Reference:", rfrnc),
          side = 3,
          line = 1.75,
          col = ref_col,
          cex = 1
        )
        
        mtext(
          text = paste("Test sound:", sgnl),
          side = 3,
          line = 0.5,
          col = test_col,
          cex = 1
        )
        
        # add y axis
        axis(side = 4, labels = FALSE)
        mtext(text = "Amplitude (PMF)",
              side = 4,
              line = 1)
        
        # blur region
        polygon(
          x = c(time.vals, rev(time.vals)),
          y = c(sgn.pmf, rev(rfrnc.pmf)),
          col = adjustcolor(blur_col, 0.2),
          border = NA
        )
        
        # add sound envelope
        lines(time.vals,
              sgn.pmf,
              col = adjustcolor(test_col, 0.7),
              lwd = 1.4)
        # add sound envelope
        lines(time.vals,
              rfrnc.pmf,
              col = adjustcolor(ref_col, 0.7),
              lwd = 1.4)
        
        # get plotting area limits
        usr <- par("usr")
        
        # and blu ratio value
        text(
          x = ((usr[1] + usr[2]) / 2) + usr[1],
          y = usr[4] * 0.95,
          paste("Blur ratio:", round(bl.rt, 2)),
          cex = 1
        )
      }
      
      # spectrum
      if (spectr) {
        # create time values for area calculation
        f.vals <-
          seq(bp[1], bp[2], length.out = length(rfrnc.pmf))
        
        # reference spectrum first
        plot(
          x = rfrnc.pmf,
          y = f.vals,
          type = "l",
          xlab = "",
          ylab = "",
          col = ref_col,
          xlim = c(min(rfrnc.pmf, sgn.pmf),
                   max(rfrnc.pmf, sgn.pmf) * 1.1),
          cex.main = 0.8,
          lwd = 1.2,
          yaxt = "n",
          xaxs = "i",
          yaxs = "i"
        )
        
        # add background color
        rect(
          par("usr")[1],
          par("usr")[3],
          par("usr")[2],
          par("usr")[4],
          col = adjustcolor("#DEF5E5FF", 0.4),
          border = NA
        )
        
        # white envelope polygon
        # add 0s at star and end so polygon doesnt twist
        rfrnc.pmf[c(1, nrow(rfrnc.pmf))] <- 0
        
        # add polygon with spectrum shape
        polygon(
          x = rfrnc.pmf,
          y = f.vals,
          col = "white",
          border = NA
        )
        # add x axis label
        mtext(text = "Power spectrum (PMF)",
              side = 1,
              line = 2.5)
        
        # add title
        mtext(
          text = paste("Sound ID:", X$sound.id[x]),
          side = 3,
          line = 3,
          cex = 1
        )
        mtext(
          text = paste("Reference:", rfrnc),
          side = 3,
          line = 1.75,
          col = ref_col,
          cex = 1
        )
        mtext(
          text = paste("Test sound:", sgnl),
          side = 3,
          line = 0.5,
          col = test_col,
          cex = 1
        )
        
        # add y axis
        axis(side = 4)
        mtext(text = "Frequency (kHz)",
              side = 4,
              line = 2)
        
        # add sound spectrum
        lines(sgn.pmf,
              f.vals,
              col = test_col,
              lwd = 1.2)
        
        # sound spectrum on top
        polygon(
          y = c(f.vals, rev(f.vals)),
          x = c(sgn.pmf, rev(rfrnc.pmf)),
          col = adjustcolor(blur_col, 0.2),
          border = NA
        )
        
        # add sound envelope
        lines(sgn.pmf,
              f.vals,
              col = adjustcolor(test_col, 0.7),
              lwd = 1.4)
        # add sound envelope
        lines(rfrnc.pmf,
              f.vals,
              col = adjustcolor(ref_col, 0.7),
              lwd = 1.4)
        
        
        # get plotting area limits
        usr <- par("usr")
        
        # and blu ratio value
        text(
          x = ((usr[1] + usr[2]) / 2),
          y = (usr[4] - usr[3]) * 0.95 + usr[3],
          paste("Spectrum blur ratio:", round(bl.rt, 2)),
          cex = 1
        )
        
        # index of reference
        rf.indx <-
          which(paste(X$sound.files, X$selec, sep = "-") == rfrnc)
        
        # freq limit of reference
        # flim <- c(X$bottom.freq[rf.indx], X$top.freq[rf.indx])
        
        # end for sound and reference
        rf.info <-
          warbleR::read_sound_file(
            X = X,
            index = rf.indx,
            header = TRUE,
            path = path
          )
        rf.dur <- rf.info$samples / rf.info$sample.rate
        
        sgnl.info <-
          warbleR::read_sound_file(
            X = X,
            index = x,
            header = TRUE,
            path = path
          )
        sgnl.dur <- sgnl.info$samples / sgnl.info$sample.rate
        
        # calculate margin for spectrogram, before and after
        mar.rf.af <-
          mar.rf.bf <- (X$end[rf.indx] - X$start[rf.indx]) / 4
        
        # start for sound and reference
        strt.sgnl <- X$start[x] - mar.rf.bf
        if (strt.sgnl < 0) {
          strt.sgnl <- 0
        }
        strt.rf <- X$start[rf.indx] - mar.rf.bf
        if (strt.rf < 0) {
          strt.rf <- 0
        }
        
        
        end.sgnl <- X$end[x] + mar.rf.af
        if (end.sgnl > sgnl.dur) {
          end.sgnl <- sgnl.dur
        }
        end.rf <- X$end[rf.indx] + mar.rf.af
        if (end.rf > rf.dur) {
          end.rf <- rf.dur
        }
        
        # extract clip reference and sound
        clp.sgnl <-
          warbleR::read_sound_file(
            X = X,
            index = x,
            from = strt.sgnl,
            to = end.sgnl,
            path = path
          )
        clp.rfnc <-
          warbleR::read_sound_file(
            X = X,
            index = rf.indx,
            from = strt.rf,
            to = end.rf,
            path = path
          )
      }
      
      # close graph
      dev.off()
    }
  }

# function to extract envelopes from wave objects
.env <-
  function(X, y, env.smooth, ovlp, wl, path, n.samples = 100) {
    # load clip
    clp <- warbleR::read_sound_file(X = X,
                                    index = which(X$.sgnl.temp == y),
                                    path = path)
    
    # define bandpass
    bp <-
      c(X$bottom.freq[X$.sgnl.temp == y], X$top.freq[X$.sgnl.temp == y])
    
    # bandpass filter
    clp <- seewave::ffilter(
      clp,
      from = bp[1] * 1000,
      ovlp = ovlp,
      to = bp[2] * 1000,
      bandpass = TRUE,
      wl = wl,
      output = "Wave"
    )
    
    # calculate envelope
    nv <-
      warbleR::envelope(x = clp@left,
                        ssmooth = env.smooth)
    
    # thin
    if (!is.null(n.samples)) {
      # reduce size of envelope
      nv <-
        stats::approx(
          x = seq_along(nv),
          y = nv,
          n = n.samples,
          method = "linear"
        )$y
    }
    
    return(nv)
  }

# function to measure envelope correlation
# y and z are the sound.files+selec names of the sounds and reference sound (model)
.env_cor <- function(X, x, envs, cor.method) {
  # if names are the same return NA
  
  # get names of sound and reference
  sgnl <- X$.sgnl.temp[x]
  rfrnc <- X$reference[x]
  
  # if reference is NA return NA
  if (is.na(rfrnc)) {
    envcor <- NA
  } else {
    # extract envelope for sound and model
    sgnl.env <- envs[[which(names(envs) == sgnl)]]
    mdl.env <- envs[[which(names(envs) == rfrnc)]]
    
    # define short and long envelope for sliding one (short) over the other (long)
    if (length(mdl.env) > length(sgnl.env)) {
      lg.env <- mdl.env
      shrt.env <- sgnl.env
    } else {
      lg.env <- sgnl.env
      shrt.env <- mdl.env
    }
    
    # get length of shortest minus 1 (1 if same length so it runs a single correlation)
    shrt.lgth <- length(shrt.env) - 1
    
    # steps for sliding one sound over the other
    stps <- length(lg.env) - shrt.lgth
    
    # calculate correlations at each step
    cors <- vapply(seq_along(stps), function(x) {
      cor(lg.env[x:(x + shrt.lgth)], shrt.env, method = cor.method)
    }, FUN.VALUE = numeric(1))
    
    # return maximum correlation
    envcor <- max(cors, na.rm = TRUE)
  }
  return(envcor)
}

# function to extract mean envelopes
mean.env <- function(y, wl, ovlp, X, path, bp) {
  # read sound clip
  clp <-
    warbleR::read_sound_file(
      X = X,
      index = which(X$.sgnl.temp == y),
      from = X$start[X$.sgnl.temp == y],
      to = X$end[X$.sgnl.temp == y],
      path = path
    )
  
  # add band-pass frequency filter
  if (!is.null(bp)) {
    # filter to bottom and top freq range
    if (bp[1] == "freq.range") {
      bp <-
        c(X$bottom.freq[X$.sgnl.temp == y], X$top.freq[X$.sgnl.temp == y])
    }
    
    clp <-
      seewave::ffilter(
        clp,
        f = clp@samp.rate,
        from = bp[1] * 1000,
        ovlp = ovlp,
        to = bp[2] * 1000,
        bandpass = TRUE,
        wl = wl,
        output = "Wave"
      )
  }
  
  sig_env <-
    mean(warbleR::envelope(x = clp@left))
  
  return(sig_env)
}

# function to measure spectrum correlation
# y and z are the sound.files+selec names of the sounds and reference sound (model)
.spctr_cor <- function(y, specs, X, cor.method) {
  # get names of sound and reference
  sgnl <- X$.sgnl.temp[y]
  rfrnc <- X$reference[y]
  
  # if reference is NA return NA
  if (is.na(rfrnc)) {
    cor.spctr <- NA
  } else {
    # extract envelope for sound and model
    sgnl.spctr <- specs[[which(names(specs) == sgnl)]]
    mdl.spctr <- specs[[which(names(specs) == rfrnc)]]
    
    if (!is.na(sgnl.spctr[1, 1]) & !is.na(mdl.spctr[1, 1])) {
      ### filter to freq range of sounds and remove freq column
      # get range as lowest bottom and highest top
      frng <-
        c(min(X$bottom.freq[X$.sgnl.temp %in% c(sgnl, rfrnc)]), max(X$top.freq[X$.sgnl.temp %in% c(sgnl, rfrnc)]))
      sgnl.spctr <-
        sgnl.spctr[sgnl.spctr[, 1] > frng[1] &
                     sgnl.spctr[, 1] < frng[2], 2]
      mdl.spctr <-
        mdl.spctr[mdl.spctr[, 1] > frng[1] &
                    mdl.spctr[, 1] < frng[2], 2]
      
      # get correlation assuming they have same length
      cor.spctr <- cor(sgnl.spctr, mdl.spctr, method = cor.method)
    } else {
      cor.spctr <- NA
    }
  }
  
  return(cor.spctr)
}

.exc_att <- function(y, X, tp, gn) {
  # get names of sound and reference
  sgnl <- X$.sgnl.temp[y]
  rfrnc <- X$reference[y]
  
  # if reference is NA return NA
  if (is.na(rfrnc)) {
    ea <- NA
  } else {
    # extract mean envelope of sounds
    sig_env_REF <- X$sig_env[X$.sgnl.temp == rfrnc]
    dist_REF <- X$distance[X$.sgnl.temp == rfrnc]
    dist_SIG <- X$distance[y]
    
    ks <- X$sig_env[y] / sig_env_REF
    
    # type Dabelsteen
    if (tp == "Dabelsteen") {
      ea <-
        (-20 * log(ks)) - (6 / (2 * (dist_SIG - dist_REF))) + gn
    }
    
    if (tp == "Darden") {
      # EA = g - 20 log(d / 10) - 20 log(k)
      ea <-
        gn - 20 * log10(dist_SIG / 10) - 20 * log(ks)
    }
  }
  return(ea)
}

# function to put together simulated sounds (synth_sounds())
.bind <- function(...) {
  suppressWarnings(rbind(...))
}

## adjust SNR
.add_noise <-
  function(x,
           mar,
           target.snr,
           precision,
           max.iterations,
           Y) {
    # extract selection as single extended selection table
    Y_x <- Y[x,]
    
    # normalize wave object
    attributes(Y_x)$wave.objects[[1]] <-
      normalize(attributes(Y_x)$wave.objects[[1]], unit = "1")
    
    # estimate current snr
    snr <-
      signal_to_noise_ratio(Y_x, mar = mar, pb = FALSE)$signal.to.noise.ratio
    
    if (snr > target.snr) {
      # reset time coordinates of sounds if lower than 0 o higher than duration
      stn <- Y$start[x] - mar
      mar1 <- mar
      
      if (stn < 0) {
        mar1 <- mar1 + stn
        stn <- 0
      }
      
      # read sound and margin
      wav <-
        warbleR::read_sound_file(
          X = Y_x,
          index = 1,
          from = 0,
          to = Inf,
          path = NULL
        )
      
      # start point for adding noise (a 1/10 of signal amplitude)
      prop_noise <- 0.3
      
      prop_noise_vector <- vector()
      snr_vector <- vector()
      
      seed <- 0
      while (all(snr > target.snr + precision |
                 snr < target.snr - precision) &
             length(prop_noise_vector) < max.iterations) {
        seed <- seed + 1
        set.seed(seed)
        noise_wav <-
          runif(length(wav@left),
                min = -1 * prop_noise,
                max = prop_noise)
        
        attributes(Y_x)$wave.objects[[1]] <- wav + noise_wav
        
        snr <-
          signal_to_noise_ratio(X = Y_x,
                                mar = mar,
                                pb = FALSE)$signal.to.noise.ratio
        
        prop_noise_vector[length(prop_noise_vector) + 1] <-
          prop_noise
        snr_vector[length(snr_vector) + 1] <- snr
        
        # increase constant to modify noise level when output snr higher than target
        if (snr > target.snr + precision) {
          prop_noise <- prop_noise * 1.3
        }
        
        # decrease constant to modify noise level when output snr lower than target
        if (snr < target.snr - precision) {
          prop_noise <- prop_noise / 1.3
        }
      }
      
      # adjust SNR using best SNR
      set.seed(which.min(abs(snr_vector - target.snr)))
      prop_noise <-
        prop_noise_vector[which.min(abs(snr_vector - target.snr))]
      attributes(Y_x)$wave.objects[[1]] <- wav +
        runif(length(wav@left),
              min = -1 * prop_noise,
              max = prop_noise)
      snr <-
        snr_vector[which.min(abs(snr_vector - target.snr))]
      modified <- TRUE
    } else {
      modified <- FALSE
    }
    seed <- NA
    
    return(list(
      wave = attributes(Y_x)$wave.objects[[1]],
      snr = snr,
      modified = modified
    ))
  }


# measure RMS of sounds referenced in X (used by signal_to_noise_ratio())
.rms <-
  function(y,
           Y,
           mar,
           sampling_rate,
           wl,
           noise.ref,
           path,
           eq.dur,
           ovlp,
           bp = bp) {
    # only calculate for non-markers and for ambient only if custom noise.ref
    if (!Y$sound.id[y] %in% c("marker", if (noise.ref != "custom")
      "ambient")) {
      if (noise.ref == "custom") {
        # read sound clip
        signal <-
          warbleR::read_sound_file(X = Y,
                                   index = y,
                                   path = path)
        
        # add band-pass frequency filter
        if (!is.null(bp)) {
          # filter to bottom and top freq range
          if (bp[1] == "freq.range") {
            bp <- c(Y$bottom.freq[y], Y$top.freq[y])
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
        sig_rms <- seewave::rms(warbleR::envelope(signal@left))
        bg_rms <- NA
      }
      
      if (noise.ref == "adjacent") {
        # set margin to half of signal duration
        if (eq.dur) {
          mar <-
            (Y$end[y] - Y$start[y])
        }
        
        # Read sound files to get sample rate and length
        r <-
          warbleR::read_sound_file(
            X = Y,
            index = y,
            header = TRUE,
            path = path
          )
        
        # reset time coordinates of sounds if lower than 0 o higher than duration
        stn <- Y$start[y] - mar
        enn <- Y$end[y] + mar
        mar1 <- mar
        
        if (stn < 0) {
          mar1 <- mar1 + stn
          stn <- 0
        }
        
        mar2 <- mar1 + Y$end[y] - Y$start[y]
        
        if (enn > r$samples / sampling_rate) {
          enn <- r$samples / sampling_rate
        }
        
        # read sound and margin
        noise_sig <-
          warbleR::read_sound_file(
            X = Y,
            index = y,
            from = stn,
            to = enn,
            path = path
          )
        
        # add band-pass frequency filter
        if (!is.null(bp)) {
          # filter to bottom and top freq range
          if (bp[1] == "freq.range") {
            bp <- c(Y$bottom.freq[y], Y$top.freq[y])
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
                        f = sampling_rate)
        
        # get RMS for signal
        sig_rms <- seewave::rms(warbleR::envelope(signal[, 1]))
        
        # cut ambient noise before sound
        noise1 <-
          seewave::cutw(noise_sig,
                        from = 0,
                        to = mar1,
                        f = sampling_rate)
        
        # get RMS for background noise
        bg_rms <- seewave::rms(warbleR::envelope(noise1[, 1]))
      }
    } else {
      sig_rms <- NA
      bg_rms <- NA
    }
    
    return(list(sig_rms = sig_rms, bg_rms = bg_rms))
  }

# measure SNR of sounds referenced in X (used by signal_to_noise_ratio())
.snr <- function(y, W, rms_list, noise.ref, type) {
  if (W$sound.id[y] != "ambient") {
    suppressWarnings({
      # sound RMS
      sig_RMS <- rms_list[[W$.y[y]]]$sig_rms
      # get reference ambient noise RMS
      if (noise.ref == "adjacent") {
        bg_RMS <- rms_list[[W$.y[y]]]$bg_rms
      } else {
        # get envelopes from ambient selections
        bg_RMS <-
          lapply(rms_list[W$.y[W$sound.files == W$sound.files[y] &
                                        W$sound.id == "ambient"]], "[", "sig_rms")
        
        # get mean RMS from combined envelopes
        bg_RMS <-
          mean(unlist(bg_RMS))
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
}

# replicated extended selection table (used by synth_sounds())
.rep_synth_sound <- function(y, sim_sounds_est, replicates) {
  rep_est_list <- lapply(y, function(x) {
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
  return(sim_sounds_est)
  
}

# simulate songs (used by sim_sounds())
.sim_song <-
  function(x,
           temp_dir,
           eg,
           frequencies,
           steps,
           am.amps,
           nharmonics,
           mar,
           sig2,
           seed,
           hrm.freqs,
           sampling.rate) {
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
        "GBM"
      } else {
        "pure.tone"
      },
      selec.table = TRUE,
      sig2 = sig2,
      steps = steps,
      file.name = paste(eg[x, ], collapse = "_"),
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
  }

# add colums to synth extended selection table (used synth_sounds())
.label_synth_est <- function(X, durations, frequencies){
  # rename sound files
  X <-
    warbleR::rename_est_waves(X = X,
                              new.sound.files = paste0("synthetic_sound_", seq_along(unique(
                                X$sound.files
                              ))))
  
  X$old.sound.file.name <- NULL
  
  # add single treatment column
  dur_label <- if (length(durations) > 1) {
    paste0("dur:", X$duration)
  } else {
    NULL
  }
  freq_label <- if (length(frequencies) > 1) {
    paste0("freq:", X$frequency)
  } else {
    NULL
  }
  freq_dur_label <- paste(dur_label, freq_label, sep = ";")
  
  X$treatment <- if (ncol(X) > 8) {
    X$treatment <-
      paste(freq_dur_label,
            apply(X[, 9:ncol(X)], 1, paste, collapse = ";"),
            sep = ";")
  } else {
    freq_dur_label
  }
  
  # add treatment column
  X$treatment <-
    gsub("^;", "", X$treatment)
  
  # add sound id column (a unique identifier for each sound)
  X$replicate <- 1
  
  for (i in 2:nrow(X)) {
    X$replicate[i] <-
      sum(X$treatment[1:i] == X$treatment[i])
  }
  
  X$sound.id <-
    paste(X$treatment, X$replicate, sep = "_")
  
  # reset row names
  rownames(X) <- seq_len(nrow(X))
  
  return(X)
}

# get same number of frequency bins for noise_profile()
.same_length_noise <- function(noise.profiles, rws){  
  
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
  
  return(noise.profiles)
  
}

# adjust wl based on hop.size if wl null
.adjust_wl <- function(wl, X, hop.size, path = NULL){
  
  # adjust wl based on hop.size
  if (is.null(wl)) {
    wl <- if (is_extended_selection_table(X))
      round(attr(X, "check.results")$sample.rate[1] * hop.size, 0) else
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
  
  return(wl)
}


# internal function to set references (used by set_reference_sounds())
.set_ref <- function(x, meth, Z) {
  # extract for single sound id and order by distance
  Y <-
    Z[Z$sound.id == Z$sound.id[Z$.sgnl.temp == x], , drop = FALSE]
  
  if (!Y$sound.id[1] %in% c("ambient", "start_marker", "end_marker")) {
    # Order by distance
    Y <- Y[order(Y$distance), ]
    
    # method 1 compare to closest distance to source
    if (meth == 1) {
      # if column transect is found select the lowest distance in that trasnect
      if (!is.null(Z$transect)) {
        W <-
          Z[Z$transect == Z$transect[Z$.sgnl.temp == x] &
              Z$sound.id == Y$sound.id[Y$.sgnl.temp == x], , drop = FALSE]
        
        # if there is another distance that is shorter in other transects for that signal, use that distance
        z <- if (min(W$distance) <= min(Y$distance)) {
          W$.sgnl.temp[which.min(W$distance)]
        } else {
          Y$.sgnl.temp[which.min(Y$distance)]
        }
      } else {
        z <- Y$.sgnl.temp[which.min(Y$distance)]
      }
    } else {
      # if method 2
      # get those from the same transect and same sound id
      W <-
        Z[Z$transect == Z$transect[Z$.sgnl.temp == x] &
            Z$sound.id == Z$sound.id[Z$.sgnl.temp == x], , drop = FALSE]
      W <- W[order(W$distance), ]
      
      # if not the first row then the previous row
      if (W$.sgnl.temp[1] != x) {
        z <- W$.sgnl.temp[which(W$.sgnl.temp == x) - 1]
      } else {
        # else the first row
        z <- x
      }
    }
    
    # set reference to NA if is the same than the current row
    if (z == x) {
      z <- NA
    }
  } else {
    z <- NA
  }
  
  return(z)
}


# measure noise profile for a single sound file
.noise_profile <-
  function(y,
           Y,
           noise.ref,
           mar,
           path,
           wl,
           PSD,
           bp,
           dB,
           norm) {
    # extract  complete sound file for custom or files in folder
    if (noise.ref == "custom") {
      noise.wv <-
        warbleR::read_sound_file(
          X = Y,
          index = y,
          from = 0,
          to = Inf,
          path = path
        )
    }
    
    if (noise.ref == "adjacent") {
      # reset time coordinates of sounds if lower than 0 o higher than duration
      stn <- Y$start[y] - mar
      
      if (stn < 0) {
        stn <- 0
      }
      
      # read ambient noise
      noise.wv <-
        warbleR::read_sound_file(
          X = Y,
          index = y,
          from = stn,
          to = Y$start[y],
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
      data.frame(sound.files = Y$sound.files[y],
                 selec = Y$selec[y],
                 mspc)
    
    # add band-pass frequency filter
    if (!is.null(bp)) {
      mspc <- mspc[mspc$freq >= bp[1] & mspc$freq <= bp[2], ]
    }
    
    return(mspc)
  }

## function to measure detection distance (used by detection_distance())

.detection_dist <-
  function(x,
           spl.cutoff,
           temp,
           rh,
           pa,
           hab.att.coef,
           max.distance,
           resolution,
           spl,
           X,
           peak_freq_list,
           ...) {
    # get names of sound and reference
    sgnl <- X$.sgnl.temp[x]
    rfrnc <- X$reference[x]
    
    # if sounds are the same or the selection is noise return NA
    # if reference is NA return NA
    if (is.na(rfrnc)) {
      detect_dist <- NA
    } else {
      # extract spectrum for sound and model
      sgnl.spl <-
        peak_freq_list[[which(names(peak_freq_list) == sgnl)]]$spl
      rfrnc.pkf <-
        peak_freq_list[[which(names(peak_freq_list) == rfrnc)]]$peakf
      
      # get detection distance
      detect_dist <-
        .detection_distance(
          spl.cutoff = spl.cutoff,
          spl = sgnl.spl,
          frequency = rfrnc.pkf * 1000,
          distance = X$distance[x],
          temp = temp,
          rh = rh,
          pa = pa,
          hab.att.coef = hab.att.coef,
          max.distance = max.distance,
          resolution = resolution
        )
    }
    return(detect_dist)
  }

# plot funciton, used by plot_degradation
.plot_degrad <-  function(x, 
                          X,
                          soundid_X,
                          flim, 
                          path,
                          dest.path,
                          img_width,
                          img_heigth,
                          res,
                          page_layout,
                          nrow,
                          ncol,
                          envelope,
                          spectrum,
                          distances,
                          wl,
                          spc_fill,
                          bg_sp_env,
                          bg_titles,
                          ovlp,
                          collevels, 
                          env.smooth,
                          palette,
                          ...) {
  
  # extract data subset for a page
  Y <- soundid_X[soundid_X$page == x,]
  
  # start graphic device
  warbleR:::img_wrlbr_int(
    filename = paste0("plot_degradation_p", x, ".jpeg"),
    path = dest.path,
    width = img_width,
    height = img_heigth,
    units = "in",
    res = res
  )
  
  # set panel layout
  invisible(close.screen(all.screens = TRUE))
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
  
  text(x = 1,
       y = 0.75,
       "Time (s)",
       cex = 1.2)
  
  
  # get combination of distances and sound id to loop over it
  grd <-
    expand.grid(distance = distances,
                sound.id.seq = unique(Y$sound.id.seq))
  
  
  # add screen number for spectrogram
  grd$screen <-
    seq(
      from = 1,
      to = nrow * ncol * (sum(c(
        envelope, spectrum
      )) + 1),
      by = sum(spectrum, envelope) + 1
    )[seq_len(nrow(grd))]
  
  grd$.sgnl.temp <- vapply(seq_len(nrow(grd)), function(o) {
    sgnl <- Y$.sgnl.temp[Y$distance == grd$distance[o] &
                           Y$sound.id.seq == grd$sound.id.seq[o]]
    if (length(sgnl) == 0) {
      sgnl <- NA
      class(sgnl) <- "character"
    }
    return(sgnl)
  }, FUN.VALUE = character(1))
  
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
      if (is.character(flim)) {
        fl <-
          c(
            min(Y$bottom.freq[Y$.sgnl.temp == sgnl]) + as.numeric(flim[1]),
            max(Y$top.freq[Y$.sgnl.temp == sgnl]) + as.numeric(flim[2])
          )
      } else {
        fl <- flim
      }
      # fix if lower than 0
      if (fl[1] < 0) {
        fl[1] <- 0
      }
      
      # fix higher if above nyquist frequency
      if (fl[2] > wave@samp.rate / 2000) {
        fl[2] <- wave@samp.rate / 2000
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
        ovlp =  ovlp,
        ...
      )
      
      # add vertical lines
      # add dotted lines
      abline(
        v = c(X$mar.start[indx], X$end[indx] - X$start[indx] + X$mar.start[indx]),
        col = "white",
        lty = 3,
        lwd = 1.5
      )
      
      if (spectrum) {
        # set screen
        i <- i + 1
        screen(i)
        par(mar = c(0, 0, 0, 0),
            new = TRUE)
        
        # get power spectrum
        spc <- spec(wave = wave, plot = FALSE)
        
        # smooth
        spc[, 2] <-
          warbleR::envelope(x = spc[, 2], ssmooth = env.smooth)
        
        # filter to flim
        spc <- spc[spc[, 1] > fl[1] & spc[, 1] < fl[2],]
        
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
        rect(
          min(spc[, 2]),
          min(spc[, 1]),
          max(spc[, 2]),
          max(spc[, 1]),
          col = "white",
          border = NA
        )
        rect(
          min(spc[, 2]),
          min(spc[, 1]),
          max(spc[, 2]),
          max(spc[, 1]),
          col = bg_sp_env,
          border = NA
        )
        
        # add 0s at star and end so polygon doesnt twist
        spc[c(1, nrow(spc)), 2] <- 0
        
        # add polygon with spectrum shape
        polygon(spc[, 2:1], col = spc_fill, border = NA)
        
        box()
        
        par(
          mar = c(0, 0, 0, 0),
          bg =  "#FFFFFF00",
          new = TRUE
        )
      }
      
      # plot frequency ticks
      if (page_layout[i, 1] <= min(page_layout[seq_along(nrow * ncol * (sum(c(
        envelope, spectrum)) + 1)), 1]) |
        is.na(prev_sgnl) & curr_dist > distances[1]) {
        at_freq <-
          pretty(seq(fl[1], fl[2], length.out = 10)[-10], n = 10)
        axis(2,
             at = at_freq,
             labels = if (envelope) {
               at_freq
             } else {
               c(at_freq[-length(at_freq)], "")
             })
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
            ssmooth = env.smooth
          )
        # punt in a matrix including time
        envlp <-
          cbind(seq(0, duration(wave), along.with = envlp), envlp)
        
        
        # set graphic parameters
        par(mar = c(0, 0, 0, 0),
            new = TRUE)
        
        
        # set white plot
        plot(
          x = envlp[, 1],
          y = envlp[, 2],
          type = "l",
          frame.plot = FALSE,
          yaxt = "n",
          xaxt = "n",
          col = spc_fill,
          xaxs = "i",
          yaxs = "i"
        )
        
        # add background color
        rect(
          0,
          min(envlp[, 2]),
          max(envlp[, 2]),
          max(envlp[, 2]),
          max(envlp[, 2]),
          col = "white",
          border = NA
        )
        rect(
          0,
          min(envlp[, 2]),
          max(envlp[, 2]),
          max(envlp[, 2]),
          col = bg_sp_env,
          border = NA
        )
        
        # add 0s at star and end so polygon doesnt twist
        envlp[c(1, nrow(envlp)), 2] <- 0
        
        # add polygon with envelope shape
        polygon(envlp,
                col = spc_fill,
                border = NA)
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
             c("", at_time[c(-1,-length(at_time))], "")
           })
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
    par(mar = c(0, 0, 0, 0),
        new = TRUE)
    
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
    par(mar = c(0, 0, 0, 0),
        new = TRUE)
    
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
  cs <- close.screen(all.screens = TRUE)
  try(grDevices::dev.off(), silent = TRUE)
}

# prepare data, used by plot_degradation
.prep_data_plot_degrad <- function(X, distances, nrow){
  X_df <- as.data.frame(X)
  
  # create data subsets (element in list) with all the copies of a sound id in transect, also including its reference if it comes from another transect
  # add text break if longer than 17 characters
  tailored_sound_id <-
    vapply(as.character(X_df$sound.id), function(x) {
      if (nchar(x) > 15) {
        paste0(substr(x, 0, 15), "\n", substr(x, 16, nchar(x)))
      } else {
        x
      }
    }, FUN.VALUE = character(1))
  
  X_df$sound.id.transect <-
    paste(tailored_sound_id, X$transect, sep = "\n")
  
  
  soundid_X_list <-
    lapply(unique(X_df$sound.id.transect), function(x) {
      Y <-
        X_df[X_df$.sgnl.temp %in% (unique(X_df$reference[X_df$sound.id.transect == x])) |
               X_df$sound.id.transect == x,]
      Y <- Y[order(Y$distance),]
      return(Y)
    })
  
  # remove those with only 1 test sound which is a reference (i.e. lowest distance)
  soundid_X_list <-
    lapply(soundid_X_list, function(x) {
      if (nrow(x) == 1 & x$distance[1] == min(distances)) {
        x <- NULL
      }
      
      return(x)
    })
  
  # add numeric id to each subset
  soundid_X_list <-
    lapply(seq_along(soundid_X_list), function(x) {
      if (!is.null(soundid_X_list[[x]])) {
        soundid_X_list[[x]]$seq_number <- x
      }
      return(soundid_X_list[[x]])
    })
  
  # put all into a single data frame
  soundid_X <- do.call(rbind, soundid_X_list)
  
  # sort by sound id as in X
  # Create a factor with desired order
  factor_order <-
    factor(soundid_X$sound.id, levels = unique(X$sound.id))
  
  # Sort
  soundid_X <- soundid_X[order(factor_order),]
  
  # add page in which will be printed each subset
  soundid_X$page <-
    as.numeric(cut(soundid_X$seq_number, breaks = nrow * 0:max(soundid_X$seq_number)))
  
  # add seq number to know which will be plotted together
  soundid_X$sound.id.seq <-
    paste(soundid_X$sound.id, soundid_X$seq_number, sep = "-")
  
  return(soundid_X)
}

# set multipanel layout, used byplot_degradation
.page_layout <- function(ncol, nrow, spectrum, envelope, heights, widths){
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
      psp <- envel <- spectr <- spec_m[y,]
      
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
        matrix(c(spectr, if (spectrum) {
          psp
        }, if (envelope) {
          envel
        }),
        ncol = 4,
        byrow = TRUE)
      
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
    c(horizontals[-length(horizontals)],
      rep(1, ncol),
      max(horizontals),
      min(horizontals))
  
  # put together in a matrix
  lab_m <- cbind(lab_lf, lab_rg, lab_btm, lab_tp)
  
  # single data frame with all panels
  page_layout <- rbind(spec_m, lab_m)
  
  # testing layout screens
  # ss <- split.screen(figs = page_layout)
  # for(i in seq_len(nrow(page_layout)))
  # {screen(i)
  #   par( mar = rep(0, 4))
  #   plot(0.5, xlim = c(0,1), ylim = c(0,1), type = "n", axes = FALSE, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
  #   box()
  #   text(x = 0.5, y = 0.5, labels = i)
  # }
  
  return(page_layout)
}

# for each sound fix start and end based on margin
.fix_margins <- function(X, i, path, margins){
  rangs <- .time_freq_range_files(X, i, path, margins)
  
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
  return(X)
}

# add top panel with master file, used by manual_realign
.master_panel <- function(Y,
                          path,
                          duration,
                          marker,
                          mar,
                          ovlp,
                          wl,
                          flim,
                          palette,
                          fast.spec,
                          grid,
                          collevels,
                          label.col,
                          srt,
                          cex,
                          ...) {
  # check duration
  # import sound data
  master_wave_info <- read_wave(index = 1,
                                Y,
                                header = TRUE,
                                path = path)
  
  # fix to if higher than sound file duration
  master_dur <-
    master_wave_info$samples / master_wave_info$sample.rate
  
  # sto if duration of spectrogram larger than master sound file
  if (master_dur < duration) {
    .stop("'duration' cannot be larger than length of master sound file")
  }
  
  ## plot master spectrogram ####
  # set start and end of clip to plot
  if (Y$start[Y$sound.id == marker] - mar >= 0) {
    from <- Y$start[Y$sound.id == marker] - mar
  } else {
    from <- 0
  }
  to <- from + duration
  
  # set margin at the end
  if (to - Y$start[Y$sound.id == marker] < mar) {
    to <- Y$start[Y$sound.id == marker] + mar
    from <- to - duration
  }
  
  # change to if larger than sound file
  if (to > master_dur) {
    to <- master_dur
    from <- to - duration
  }
  
  # set margins from the start of marker to be used by spectrograms of test sounds
  left_mar <- Y$start[Y$sound.id == marker] - from
  right_mar <- to - Y$start[Y$sound.id == marker]
  
  # import sound data
  master_wave <- read_wave(
    Y,
    index = 1,
    from = from,
    to = to,
    path = path
  )
  
  # frequency label
  screen(3)
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
    cex = 1.4
  )
  
  # top box
  screen(6)
  par(mar = c(0, 0, 0, 0),
      new = TRUE)
  
  plot(
    1,
    frame.plot = FALSE,
    type = "n",
    yaxt = "n",
    xaxt = "n"
  )
  
  # add background color
  rect(0, 0, 2, 2, col = "white", border = NA)
  rect(0, 0, 2, 2, col = "#366A9FFF", border = NA)
  
  top_lab <-
    vapply(as.character(Y$sound.files[1]), function(x) {
      if (nchar(x) > 20)
        paste0(substr(x, 0, 20), "\n", substr(x, 21, nchar(x))) else
          x
    }, FUN.VALUE = character(1))
  
  text(
    x = 1,
    y = 1,
    top_lab,
    cex = 1.2,
    col = "white",
    font = 2,
    srt = 270
  )
  box()
  
  # top spectrogram
  screen(1)
  par(mar = c(0, 2, 0, 0))
  
  # set flim (mostly for sound.id labels below)
  if (is.null(flim)) {
    flim <-
      c(0, master_wave@samp.rate / 2000.1)
  } # use 2000.1 to avoid errors at the highest of nyquist frequency
  
  # plot spectrogram
  warbleR:::spectro_wrblr_int2(
    wave = master_wave,
    collevels = collevels,
    ovlp = ovlp,
    wl = wl,
    flim = flim,
    palette = palette,
    grid = FALSE,
    axisX = FALSE,
    axisY = FALSE,
    flab = "",
    tlab = "",
    fast.spec = fast.spec,
    ...
  )
  
  # plot grid
  if (grid > 0)
    abline(v = seq(grid, duration(master_wave), by = grid), lty = 4)
  
  # add Freq axis
  at_freq <-
    pretty(seq(0, master_wave@samp.rate / 2000, length.out = 6)[-6], n = 6)
  axis(2,
       at = at_freq,
       labels =
         c("", at_freq[-1]))
  
  
  # add dotted lines
  abline(
    v = c(Y$start - from, Y$end - from),
    col = label.col,
    lty = 3,
    lwd = 1.5
  )
  
  # add sound.id labels
  # position of sound.id labels in the freq axis
  y_pos <- flim[2] - 2 * ((flim[2] - flim[1]) / 12)
  
  # plot sound id labels
  text(
    labels = Y$sound.id,
    x = ((Y$end + Y$start) / 2) - from,
    y = y_pos,
    pos = 3,
    col = label.col,
    srt = srt,
    cex = cex
  )
  
  output <-
    list(
      left_mar = left_mar,
      right_mar = right_mar,
      left_mar = left_mar,
      master_wave = master_wave
    )
  
  return(output)
}


# add "stop_by_user" and "last file" messages, used by manual_realign
.stop_by_user <- function(i, xy, xs, grYs, rerec_files) {
  if (xy$x > min(xs) &
      xy$x < max(xs) &
      xy$y > min(grYs$stop) &
      xy$y < max(grYs$stop) & i != length(rerec_files)) {
    text(
      x = mean(par("usr")[1:2]),
      y = mean(par("usr")[3:4]),
      labels = "Stopped by user",
      cex = 1.2,
      col = "black",
      font = 2
    )
  }
  
  if (xy$x > min(xs) &
      xy$x < max(xs) &
      xy$y > min(grYs$`next`) &
      xy$y < max(grYs$`next`) |
      xy$x > min(xs) &
      xy$x < max(xs) &
      xy$y > min(grYs$stop) &
      xy$y < max(grYs$stop) & i == length(rerec_files)) {
    text(
      x = mean(par("usr")[1:2]),
      y = mean(par("usr")[3:4]),
      labels = "All sound files were re-aligned",
      cex = 1.2,
      col = "black",
      font = 2
    )
  }
}

# add bottom responsive panel to plot, used by manual_realign
.responsive_panel <-
  function(X,
           rerec_files,
           i,
           marker,
           step_sum,
           path,
           omp,
           flim,
           collevels,
           ovlp,
           wl,
           palette,
           fast.spec,
           grid,
           label.col,
           border_col,
           fill_col,
           ...) {
    # get data subset for 1 file
    W <- X[X$sound.files == rerec_files[i], ]
    
    # bottom box
    screen(5)
    par(mar = c(0, 0, 0, 0),
        new = TRUE)
    
    plot(
      1,
      frame.plot = FALSE,
      type = "n",
      yaxt = "n",
      xaxt = "n"
    )
    
    # add background color
    rect(0, 0, 2, 2, col = "white", border = NA)
    rect(0, 0, 2, 2, col = "#366A9FFF", border = NA)
    
      bottom_lab <-
      vapply(as.character(W$sound.files[1]), function(x) {
        if (nchar(x) > 15)
          paste0(substr(x, 0, 15), "\n", substr(x, 16, nchar(x))) else
            x
      }, FUN.VALUE = character(1))
    
    
    # spectrogram title vertical boxes
    text(
      x = 1,
      y = 1,
      bottom_lab,
      cex = 1.2,
      col = "white",
      font = 2,
      srt = 270
    )
    box()
    
    # plot spectrogram
    screen(2)
    par(mar = c(0, 2, 0, 0))
    
    from_stp <-
      W$start[W$sound.id == marker] - omp$left_mar + step_sum
    to_stp <- W$start[W$sound.id == marker] + omp$right_mar + step_sum
    
    # import sound data
    wave_info <- read_wave(index = 1,
                           W,
                           header = TRUE,
                           path = path)
    
    # fix to if higher than sound file duration
    wave_dur <- wave_info$samples / wave_info$sample.rate
    
    # read wave
    wave <- read_wave(
      W,
      index = 1,
      from = if (from_stp < 0)
        0 else
          from_stp,
      to = if (to_stp > wave_dur)
        wave_dur else
          to_stp,
      path = path
    )
    
    if (from_stp < 0) {
      wave <-
        seewave::pastew(
          wave1 = wave,
          wave2 = tuneR::silence(
            duration = abs(from_stp),
            xunit = "time",
            samp.rate = wave@samp.rate
          ),
          output = "Wave"
        )
    }
    
    if (to_stp > wave_dur) {
      wave <-
        seewave::pastew(
          wave2 = wave,
          wave1 = tuneR::silence(
            duration = to_stp - wave_dur,
            xunit = "time",
            samp.rate = wave@samp.rate
          ),
          output = "Wave"
        )
    }
    
    # set flim (mostly for sound.id labels below)
    if (is.null(flim)) {
      flim <-
        c(0, wave@samp.rate / 2000.1)
    } # use 2000.1 to avoid errors at the highest of nyquist frequency


    
    # plot spectrogram
    warbleR:::spectro_wrblr_int2(
      wave = wave,
      collevels = collevels,
      ovlp = ovlp,
      wl = wl,
      flim = flim,
      palette = palette,
      grid = FALSE,
      axisX = FALSE,
      axisY = FALSE,
      flab = "",
      tlab = "",
      fast.spec = fast.spec,
      ...
    )
      
    # add white rectangle at the begining if spectrogram shorter than view range
    if (from_stp < 0)
      rect(
        par("usr")[1],
        par("usr")[3],
        abs(from_stp),
        par("usr")[4],
        col = "white",
        border = NA
      )
  
    # add white rectangle at the end if spectrogram shorter than view range
  if (to_stp > wave_dur)
      rect(
        par("usr")[2] - to_stp + wave_dur,
        par("usr")[3],
        par("usr")[2],
        par("usr")[4],
        col = "white",
        border = NA
      )
    
    # plot grid
    if (grid > 0)
      abline(v = seq(grid, duration(omp$master_wave), by = grid), lty = 4)
    
    # progress bar
    prct <- (i / length(rerec_files))
    y <- grconvertY(y = 0.99, from = "npc", to = "user")
    lines(
      x = c(0, grconvertX(
        x = prct - 0.02,
        from = "npc",
        to = "user"
      )),
      y = rep(y, 2),
      lwd = 7,
      col = adjustcolor("#E37222", alpha.f = 0.6),
      xpd = TRUE
    )
    text(
      x = grconvertX(x = prct, from = "npc", to = "user"),
      y = y,
      xpd = TRUE,
      labels = paste0(floor(prct * 100), "%"),
      col = "#E37222",
      cex = 0.8
    )
    
    # add dotted lines on the position of test sounds
    abline(
      v = c(W$start - from_stp + step_sum, W$end - from_stp + step_sum),
      col = label.col,
      lty = 3,
      lwd = 1.5
    )
    
    # add Freq axis
    at_freq <-
      pretty(seq(0, wave@samp.rate / 2000, length.out = 6)[-6], n = 6)
    axis(2,
         at = at_freq,
         labels = c(at_freq[-length(at_freq)], ""))
    
    # negative fake coordinates to get started
    xy <- list(x = -1000, y = -1000)
    
    # plot time ticks
    if (i == 1 & step_sum == 0 & xy$x == -1000 & xy$y == -1000) {
      at_time <-
        pretty(seq(0, duration(wave), length.out = 7)[-7], n = 7)
      axis(1,
           at = at_time,
           labels = at_time)
      
      mtext(
        text = "Time (s)",
        side = 1,
        line = 2.2,
        cex = 1.4
      )
    }
    
    # add buttons
    xs <-
      grconvertX(x = c(0.92, 0.92, 0.99, 0.99),
                 from = "npc",
                 to = "user")
    labels <-
      c("stop",
        "long left",
        "long right",
        "short left",
        "short right",
        "next",
        "reset")
    
    # relative position of buttons in y axis
    cpy <- seq(0.93, 0.07, length.out = length(labels))
    
    mrg <- (cpy[1] - cpy[2]) / 3
    
    # mid position ofbuttons in c(0, 1) range
    ys <- c(-mrg, mrg, mrg, -mrg)
    
    grYs <- lapply(seq_len(length(labels)), function(x) {
      grY <- grconvertY(y = cpy[x] - ys,
                        from = "npc",
                        to = "user")
      polygon(
        x = xs,
        y = grY,
        border = border_col,
        col = fill_col,
        lwd = 2
      )
      # "#440154FF" "#482878FF" "#3E4A89FF" "#31688EFF" "#26828EFF" "#1F9E89FF"
      # [7] "#35B779FF" "#6DCD59FF" "#B4DE2CFF" border_col
      # plot symbols
      if (labels[x] == "stop") {
        points(
          x = mean(xs),
          y = mean(grY),
          pch = 15,
          cex = 1.5,
          col = border_col
        )
      }
      
      if (labels[x] == "long right") {
        text(
          x = mean(xs),
          y = mean(grY),
          labels = ">>",
          cex = 1.2,
          font = 2,
          col = border_col
        )
      }
      
      if (labels[x] == "long left") {
        text(
          x = mean(xs),
          y = mean(grY),
          labels = "<<",
          cex = 1.2,
          font = 2,
          col = border_col
        )
      }
      
      if (labels[x] == "short right") {
        text(
          x = mean(xs),
          y = mean(grY),
          labels = ">",
          cex = 1.2,
          font = 2,
          col = border_col
        )
      }
      
      if (labels[x] == "short left") {
        text(
          x = mean(xs),
          y = mean(grY),
          labels = "<",
          cex = 1.2,
          font = 2,
          col = border_col
        )
      }
      
      if (labels[x] == "next") {
        text(
          x = mean(xs),
          y = mean(grY),
          labels = "next",
          cex = 1.2,
          font = 2,
          col = border_col
        )
      }
      
      if (labels[x] == "reset") {
        text(
          x = mean(xs),
          y = mean(grY),
          labels = "reset",
          cex = 1.2,
          font = 2,
          col = border_col
        )
      }
      
      return(grY)
    })
    
    names(grYs) <- labels
    
    output <- list(grYs = grYs, xs = xs, xy = xy)
    
    return(output)
  }

# cycle while no next, used by manual_realign
.not_next <-
  function(xy,
           rp,
           step_sum,
           step_sum_vector,
           step.lengths,
           i) {
    while (!(
      xy$x > min(rp$xs) &
      xy$x < max(rp$xs) &
      xy$y > min(rp$grYs$stop) &
      xy$y < max(rp$grYs$stop)
    ) &
    !(
      xy$x > min(rp$xs) &
      xy$x < max(rp$xs) &
      xy$y > min(rp$grYs$`next`) &
      xy$y < max(rp$grYs$`next`)
    )) {
      xy <- locator(n = 1, type = "n")
      
      # if reset
      if (xy$x > min(rp$xs) &
          xy$x < max(rp$xs) &
          xy$y > min(rp$grYs$reset) & xy$y < max(rp$grYs$reset)) {
        step_sum <- 0
        step_sum_vector[i] <- step_sum
        break
      }
      # if next
      if (xy$x > min(rp$xs) &
          xy$x < max(rp$xs) &
          xy$y > min(rp$grYs$`next`) &
          xy$y < max(rp$grYs$`next`)) {
        i <- i + 1
        step_sum <- 0
        step_sum_vector[i] <- step_sum
        break
      }
      
      # if long left
      if (xy$x > min(rp$xs) &
          xy$x < max(rp$xs) &
          xy$y > min(rp$grYs$`long left`) &
          xy$y < max(rp$grYs$`long left`)) {
        step_sum <- step_sum + max(step.lengths / 1000)
        step_sum_vector[i] <- step_sum
        break
      }
      
      # if long right
      if (xy$x > min(rp$xs) &
          xy$x < max(rp$xs) &
          xy$y > min(rp$grYs$`long right`) &
          xy$y < max(rp$grYs$`long right`)) {
        step_sum <- step_sum - max(step.lengths / 1000)
        step_sum_vector[i] <- step_sum
        break
      }
      
      # if short left
      if (xy$x > min(rp$xs) &
          xy$x < max(rp$xs) &
          xy$y > min(rp$grYs$`short left`) &
          xy$y < max(rp$grYs$`short left`)) {
        step_sum <- step_sum + min(step.lengths / 1000)
        step_sum_vector[i] <- step_sum
        break
      }
      
      # if short right
      if (xy$x > min(rp$xs) &
          xy$x < max(rp$xs) &
          xy$y > min(rp$grYs$`short right`) &
          xy$y < max(rp$grYs$`short right`)) {
        step_sum <- step_sum - min(step.lengths / 1000)
        step_sum_vector[i] <- step_sum
       
        break
      }
      
      # if stop
      if (xy$x > min(rp$xs) &
          xy$x < max(rp$xs) &
          xy$y > min(rp$grYs$stop) & xy$y < max(rp$grYs$stop)) {
        break
      }
    }
    
    output <-
      list(
        step_sum = step_sum,
        step_sum_vector = step_sum_vector,
        i = i,
        xy = xy
      )
  }

# adjust annotation time coordiantes, used by manual_realign
.adjust_coors <- function(rerec_files, X, step_sum_vector, arg_call){
  
  for (u in rerec_files) {
    X$start[X$sound.files == u] <-
      X$start[X$sound.files == u] + step_sum_vector[u]
    X$end[X$sound.files == u] <-
      X$end[X$sound.files == u] + step_sum_vector[u]
    
    # fix internally for extended selection tables
    if (is_selection_table(X) | is_extended_selection_table(X)) {
      attributes(X)$check.res$start[attributes(X)$check.res$sound.files == u] <-
        attributes(X)$check.res$start[X$sound.files == u] + step_sum_vector[u]
      attributes(X)$check.res$end[attributes(X)$check.res$sound.files == u] <-
        attributes(X)$check.res$end[X$sound.files == u] + step_sum_vector[u]
      
    }
  }
  
  # fix call attribute
  attributes(X)$call <- arg_call
  
  return(X)
}

### CHECK ARGUMENTS  ###
.report_assertions <- function(collection) {
  checkmate::assertClass(collection, "AssertCollection")
  if (!collection$isEmpty()) {
    msgs <- collection$getMessages()
    
    # modfied to get more informative message
    msgs <-
      gsub(
        pattern = "Variable 'cores': Element 1 is not <=",
        replacement = "The maximum number of cores available on this computer is",
        x = msgs,
        fixed = TRUE
      )
    
    msgs <-
      gsub(
        pattern = "Variable 'X$distances':",
        replacement = "",
        x = msgs,
        fixed = TRUE
      )
    
    msgs <-
      gsub(
        pattern = "Variable 'X$sound.files': ",
        replacement = "",
        x = msgs,
        fixed = TRUE
      )
    
    msgs <-
      gsub(
        pattern = "Variable 'names(X)': Names",
        replacement = "Columns in 'X':",
        x = msgs,
        fixed = TRUE
      )
    
    msgs <-
      gsub(
        pattern = "Variable 'transect column': Names must include the elements {'transect'}, but is missing elements {'transect'}.",
        replacement = "Column 'transect' in 'X' is required when 'method = 2'",
        x = msgs,
        fixed = TRUE
      )
    
    
    msgs <-
      gsub(
        pattern = "the elements ",
        replacement = "",
        x = msgs,
        fixed = TRUE
      )
    
    msgs <-
      gsub(
        pattern = "elements ",
        replacement = "",
        x = msgs,
        fixed = TRUE
      )
    
    msgs <-
      gsub(
        pattern = "Variable",
        replacement = "Argument",
        x = msgs,
        fixed = TRUE
      )
    
    msgs <-
      gsub(
        pattern = "must be disjunct from {'reference'}, but has {'reference'}.",
        replacement = "Cannot include a column 'reference'. It must be removed.",
        x = msgs,
        fixed = TRUE
      )
    
    msgs <-
      gsub(
        pattern = "Argument 'X$distance'",
        replacement = "Column 'distance' in 'X'",
        x = msgs,
        fixed = TRUE
      )
    
    msgs <-
      gsub(
        pattern = "Argument 'transect column': Names must include the elements {'transect'}, but is missing elements {'transect'}",
        replacement = "Argument 'X': requires a 'transect' column when 'method = 2'",
        x = msgs,
        fixed = TRUE
      )
    
    msgs <-
      gsub(
        pattern = "but is missing {'reference'}",
        replacement = "but is missing {'reference'} (Did you forget to run 'set_reference_sounds()' first?)",
        x = msgs,
        fixed = TRUE
      )
    
    msgs <-
      gsub(
        pattern = "Argument 'transect values'",
        replacement = "Column 'transect' in 'X'",
        x = msgs,
        fixed = TRUE
      )
    
    
    context <- "\n %i argument check(s) failed:"
    err <- c(sprintf(context, length(msgs)), strwrap(msgs,
                                                     prefix = " * "))
    stop(simpleError(paste0(err, collapse = "\n"), call = sys.call(1L)))
  }
  invisible(TRUE)
}

# custom assert functions to check arguments
# check no duplicated selection labels
.check_unique_sels <- function(x, fun) {
  if (anyDuplicated(paste(x$sound.files, x$selec)) > 0) {
    "Duplicated 'selec' labels within at least one sound file"
  } else {
    TRUE
  }
}

.assert_unique_sels <-
  checkmate::makeAssertionFunction(.check_unique_sels)

# check if X has any ambient reference in sound id column
.check_ambient_ref <- function(x, noise.ref) {
    if (noise.ref == "custom" & !any(x$sound.id == "ambient")) {
      "'noise.ref = custom' but no 'ambient' label found in 'sound.id' column"
  } else {
    TRUE
  }
}

.assert_ambient_ref <-
  checkmate::makeAssertionFunction(.check_ambient_ref)


# check sound files exist
.check_sound_files_found <- function(files, fun) {
  if (all(!file.exists(files))) {
    if (!fun != "noise_profile") {
      "Not a single sound files in 'X$sound.files' can be found. Use 'path' to set the directory where the sound files are found"
    } else {
      "Not a single sound files in 'files' can be found. Use 'path' to set the directory where the sound files are found"
    }
  } else {
    if (!all(file.exists(files))) {
      if (!fun != "noise_profile") {
        "Some sound files in 'X$sound.files' cannot be found. Make sure all sound files referenced in 'X' are in the 'path' supplied (or current working directory if 'path' was not supplied)"
      } else {
        "Some sound files in 'files' cannot be found. Make sure all sound files are in the 'path' supplied (or current working directory if 'path' was not supplied)"
      }
    } else {
      TRUE
    }
  }
}

.assert_sound_files_found <-
  checkmate::makeAssertionFunction(.check_sound_files_found)


# check unique sound.id
.check_unique_sound.id <- function(x, fun) {
  if (!is.null(x$sound.id) & is.null(x$distance)) {
    if (anyDuplicated(paste0(x$sound.files, x$sound.id)) > 0) {
      "Duplicated 'sound.id' labels are not allowed within a sound file"
    }
  }
  if (!is.null(x$sound.id) & !is.null(x$distance)) {
    if (anyDuplicated(paste0(x$sound.files, x$sound.id, x$distance)) > 0) {
      "Duplicated 'sound.id' labels are not allowed within a sound file or sound file/distance combination"
    }
  }
  if (is.null(x$sound.id) & is.null(x$distance)) {
    TRUE
  }
}

.assert_unique_sound.id <-
  checkmate::makeAssertionFunction(.check_unique_sound.id)

# check than more than 1 distance is found
.check_several_distances <- function(x, fun) {
  if (!is.null(x$distance)) {
    if (length(unique(x$distance)) == 1) {
      "Column 'distance'in 'X' must include more than 1 distance"
    }
  } else {
    TRUE
  }
}

.assert_several_distances <-
  checkmate::makeAssertionFunction(.check_several_distances)

# check if argument has been deprecated
.check_deprecated <- function(x) {
  if (!is.null(x)) {
    "has been deprecated"
  } else {
    TRUE
  }
}

.assert_deprecated <-
  checkmate::makeAssertionFunction(.check_deprecated)

# check if X is an extended selection table
.check_extended_selection_table <- function(x) {
  if (!is.null(x)) {
    if (!is_extended_selection_table(x)) {
      "'X' must be of class 'extended_selection_table'"
    }
  } else {
    TRUE
  }
}

.assert_extended_selection_table <-
  checkmate::makeAssertionFunction(.check_extended_selection_table)

.check_est_by_element <- function(x) {
  if (!is.null(x)) {
    if (is_extended_selection_table(x)) {
      if (attr(x, "by.song")[[1]]) {
        "Extended selection table 'X' must be created 'by element', not 'by song'. Use warbleR::by_element_est(X) to convert it to the right format"
      }
    }
  } else {
    TRUE
  }
}

.assert_est_by_element <-
  checkmate::makeAssertionFunction(.check_est_by_element)


.check_est_by_song <- function(x) {
  if (!is.null(x)) {
    if (is_extended_selection_table(x)) {
      if (!attr(x, "by.song")[[1]]) {
        "Extended selection table must be created 'by song', not 'by element' to be used in manual_realign()"
      }
    }
  } else {
    TRUE
  }
}

.assert_est_by_song <-
  checkmate::makeAssertionFunction(.check_est_by_song)

.check_no_margin <- function(x) {
  if (!is.null(x)) {
    if (is_extended_selection_table(x)) {
      if (any(x$start == 0))
        "Some annotations have no margin in which to measure background noise (X$start == 0)"
      
    }
  }
  else {
    TRUE
  }
}

.assert_no_margin <-
  checkmate::makeAssertionFunction(.check_no_margin)

# check than more than 1 distance is found
.check_several_distances <- function(x, fun) {
  if (!is.null(x$distance)) {
    if (length(unique(x$distance)) == 1) {
      "Column 'distance'in 'X' must include more than 1 distance"
    }
  } else {
    TRUE
  }
}

.assert_several_distances <-
  checkmate::makeAssertionFunction(.check_several_distances)

## function to check arguments
.check_arguments <- function(fun, args) {
  # make function name a character
  fun <- as.character(fun)
  
  # create object to store check results
  check_collection <- checkmate::makeAssertCollection()
  
  # run functions to check arguments
  check_collection <- .check_X(args, check_collection, fun)
  
  check_collection <- .check_envelopes(args, check_collection)
  
  check_collection <- .check_spectra(args, check_collection)
  
  check_collection <- .check_sampling.rate(args, check_collection)
  
  check_collection <- .check_shuffle(args, check_collection)
  
  check_collection <- .check_seed(args, check_collection)
  
  check_collection <- .check_srt(args, check_collection)
  
  check_collection <- .check_sig2(args, check_collection)
  
  check_collection <- .check_fm(args, check_collection)
  
  check_collection <- .check_am(args, check_collection)
  
  check_collection <- .check_frequencies(args, check_collection)
  
  check_collection <- .check_durations(args, check_collection)
  
  check_collection <- .check_am.amps(args, check_collection)
  
  check_collection <- .check_hrm.freqs(args, check_collection)
  
  check_collection <- .check_hrm.freqs(args, check_collection)
  
  check_collection <- .check_noise.ref(args, check_collection)
  
  check_collection <- .check_label(args, check_collection)
  
  check_collection <- .check_fast.spec(args, check_collection)
  
  check_collection <- .check_width(args, check_collection, fun)
  
  check_collection <- .check_height(args, check_collection, fun)
  
  check_collection <- .check_env.smooth(args, check_collection)
  
  check_collection <- .check_res(args, check_collection)
  
  check_collection <- .check_dest.path(args, check_collection)
  
  check_collection <- .check_flim(args, check_collection)
  
  check_collection <- .check_mar(args, check_collection, fun)
  
  check_collection <- .check_duration(args, check_collection)
  
  check_collection <- .check_collevels(args, check_collection)
  
  check_collection <- .check_palette(args, check_collection)
  
  check_collection <- .check_PSD(args, check_collection)
  
  check_collection <- .check_norm(args, check_collection)
  
  check_collection <- .check_dB(args, check_collection)
  
  check_collection <- .check_averaged(args, check_collection)
  
  check_collection <- .check_frequency(args, check_collection)
  
  check_collection <- .check_temp(args, check_collection)
  
  check_collection <- .check_rh(args, check_collection)
  
  check_collection <- .check_pa(args, check_collection)
  
  check_collection <- .check_dist(args, check_collection)
  
  check_collection <- .check_dist0(args, check_collection)
  
  check_collection <- .check_hab.att.coef(args, check_collection)
  
  check_collection <- .check_spl.cutoff(args, check_collection)
  
  check_collection <- .check_files(args, check_collection, fun)
  
  check_collection <- .check_output(args, check_collection)
  
  check_collection <- .check_marker(args, check_collection, fun)
  
  check_collection <- .check_parallel(args, check_collection)
  
  check_collection <- .check_template.rows(args, check_collection)
  
  check_collection <- .check_cores(args, check_collection)
  
  check_collection <- .check_pb(args, check_collection)
  
  check_collection <- .check_path(args, check_collection)
  
  check_collection <- .check_n.samples(args, check_collection)
  
  check_collection <- .check_hop.size(args, check_collection)
  
  check_collection <- .check_wl(args, check_collection)
  
  check_collection <- .check_bp(args, check_collection)
  
  check_collection <- .check_img(args, check_collection)
  
  check_collection <- .check_col(args, check_collection)
  
  check_collection <- .check_ovlp(args, check_collection)
  
  check_collection <- .check_only.sels(args, check_collection)
  
  check_collection <- .check_method(args, check_collection)
  
  check_collection <- .check_cor.method(args, check_collection)
  
  check_collection <- .check_wn(args, check_collection)

  check_collection <- .check_Y(args, check_collection, fun)
  
  return(check_collection)
}

## single functions to check individual arguments
.check_X <- function(args, check_collection, fun) {
  ### check arguments
  if (any(names(args) == "X")) {
    
    if (!is_extended_selection_table(args$X)) {
      if (!is.null(args$path)) {
        files <- file.path(args$path, unique(args$X$sound.files))
      } else {
        if (!is.null(getOption("sound.files.path"))) {
          files <-
            file.path(getOption("sound.files.path"),
                      unique(args$X$sound.files))
        } else {
          files <- unique(args$X$sound.files)
        }
      }
      .assert_sound_files_found(
        files = files,
        fun = fun,
        add = check_collection,
        .var.name = "X$sound.files"
      )
    }
    
    
    if (fun != "noise_profile") {
      checkmate::assert_data_frame(
        x = args$X,
        any.missing = TRUE,
        min.rows = if (!fun %in% c("tail_to_signal_ratio",
                                   "signal_to_noise_ratio",
                                   "add_noise")) {
          2
        } else {
          1
        },
        add = check_collection,
        .var.name = "X"
      )
      
      if (fun != "add_noise") {
        checkmate::assert_multi_class(
          x = args$X,
          classes = c(
            "data.frame",
            "selection.table",
            "extended.selection.table"
          ),
          add = check_collection,
          .var.name = "X"
        )
        
      if   (fun == "manual_realign"){
        .assert_est_by_song(x = args$X, 
                            add = check_collection,
                            .var.name = "X")
        
      } 
        
      } else {
        .assert_extended_selection_table(x = args$X,
                                         add = check_collection,
                                         .var.name = "X")
        
        .assert_est_by_element(x = args$X, add = check_collection,
                                         .var.name = "X")
      }
      
      # default columns
      cols <- c("sound.files", "selec", "start", "end", "sound.id")
      
      # overwrite for other functions
      if (fun == "master_sound_file") {
        cols <- c("sound.files", "selec", "start", "end")
      }
      # functions that compare by distance
      if (fun %in% c(
        "blur_ratio",
        "plot_blur_ratio",
        "detection_distance",
        "envelope_correlation",
        "excess_attenuation",
        "spcc",
        "spectrum_blur_ratio",
        "spectrum_correlation",
        "plot_degradation"
      )) {
        cols <-
          c("sound.files",
            "selec",
            "start",
            "end",
            "sound.id",
            "distance",
            "reference")
        
        checkmate::assert_numeric(
          x = args$X$distance,
          any.missing = FALSE,
          all.missing = FALSE,
          lower = 0,
          null.ok = TRUE,
          add = check_collection,
          .var.name = "X$distance"
        )
        
        .assert_several_distances(
          x = args$X,
          fun = fun,
          add = check_collection,
          .var.name = "X$distances"
        )
      }
      
      
      if (fun == "set_reference_sounds") {
        cols <- c("sound.files",
                  "selec",
                  "start",
                  "end",
                  "sound.id",
                  "distance")
        checkmate::assert_numeric(
          x = args$X$distance,
          any.missing = FALSE,
          all.missing = FALSE,
          lower = 0,
          null.ok = TRUE,
          add = check_collection,
          .var.name = "X$distance"
        )
      }
   
      checkmate::assert_names(
        x = names(args$X),
        type = "unique",
        must.include = cols,
        add = check_collection,
        disjunct.from = if (fun == "set_reference_sounds") {
          "reference"
        } else {
          NULL
        },
        .var.name = "names(X)"
      )
      try(checkmate::assert_data_frame(
        x = args$X[, cols],
        any.missing = TRUE,
        add = check_collection,
        .var.name = "X"
      ),
      silent = TRUE)
      
      .assert_unique_sels(
        x = args$X,
        fun = fun,
        add = check_collection,
        .var.name = "X"
      )
      
      .assert_unique_sound.id(
        x = args$X,
        fun = fun,
        add = check_collection,
        .var.name = "X"
      )
      
      if (!fun %in% c(
        "plot_aligned_sounds",
        "degrad_catalog",
        "detection_distance",
        "signal_to_noise_ratio",
        "tail_to_signal_ratio",
        "add_noise"
      )) {
        if (is_extended_selection_table(args$X)) {
          if (length(unique(attr(args$X, "check.results")$sample.rate)) > 1) {
            .stop(
              "all wave objects in the extended selection table must have the same sampling rate (they can be homogenized using warbleR::resample_est())"
            )
          }
        } else {
          if (!fun %in% c("find_markers",
                          "align_test_files",
                          "set_reference_sounds",
                          "add_noise")) {
            .warning("assuming all sound files have the same sampling rate")
          }
          
          if (fun == "add_noise") {
            if (attr(args$X, "by.song")[[1]]) {
              .stop(
                "Extended selection table 'X' must be created 'by element', not 'by song'. Use warbleR::by_element_est(X) to convert it to the right format."
              )
            }
            if (!is.null(args$X$start)) {
              if (any(args$X$start == 0))
                .stop(
                  "Some annotations have no margin in which to measure background noise (X$start == 0)"
                )
            }
          }
        }
      }
    } else {
      checkmate::assert_data_frame(
        x = args$X,
        any.missing = TRUE,
        add = check_collection,
        .var.name = "X",
        null.ok = TRUE
      )
      
      .assert_ambient_ref(
        x = args$X, 
        noise.ref = if(!is.null(args$noise.ref)) args$noise.ref  else "adjacent",
        add = check_collection,
        .var.name = "X"
      )
      
      
    }
  }
  return(check_collection)
}

.check_envelopes <- function(args, check_collection) {
  if (any(names(args) == "envelopes")) {
    checkmate::assert_logical(
      x = args$envelopes,
      len = 1,
      null.ok = FALSE,
      add = check_collection,
      .var.name = "envelopes",
      any.missing = FALSE
    )
  }
  return(check_collection)
}

.check_spectra <- function(args, check_collection) {
  if (any(names(args) == "spectra")) {
    checkmate::assert_logical(
      x = args$spectra,
      len = 1,
      null.ok = FALSE,
      add = check_collection,
      .var.name = "spectra",
      any.missing = FALSE
    )
  }
  return(check_collection)
}

.check_sampling.rate <- function(args, check_collection) {
  if (any(names(args) == "sampling.rate")) {
    checkmate::assert_number(
      x = args$sampling.rate,
      lower = 1,
      add = check_collection,
      .var.name = "sampling.rate",
      null.ok = FALSE,
      na.ok = FALSE
    )
  }
  return(check_collection)
}

.check_shuffle <- function(args, check_collection) {
  if (any(names(args) == "shuffle")) {
    checkmate::assert_logical(
      x = args$shuffle,
      len = 1,
      null.ok = FALSE,
      add = check_collection,
      .var.name = "shuffle"
    )
  }
  return(check_collection)
}

.check_seed <- function(args, check_collection) {
  if (any(names(args) == "seed")) {
    checkmate::assert_number(
      x = args$seed,
      add = check_collection,
      .var.name = "seed",
      null.ok = TRUE
    )
  }
  return(check_collection)
}

.check_srt <- function(args, check_collection) {
  if (any(names(args) == "srt")) {
    checkmate::assert_number(
      x = args$srt,
      add = check_collection,
      .var.name = "srt",
      null.ok = FALSE
    )
  }
  return(check_collection)
}

.check_sig2 <- function(args, check_collection) {
  if (any(names(args) == "sig2")) {
    checkmate::assert_number(
      x = args$sig2,
      add = check_collection,
      lower = 0.00001,
      .var.name = "sig2",
      null.ok = TRUE
    )
  }
  return(check_collection)
}

.check_fm <- function(args, check_collection) {
  if (any(names(args) == "fm")) {
    checkmate::assert_logical(
      x = args$fm,
      len = 1,
      null.ok = FALSE,
      add = check_collection,
      .var.name = "fm"
    )
  }
  return(check_collection)
}

.check_am <- function(args, check_collection) {
  if (any(names(args) == "am")) {
    checkmate::assert_logical(
      x = args$am,
      len = 1,
      null.ok = FALSE,
      add = check_collection,
      .var.name = "am"
    )
  }
  return(check_collection)
}

.check_frequencies <- function(args, check_collection) {
  if (any(names(args) == "frequencies")) {
    checkmate::assert_numeric(
      x = args$frequencies,
      any.missing = FALSE,
      all.missing = FALSE,
      unique = TRUE,
      lower = 0.001,
      add = check_collection,
      .var.name = "frequencies"
    )
  }
  return(check_collection)
}

.check_durations <- function(args, check_collection) {
  if (any(names(args) == "durations")) {
    checkmate::assert_numeric(
      x = args$durations,
      any.missing = FALSE,
      all.missing = FALSE,
      unique = TRUE,
      lower = 0.001,
      add = check_collection,
      .var.name = "durations"
    )
  }
  return(check_collection)
}

.check_am.amps <- function(args, check_collection) {
  if (any(names(args) == "am.amps")) {
    checkmate::assert_numeric(
      x = args$am.amps,
      any.missing = FALSE,
      all.missing = FALSE,
      unique = TRUE,
      null.ok = TRUE,
      lower = 0.0001,
      add = check_collection,
      .var.name = "am.amps"
    )
  }
  return(check_collection)
}

.check_hrm.freqs <- function(args, check_collection) {
  if (any(names(args) == "hrm.freqs")) {
    checkmate::assert_numeric(
      x = args$hrm.freqs,
      any.missing = FALSE,
      all.missing = FALSE,
      unique = TRUE,
      null.ok = TRUE,
      lower = 0.0001,
      add = check_collection,
      .var.name = "hrm.freqs"
    )
  }
  return(check_collection)
}

.check_noise.ref <- function(args, check_collection) {
  if (any(names(args) == "noise.ref")) {
    checkmate::assert_character(
      x = args$noise.ref,
      null.ok = FALSE,
      add = check_collection,
      .var.name = "noise.ref",
      len = 1,
      any.missing = FALSE,
      all.missing = FALSE,
      ignore.case = FALSE
    )
    checkmate::assert_choice(
      x = args$noise.ref,
      choices = c("adjacent", "custom"),
      add = check_collection,
      .var.name = "noise.ref"
    )
  }
  return(check_collection)
}

.check_label <- function(args, check_collection) {
  if (any(names(args) == "label")) {
    checkmate::assert_logical(
      x = args$label,
      len = 1,
      null.ok = FALSE,
      add = check_collection,
      .var.name = "label"
    )
  }
  return(check_collection)
}

.check_fast.spec <- function(args, check_collection) {
  if (any(names(args) == "fast.spec")) {
    checkmate::assert_logical(
      x = args$fast.spec,
      len = 1,
      null.ok = FALSE,
      add = check_collection,
      .var.name = "fast.spec"
    )
  }
  return(check_collection)
}

.check_width <- function(args, check_collection, fun) {
  if (any(names(args) == "width")) {
    if (fun == "plot_degradation") {
      checkmate::assert_numeric(
        x = args$width,
        lower = 0.0001,
        finite = TRUE,
        any.missing = FALSE,
        all.missing = FALSE,
        len = 2,
        add = check_collection,
        .var.name = "width"
      )
    } else {
      checkmate::assert_number(
        x = args$width,
        lower = 0.1,
        add = check_collection,
        .var.name = "width",
        null.ok = FALSE,
        na.ok = FALSE
      )
    }
  }
  return(check_collection)
}

.check_height <- function(args, check_collection, fun) {
  if (any(names(args) == "height")) {
    if (fun == "plot_degradation") {
      checkmate::assert_numeric(
        x = args$height,
        lower = 0.0001,
        finite = TRUE,
        any.missing = FALSE,
        all.missing = FALSE,
        len = 2,
        add = check_collection,
        .var.name = "height"
      )
    } else {
      checkmate::assert_number(
        x = args$height,
        lower = 0.1,
        add = check_collection,
        .var.name = "height",
        null.ok = FALSE,
        na.ok = FALSE
      )
    }
  }
  return(check_collection)
}

.check_env.smooth <- function(args, check_collection) {
  if (any(names(args) == "env.smooth")) {
    checkmate::assert_number(
      x = args$env.smooth,
      lower = 1,
      add = check_collection,
      .var.name = "env.smooth",
      null.ok = FALSE,
      na.ok = FALSE
    )
  }
  return(check_collection)
}

.check_res <- function(args, check_collection) {
  if (any(names(args) == "res")) {
    checkmate::assert_number(
      x = args$res,
      lower = 1,
      add = check_collection,
      .var.name = "res",
      null.ok = FALSE,
      na.ok = FALSE
    )
  }
  return(check_collection)
}

.check_dest.path <- function(args, check_collection) {
  if (any(names(args) == "dest.path")) {
    checkmate::assert_directory(
      x = args$dest.path,
      access = "r",
      add = check_collection,
      .var.name = "dest.path"
    )
  }
  return(check_collection)
}

.check_flim <- function(args, check_collection) {
  if (any(names(args) == "flim")) {
    checkmate::assert_vector(
      x = args$flim,
      any.missing = FALSE,
      all.missing = FALSE,
      null.ok = FALSE,
      len = 2,
      add = check_collection,
      .var.name = "flim"
    )
  }
  return(check_collection)
}

.check_mar <- function(args, check_collection, fun) {
  if (any(names(args) == "mar")) {
    if (fun == "plot_degradation") {
      checkmate::assert_numeric(
        x = args$mar,
        lower = 0.0001,
        finite = TRUE,
        any.missing = FALSE,
        all.missing = FALSE,
        len = 2,
        add = check_collection,
        .var.name = "mar"
      )
    } else {
      checkmate::assert_number(
        x = args$mar,
        lower = 0.00001,
        add = check_collection,
        .var.name = "mar",
        null.ok = FALSE,
        na.ok = FALSE
      )
    }
  }
  return(check_collection)
}

.check_duration <- function(args, check_collection) {
  if (any(names(args) == "duration")) {
    checkmate::assert_number(
      x = args$duration,
      lower = 0.00001,
      add = check_collection,
      .var.name = "duration",
      null.ok = FALSE,
      na.ok = FALSE
    )
  }
  return(check_collection)
}

.check_collevels <- function(args, check_collection) {
  if (any(names(args) == "collevels")) {
    checkmate::assert_numeric(
      x = args$collevels,
      any.missing = FALSE,
      all.missing = FALSE,
      unique = TRUE,
      upper = 0,
      add = check_collection,
      .var.name = "collevels"
    )
  }
  return(check_collection)
}

.check_palette <- function(args, check_collection) {
  if (any(names(args) == "palette")) {
    checkmate::assert_function(
      x = args$palette,
      null.ok = FALSE,
      add = check_collection,
      .var.name = "palette"
    )
  }
  return(check_collection)
}

.check_PSD <- function(args, check_collection) {
  if (any(names(args) == "PSD")) {
    checkmate::assert_logical(
      x = args$PSD,
      len = 1,
      null.ok = FALSE,
      add = check_collection,
      .var.name = "PSD"
    )
  }
  return(check_collection)
}

.check_norm <- function(args, check_collection) {
  if (any(names(args) == "norm")) {
    checkmate::assert_logical(
      x = args$norm,
      len = 1,
      null.ok = FALSE,
      add = check_collection,
      .var.name = "norm"
    )
  }
  return(check_collection)
}

.check_dB <- function(args, check_collection) {
  if (any(names(args) == "dB")) {
    checkmate::assert_character(
      x = args$dB,
      null.ok = FALSE,
      add = check_collection,
      .var.name = "dB",
      len = 1,
      any.missing = FALSE,
      all.missing = FALSE,
      ignore.case = FALSE
    )
    checkmate::assert_choice(
      x = args$dB,
      choices = c("max0", "A", "B", "C", "D", "ITU"),
      add = check_collection,
      .var.name = "dB"
    )
  }
  return(check_collection)
}

.check_averaged <- function(args, check_collection) {
  if (any(names(args) == "averaged")) {
    checkmate::assert_logical(
      x = args$averaged,
      len = 1,
      null.ok = FALSE,
      add = check_collection,
      .var.name = "averaged"
    )
  }
  return(check_collection)
}

.check_frequency <- function(args, check_collection) {
  if (any(names(args) == "frequency")) {
    checkmate::assert_number(
      x = args$frequency,
      lower = 0.0001,
      add = check_collection,
      .var.name = "frequency",
      null.ok = FALSE,
      na.ok = FALSE
    )
  }
  return(check_collection)
}

.check_temp <- function(args, check_collection) {
  if (any(names(args) == "temp")) {
    # lowest is absolute 0 and upper is plank number
    checkmate::assert_number(
      x = args$temp,
      lower = -273.15,
      upper = 1.416808e+32,
      add = check_collection,
      .var.name = "temp",
      null.ok = FALSE,
      na.ok = FALSE
    )
  }
  return(check_collection)
}

.check_rh <- function(args, check_collection) {
  if (any(names(args) == "rh")) {
    checkmate::assert_number(
      x = args$rh,
      lower = 0,
      upper = 100,
      add = check_collection,
      .var.name = "rh",
      null.ok = FALSE,
      na.ok = FALSE
    )
  }
  return(check_collection)
}

.check_pa <- function(args, check_collection) {
  if (any(names(args) == "pa")) {
    checkmate::assert_number(
      x = args$pa,
      lower = 0,
      add = check_collection,
      .var.name = "pa",
      null.ok = FALSE,
      na.ok = FALSE
    )
  }
  return(check_collection)
}

.check_dist <- function(args, check_collection) {
  if (any(names(args) == "dist")) {
    checkmate::assert_number(
      x = args$dist,
      lower = -273.15,
      add = check_collection,
      .var.name = "dist",
      null.ok = FALSE,
      na.ok = FALSE
    )
  }
  return(check_collection)
}

.check_dist0 <- function(args, check_collection) {
  if (any(names(args) == "dist0")) {
    checkmate::assert_number(
      x = args$dist0,
      lower = -273.15,
      add = check_collection,
      .var.name = "dist0",
      null.ok = FALSE,
      na.ok = FALSE
    )
  }
  return(check_collection)
}

.check_hab.att.coef <- function(args, check_collection) {
  if (any(names(args) == "hab.att.coef")) {
    checkmate::assert_number(
      x = args$hab.att.coef,
      add = check_collection,
      .var.name = "hab.att.coef",
      null.ok = FALSE,
      na.ok = FALSE
    )
  }
  return(check_collection)
}

.check_spl.cutoff <- function(args, check_collection) {
  if (any(names(args) == "spl.cutoff")) {
    checkmate::assert_number(
      x = args$spl.cutoff,
      lower = 0,
      add = check_collection,
      .var.name = "spl.cutoff",
      null.ok = TRUE,
      na.ok = FALSE
    )
  }
  return(check_collection)
}

.check_files <- function(args, check_collection, fun) {
  if (any(names(args) == "files")) {
    if (!is.null(args$path)) {
      files <- file.path(args$path, unique(args$files))
    } else {
      if (!is.null(getOption("sound.files.path"))) {
        files <-
          file.path(getOption("sound.files.path"), unique(args$files))
      } else {
        files <- unique(args$files)
      }
    }
    .assert_sound_files_found(
      files = files,
      fun = fun,
      add = check_collection,
      .var.name = "files"
    )
  }
  return(check_collection)
}

.check_output <- function(args, check_collection) {
  if (any(names(args) == "output")) {
    .assert_deprecated(x = args$output,
                       add = check_collection,
                       .var.name = "output")
  }
  return(check_collection)
}

.check_marker <- function(args, check_collection, fun) {
  if (fun != "manual_realign" & any(names(args) == "marker")) {
    .assert_deprecated(x = args$marker,
                       add = check_collection,
                       .var.name = "marker")
  }
  
  if (fun == "manual_realign" & any(names(args) == "marker")) {
    .assert_deprecated(x = args$marker,
                       add = check_collection,
                       .var.name = "marker")
  }
  return(check_collection)
}

.check_parallel <- function(args, check_collection) {
  if (any(names(args) == "parallel")) {
    .assert_deprecated(x = args$parallel,
                       add = check_collection,
                       .var.name = "parallel")
  }
  return(check_collection)
}

.check_template.rows <- function(args, check_collection) {
  if (any(names(args) == "template.rows")) {
    .assert_deprecated(x = args$template.rows,
                       add = check_collection,
                       .var.name = "template.rows")
  }
  return(check_collection)
}

.check_cores <- function(args, check_collection) {
  if (any(names(args) == "cores")) {
    checkmate::assert_integerish(
      args$cores,
      add = check_collection,
      lower = 1,
      upper = parallel::detectCores(),
      .var.name = "cores",
      any.missing = FALSE,
      all.missing = FALSE,
      null.ok = FALSE
    )
  }
  return(check_collection)
}

.check_pb <- function(args, check_collection) {
  if (any(names(args) == "pb")) {
    checkmate::assert_logical(
      x = args$pb,
      len = 1,
      add = check_collection,
      .var.name = "pb"
    )
  }
  return(check_collection)
}

.check_path <- function(args, check_collection) {
  if (any(names(args) == "path")) {
    checkmate::assert_directory(
      x = args$path,
      access = "r",
      add = check_collection,
      .var.name = "path"
    )
  }
  return(check_collection)
}

.check_n.samples <- function(args, check_collection) {
  if (any(names(args) == "n.samples")) {
    checkmate::assert_number(
      x = args$n.samples,
      lower = 10,
      null.ok = TRUE,
      add = check_collection,
      .var.name = "n.samples"
    )
  }
  return(check_collection)
}

.check_hop.size <- function(args, check_collection) {
  if (any(names(args) == "hop.size")) {
    checkmate::assert_number(
      x = args$hop.size,
      lower = 0.0001,
      add = check_collection,
      .var.name = "hop.size"
    )
  }
  return(check_collection)
}

.check_wl <- function(args, check_collection) {
  if (any(names(args) == "wl")) {
    checkmate::assert_number(
      x = args$wl,
      lower = 2,
      add = check_collection,
      .var.name = "wl",
      na.ok = FALSE
    )
  }
  return(check_collection)
}

.check_bp <- function(args, check_collection) {
  if (any(names(args) == "bp")) {
    if (is.numeric(args$bp)) {
      checkmate::assert_numeric(
        x = args$bp,
        any.missing = FALSE,
        all.missing = FALSE,
        len = 2,
        unique = TRUE,
        lower = 0,
        add = check_collection,
        .var.name = "bp"
      )
    } else {
      checkmate::assert_character(
        x = args$bp,
        null.ok = FALSE,
        add = check_collection,
        .var.name = "bp",
        len = 1,
        pattern = "^freq.range$",
        any.missing = FALSE,
        all.missing = FALSE,
        ignore.case = FALSE
      )
    }
  }
  return(check_collection)
}

.check_img <- function(args, check_collection) {
  if (any(names(args) == "img")) {
    checkmate::assert_logical(
      x = args$img,
      len = 1,
      add = check_collection,
      .var.name = "img"
    )
  }
  return(check_collection)
}

.check_col <- function(args, check_collection) {
  if (any(names(args) == "col")) {
    checkmate::assert_character(
      x = args$col,
      add = check_collection,
      .var.name = "col",
      any.missing = FALSE,
      all.missing = FALSE
    )
  }
  return(check_collection)
}

.check_ovlp <- function(args, check_collection) {
  if (any(names(args) == "ovlp")) {
    checkmate::assert_numeric(
      x = args$ovlp,
      any.missing = FALSE,
      all.missing = FALSE,
      unique = TRUE,
      lower = 0,
      upper = 99.9,
      len = 1,
      add = check_collection,
      .var.name = "ovlp"
    )
  }
  return(check_collection)
}

.check_only.sels <- function(args, check_collection) {
  if (any(names(args) == "only.sels")) {
    checkmate::assert_logical(
      x = args$only.sels,
      len = 1,
      add = check_collection,
      .var.name = "only.sels"
    )
  }
  return(check_collection)
}

.check_method <- function(args, check_collection) {
  if (any(names(args) == "method")) {
    checkmate::assert_choice(
      x = args$method,
      null.ok = FALSE,
      add = check_collection,
      .var.name = "method",
      choices = c(1, 2)
    )
    
    if (args$method == 2) {
      checkmate::assert_names(
        x = names(args$X),
        type = "unique",
        must.include = "transect",
        add = check_collection,
        .var.name = "transect column"
      )
      
      if (!is.null(args$X$transect)) {
        checkmate::assert_vector(
          x = args$X$transect,
          any.missing = FALSE,
          all.missing = FALSE,
          null.ok = FALSE,
          add = check_collection,
          .var.name = "transect values"
        )
      }
    }
  }
  return(check_collection)
}

.check_cor.method <- function(args, check_collection) {
  if (any(names(args) == "cor.method")) {
    checkmate::assert_character(
      x = args$cor.method,
      null.ok = FALSE,
      add = check_collection,
      .var.name = "cor.method",
      len = 1,
      any.missing = FALSE,
      all.missing = FALSE,
      ignore.case = FALSE
    )
    checkmate::assert_choice(
      x = args$cor.method,
      choices = c("pearson", "kendall", "spearman"),
      add = check_collection,
      .var.name = "cor.method"
    )
  }
  return(check_collection)
}

.check_wn <- function(args, check_collection) {
  if (any(names(args) == "wn")) {
    checkmate::assert_character(
      x = args$wn,
      null.ok = FALSE,
      add = check_collection,
      .var.name = "wn",
      len = 1,
      any.missing = FALSE,
      all.missing = FALSE,
      ignore.case = FALSE
    )
    checkmate::assert_choice(
      x = args$wn,
      choices = c(
        "bartlett",
        "blackman",
        "flattop",
        "hamming",
        "hanning",
        "rectangle"
      ),
      add = check_collection,
      .var.name = "wn"
    )
  }
  return(check_collection)
}


.check_Y <- function(args, check_collection, fun) {
  if (any(names(args) == "Y")){
  if (fun == "manual_realign") {
    .assert_est_by_song(x = args$Y,
                        add = check_collection,
                        .var.name = "Y")
    
   }
    }
  return(check_collection)
}
