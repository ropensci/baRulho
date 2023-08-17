# internal function not to be called by users
# stop function that doesn't print call
stop2 <- function(...) {
  stop(..., call. = FALSE)
}

# internal baRulho function, not to be called by users. It prepares X for comparing sounds
# @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})
# last modification on sep-2022 (MAS)

prep_X_bRlo_int <- function(X, method, cores, pb) {
  # set pb options
  on.exit(pbapply::pboptions(type = .Options$pboptions$type), add = TRUE)

  # add sound file selec colums to X (weird column name so it does not overwrite user columns)
  X$.sgnl.temp <- paste(X$sound.files, X$selec, sep = "-")

  # make it a regular data frame (no est)
  X2 <- as.data.frame(X)

  # set clusters for windows OS
  if (Sys.info()[1] == "Windows" & cores > 1) {
    cl <- parallel::makePSOCKcluster(getOption("cl.cores", cores))
  } else {
    cl <- cores
  }

  # set pb options
  pbapply::pboptions(type = ifelse(as.logical(pb), "timer", "none"))

  # add second column with names of the reference sounds to be compared against
  X$reference <-
    pbapply::pbsapply(X2$.sgnl.temp, cl = cl, function(x, meth = method) {
      # extract for single sound id and order by distance
      Y <-
        X2[X2$sound.id == X$sound.id[X2$.sgnl.temp == x], , drop = FALSE]
      Y <- Y[order(Y$distance), ]

      # method 1 compare to closest distance to source
      if (meth == 1) {
        # if column transect if found select the lowest distance in that trasnect
        if (!is.null(X2$transect)) {
          W <- X2[X2$transect == X2$transect[X2$.sgnl.temp == x] & X2$sound.id == X$sound.id[X2$.sgnl.temp == x], , drop = FALSE]

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
        W <- X2[X2$transect == X2$transect[X2$.sgnl.temp == x] & X2$sound.id == X$sound.id[X2$.sgnl.temp == x], , drop = FALSE]
        W <- W[order(W$distance), ]

        # if not the first row then the previous row
        if (W$.sgnl.temp[1] != x) {
          z <- W$.sgnl.temp[which(W$.sgnl.temp == x) - 1]
        } else {
          # else the first row
          z <- x
        }
      }

      return(z)
    })

  return(X)
}


# internal baRulho function, not to be called by users. It is a modified version of warbleR::find_peaks
# that allows to define internally if progress bar would be used (pbapply::pblapply uses pboptions to do this)
# Find cross-correlation peaks

find_peaks_bRlh_int <-
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
            stop2("'xc.output' must be and object of class 'xcorr.output'")
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
          regexpr("\\-[^\\-]*$", peaks$sound.files) - 1
        )

      #### add start and end
      # add template column to selection table in xc.output
      Y <- xc.output$org.selection.table
      Y$template <- paste(Y$sound.files, Y$selec, sep = "-")

      # Y <- Y[Y$template %in% comp_mat[, 1], ]

      # add start as time - half duration of template
      peaks$start <- sapply(seq_len(nrow(peaks)), function(i) {
        peaks$time[i] -
          ((Y$end[Y$template == peaks$template[i]] -
            Y$start[Y$template == peaks$template[i]]) / 2)
      })

      # add end as time + half duration of template
      peaks$end <- sapply(seq_len(nrow(peaks)), function(i) {
        peaks$time[i] +
          ((Y$end[Y$template == peaks$template[i]] -
            Y$start[Y$template == peaks$template[i]]) / 2)
      })

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
      "Araya-Salas, M. (2020), baRulho: quantifying habitat-induced degradation of (animal) acoustic signals in R. R package version 1.0.0"
    )

    invisible(TRUE)
  }

report_assertions2 <- function(collection) {
  checkmate::assertClass(collection, "AssertCollection")
  if (!collection$isEmpty()) {
    msgs <- collection$getMessages()

    # modfied to get more informative message
    msgs <-
      gsub(
        pattern = "Variable 'names(X)'",
        replacement = "Column names in 'X'",
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
        pattern = "Argument 'transect values'",
        replacement = "Column 'transect' in 'X'",
        x = msgs,
        fixed = TRUE
      )


    context <- "\n %i argument check(s) failed:"
    err <- c(sprintf(context, length(msgs)), strwrap(msgs,
      prefix = " * "
    ))
    stop(simpleError(paste0(err, collapse = "\n"), call = sys.call(1L)))
  }
  invisible(TRUE)
}

# custom assert functions to check arguments
# check no duplicated selection labels
check_unique_sels <- function(x, fun) {
  if (anyDuplicated(paste(x$sound.files, x$selec)) > 0) {
    "Duplicated 'selec' labels within at least one sound file"
  } else {
    TRUE
  }
}

assert_unique_sels <-
  checkmate::makeAssertionFunction(check_unique_sels)

# duplicated sound.id within  sound file and distance
check_unique_sound.id <- function(x, fun) {
  if (anyDuplicated(paste0(x$sound.files, x$sound.id)) > 0) {
    "Duplicated 'selec' labels within at least one sound file"
  } else {
    TRUE
  }
}

assert_unique_sels <-
  checkmate::makeAssertionFunction(check_unique_sels)

# check unique sound.id
check_unique_sound.id <- function(x, fun) {
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

assert_unique_sound.id <-
  checkmate::makeAssertionFunction(check_unique_sound.id)


# check if argument has been deprecated
check_deprecated <- function(x) {
  if (!is.null(x)) {
    "has been deprecated"
  } else {
    TRUE
  }
}

assert_deprecated <-
  checkmate::makeAssertionFunction(check_deprecated)


## function to check arguments
check_arguments <- function(fun, args) {
  # make function name a character
  fun <- as.character(fun)

  # create object to store check results
  check_collection <- checkmate::makeAssertCollection()

  ### check arguments
  if (any(names(args) == "X")) {
    if (fun != "noise_profile") {
      checkmate::assert_data_frame(
        x = args$X,
        any.missing = TRUE,
        min.rows = if (fun == "master_sound_file") {
          2
        } else {
          NULL
        },
        add = check_collection,
        .var.name = "X"
      )

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

      if (fun != "master_sound_file") {
        # functions that compare by distance
        if (fun %in% c(
          "blur_ratio",
          "detection_distance",
          "envelope_correlation",
          "excess_attenuation",
          "spcc",
          "spectrum_blur_ratio",
          "spectrum_correlation"
        )) {
          cols <-
            c(
              "sound.files",
              "selec",
              "start",
              "end",
              "sound.id",
              "distance"
            )

          checkmate::assert_numeric(
            x = args$X$distance,
            any.missing = FALSE,
            all.missing = FALSE,
            lower = 0,
            add = check_collection,
            .var.name = "X$distance"
          )
        } else {
          cols <-
            c("sound.files", "selec", "start", "end", "sound.id")
        }
      } else {
        cols <- c("sound.files", "selec", "start", "end")
      }

      checkmate::assert_names(
        x = names(args$X),
        type = "unique",
        must.include = cols,
        add = check_collection,
        .var.name = "names(X)"
      )
      try(
        checkmate::assert_data_frame(
          x = args$X[, cols],
          any.missing = TRUE,
          add = check_collection,
          .var.name = "X"
        ),
        silent = TRUE
      )

      assert_unique_sels(
        x = args$X,
        fun = fun,
        add = check_collection,
        .var.name = "X"
      )

      assert_unique_sound.id(
        x = args$X,
        fun = fun,
        add = check_collection,
        .var.name = "X"
      )

      if (!fun %in% c(
        "plot_align_sounds",
        "degrad_catalog",
        "detection_distance",
        "signal_to_noise_ratio",
        "tail_to_signal_ratio"
      )) {
        if (is_extended_selection_table(args$X)) {
          if (length(unique(attr(args$X, "check.results")$sample.rate)) > 1) {
            stop2(
              "all wave objects in the extended selection table must have the same sampling rate (they can be homogenized using warbleR::resample_est())"
            )
          }
        } else {
          if (!fun %in% c("find_markers", "align_test_files")) {
            warning2("assuming all sound files have the same sampling rate")
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
    }
  }

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

  if (any(names(args) == "shuffle")) {
    checkmate::assert_logical(
      x = args$shuffle,
      len = 1,
      null.ok = FALSE,
      add = check_collection,
      .var.name = "shuffle"
    )
  }

  if (any(names(args) == "seed")) {
    checkmate::assert_number(
      x = args$seed,
      add = check_collection,
      .var.name = "seed",
      null.ok = TRUE
    )
  }


  if (any(names(args) == "srt")) {
    checkmate::assert_number(
      x = args$srt,
      add = check_collection,
      .var.name = "srt",
      null.ok = FALSE
    )
  }

  if (any(names(args) == "sig2")) {
    checkmate::assert_number(
      x = args$sig2,
      add = check_collection,
      lower = 0.00001,
      .var.name = "sig2",
      null.ok = TRUE
    )
  }

  if (any(names(args) == "fm")) {
    checkmate::assert_logical(
      x = args$fm,
      len = 1,
      null.ok = FALSE,
      add = check_collection,
      .var.name = "fm"
    )
  }

  if (any(names(args) == "am")) {
    checkmate::assert_logical(
      x = args$am,
      len = 1,
      null.ok = FALSE,
      add = check_collection,
      .var.name = "am"
    )
  }

  if (any(names(args) == "frequencies")) {
    checkmate::assert_numeric(
      x = args$frequencies,
      any.missing = FALSE,
      all.missing = FALSE,
      unique = TRUE,
      lower = 0001,
      add = check_collection,
      .var.name = "frequencies"
    )
  }

  if (any(names(args) == "durations")) {
    checkmate::assert_numeric(
      x = args$durations,
      any.missing = FALSE,
      all.missing = FALSE,
      unique = TRUE,
      lower = 0001,
      add = check_collection,
      .var.name = "durations"
    )
  }

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

  if (any(names(args) == "label")) {
    checkmate::assert_logical(
      x = args$label,
      len = 1,
      null.ok = FALSE,
      add = check_collection,
      .var.name = "label"
    )
  }

  if (any(names(args) == "fast.spec")) {
    checkmate::assert_logical(
      x = args$fast.spec,
      len = 1,
      null.ok = FALSE,
      add = check_collection,
      .var.name = "fast.spec"
    )
  }

  if (any(names(args) == "width")) {
    checkmate::assert_number(
      x = args$width,
      lower = 0.1,
      add = check_collection,
      .var.name = "width",
      null.ok = FALSE,
      na.ok = FALSE
    )
  }

  if (any(names(args) == "height")) {
    checkmate::assert_number(
      x = args$height,
      lower = 0.1,
      add = check_collection,
      .var.name = "height",
      null.ok = FALSE,
      na.ok = FALSE
    )
  }

  if (any(names(args) == "ssmooth")) {
    checkmate::assert_number(
      x = args$ssmooth,
      lower = 1,
      add = check_collection,
      .var.name = "ssmooth",
      null.ok = FALSE,
      na.ok = FALSE
    )
  }

  if (any(names(args) == "msmooth")) {
    checkmate::assert_numeric(
      x = args$msmooth,
      lower = 0,
      finite = TRUE,
      any.missing = FALSE,
      all.missing = FALSE,
      len = 2,
      add = check_collection,
      .var.name = "msmooth"
    )
  }

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

  if (any(names(args) == "dest.path")) {
    checkmate::assert_directory(
      x = args$dest.path,
      access = "r",
      add = check_collection,
      .var.name = "dest.path"
    )
  }

  if (any(names(args) == "flim")) {
    checkmate::assert_numeric(
      x = args$flim,
      lower = 0,
      finite = TRUE,
      any.missing = FALSE,
      all.missing = FALSE,
      len = 2,
      add = check_collection,
      .var.name = "flim"
    )
  }

  if (any(names(args) == "mar")) {
    checkmate::assert_number(
      x = args$mar,
      lower = 0.00001,
      add = check_collection,
      .var.name = "mar",
      null.ok = FALSE,
      na.ok = FALSE
    )
  }

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

  if (any(names(args) == "palette")) {
    checkmate::assert_function(
      x = args$palette,
      null.ok = FALSE,
      add = check_collection,
      .var.name = "palette"
    )
  }

  if (any(names(args) == "PSD")) {
    checkmate::assert_logical(
      x = args$PSD,
      len = 1,
      null.ok = FALSE,
      add = check_collection,
      .var.name = "PSD"
    )
  }

  if (any(names(args) == "norm")) {
    checkmate::assert_logical(
      x = args$norm,
      len = 1,
      null.ok = FALSE,
      add = check_collection,
      .var.name = "norm"
    )
  }

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

  if (any(names(args) == "averaged")) {
    checkmate::assert_logical(
      x = args$averaged,
      len = 1,
      null.ok = FALSE,
      add = check_collection,
      .var.name = "averaged"
    )
  }


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

  if (any(names(args) == "hab.att.coef")) {
    checkmate::assert_number(
      x = args$hab.att.coef,
      add = check_collection,
      .var.name = "hab.att.coef",
      null.ok = FALSE,
      na.ok = FALSE
    )
  }

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
  # if (any(names(args) == "files"))
  #   if (!is.null(args$files))
  #   checkmate::assert_file_exists(x = args$files, access = "r", extension = c("mp3", "wav", "wac", "flac"), add = check_collection, .var.name = "files")

  if (any(names(args) == "output")) {
    assert_deprecated(
      x = args$output,
      add = check_collection,
      .var.name = "output"
    )
  }

  if (any(names(args) == "marker")) {
    assert_deprecated(
      x = args$marker,
      add = check_collection,
      .var.name = "marker"
    )
  }

  if (any(names(args) == "parallel")) {
    assert_deprecated(
      x = args$parallel,
      add = check_collection,
      .var.name = "parallel"
    )
  }

  if (any(names(args) == "template.rows")) {
    assert_deprecated(
      x = args$template.rows,
      add = check_collection,
      .var.name = "template.rows"
    )
  }

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

  if (any(names(args) == "pb")) {
    checkmate::assert_logical(
      x = args$pb,
      len = 1,
      add = check_collection,
      .var.name = "pb"
    )
  }

  if (any(names(args) == "path")) {
    checkmate::assert_directory(
      x = args$path,
      access = "r",
      add = check_collection,
      .var.name = "path"
    )
  }

  if (any(names(args) == "hop.size")) {
    checkmate::assert_number(
      x = args$hop.size,
      lower = 0.0001,
      add = check_collection,
      .var.name = "hop.size"
    )
  }

  if (any(names(args) == "wl")) {
    checkmate::assert_number(
      x = args$wl,
      lower = 2,
      add = check_collection,
      .var.name = "wl",
      na.ok = FALSE
    )
  }

  if (any(names(args) == "bp")) {
    checkmate::assert_numeric(
      x = args$bp,
      any.missing = FALSE,
      all.missing = FALSE,
      len = 2,
      unique = TRUE,
      lower = 0,
      add = check_collection,
      .var.name = "bp",
      na.ok = FALSE
    )
  }

  if (any(names(args) == "img")) {
    checkmate::assert_logical(
      x = args$img,
      len = 1,
      add = check_collection,
      .var.name = "img"
    )
  }

  if (any(names(args) == "col")) {
    checkmate::assert_character(
      x = args$col,
      add = check_collection,
      .var.name = "col",
      any.missing = FALSE,
      all.missing = FALSE
    )
  }


  if (any(names(args) == "palette")) {
    checkmate::assert_function(
      x = args$palette,
      null.ok = FALSE,
      add = check_collection,
      .var.name = "palette"
    )
  }

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


  if (any(names(args) == "only.sels")) {
    checkmate::assert_logical(
      x = args$only.sels,
      len = 1,
      add = check_collection,
      .var.name = "only.sels"
    )
  }

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
warning2 <- function(...) {
  warning(..., call. = FALSE)
}

# message function that changes colors
message2 <- function(x, color = "black") {
  message(colortext(x, as = color))
}


colortext <-
  function(text,
           as = c(
             "red",
             "blue",
             "green",
             "magenta",
             "cyan",
             "orange",
             "black",
             "silver"
           )) {
    if (has_color()) {
      unclass(cli::make_ansi_style(baRulho_style(as))(text))
    } else {
      text
    }
  }

has_color <- function() {
  cli::num_ansi_colors() > 1
}

baRulho_style <-
  function(color = c(
             "red",
             "blue",
             "green",
             "magenta",
             "cyan",
             "orange",
             "black",
             "silver"
           )) {
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
detection_distance_ohn_int <-
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
# y and z are the sound.files+selec names of the sounds and reference sound (model)
# envelope mismatch ratio
blur_FUN <-
  function(X,
           envs,
           path,
           img,
           dest.path,
           x,
           res,
           ovlp,
           wl,
           collevels,
           palette,
           ...) {
    # get names of sound and reference
    sgnl <- X$.sgnl.temp[x]
    rfrnc <- X$reference[x]

    # if sounds are the same or the selection is noise return NA
    if (sgnl == rfrnc |
      any(c(X$sound.id[X$.sgnl.temp == sgnl], X$sound.id[X$reference == rfrnc]) == "ambient")) {
      out <- NA
    } else {
      # extract envelope for sound and model
      sgnl.env <- envs[[which(names(envs) == sgnl)]]
      rfrnc.env <- envs[[which(names(envs) == rfrnc)]]

      # make them the same length as the shortest one
      if (length(sgnl.env) > length(rfrnc.env)) {
        sgnl.env <- sgnl.env[1:length(rfrnc.env)]
      }
      if (length(rfrnc.env) > length(sgnl.env)) {
        rfrnc.env <- rfrnc.env[1:length(sgnl.env)]
      }

      # duration (any works as they all must have the same sampling rate)
      # dur <- length(sgnl.env) / (attr(X, "check.results")$sample.rate[1] * 1000)
      dur <-
        length(sgnl.env) / warbleR::read_sound_file(
          X = X,
          index = x,
          path = path,
          header = TRUE
        )$sample.rate

      # convert envelopes to PMF (probability mass function)
      rfrnc.pmf <- rfrnc.env / sum(rfrnc.env)
      sgn.pmf <- sgnl.env / sum(sgnl.env)

      # get blur ratio as half the sum of absolute differences between envelope PMFs
      bl.rt <- sum(abs(rfrnc.pmf - sgn.pmf)) / 2

      # plot
      if (img) {
        warbleR:::img_wrlbr_int(
          filename = paste0(
            "blur_ratio_",
            X$sound.id[x],
            "-",
            rfrnc,
            "-",
            sgnl,
            ".jpeg"
          ),
          path = dest.path,
          width = 10.16 * 1.5,
          height = 10.16,
          units = "cm",
          res = res
        )

        # time values for plots
        time.vals <- seq(0, dur, length.out = length(sgnl.env))

        # difference between envelopes
        env.diff <- rfrnc.pmf - sgn.pmf

        # matrix for layout
        ly.mat <- matrix(
          c(
            0, 0.3, 0, 0.5, # bottom left spectrogram
            0, 0.3, 0.5, 1, # top left spectrogram
            0.2, 1, 0, 1
          ),
          # right pannel envelopes
          nrow = 3,
          byrow = TRUE
        )

        # save par settings
        oldpar <- par(no.readonly = TRUE)
        on.exit(par(oldpar))

        # close if open any screen
        invisible(close.screen(all.screens = TRUE))

        # split screen
        split.screen(ly.mat)

        # plot envelopes
        screen(3)

        # set image margins
        par(mar = c(4, 4, 4, 3))

        # reference envelope first
        plot(
          time.vals,
          rfrnc.pmf,
          type = "l",
          xlab = "",
          ylab = "",
          col = "#31688E",
          ylim = c(min(rfrnc.pmf, sgn.pmf), max(rfrnc.pmf, sgn.pmf) * 1.1),
          cex.main = 0.8,
          lwd = 1.2,
          yaxt = "n"
        )

        # add x axis label
        mtext(
          text = "Time (s)",
          side = 1,
          line = 2.5
        )

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
          col = "#31688E",
          cex = 1
        )
        mtext(
          text = paste("Test sound:", sgnl),
          side = 3,
          line = 0.5,
          col = "#B4DE2C",
          cex = 1
        )

        # add y axis
        axis(side = 4, labels = FALSE)
        mtext(
          text = "Amplitude (PMF)",
          side = 4,
          line = 1.5
        )

        # add sound envelope
        lines(time.vals,
          sgn.pmf,
          col = "#B4DE2CFF",
          lwd = 1.2
        )

        # sound envelope on top
        polygon(
          x = c(time.vals, rev(time.vals)),
          y = c(sgn.pmf, rev(rfrnc.pmf)),
          col = "#FDE72533",
          border = NA
        )

        # get plotting area limits
        usr <- par("usr")

        # and blu ratio value
        text(
          x = ((usr[1] + usr[2]) / 2) + usr[1],
          y = usr[4] * 0.9,
          paste("Blur ratio:", round(bl.rt, 2)),
          cex = 1
        )

        # index of reference
        rf.indx <-
          which(paste(X$sound.files, X$selec, sep = "-") == rfrnc)

        # freq limit of reference
        flim <- c(X$bottom.freq[rf.indx], X$top.freq[rf.indx])

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

        ## plot spectros
        # sound at bottom left
        screen(1)
        par(mar = c(0.3, 0.3, 0.15, 0.3))

        warbleR:::spectro_wrblr_int2(
          wave = clp.sgnl,
          f = clp.sgnl@samp.rate,
          flim = flim,
          axisX = FALSE,
          axisY = FALSE,
          tlab = NULL,
          flab = NULL,
          main = NULL,
          grid = FALSE,
          rm.zero = TRUE,
          cexaxis = 1.2,
          add = TRUE,
          ovlp = ovlp,
          wl = wl,
          collevels = collevels,
          palette = palette
        )

        # lines showing position of sound
        abline(
          v = c(mar.rf.bf, X$end[x] - X$start[x] + mar.rf.bf),
          col = "#B4DE2CFF",
          lty = 2
        )

        # add box with sound color
        box(col = "#B4DE2CFF", lwd = 3)

        # reference at top left
        screen(2)
        par(mar = c(0.15, 0.3, 0.3, 0.3))

        warbleR:::spectro_wrblr_int2(
          wave = clp.rfnc,
          f = clp.rfnc@samp.rate,
          flim = flim,
          axisX = FALSE,
          axisY = FALSE,
          tlab = NULL,
          flab = NULL,
          main = NULL,
          grid = FALSE,
          rm.zero = TRUE,
          cexaxis = 1.2,
          add = TRUE,
          ovlp = ovlp,
          wl = wl,
          collevels = collevels,
          palette = palette
        )

        # lines showing position of sound
        abline(
          v = c(mar.rf.bf, X$end[rf.indx] - X$start[rf.indx] + mar.rf.bf),
          col = "#31688ECC",
          lty = 2
        )

        # add box with reference color
        box(col = "#31688E", lwd = 3)

        # close graph
        dev.off()
      }

      # return maximum correlation
      return(bl.rt)
    }
    return(out)
  }

# function to measure envelope correlation
# y and z are the sound.files+selec names of the sounds and reference sound (model)
env_cor_FUN <- function(Y, y, z, envs, cor.method) {
  # if names are the same return NA
  if (y == z | Y$sound.id[Y$.sgnl.temp == y] == "ambient") {
    out <- NA
  } else {
    # extract envelope for sound and model
    sgnl.env <- envs[[which(names(envs) == y)]]
    mdl.env <- envs[[which(names(envs) == z)]]

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
    cors <- sapply(1:stps, function(x) {
      cor(lg.env[x:(x + shrt.lgth)], shrt.env, method = cor.method)
    })

    # return maximum correlation
    out <- max(cors, na.rm = TRUE)
  }

  return(out)
}

# function to extract mean envelopes
meanenv_FUN <- function(y, wl, ovlp, X, path, bp) {
  if (X$sound.id[y] == "ambient") {
    sig_env <- NA
  } else {
    # read sound clip
    clp <-
      warbleR::read_sound_file(
        X = X,
        index = y,
        from = 0,
        to = X$end[y],
        path = path
      )

    if (X$sound.id[y] != "ambient") {
      noise_clp <-
        warbleR::read_sound_file(
          X = X,
          index = y,
          from = 0,
          to = X$start[y] - 0.001,
          path = path
        )
    }

    # add band-pass frequency filter
    if (!is.null(bp)) {
      # filter to bottom and top freq range
      if (bp == "freq.range") {
        bp <- c(X$bottom.freq[y], X$top.freq[y])
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

      if (X$sound.id[y] != "ambient") {
        noise_clp <-
          seewave::ffilter(
            noise_clp,
            f = noise_clp@samp.rate,
            from = bp[1] * 1000,
            ovlp = ovlp,
            to = bp[2] * 1000,
            bandpass = TRUE,
            wl = wl,
            output = "Wave"
          )
      }
    }

    # get mean envelopes
    sig_env <-
      mean(seewave::env(
        clp,
        f = clp@samp.rate,
        envt = "hil",
        plot = FALSE
      ))
  }

  return(data.frame((X[y, , drop = FALSE]), sig_env))
}

# function to measure spectrum correlation
# y and z are the sound.files+selec names of the sounds and reference sound (model)
spctr_cor_FUN <- function(y, z, spcs, X, cor.method) {
  # if names are the same return NA
  if (y == z | X$sound.id[X$.sgnl.temp == y] == "ambient") {
    cor.spctr <- NA
  } else {
    # extract envelope for sound and model
    sgnl.spctr <- spcs[[which(names(spcs) == y)]]
    mdl.spctr <- spcs[[which(names(spcs) == z)]]


    ### filter to freq range of sounds and remove freq column
    # get range as lowest bottom and highest top
    frng <-
      c(min(X$bottom.freq[X$.sgnl.temp %in% c(y, z)]), max(X$top.freq[X$.sgnl.temp %in% c(y, z)]))
    sgnl.spctr <-
      sgnl.spctr[sgnl.spctr[, 1] > frng[1] &
        sgnl.spctr[, 1] < frng[2], 2]
    mdl.spctr <-
      mdl.spctr[mdl.spctr[, 1] > frng[1] &
        mdl.spctr[, 1] < frng[2], 2]

    # get correlation assuming they have same length
    cor.spctr <- cor(sgnl.spctr, mdl.spctr, method = cor.method)
  }

  return(cor.spctr)
}

exc_att_FUN <- function(y, X, meth, tp, gn) {
  # print(y)
  # get names of sound and reference
  sgnl <- X$.sgnl.temp[y]
  rfrnc <- X$reference[y]

  if (X$sound.id[X$.sgnl.temp == sgnl] == "ambient" | sgnl == rfrnc) {
    ea <- NA
  } else {
    # method 1 compare to closest distance to source
    # if (meth == 1) {
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
rbind2 <- function(...) {
  suppressWarnings(rbind(...))
}

# ## internal function to subtract SPL from background noise
# # sound = sound SPL
# # noise = noise SPL
# lessdB <- function(sound.noise, noise){
#
#   puttative_SPLs <- seq(0.01, sound.noise, by = 0.01)
#
#   sum_SPLs <-  20 * log10((10^(puttative_SPLs/20)) + (10^(noise/20)))
#
#   sound_SPL <- puttative_SPLs[which.min(abs(sum_SPLs - sound.noise))]
#
#   return(sound_SPL)
#   }
