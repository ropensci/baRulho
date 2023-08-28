# internal function not to be called by users
# stop function that doesn't print call
stop2 <- function(...) {
  stop(..., call. = FALSE)
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
        pattern = "Variable 'X$distances':",
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
        replacement = "but is missing {'reference'} (Did you forget to run 'set_reference_sounds() first?)",
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

# duplicated sound.id within sound file and distance
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

# check than more than 1 distance is found
check_several_distances <- function(x, fun) {
  if (!is.null(x$distance)) {
    if (length(unique(x$distance)) == 1) {
      "Column 'distance'in 'X' must include more than 1 distance"
    }
  } else {
    TRUE
  }
}

assert_several_distances <-
  checkmate::makeAssertionFunction(check_several_distances)

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
        min.rows = if (!fun %in% c("tail_to_signal_ratio", "signal_to_noise_ratio")) {
          2
        } else {
          1
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
          c(
            "sound.files",
            "selec",
            "start",
            "end",
            "sound.id",
            "distance",
            "reference"
          )

        checkmate::assert_numeric(
          x = args$X$distance,
          any.missing = FALSE,
          all.missing = FALSE,
          lower = 0,
          null.ok = TRUE,
          add = check_collection,
          .var.name = "X$distance"
        )

        assert_several_distances(
          x = args$X,
          fun = fun,
          add = check_collection,
          .var.name = "X$distances"
        )
      }


      if (fun == "set_reference_sounds") {
        cols <- c(
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
        disjunct.from = if (fun == "set_reference_sounds") "reference" else NULL,
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
        "plot_aligned_sounds",
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
          if (!fun %in% c("find_markers", "align_test_files", "set_reference_sounds")) {
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

spctr_FUN <- function(y, spec.smooth, wl, X, path, meanspc = FALSE, ovlp) {
  # load clip
  clp <- warbleR::read_sound_file(
    X = X,
    index = which(X$.sgnl.temp == y),
    path = path
  )

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
    warbleR::envelope(
      x = clp.spc[, 2],
      ssmooth = spec.smooth
    )

  return(clp.spc)
}

## function to measure blur ratio
blur_FUN <-
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
        sgnl.env <- sgnl.env[1:length(rfrnc.env)]
      }
      if (length(rfrnc.env) > length(sgnl.env)) {
        rfrnc.env <- rfrnc.env[1:length(sgnl.env)]
      }

      # duration (any sampling rate works as they all must have the same sampling rate)
      # dur <- length(sgnl.env) / (attr(X, "check.results")$sample.rate[1] * 1000)
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
blur_sp_FUN <-
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

      # make them the same number of rows
      sgnl.spc <- sgnl.spc[1:(min(c(nrow(sgnl.spc), nrow(rfrnc.spc)))), ]
      rfrnc.spc <- rfrnc.spc[1:(min(c(nrow(sgnl.spc), nrow(rfrnc.spc)))), ]

      # make test the same frequency range as reference
      bp <- c(X$bottom.freq[X$.sgnl.temp == rfrnc], X$top.freq[X$.sgnl.temp == rfrnc])

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
    }
    return(sp.bl.rt)
  }

## function to plot  blur ratios
plot_blur_FUN <-
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
      sgnl.energy <- energy_vectors[[which(names(energy_vectors) == sgnl)]]
      rfrnc.energy <- energy_vectors[[which(names(energy_vectors) == rfrnc)]]

      # make them the same length as the shortest one
      if (length(sgnl.energy) > length(rfrnc.energy)) {
        sgnl.energy <- sgnl.energy[1:length(rfrnc.energy)]
      }
      if (length(rfrnc.energy) > length(sgnl.energy)) {
        rfrnc.energy <- rfrnc.energy[1:length(sgnl.energy)]
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
        bp <- c(X$bottom.freq[X$.sgnl.temp == rfrnc], X$top.freq[X$.sgnl.temp == rfrnc])

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

      img_name <- paste0(
        if (spectr) "spectrum_blur_ratio_" else "blur_ratio_",
        X$sound.id[x],
        "-",
        rfrnc,
        "-",
        sgnl,
        ".jpeg"
      )

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
          0.06, 0.4, 0, 0.562, # bottom left spectrogram
          0.06, 0.4, 0.562, 1, # top left spectrogram
          0.4, 1, 0, 1,  # right pannel with blur ratio
          0, 0.06, 0.1, 1
        ),
        nrow = 4,
        byrow = TRUE
      )

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
          c(
            flm[1] + as.numeric(flim[1]),
            flm[2] + as.numeric(flim[2])
          )
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
             c(at_freq[-length(at_freq)], "")
      )
      
      # plot time ticks
      at_time <-
        pretty(seq(0, duration(clp.sgnl), length.out = 10)[-10], n = 4)
      axis(1,
           at = at_time,
           labels =c(at_time[c(-length(at_time))], "")
      )
      
      # add x axis label
      mtext(
        text = "Time (s)",
        side = 1,
        line = 2
      )

      
      
      # lines showing position of sound
      abline(
        v = c(mar.rf.bf, X$end[x] - X$start[x] + mar.rf.bf),
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
        v = c(mar.rf.bf, X$end[rf.indx] - X$start[rf.indx] + mar.rf.bf),
        col = "white",
        lty = 2
      )

      at_freq <-
        pretty(seq(0, clp.rfnc@samp.rate / 2000, length.out = 10)[-10], n = 10)
      axis(2,
           at = at_freq,
           labels = 
             c(at_freq[-length(at_freq)], "")
      )
      
      # plot envelopes
      screen(3)

      # set image margins
      par(mar = c(4, 1, 4, if (spectr) 4.2 else 3.2))

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
          xaxs="i",
          yaxs="i"
        )

        # add background color
        rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col =adjustcolor("#DEF5E5FF", 0.4), border = NA)
        
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
        mtext(
          text = "Amplitude (PMF)",
          side = 4,
          line = 1
        )

        
        # sound envelope on top
        polygon(
          x = c(time.vals, rev(time.vals)),
          y = c(sgn.pmf, rev(rfrnc.pmf)),
          col = adjustcolor(blur_col, 0.2),
          border = NA
        )
        
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
              lwd = 1.4
        )
        # add sound envelope
        lines(time.vals,
              rfrnc.pmf,
              col = adjustcolor(ref_col, 0.7),
              lwd = 1.4
        )
        
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
          xlim = c(
            min(rfrnc.pmf, sgn.pmf),
            max(rfrnc.pmf, sgn.pmf) * 1.1
          ),
          cex.main = 0.8,
          lwd = 1.2,
          yaxt = "n",
          xaxs="i",
          yaxs="i"
        )
        
        # add background color
        rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col ="#D3D3D380", border = NA)
        
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
        mtext(
          text = "Power spectrum (PMF)",
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
        mtext(
          text = "Frequency (kHz)",
          side = 4,
          line = 1.8
        )

        # add sound spectrum
        lines(sgn.pmf,
          f.vals,
          col = test_col,
          lwd = 1.2
        )

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
              lwd = 1.4
        )
        # add sound envelope
        lines(rfrnc.pmf,
              f.vals,
              col = adjustcolor(ref_col, 0.7),
              lwd = 1.4
        )
        
        
        # get plotting area limits
        usr <- par("usr")

        # and blu ratio value
        text(
          x = ((usr[1] + usr[2]) / 2),
          y = (usr[4] - usr[3]) * 0.9 + usr[3],
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
env_FUN <- function(X, y, env.smooth, ovlp, wl, path) {
  # load clip
  clp <- warbleR::read_sound_file(
    X = X,
    index = which(X$.sgnl.temp == y),
    path = path
  )

  # define bandpass
  bp <- c(X$bottom.freq[X$.sgnl.temp == y], X$top.freq[X$.sgnl.temp == y])

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
    warbleR::envelope(
      x = clp@left,
      ssmooth = env.smooth
    )

  return(nv)
}

# function to measure envelope correlation
# y and z are the sound.files+selec names of the sounds and reference sound (model)
env_cor_FUN <- function(X, x, envs, cor.method) {
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
    cors <- sapply(1:stps, function(x) {
      cor(lg.env[x:(x + shrt.lgth)], shrt.env, method = cor.method)
    })

    # return maximum correlation
    envcor <- max(cors, na.rm = TRUE)
  }
  return(envcor)
}

# function to extract mean envelopes
meanenv_FUN <- function(y, wl, ovlp, X, path, bp) {
  # read sound clip
  clp <-
    warbleR::read_sound_file(
      X = X,
      index = which(X$.sgnl.temp == y),
      from = 0,
      to = X$end[X$.sgnl.temp == y],
      path = path
    )

  noise_clp <-
    warbleR::read_sound_file(
      X = X,
      index = which(X$.sgnl.temp == y),
      from = 0,
      to = X$start[X$.sgnl.temp == y] - 0.001,
      path = path
    )

  # add band-pass frequency filter
  if (!is.null(bp)) {
    # filter to bottom and top freq range
    if (bp[1] == "freq.range") {
      bp <- c(X$bottom.freq[X$.sgnl.temp == y], X$top.freq[X$.sgnl.temp == y])
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

  # get mean envelopes
  sig_env <-
    mean(seewave::env(
      clp,
      f = clp@samp.rate,
      envt = "hil",
      plot = FALSE
    ))

  return(sig_env)
}

# function to measure spectrum correlation
# y and z are the sound.files+selec names of the sounds and reference sound (model)
spctr_cor_FUN <- function(y, specs, X, cor.method) {
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
  }

  return(cor.spctr)
}

exc_att_FUN <- function(y, X, tp, gn) {
  
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
rbind2 <- function(...) {
  suppressWarnings(rbind(...))
}
