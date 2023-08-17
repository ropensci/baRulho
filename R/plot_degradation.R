#' #' Measure spectrographic cross-correlation as a measure of sound distortion
#' #'
#' #' \code{plot_degradation} measures spectrographic cross-correlation as a measure of sound distortion in sounds referenced in an extended selection table.
#' #' @usage plot_degradation(X, parallel = NULL, cores = getOption("mc.cores", 1), pb = getOption("pb", TRUE),
#' #' method = getOption("method", 1), cor.method = "pearson", output = "est",
#' #' hop.size = getOption("hop.size", 11.6), wl = getOption("wl", NULL),
#' #' ovlp = getOption("ovlp", 90), wn = 'hanning', path = getOption("sound.files.path", "."))
#' #' @param X Object of class 'data.frame', 'selection_table' or 'extended_selection_table' (the last 2 classes are created by the function \code{\link[warbleR]{selection_table}} from the warbleR package) with the reference to the sounds in the master sound file. Must contain the following columns: 1) "sound.files": name of the .wav files, 2) "selec": unique selection identifier (within a sound file), 3) "start": start time and 4) "end": end time of selections, 5)  "bottom.freq": low frequency for bandpass, 6) "top.freq": high frequency for bandpass and 7) "sound.id": ID of sounds used to identify counterparts across distances. Each sound must have a unique ID within a distance.
#'
#' #'
#' #' @return One ore more image files with a multipannel figure. Each row includes  all the copies of a sound id for a given transect, also including its reference if it comes from another transect, ordered in increasing 
#' #' with the spectrogram cross-correlation coefficients.
#' #' @export
#' #' @name plot_degradation
#' #' @details
#' #' @examples
#' #' {
#' #'  #'   # load example data
#'   data("degradation_est")
#'
#'   # create subset of data with only re-recorded files
#'   rerecorded_est <- degradation_est[degradation_est$sound.files != "master.wav", ]
#' #'
#' #'   # method 1
#' #'   plot_degradation(X = rerecorded_est, method = getOption("method", 1))
#' #'
#' #'   # method 2
#' #'   plot_degradation(X = rerecorded_est, method = 2)
#' #' }
#' #'
#' #' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})
#' #' @seealso \code{\link{blur_ratio}}, \code{\link{realign_test_sounds}}, \code{\link[warbleR]{cross_correlation}}
#' #' @references {
#' #' Araya-Salas, M. (2020). baRulho: baRulho: quantifying habitat-induced degradation of (animal) acoustic signals in R. R package version 1.0.2
#' #' }
#'
  plot_degradation <- function(X, nrow = 4, mar = 0.05, hop.size = getOption("hop.size", 1), wl = getOption("wl", NULL), ovlp = getOption("ovlp", 0), path = getOption("sound.files.path", "."),  dest.path = getOption("dest.path", "."), cores = getOption("mc.cores", 1), pb = getOption("pb", TRUE), collevels = seq(-120, 0, 5), pal = viridis, box = FALSE, ...) {
 
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
  if (!(wl %% 2) == 0)
    wl <- wl + 1
  
  # set clusters for windows OS
  if (Sys.info()[1] == "Windows" & cores > 1) {
    cl <- parallel::makePSOCKcluster(getOption("cl.cores", cores))
  } else {
    cl <- cores
  }
  
  X <- X[X$sound.id != "ambient", ]
  
  # stop if more than 1 sample per distance but no transect info
  if (any(table(X$sound.id, X$distance)) > 1 & is.null(X$transect)) {
    stop2("There are more than 1 test sound per sound.id/distance combination but no 'transect' column to group by transect")
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
    transects <-unique(X$transect)
  }
  
  # add sound file selec column and names to X (weird column name so it does not overwrite user columns)
  if (pb) {
    write(file = "",
          x = paste0("Preparing data for analysis (step 1 out of 3):"))
  }
  
  X <- prep_X_bRlo_int(X,
                       method = 1,
                       cores = cores,
                       pb = pb)
  
  X2 <- as.data.frame(X)
  
  # create data subsets (element in list) with all the copies of a sound id in transect, also including its reference if it comes from another transect
  X2$sound.id.transect <- paste(X2$sound.id, X$transect, sep = "-")
  
  soundid_X_list <- lapply(unique(X2$sound.id.transect), function(x) {
    Y <- X2[X2$.sgnl.temp %in% unique(X2$reference[X2$sound.id.transect == x]) | X2$sound.id.transect == x, ]
  Y <- Y[order(Y$distance), ]
  return(Y)
  })

  # add numeric id to each subset 
  soundid_X_list <- lapply(seq_along(soundid_X_list), function(x) {
    soundid_X_list[[x]]$seq_number <- x 
    return(soundid_X_list[[x]])
  })
  
  # get unique distances and maximum number of columns (distances) for any sound id
  distances <- sort(unique(X$distance))
  ncol <- length(distances)

  # put all into a single data frame
  soundid_X <- do.call(rbind, soundid_X_list)
    
  # add page in which will be printed each subset
  soundid_X$page <- as.numeric(cut(soundid_X$seq_number, breaks = nrow * 0:max(soundid_X$seq_number)))
  
  # set basic layout for all pages
  # spectrogram pannels
  # vertical lines
  verticals <- seq(from = 0.06, to = 0.95, length.out = ncol + 1)
  
  # left margin
  spec_lf <- rep(verticals[-length(verticals)], times = nrow)
  
  # right margin
  spec_rg <- rep(verticals[-1], times = nrow)
  
  # horizontal lines
  horizontals <- seq(from = 0.95, to = 0.075, length.out = nrow + 1)
  
  # bottom margin
  spec_btm <- rep(horizontals[-1], each = ncol)
  
  # top margin
  spec_tp <- rep(horizontals[-length(horizontals)], each = ncol)
  
  # put together in a matrix
  spec_m <- cbind(spec_lf, spec_rg, spec_btm, spec_tp)
  
  ## label panels
  # left margin
  lab_lf <- c(rep(max(verticals), each = nrow), verticals[-length(verticals)], 0, 0)
  
  # right margin
  lab_rg <- c(rep(1, each = nrow), verticals[-1], min(verticals), 1)
  
  # bottom margin
  lab_btm <- c(horizontals[-1], rep(max(horizontals), ncol), min(horizontals), 0)
  
  # top margin
  lab_tp <- c(horizontals[-length(horizontals)], rep(1, ncol), max(horizontals), min(horizontals))
  
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
  
  out <- warbleR:::pblapply_wrblr_int(unique(soundid_X$page), function(x) {
    # extract data subset for a page
    Y <- soundid_X[soundid_X$page == x, ]
    
      # start graphic device
      warbleR:::img_wrlbr_int(
      filename = paste0(
        "plot_degradation_p", x, ".jpeg"
      ),
      path = dest.path,
      # width = 1000,
      # height = 1000,
      units = "in",
      res = 120
    )
    
      # set panel layout 
      ss <- split.screen(figs = page_layout)

      # get combination of distances and sound id to loop over it
      grd <- expand.grid(distance = distances, sound.id =  unique(Y$sound.id))
      
      # starting screen
      scrn <- 1
      
      for (i in seq_len(nrow(grd))) {
        
        # get id of signal
        sgnl <- Y$.sgnl.temp[Y$distance == grd$distance[i] & Y$sound.id == grd$sound.id[i]]      
        
        # get row index for signals in X
        indx <- which(X$.sgnl.temp == sgnl)
        
        # read wave
        wave <-
          read_wave(X = X,
                    index = indx,
                    from = if (X$start[indx] - mar < 0) 0 else X$start[indx] - mar,
                    to = X$end[indx] + mar,
                    path = path)
        screen(i)
        par(mar = c(0, 0, 0, 0))

          
        warbleR:::spectro_wrblr_int2(
          wave = wave,
          collevels = collevels,
          ovlp = ovlp,
          wl = wl,
          flim = c(0.1, 10.6),
          palette = viridis,
          axisX = FALSE,
          axisY = FALSE,
          grid = FALSE
        )

        # # add frequency axis
        # if (grepl("open", aligned_tests$sound.files[aligned_tests$sound.files == files[i]])) {
        #   axis(2, at = c(seq(2, 10, 2)))
        # }
        # 
        # # add time axis
        # if (grepl("100m", aligned_tests$sound.files[aligned_tests$sound.files == files[i]])) {
        #   axis(1)
        # }
        # 
        # 
        # lns <-
        #   c(
        #     aligned_tests$start[aligned_tests$sound.files == files[i]] - min(aligned_tests$start[aligned_tests$sound.files == files[i]]),
        #     aligned_tests$end[aligned_tests$sound.files == files[i]] - min(aligned_tests$start[aligned_tests$sound.files == files[i]])
        #   ) + 1.2
        # 
        # lns <- c(lns, min(lns) - 0.1, min(lns) - 1.05)
        # 
        # abline(
        #   v = lns,
        #   col = "white",
        #   lty = 3,
        #   lwd = 1.2
        # )
        # abline(
        #   v = lns,
        #   col = "white",
        #   lty = 3,
        #   lwd = 1.2
        # )
        scrn <- scrn + 1
      }

      # frequency label
      screen(((nrow + 1) * (ncol + 1)))
      par(mar = c(0, 0, 0, 0), new = TRUE)
      plot(
        1,
        frame.plot = FALSE,
        type = "n",
        yaxt = "n",
        xaxt = "n"
      )
      text(
        x = 0.9,
        y = 1,
        "Frequency (kHz)",
        srt = 90,
        cex = 1.2
      )
      
      # time label
      screen(((nrow + 1) * (ncol + 1)) + 1)
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
      
      
      dist_labs <- paste(distances, "m")

      par(
        mar = c(0, 0, 0, 0),
        bg = "black",
        new = TRUE
      )
      
      # add vertical labels
      for (i in seq_len(ncol) + (nrow * (ncol + 1))) {
        screen(i)
        # par(mar = c(0, 0, 0, 0))
        par(
          mar = c(0, 0, 0, 0),
          bg = "black",
          new = TRUE
        )
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
          dist_labs[i - (nrow * (ncol + 1))],
          cex = 1.2,
          col = "white",
          font = 2
        )
        box()
      # }

      # soundid_transect <- c("Open habitat", "Closed habitat")
      # for (i in seq_len(nrow) + (ncol * nrow)) {
      #   screen(i)
      #   par(
      #     mar = c(0, 0, 0, 0),
      #     bg = "black",
      #     new = TRUE
      #   )
      #   plot(
      #     1,
      #     frame.plot = FALSE,
      #     type = "n",
      #     yaxt = "n",
      #     xaxt = "n"
      #   )
      #   text(
      #     x = 1,
      #     y = 1,
      #     hlabs[i - 12],
      #     font = 2,
      #     cex = 1.6,
      #     col = "white"
      #   )
      #   box()
      # }

      dev.off()
    }
  })
}
