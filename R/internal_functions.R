# internal function not to be called by users
# stop function that doesn't print call
stop2 <- function (...)
{
  stop(..., call. = FALSE)
}

# internal baRulho function, not to be called by users. It prepares X for comparing sounds
# @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})
# last modification on sep-2022 (MAS)

prep_X_bRlo_int <- function(X, method = getOption("method", 1), cores = getOption("mc.cores", 1), pb = getOption("pb", TRUE)) {
  
  # set pb options 
  on.exit(pbapply::pboptions(type = .Options$pboptions$type), add = TRUE)
  
  # add sound file selec colums to X (weird column name so it does not overwrite user columns)
  X$TEMP....sgnl <- paste(X$sound.files, X$selec, sep = "-")
  
  # make it a regular data frame (no est)
  X2 <- as.data.frame(X) 
  
  # set clusters for windows OS
  if (Sys.info()[1] == "Windows" & cores > 1)
    cl <- parallel::makePSOCKcluster(getOption("cl.cores", cores)) else cl <- cores
  
  # set pb options 
  pbapply::pboptions(type = ifelse(as.logical(pb), "timer", "none"))
  
  # add second column with names of the reference sounds to be compared against
  X$reference <- pbapply::pbsapply(1:nrow(X), cl = cl, function(x, meth = method){
    
    # extract for single sound and order by distance
    Y <- X2[X2$sound.id == X$sound.id[X2$TEMP....sgnl == X2$TEMP....sgnl[x]], , drop = FALSE]
    Y <- Y[order(Y$distance), ]
    
    # method 1 compare to closest distance to source
    if (meth == 1) z <- Y$TEMP....sgnl[which.min(Y$distance)] else # if method 2
      # if not the first row then the previous row
      if (Y$TEMP....sgnl[1] != X$TEMP....sgnl[x]) z <- X$TEMP....sgnl[x - 1] else # else the first row
        z <- Y$TEMP....sgnl[1] 
    
    return(z)
  })
  
  return(X)
}


### functions wrappers over warbleR's functions
# img_wrlbr_int and were copied from warbleR 1.25
# author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})

# copied from warbleR::img_wrlbr_int()
img_bRlo_int <- function(filename, path = getOption("sound.files.path", "."), res = 160, units = "in", width = 8.5, height = 11, horizontal = FALSE){
  
  if (horizontal & missing(width)) {
    width <- 11
    height <- 8.5
  }
  
  # add path to filename
  flnm <- file.path(path, filename)
  
  # jpeg 
  if (grepl("jpeg$|jpg$", filename)) 
    jpeg(filename = flnm, res = res, units = units, width = width, height = height) else # or tiff
      tiff(filename = flnm, res = res, units = units, width = width, height = height) 
}

# copied from warbleR::spectro_wrblr_int2()
spectro_bRlo_int2 <- function(wave, f, wl = 512, wn = "hanning", zp = 0, ovlp = 0, 
                              complex = FALSE, norm = TRUE, fftw = FALSE, dB = "max0", 
                              dBref = NULL, plot = TRUE, grid = TRUE, 
                              cont = FALSE, collevels = NULL, palette = spectro.colors, 
                              contlevels = NULL, colcont = "black", colbg = "white", colgrid = "gray", 
                              colaxis = "black", collab = "black", cexlab = 1, cexaxis = 1, 
                              tlab = "Time (s)", flab = "Frequency (kHz)", alab = "Amplitude", 
                              scalelab = "Amplitude\n(dB)", main = NULL, scalefontlab = 1, 
                              scalecexlab = 0.75, axisX = TRUE, axisY = TRUE, tlim = NULL, 
                              trel = TRUE, flim = NULL, flimd = NULL, widths = c(6, 1), 
                              heights = c(3, 1), oma = rep(0, 4), listen = FALSE, fast.spec = FALSE, 
                              rm.zero = FALSE, amp.cutoff = NULL, X = NULL, palette.2 = reverse.topo.colors, bx = TRUE, add = FALSE, collev.min = NULL) 
{
  if (!is.null(dB) && all(dB != c("max0", "A", "B", "C", "D"))) 
    stop2("'dB' has to be one of the following character strings: 'max0', 'A', 'B', 'C' or 'D'")
  sel.tab <- X
  
  if (is.list(palette)) palette <- unlist(palette[[1]])
  if (is.null(palette)) palette <- spectro.colors  
  if (!is.function(palette)) palette <- get(palette)
  
  if (is.null(collevels) & !is.null(collev.min))
    collevels <- seq(collev.min, 0, 1)
  
  if (!is.null(sel.tab)) fast.spec <- TRUE 
  
  if (complex & norm) {
    norm <- FALSE
    warning("\n'norm' was turned to 'FALSE'")
  }
  if (complex & !is.null(dB)) {
    dB <- NULL
    warning("\n'dB' was turned to 'NULL'")
  }
  input <- seewave::inputw(wave = wave, f = f)
  
  wave <- input$w
  
  f <- input$f
  rm(input)
  if (!is.null(tlim)) 
    wave <- cutw(wave, f = f, from = tlim[1], to = tlim[2])
  if (!is.null(flimd)) {
    mag <- round((f/2000)/(flimd[2] - flimd[1]))
    wl <- wl * mag
    if (ovlp == 0) 
      ovlp <- 100
    ovlp <- 100 - round(ovlp/mag)
    flim <- flimd
  }
  n <- nrow(wave)
  step <- seq(1, n - wl, wl - (ovlp * wl/100))
  
  # to fix function name change in after version 2.0.5
  # if (exists("stdft")) stft <- stdft
  z <- stft_bRlo_int(wave = wave, f = f, wl = wl, zp = zp, step = step, 
                     wn = wn, fftw = fftw, scale = norm, complex = complex)
  if (!is.null(tlim) && trel) {
    X <- seq(tlim[1], tlim[2], length.out = length(step))
  }  else {
    X <- seq(0, n/f, length.out = length(step))
  }
  if (is.null(flim)) {
    Y <- seq(0, (f/2) - (f/wl), length.out = nrow(z))/1000
  } else {
    fl1 <- flim[1] * nrow(z) * 2000/f
    fl2 <- flim[2] * nrow(z) * 2000/f
    z <- z[(fl1:fl2) + 1, ]
    Y <- seq(flim[1], flim[2], length.out = nrow(z))
  }
  if (!is.null(dB)) {
    if (is.null(dBref)) {
      z <- 20 * log10(z)
    } else {
      z <- 20 * log10(z/dBref)
    }
    if (dB != "max0") {
      if (dB == "A") 
        z <- dBweight(Y * 1000, dBref = z)$A
      if (dB == "B") 
        z <- dBweight(Y * 1000, dBref = z)$B
      if (dB == "C") 
        z <- dBweight(Y * 1000, dBref = z)$C
      if (dB == "D") 
        z <- dBweight(Y * 1000, dBref = z)$D
    }
  }
  Z <- t(z)
  
  maxz <- round(max(z, na.rm = TRUE))
  if (!is.null(dB)) {
    if (is.null(collevels)) 
      collevels <- seq(maxz - 30, maxz, by = 1)
    if (is.null(contlevels)) 
      contlevels <- seq(maxz - 30, maxz, by = 10)
  } else {
    if (is.null(collevels)) 
      collevels <- seq(0, maxz, length = 30)
    if (is.null(contlevels)) 
      contlevels <- seq(0, maxz, length = 3)
  }
  Zlim <- range(Z, finite = TRUE, na.rm = TRUE)
  
  if (!is.null(amp.cutoff)) Z[Z >= (diff(range(Z)) * amp.cutoff) + min(Z)] <- 0 
  
  if (!fast.spec)
    filled_contour_bRlo_int(x = X, y = Y, z = Z, bg.col = colbg, levels = collevels, 
                            nlevels = 20, plot.title = title(main = main, 
                                                             xlab = tlab, ylab = flab), color.palette = palette, 
                            axisX = FALSE, axisY = axisY, col.lab = collab, 
                            colaxis = colaxis, add = add) else {
                              image(x = X, y = Y, z = Z, col = palette(30), xlab = tlab, ylab = flab, axes = FALSE)
                              if (!is.null(sel.tab))    
                                out <- lapply(1:nrow(sel.tab), function(i)
                                  image(x = X[X > sel.tab$start[i] & X < sel.tab$end[i]], y = Y[Y > sel.tab$bottom.freq[i] & Y < sel.tab$top.freq[i]], z = Z[X > sel.tab$start[i] & X < sel.tab$end[i], Y > sel.tab$bottom.freq[i] & Y < sel.tab$top.freq[i]], col = palette.2(30), xlab = tlab, ylab = flab, axes = FALSE, xlim = range(X), add = TRUE)      
                                )
                              
                              
                              if (axisY) axis(2, at = pretty(Y), labels = pretty(Y), cex.axis = cexlab)
                              if (bx)  box()
                              if (!is.null(main)) title(main)       
                            }
  
  if (axisX) {
    if (rm.zero)
      axis(1, at = pretty(X)[-1], labels = pretty(X)[-1], cex.axis = cexaxis)  else
        axis(1, at = pretty(X), labels = pretty(X), cex.axis = cexaxis) 
  }
  
  if (grid) 
    grid(nx = NA, ny = NULL, col = colgrid)
  
}

####

stft_bRlo_int <- function (wave, f, wl, zp, step, wn, scale = TRUE, norm = FALSE,                           correction = "none", fftw = FALSE, complex = FALSE) 
{
  if (zp < 0) 
    stop2("zero-padding cannot be negative")
  W <- ftwindow(wl = wl, wn = wn, correction = correction)
  if (fftw) {
    p <- fftw::planFFT(wl + zp)
    z <- apply(as.matrix(step), 1, function(x) fftw::FFT(c(wave[x:(wl + 
                                                                     x - 1), ] * W, rep(0, zp)), plan = p))
  }
  else {
    z <- apply(as.matrix(step), 1, function(x) stats::fft(c(wave[x:(wl + 
                                                                      x - 1), ] * W, rep(0, zp))))
  }
  z <- z[1:((wl + zp)%/%2), , drop = FALSE]
  z <- z/(wl + zp)
  if (complex == FALSE) {
    z <- 2 * Mod(z)
    if (scale) {
      if (norm) {
        z <- z/apply(X = z, MARGIN = 2, FUN = max)
      }
      else {
        z <- z/max(z)
      }
    }
  }
  return(z)
}


#####

filled_contour_bRlo_int <- function (x = seq(0, 1, len = nrow(z)), y = seq(0, 1, len = ncol(z)), 
                                     z, xlim = range(x, finite = TRUE), ylim = range(y, finite = TRUE), col.lab, colaxis,
                                     zlim = range(z, finite = TRUE), levels = pretty(zlim, nlevels), add = FALSE,
                                     nlevels = 20, color.palette = cm.colors, col = color.palette(length(levels) - 
                                                                                                    1), plot.title, plot.axes, key.title, asp = NA, xaxs = "i", 
                                     yaxs = "i", las = 1, axisX = TRUE, axisY = TRUE, bg.col = "white") 
{
  if (missing(z)) {
    if (!missing(x)) {
      if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      } else {
        z <- x
        x <- seq(0, 1, len = nrow(z))
      }
    } else stop2("no 'z' matrix specified")
  } else if (is.list(x)) {
    y <- x$y
    x <- x$x
  }
  if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
    stop2("increasing 'x' and 'y' values expected")
  if (!add) plot.new()
  plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
  usr <- par("usr")
  rect(xleft = usr[1], xright = usr[2], ybottom = usr[3], usr[4], col = bg.col, border = bg.col)
  
  
  if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1) 
    stop2("no proper 'z' matrix specified")
  if (!is.double(z)) 
    storage.mode(z) <- "double"
  .filled.contour(as.double(x), as.double(y), z, as.double(levels), 
                  col = col)
  if (missing(plot.axes)) {
    if (axisX) {
      title(main = "", xlab = "", ylab = "")
      axis(1)
    }
    if (axisY) {
      title(main = "", xlab = "", ylab = "")
      axis(2)
    }
  }
  else plot.axes
  box()
  if (missing(plot.title)) 
    title()
  else plot.title
  invisible()
}



# internal baRulho function, not to be called by users. It is a modified version of warbleR::find_peaks
# that allows to define internally if progress bar would be used (pbapply::pblapply uses pboptions to do this) 
# Find cross-correlation peaks

# last modification on jan-014-2021 (MAS)
find_peaks_bRlh_int <- function(xc.output, cores = getOption("mc.cores", 1), cutoff = 0.4, pb = getOption("pb", TRUE), max.peak = FALSE, output = "data.frame") 
{
  
  # set clusters for windows OS and no soz
  if (Sys.info()[1] == "Windows" & cores > 1)
    cl <- parallel::makePSOCKcluster(getOption("cl.cores", cores)) else cl <- cores
    
    # loop over scores of each dyad
    pks <- warbleR:::pblapply_wrblr_int(pbar = pb, X = unique(xc.output$scores$dyad), cl = cl, FUN = function(i) {
      
      # extract data for a dyad
      dat <- xc.output$scores[xc.output$scores$dyad == i, ]
      
      # check xc.output being a autodetec.output object
      if (!(is(xc.output, "xcorr.output") | is(xc.output, "xc.output"))) 
        stop2("'xc.output' must be and object of class 'xcorr.output'")
      
      ## get peaks as the ones higher than previous and following scores  
      pks <- dat[c(FALSE, diff(dat$score) > 0) & c(rev(diff(rev(dat$score)) > 0), FALSE) & dat$score > cutoff, , drop = FALSE]
      
      # get the single highest peak
      if (max.peak)
        pks <- dat[which.max(dat$score), , drop = FALSE]
      
      return(pks)
    })
    
    # put results in a data frame
    peaks <- do.call(rbind, pks)
    
    # relabel rows
    if (nrow(peaks) > 0)
    {  rownames(peaks) <- 1:nrow(peaks)
    
    # remove dyad column
    peaks$dyad <- NULL
    
    #### name as in a warbleR selection table
    # remove selec info at the end
    peaks$sound.files <- substr(peaks$sound.files, start = 0, regexpr("\\-[^\\-]*$", peaks$sound.files) - 1)
    
    #### add start and end
    # add template column to selection table in xc.output
    Y <- xc.output$org.selection.table
    Y$template <- paste(Y$sound.files, Y$selec, sep = "-")
    
    # Y <- Y[Y$template %in% comp_mat[, 1], ]
    
    # add start as time - half duration of template
    peaks$start <- sapply(1:nrow(peaks), function(i){
      
      peaks$time[i] - 
        ((Y$end[Y$template == peaks$template[i]] - 
            Y$start[Y$template == peaks$template[i]])  / 2)
      
    })
    
    # add end as time + half duration of template
    peaks$end <- sapply(1:nrow(peaks), function(i){
      
      peaks$time[i] + 
        ((Y$end[Y$template == peaks$template[i]] - 
            Y$start[Y$template == peaks$template[i]]) / 2)
      
    })
    
    # add selec labels
    peaks$selec <- 1
    
    if (nrow(peaks) > 1)
      for(i in 2:nrow(peaks)) 
        if (peaks$sound.files[i] == peaks$sound.files[i - 1])
          peaks$selec[i] <- peaks$selec[i - 1] + 1
    
    # sort columns in a intuitive order
    peaks <- warbleR::sort_colms(peaks)
    
    # output results
    if (output == "data.frame") return(peaks) else{
      
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

#
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
