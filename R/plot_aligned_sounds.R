#' Plot spectrograms to check test sound files alignment
#'
#' \code{plot_aligned_sounds} plots spectrograms to visually inspect alignment precision on test sound files.
#' @usage plot_aligned_sounds(X, hop.size = getOption("hop.size", 11.6),
#' wl = getOption("wl", NULL), ovlp = getOption("ovlp", 50),
#' path = getOption("sound.files.path", "."), cores = getOption("mc.cores", 1),
#' pb = getOption("pb", TRUE), collevels = seq(-120, 0, 5),
#' palette = viridis::viridis, duration = 2,
#' mar = 0.2, dest.path = getOption("dest.path", "."), flim = NULL, col = "white", width = 7,
#' height = 4, res = 100, label = TRUE, fast.spec = FALSE, srt = 0, cex = 1, ...)
#' @param X Object of class 'data.frame', 'selection_table' or 'extended_selection_table' (the last 2 classes are created by the function \code{\link[warbleR]{selection_table}} from the warbleR package) with the reference to the sounds in the master sound file. Must contain the following columns: 1) "sound.files": name of the .wav files, 2) "selec": unique selection identifier (within a sound file), 3) "start": start time and 4) "end": end time of selections, 5)  "bottom.freq": low frequency for bandpass, 6) "top.freq": high frequency for bandpass and 7) "sound.id": ID of sounds used to identify counterparts across distances. Each sound must have a unique ID within a distance.
#' @param hop.size A numeric vector of length 1 specifying the time window duration (in ms). Default is 1 1. 6ms, which is equivalent to 512 wl for a 44.1 kHz sampling rate. Ignored if 'wl' is supplied.
#' @param wl A numeric vector of length 1 specifying the window length of the spectrogram, default
#' is NULL. Ignored if \code{bp = NULL}. If supplied, 'hop.size' is ignored.
#' @param ovlp Numeric vector of length 1 specifying the percent overlap between two
#'   consecutive windows, as in \code{\link[seewave]{spectro}}. Default is 0.
#' @param path Character string containing the directory path where the sound files are found. Only needed when 'X' is not an extended selection table.
#' @param cores Numeric vector of length 1. Controls whether parallel computing is applied by specifying the number of cores to be used. Default is 1 (i.e. no parallel computing).
#' @param pb Logical argument to control if progress bar is shown. Default is \code{TRUE}.
#' @param collevels A numeric vector of length 3. Specifies levels to partition the
#'   amplitude range of the spectrogram (in dB). The more levels the higher the
#'   resolution of the spectrogram. Default is seq(-40, 0, 1). seq(-115, 0, 1) will produces spectrograms
#'   similar to other acoustic analysis software packages.
#' @param palette Color palette function for spectrogram. Default is  \code{\link[viridis]{viridis}}. See
#' \code{\link[seewave]{spectro}} for more palettes. Palettes as \code{\link[monitoR:specCols]{gray.2}} may work better when \code{fast.spec = TRUE}.
#' @param duration A numeric vector of length 1. Specifies the overall duration of the clip that will be plotted. Notice that only the initial part of the test files are plotted as this is enough to tell the precision of the alignment.
#' @param mar numeric vector of length 1. Specifies the margins adjacent to the start of the first annotation to be included in the plot.
#' @param dest.path Character string containing the directory path where the image files will be save. If not supplied files will be save in the current working directory.
#' @param flim A numeric vector of length 2 indicating the highest and lowest frequency limits (kHz) of the spectrogram, as in \code{\link[seewave]{spectro}}. Default is \code{NULL} which will plot spectrograms in the full frequency range (0 - nyquist frequency).
#' @param col Character string controlling the color of lines and sound ID labels.
#' @param width Numeric vector of length 1. Single value (in inches) indicating the width of the output image files. Default is 7.
#' @param height Numeric vector of length 1. Single value (in inches) indicating the height of the output image files. Default is 4.
#' @param res Numeric argument of length 1. Controls image resolution. Default is 100 (faster) although 300 - 400 is recommended for publication/presentation quality.
#' @param label Logical to control if labels (from 'sound.id' column in 'X') are plotted. Default is  \code{TRUE}.
#' @param fast.spec Logical. If \code{TRUE} then image function is used internally to create spectrograms, which substantially
#' increases performance (much faster), although some options become unavailable, as collevels (amplitude scale). Default is \code{FALSE}.
#' @param srt Numeric argument of length 1. The rotation (in degrees) of the sound id labels. Default is 0.
#' @param cex Numeric argument of length 1controlling the size of sound id text labels. Default is 1.
#' @param ... Additional arguments to be passed to the internal spectrogram
#' creating function for customizing graphical output. The function is a modified
#' version of \code{\link[seewave]{spectro}}, so it takes the same arguments.
#' @return Image files in jpeg format with spectrograms in the working directory, one for each sound file in 'X'.
#' @export
#' @name plot_aligned_sounds
#' @details This functions aims to simplify the evaluation of the alignment of test sound files from  \code{\link{align_test_files}}. The function creates a single spectrogram for each sound file (saved at 'dest.path'). Spectrograms include the first few seconds of the sound files (controlled by 'duration') which is usually enough to tell the precision of the alignment. The plots include vertical lines denoting the start and end of each sound as well as the sound ID ('sound.id' column in 'X'). Note that no plot is created in the R graphic device.
#' @family test sound alignment
#' @seealso \code{\link{manual_realign}};  \code{\link{auto_realign}}; \code{\link{find_markers}}; \code{\link{align_test_files}}
#' @examples {
#'   # load example data
#'   data("test_sounds_est")
#'
#'   # plot (look into temporary working directory `tempdir()`)
#'   plot_aligned_sounds(X = test_sounds_est, dest.path = tempdir(), duration = 3, ovlp = 0)
#' }
#'
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})
#' @references {
#' Araya-Salas, M. (2020). baRulho: baRulho: quantifying degradation of (animal) acoustic signals in R. R package version 1.0.2
#' }


plot_aligned_sounds <-
  function(X,
           hop.size = getOption("hop.size", 11.6),
           wl = getOption("wl", NULL),
           ovlp = getOption("ovlp", 50),
           path = getOption("sound.files.path", "."),
           cores = getOption("mc.cores", 1),
           pb = getOption("pb", TRUE),
           collevels = seq(-120, 0, 5),
           palette = viridis::viridis,
           duration = 2,
           mar = 0.2,
           dest.path = getOption("dest.path", "."),
           flim = NULL,
           col = "white",
           width = 7,
           height = 4,
           res = 100,
           label = TRUE,
           fast.spec = FALSE,
           srt = 0,
           cex = 1,
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
    
    
    # get sampling rate
    sampling_rate <-
      warbleR::read_sound_file(
        X = X,
        index = 1,
        path = path,
        header = TRUE
      )$sample.rate
    
    # adjust wl based on hope.size
    if (is.null(wl)) {
      wl <- round(sampling_rate * hop.size / 1000,
                  0)
    }
    
    # make wl even if odd
    if (!(wl %% 2) == 0) {
      wl <- wl + 1
    }
    
    # split data set by sound file
    X_by_sound_file <- split(X, X$sound.files)
    
    ### run loop over sound files
    # set clusters for windows OS
    if (Sys.info()[1] == "Windows" & cores > 1) {
      cl <- parallel::makePSOCKcluster(getOption("cl.cores", cores))
    } else {
      cl <- cores
    }
    
    # run loop
    out <-
      warbleR:::pblapply_wrblr_int(
        pbar = pb,
        X = X_by_sound_file,
        cl = cl,
        FUN = function(Y) {
          warbleR:::img_wrlbr_int(
            filename = paste0("plot_align_", gsub(".wav", "", Y$sound.files[1]), ".jpeg"),
            path = dest.path,
            width = width,
            height = height,
            units = "in",
            res = res
          )
          
          # set start and end of clip to plot
          from <- if (min(Y$start) - mar >= 0) {
            min(Y$start) - mar
          } else {
            0
          }
          to <- from + duration
          
          # import sound data
          wave <- read_wave(
            Y,
            index = 1,
            from = from,
            to = to,
            path = path
          )
          
          
          # adjust spectrogram margins
          prev_mar <- par("mar")
          par(mar = c(5.1, 4, 3, 1))
          
          # reset graphic parameterswhen function is don
          on.exit(par(mar = prev_mar))
          
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
            main = Y$sound.files[1],
            fast.spec = fast.spec,
            ...
          )
          
          # add dotted lines
          abline(
            v = c(Y$start - from, Y$end - from),
            col = col,
            lty = 3,
            lwd = 1.5
          )
          
          # add sound.id labels
          if (label) {
            # position of sound.id labels in the freq axis
            y_pos <- flim[2] - 2 * ((flim[2] - flim[1]) / 12)
            
            # plot sound id labels
            text(
              labels = Y$sound.id,
              x = ((Y$end + Y$start) / 2) - from,
              y = y_pos,
              pos = 3,
              col = col,
              srt = srt,
              cex = cex
            )
          }
          
          dev.off()
        }
      )
  }
