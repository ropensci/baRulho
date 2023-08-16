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
#' #' @return A data frame, or extended selection table similar to input data (depending on argument 'output'), but includes a new column (cross.correlation)
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
#' plot_degradation <- function(X, hop.size = getOption("hop.size", 1), wl = getOption("wl", NULL), ovlp = getOption("ovlp", 0), path = getOption("sound.files.path", "."), cores = getOption("mc.cores", 1), pb = getOption("pb", TRUE), collevels = seq(-120, 0, 5), pal = viridis, box = FALSE, ...) {
#'   X <- X[X$sound.id != "ambient", ]
#'
#'   # stop if more than 1 sample per distance but no transect info
#'   if (any(table(X$sound.id, X$distance)) > 1 & is.null(X$transect)) {
#'     stop2("There are more than 1 test sound per sound.id/distance combination but no 'transect' column to group by transect")
#'   }
#'
#'   # stop if there are more than 1 sample for distance sound id and transect combination (only 1 samp)
#'   if (!is.null(X$transect)) {
#'     if (any(table(X$sound.id, X$distance, X$transect)) > 1) {
#'       stop2("There are more than 1 test sound per sound.id/distance/transect combination")
#'     }
#'   }
#'
#'   if (is.null(X$transect)) {
#'     transects <- 1
#'     X$transect <- 1
#'   } else {
#'     unique(X$transect)
#'   }
#'
#'   X_by_transect <- split(X, f = X$transect)
#'
#'
#'
#'
#'
#'   out <- warbleR:::pblapply_wrblr_int(X_by_transect, function(Y) {
#'     # warbleR::catalog(X = Y, nrow = length(unique(Y$sound.id)), ncol = length(unique(Y$distance)), same.time.scale = FALSE, path = path, flim = c(0, 10), parallel = cores, collevels = collevels, pal = pal, box = box, spec.mar = 0.01, ovlp = ovlp, pb = FALSE, rm.axes = TRUE, img.prefix = paste0("transect-", Y$transect[1]), title = paste("Transect", Y$transect[1]), group.tag = "sound.id", tag.pal = list(magma), ...)
#'     lf <- rep(c(0.06, 0.5), each = 4)
#'     rg <- rep(c(0.5, 0.95), each = 4)
#'     horiz <- seq(0.95, 0.075, length.out = 5)
#'     btm <- rep(horiz[-1], 2)
#'     tp <- rep(horiz[-length(horiz)], 2)
#'
#'     m <- cbind(lf, rg, btm, tp)
#'
#'     lf <- c(rep(0.95, each = 4), 0.06, 0.5, 0, 0)
#'     rg <- c(rep(1, each = 4), 0.5, 0.95, 0.05, 1)
#'     horiz <- seq(0.95, 0.075, length.out = 5)
#'     btm <- c(horiz[-1], 0.95, 0.95, 0.075, 0)
#'     tp <- c(horiz[-length(horiz)], 1, 1, 0.95, 0.075)
#'
#'     m2 <- cbind(lf, rg, btm, tp)
#'     m <- rbind(m, m2)
#'
#'
#'     for (e in paste0("trnsc", 1:3)) {
#'       png(filename = paste0("./output/spectrograms_by_habitat_and_distance_", e, ".png"), res = 300, width = 4000, height = 3000)
#'
#'       ss <- split.screen(figs = m)
#'
#'       # # testing layout screens
#'       # for(i in 1:nrow(m))
#'       # {screen(i)
#'       #   par( mar = rep(0, 4))
#'       #   plot(0.5, xlim = c(0,1), ylim = c(0,1), type = "n", axes = FALSE, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
#'       #   box()
#'       #   text(x = 0.5, y = 0.5, labels = i)
#'       # }
#'       # close.screen(all.screens = T)
#'
#'
#'       ovlp <- 1
#'       colls <- seq(-110, 0, 5)
#'       wl <- 512
#'
#'       lab_bg <- viridis(10, alpha = 0.25)[8]
#'
#'       lab_bg <- "#3A3B39"
#'       files <-
#'         c(
#'           "_10m_open.wav",
#'           "_30m_open.wav",
#'           "_65m_open.wav",
#'           "_100m_open.wav",
#'           "_10m_closed.wav",
#'           "_30m_closed.wav",
#'           "_65m_closed.wav",
#'           "_100m_closed.wav"
#'         )
#'
#'       files <- paste0(e, files)
#'
#'       titles <- gsub(paste0(e, "_|_|.wav"), " ", files)
#'
#'
#'       # frequency label
#'       screen(15)
#'       par(mar = c(0, 0, 0, 0), new = TRUE)
#'       plot(
#'         1,
#'         frame.plot = FALSE,
#'         type = "n",
#'         yaxt = "n",
#'         xaxt = "n"
#'       )
#'       text(
#'         x = 0.9,
#'         y = 1,
#'         "Frequency (kHz)",
#'         srt = 90,
#'         cex = 1.6
#'       )
#'
#'       # time label
#'       screen(16)
#'       par(mar = c(0, 0, 0, 0), new = TRUE)
#'       plot(
#'         1,
#'         frame.plot = FALSE,
#'         type = "n",
#'         yaxt = "n",
#'         xaxt = "n"
#'       )
#'       text(
#'         x = 1,
#'         y = 0.75,
#'         "Time (s)",
#'         cex = 1.6
#'       )
#'
#'       for (i in seq_along(files)) {
#'         # print(i)
#'         screen(i)
#'         par(mar = c(0, 0, 0, 0))
#'
#'         warbleR:::spectro_wrblr_int2(
#'           wave = readWave(
#'             file.path(path_to_files, files[i]),
#'             from = min(aligned_tests$start[aligned_tests$sound.files == files[i]]) - 1.2,
#'             to = min(aligned_tests$start[aligned_tests$sound.files == files[i]]) + 4,
#'             units = "seconds"
#'           ),
#'           collevels = colls,
#'           ovlp = ovlp,
#'           wl = wl,
#'           flim = c(0.1, 10.6),
#'           palette = viridis,
#'           axisX = FALSE,
#'           axisY = FALSE,
#'           grid = FALSE
#'         )
#'
#'         # add frequency axis
#'         if (grepl("open", aligned_tests$sound.files[aligned_tests$sound.files == files[i]])) {
#'           axis(2, at = c(seq(2, 10, 2)))
#'         }
#'
#'         # add time axis
#'         if (grepl("100m", aligned_tests$sound.files[aligned_tests$sound.files == files[i]])) {
#'           axis(1)
#'         }
#'
#'
#'         lns <-
#'           c(
#'             aligned_tests$start[aligned_tests$sound.files == files[i]] - min(aligned_tests$start[aligned_tests$sound.files == files[i]]),
#'             aligned_tests$end[aligned_tests$sound.files == files[i]] - min(aligned_tests$start[aligned_tests$sound.files == files[i]])
#'           ) + 1.2
#'
#'         lns <- c(lns, min(lns) - 0.1, min(lns) - 1.05)
#'
#'         abline(
#'           v = lns,
#'           col = "white",
#'           lty = 3,
#'           lwd = 1.2
#'         )
#'         abline(
#'           v = lns,
#'           col = "white",
#'           lty = 3,
#'           lwd = 1.2
#'         )
#'       }
#'
#'       vlabs <- paste(c(10, 30, 65, 100), "m")
#'
#'       par(
#'         mar = c(0, 0, 0, 0),
#'         bg = lab_bg,
#'         new = TRUE
#'       )
#'       # add vertical labels
#'       for (i in 9:12) {
#'         screen(i)
#'         # par(mar = c(0, 0, 0, 0))
#'         par(
#'           mar = c(0, 0, 0, 0),
#'           bg = lab_bg,
#'           new = TRUE
#'         )
#'         plot(
#'           1,
#'           frame.plot = FALSE,
#'           type = "n",
#'           yaxt = "n",
#'           xaxt = "n"
#'         )
#'         text(
#'           x = 1,
#'           y = 1,
#'           vlabs[i - 8],
#'           srt = 270,
#'           cex = 1.6,
#'           col = "white",
#'           font = 2
#'         )
#'         box()
#'       }
#'
#'       hlabs <- c("Open habitat", "Closed habitat")
#'       for (i in 13:14) {
#'         screen(i)
#'         par(
#'           mar = c(0, 0, 0, 0),
#'           bg = lab_bg,
#'           new = TRUE
#'         )
#'         plot(
#'           1,
#'           frame.plot = FALSE,
#'           type = "n",
#'           yaxt = "n",
#'           xaxt = "n"
#'         )
#'         text(
#'           x = 1,
#'           y = 1,
#'           hlabs[i - 12],
#'           font = 2,
#'           cex = 1.6,
#'           col = "white"
#'         )
#'         box()
#'       }
#'
#'       dev.off()
#'     }
#'
#'
#'
#'     # split.screen(figs = m)
#'     # if (nrow(X) < nrow * ncol)
#'     #
#'   })
#' }
