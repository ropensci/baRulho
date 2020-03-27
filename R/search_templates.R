#' Search acoustic templates on test sound files
#' 
#' \code{search_templates} searches acoustic templates on test (re-recorded) sound files.
#' @usage search_templates(X, template.rows, test.files, path = NULL, pb = TRUE, ...)
#' @param X object of class 'data.frame', 'selection_table' or 'extended_selection_table' (the last 2 classes are created by the function \code{\link[warbleR]{selection_table}} from the warbleR package). Required.
#' @param template.rows numeric vector with the index of the rows from 'X' to be used as templates. If only 1 is supplied the same template will be run over all 'test.files'. Otherwise, 'template.rows' must be the same length as 'test.files'. Required.
#' @param test.files Character vector of length 1 with the name(s) of the test (re-recorded) file(s) in which to search for the template(s) (see argument 'template.rows').
#' @param path Character string containing the directory path where test (re-recorded) sound files are found.
#' @param pb Logical argument to control if progress bar is shown. Default is \code{TRUE}.
#' @param ...	Additional arguments to be passed to \code{\link[warbleR]{xcorr}} for setting cross-correlation parameters (e.g. 'wl', 'ovlp', etc).
#' @return A data frame with the time, start, end, test file names, template name, maximum cross-correlation score and the time where it was detected.
#' @export
#' @name search_templates
#' @details The function takes the output of \code{\link[warbleR]{find_peaks}} ('Y') and aligns signals found in re-recorded sound files according to a master sound file referenced in 'X'. The function outputs a 'extended selection table'.
#' @seealso \code{\link{spcc_align}}; \code{\link{align_test_files}}
#' @examples
#' \dontrun{
#' # use a temporary directory
#' td <- tempdir()  
#' 
#' # load example data from warbleR
#' data(list = c("Phae.long1", "Phae.long2", "Phae.long3", "Phae.long4", 
#' "lbh_selec_table"))
#' 
#' # save sound files to temporary folder
#' writeWave(Phae.long1, file.path(td, "Phae.long1.wav"))
#' writeWave(Phae.long2, file.path(td, "Phae.long2.wav"))
#' writeWave(Phae.long3, file.path(td, "Phae.long3.wav"))
#' writeWave(Phae.long4, file.path(td, "Phae.long4.wav"))
#' 
#' # make an extended selection table
#' est <- selection_table(X = lbh_selec_table, extended = TRUE, confirm.extended = FALSE, 
#' path = td, pb = FALSE)
#' 
#' # create master sound file
#' master.sf <- master_sound_file(X = est, file.name = "example_master", 
#' dest.path = td, gap.duration = 0.3)
#' 
#' # read master
#' exmp.master <- readWave(file.path(td, "example_master.wav"))
#' 
#' # add 1 s silence
#' exmp.test1 <- addsilw(wave = exmp.master, at = "start", d = 1,
#' output = "Wave", f = exmp.master@samp.rate)
#' 
#' exmp.test2 <- addsilw(wave = exmp.master, at = "start", d = 2, 
#' output = "Wave", f = exmp.master@samp.rate)
#' 
#' # create noise
#' ns <- noisew(f = exmp.master@samp.rate, d = duration(exmp.test2) + 1, 
#' output = "Wave")
#' 
#' # make noise exactly the same length and add noise to 2 examples
#' exmp.test1@left <- exmp.test1@left + (ns@left[1:length(exmp.test1@left)] * 500)
#' exmp.test2@left <- exmp.test2@left + (ns@left[1:length(exmp.test2@left)] * 500)
#' 
#' exmp.test1 <- normalize(exmp.test1, unit = "16")
#' exmp.test2 <- normalize(exmp.test2, unit = "16")
#' 
#' # save examples
#' writeWave(object = exmp.test1, filename = file.path(td, "example_test1.wav"), extensible = FALSE)
#' writeWave(object = exmp.test2, filename = file.path(td, "example_test2.wav"), extensible = FALSE)
#' 
#' # search using start marker as template
#' search_templates(X = master.sf, 
#' template.rows = which(master.sf$orig.sound.file == "start_marker"), 
#' test.files = c("example_test1.wav", "example_test2.wav"), path = td, pb = FALSE)
#' 
#' # search using start marker as template
#' search_templates(X = master.sf, template.rows = which(master.sf$orig.sound.file == "start_marker"), 
#' test.files = c("example_test1.wav", "example_test2.wav", "example_test1.wav"), 
#' path = td, pb = FALSE)
#' 
#' # search using both start and end markers as template
#' search_templates(X = master.sf, 
#' template.rows = which(master.sf$orig.sound.file == "start_marker" | 
#' master.sf$orig.sound.file == "end_marker"), 
#' test.files = c("example_test1.wav", "example_test1.wav"), 
#' path = td, pb = FALSE)
#' }
#' @author Marcelo Araya-Salas (\email{marceloa27@@gmail.com})
#' @references {
#' Araya-Salas, M. (2020). baRulho: baRulho: quantifying habitat-induced degradation of (animal) acoustic signals in R. R package version 1.0.2
#' }
#last modification on dec-26-2019 (MAS)

search_templates <- function(X, template.rows, test.files, path = NULL,  pb = TRUE, ...){
  
 if (length(template.rows) != 1 & length(test.files) != length(template.rows)) stop("'template.rows' must be 1 or the same length than 'test.files'")
  
  # create a matrix that contains the selection/files to be cross-correlated
  if (length(template.rows) == 1) # in only 1 template, repeate it
  comp_mat <- matrix(c(rep(paste(X$sound.files[template.rows], X$selec[template.rows], sep = "-"), length(test.files)), test.files), ncol = 2) else 
    comp_mat <- matrix(c(paste(X$sound.files[template.rows], X$selec[template.rows], sep = "-"), test.files), ncol = 2) 
  
  # run cross correlation
  xc <- xcorr(X, compare.matrix = comp_mat, path = path, output = "list", pb = pb, ...)
  
  # find peaks
  pks <- find_peaks(xc.output = xc, max.peak = TRUE, path = path, pb = pb, cutoff = 0)
  
  # rename sound file column
  names(pks)[names(pks) == "sound.files"] <- "test.files"
  return(pks)

  }
