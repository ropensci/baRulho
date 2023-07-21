#' Search acoustic templates on test sound files
#' 
#' \code{search_templates} searches acoustic templates on test (re-recorded) sound files.
#' @usage search_templates(X, template.rows, test.files = NULL, path = NULL, pb = TRUE, ...)
#' @param X Object of class 'data.frame', 'selection_table' or 'extended_selection_table' (the last 2 classes are created by the function \code{\link[warbleR]{selection_table}} from the warbleR package) with the reference to the sounds in the master sound file. Must contain the following columns: 1) "sound.files": name of the .wav files, 2) "selec": unique selection identifier (within a sound file), 3) "start": start time and 4) "end": end time of selections. Columns for 'top.freq', 'bottom.freq' and 'channel' are optional. Required. 
#' @param template.rows Numeric vector with the index of the rows from 'X' to be used as templates. If only 1 is supplied the same template will be run over all 'test.files'. Otherwise, 'template.rows' must be the same length as 'test.files'. Required.
#' @param test.files Character vector of length 1 with the name(s) of the test (re-recorded) file(s) in which to search for the template(s) (see argument 'template.rows'). If not supplied all sound files in 'path' are used instead.
#' @param path Character string containing the directory path where test (re-recorded) sound files are found.
#' @param pb Logical argument to control if progress bar is shown. Default is \code{TRUE}.
#' @param ...	Additional arguments to be passed to \code{\link[warbleR]{cross_correlation}} for setting cross-correlation parameters (e.g. 'wl', 'ovlp', etc).
#' @return A data frame with the time, start, end, test file names, template name, maximum cross-correlation score and the time where it was detected.
#' @export
#' @name search_templates
#' @details The function takes a master sound file's reference data ('X') and finds the position of acoustics templates (included as selections in 'X') in the re-recorded sound files. This is used to align signals found in re-recorded sound files according to a master sound file referenced in 'X'. \strong{Make sure the master sound file (that refered to in 'X') is found in the same folder than the re-recorded sound files}.Take a look at the package vignette for information on how to incorporate this function into a sound degradation analysis workflow.
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
#' exmp.test1 <- tuneR::normalize(exmp.test1, unit = "16")
#' exmp.test2 <- tuneR::normalize(exmp.test2, unit = "16")
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
#' # search using end marker as template
#' search_templates(X = master.sf, template.rows = which(master.sf$orig.sound.file == "end_marker"), 
#' test.files = c("example_test1.wav", "example_test2.wav"), 
#' path = td, pb = FALSE)
#' 
#' # search using both start and end markers as template
#' search_templates(X = master.sf, 
#' template.rows = which(master.sf$orig.sound.file == "start_marker" | 
#' master.sf$orig.sound.file == "end_marker"), 
#' test.files = c("example_test1.wav", "example_test2.wav"), 
#' path = td, pb = FALSE)
#' }
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})
#' @references {
#' Araya-Salas, M. (2020). baRulho: baRulho: quantifying habitat-induced degradation of (animal) acoustic signals in R. R package version 1.0.2
#' }
#last modification on dec-26-2019 (MAS)

search_templates <- function(X, template.rows, test.files = NULL, path = NULL, pb = TRUE, ...){
  
  # get sound files in path
  files_in_path <- list.files(path = path, pattern = "\\.wav$", ignore.case = TRUE)
  # check if there are files
  if (length(files_in_path) == 0) stop2("No .wav files found in 'path'")
  
  # check for master sound file
  if (!any(files_in_path %in% unique(X$sound.files))) stop2("sound file referenced in 'X' not found in 'path' (make sure you put the master sound file in the same folder than the re-recorded files)")
  
  # remove sound files
  if (!is.null(test.files)) {
    if (!all(test.files %in% files_in_path)) stop2("Not all 'test.files' were found in 'path'")
  } else test.files <- files_in_path[!files_in_path %in% unique(X$sound.files)] # remove master sound file
  
 if (length(template.rows) != 1 & length(test.files) != length(template.rows)) stop2("'template.rows' must be 1 or the same length than 'test.files'")
  
  # create a matrix that contains the selection/files to be cross-correlated
  if (length(template.rows) == 1) # in only 1 template, repeate it
  comp_mat <- matrix(c(rep(paste(X$sound.files[template.rows], X$selec[template.rows], sep = "-"), length(test.files)), test.files), ncol = 2) else 
    comp_mat <- matrix(c(rep(paste(X$sound.files[template.rows], X$selec[template.rows], sep = "-"), length(test.files)), rep((test.files), each = length(template.rows))), ncol = 2) 
  
  # patch to work on warbleR versions before 1.1.26
  wi <- warbleR::info_sound_files(path = path, parallel = 1, pb = FALSE, skip.error = TRUE, files = c(test.files, unique(X$sound.files))) 
  
  wi <- wi[wi$sound.files %in% unique(c(X$sound.files, test.files)), ]
  
  if (length(unique(wi$sample.rate)) > 1) 
    stop2("Not all sound files share the same sampling rate (check wave properties with warbleR::wav_info())")
  
  # steps for warbleR message
  if (pb){
  options("int_warbleR_steps" = c(current = 0, total = 2))
  
  on.exit(options("int_warbleR_steps" = c(current = 0, total = 0)), add = TRUE)
  }
  
  # run cross correlation
  xc <- warbleR::cross_correlation(X, compare.matrix = comp_mat, path = path, output = "list", pb = pb, ...)
  
  # find peaks
  if (pb) 
    write(file = "", x = paste0("running peak detection (step 2 out of 2):"))
  pks <- find_peaks_bRlh_int(xc.output = xc, max.peak = if(length(template.rows) == 1) TRUE else FALSE, pb = pb, cutoff = 0)
  
  if(length(template.rows) > 1)
  pks <- pks[stats::ave(x = -pks$score, as.factor(pks$sound.files), as.factor(pks$template), FUN = rank) <= 1, ]
  
  # rename sound file column
  try(names(pks)[names(pks) == "sound.files"] <- "test.files")
  return(pks)

  }
