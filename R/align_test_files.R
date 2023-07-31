#' Align test sound files
#' 
#' \code{align_test_files} aligns test (re-recorded) sound files.
#' @usage align_test_files(X, Y, output = "est", path = NULL, 
#' by.song = TRUE, marker = NULL, cores = 1, pb = TRUE, ...)
#' @param X object of class 'data.frame', 'selection_table' or 'extended_selection_table' (the last 2 classes are created by the function \code{\link[warbleR]{selection_table}} from the warbleR package). This should be the same data than that was used for finding the position of markers in \code{\link{find_markers}}. It should also contain a 'sound.id' column that will be used to label re-recorded sounds according to their counterpart in the master sound file.
#' @param Y object of class 'data.frame' with the output of \code{\link{find_markers}}. This object contains the position of markers in the re-recorded sound files. If more than one marker is supplied for a sound file only the one with the highest correlation score ('scores' column in 'X') is used.
#' @param output Character vector of length 1 to determine if an extended selection table ('est', default) or a data.frame ("data.frame").
#' @param path Character string containing the directory path where test (re-recorded) sound files are found. 
#' @param by.song Logical argument to indicate if the extended selection table should be created by song (see 'by.song' \code{\link[warbleR]{selection_table}} argument). Default is \code{TRUE}.
#' @param marker Character string to define whether a "start" or "end" marker would be used for aligning re-recorded sound files. Default is \code{NULL}. DEPRECATED.
#' @param cores Numeric vector of length 1. Controls whether parallel computing is applied by specifying the number of cores to be used. Default is 1 (i.e. no parallel computing).
#' @param pb Logical argument to control if progress bar is shown. Default is \code{TRUE}.
#' @param ...	Additional arguments to be passed to \code{\link[warbleR]{selection_table}} for customizing extended selection table.
#' @return An extended selection table with the aligned sounds from test (re-recorded) sound files.
#' @export
#' @name align_test_files
#' @details The function aligns sounds found in re-recorded sound files (referenced in 'Y') according to a master sound file (referenced in 'X'). The function outputs an 'extended selection table'by default.
#' @seealso \code{\link{realign_test_sounds}}; \code{\link{find_markers}}
#' @examples
#' \dontrun{
#' # set temporary directory
#' td <- tempdir()  
#' 
#' # load example data from warbleR
#' data(list = c("Phae.long1", "Phae.long2", "Phae.long3", "Phae.long4", "lbh_selec_table"))
#' 
#' # save sound files to temporary folder
#' writeWave(Phae.long1, file.path(td, "Phae.long1.wav"))
#' writeWave(Phae.long2, file.path(td, "Phae.long2.wav"))
#' writeWave(Phae.long3, file.path(td, "Phae.long3.wav"))
#' writeWave(Phae.long4, file.path(td, "Phae.long4.wav"))
#' 
#' # make an extended selection table
#' est <- warbleR::selection_table(X = lbh_selec_table, extended = TRUE, 
#' confirm.extended = FALSE, path = td, pb = FALSE)
#' 
#' # create master sound file
#' master.sf <- master_sound_file(X = est, file.name = "example_master", 
#' dest.path = td, gap.duration = 0.3)
#' 
#' # read master
#' exmp.master <- readWave(file.path(td, "example_master.wav"))
#' 
#' # add 1 s silence
#' exmp.test1 <- addsilw(wave = exmp.master, at = "start", d = 1, output = "Wave", 
#' f = exmp.master@samp.rate)
#' 
#' exmp.test2 <- addsilw(wave = exmp.master, at = "start", d = 2, output = "Wave", 
#' f = exmp.master@samp.rate)
#' 
#' # create noise
#' ns <- noisew(f = exmp.master@samp.rate, d = duration(exmp.test2) + 1, 
#' output = "Wave")
#' 
#' # make noise exactly the same length and add noise to 2 examples
#' exmp.test1@left <- exmp.test1@left + (ns@left[1:length(exmp.test1@left)] * 500)
#' exmp.test2@left <- exmp.test2@left + (ns@left[1:length(exmp.test2@left)] * 500)
#' 
#' # normalize before saving
#' exmp.test1 <- normalize(exmp.test1, unit = "16")
#' exmp.test2 <- normalize(exmp.test2, unit = "16")
#' 
#' # save examples
#' writeWave(object = exmp.test1, filename = file.path(td, "example_test1.wav"), 
#' extensible = FALSE)
#' 
#' writeWave(object = exmp.test2, filename = file.path(td, "example_test2.wav"), 
#' extensible = FALSE)
#' 
#' # find  tempaltes
#' found.templts <- find_markers(X = master.sf, 
#' template.rows = which(master.sf$orig.sound.file == "start_marker"), 
#' test.files = c("example_test1.wav", "example_test2.wav"), path = td, pb = FALSE)
#' 
#' # align sounds and output extended selection table
#' alg.tests <- align_test_files(X =  master.sf, Y = found.templts, path = td, pb = FALSE)
#' }
#' 
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr}) 
#' @references {
#' Araya-Salas, M. (2020). baRulho: baRulho: quantifying habitat-induced degradation of (animal) acoustic signals in R. R package version 1.0.2
#' }
#last modification on dec-26-2019 (MAS)

align_test_files <- function(X, Y, output = "est", path = NULL, by.song = TRUE, marker = NULL, cores = 1, pb = TRUE, ...){
  
  # deprecated message
  if (!is.null(marker))
    stop2("'marker' has been deprecated")

  # If cores is not numeric
  if (!is.numeric(cores)) stop2("'cores' must be a numeric vector of length 1") 
  if (any(!(cores %% 1 == 0), cores < 1)) stop2("'cores' should be a positive integer")
  
  #check output
  if (!any(output %in% c("est", "data.frame"))) stop2("'output' must be either 'est' or 'data.frame'")  
  
  # check sound.id column 
  if (is.null(X$sound.id)) stop2("'X' must contain a 'sound.id' column")
  
  
  # if more than one marker per test files then keep only the marker with the highest score
  if (any(table(Y$test.files) > 1))
    Y <- Y[stats::ave(x = -Y$scores, as.factor(Y$test.files), FUN = rank) <= 1, ]
  
  # set clusters for windows OS
  if (Sys.info()[1] == "Windows" & cores > 1)
    cl <- parallel::makePSOCKcluster(getOption("cl.cores", cores)) else cl <- cores
    
  # align each file
 out <- warbleR:::pblapply_wrblr_int(pbar = pb, X = 1:nrow(Y), cl = cl, FUN = function(y) {

   # compute start and end as the difference in relation to the template position in the master sound file 
   start <- X$start + (Y$start[y] - X$start[X$sound.id == Y$template[y]])
   end <- X$end + (Y$start[y] - X$start[X$sound.id == Y$template[y]])
   
    # make data frame
    W <- data.frame(sound.files = Y$test.files[y], selec = 1:length(start), start, end, bottom.freq = X$bottom.freq, top.freq = X$top.freq, sound.id = X$sound.id, template = Y$template[y])
    
    return(W)
  })
  
  # put data frames togheter
  sync.sls <- do.call(rbind, out)
  
  # check if any selection exceeds length of recordings
  #if(exists("wav_dur"))
    wvdr <- wavdur(path = path, files = unique(sync.sls$sound.files)) #else
  #wvdr <- warbleR::duration_wavs(path = path)
  
  # add duration to data frame
  sync.sls <- merge(sync.sls, wvdr)
  
  # start empty vector to add name of problematic files
  problematic_files <- character()
  
  if (any(sync.sls$end > sync.sls$duration)) {
    write(file = "", x = paste(sum(sync.sls$end > sync.sls$duration), "selection(s) exceeded sound file length and were removed (run .Options$baRulho to see which test files were involved)"))
    
    problematic_files <- append(problematic_files, unique(sync.sls$sound.files[sync.sls$end > sync.sls$duration]))
    
    # remove exceeding selections
    sync.sls <- sync.sls[!sync.sls$end > sync.sls$duration, ]
  }

  if (any(sync.sls$start < 0)) {
    write(file = "", x = paste(sum(sync.sls$start < 0), "selection(s) were absent at the start of the files (negative start values) and were removed"))
    
    problematic_files <- append(problematic_files, unique(sync.sls$sound.files[sync.sls$start >= 0]))
    
    # remove exceeding selections
    sync.sls <- sync.sls[sync.sls$start >= 0, ]
  }
  
  on.exit(options(baRulho = list(files_to_check_align_test_files = unique(problematic_files))))
  
  # remove duration column and template
  sync.sls$duration <- NULL
  # sync.sls$sound.id <- sync.sls$template
  # sync.sls$template <- NULL
    
  if (output == "est")
  {
    if (by.song) # if by song add a numeric column to represent sound files
    {
      sync.sls$song <- as.numeric(as.factor(sync.sls$sound.files))
      by.song <- "song"
      } else 
        by.song <- NULL # rewrite by song as null
        
    sync.sls <- selection_table(sync.sls, extended = TRUE, confirm.extended = FALSE, path = path, by.song = by.song, pb = pb, verbose = pb, ...)
    
    # fix call attribute
    attributes(sync.sls)$call <- base::match.call()
    
    }
  
  return(sync.sls)
}
