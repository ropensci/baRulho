#' Align test sound files
#' 
#' \code{align_test_files} aligns test (re-recorded) sound files.
#' @usage align_test_files(X, Y, output = "est", marker = "start_marker", path = NULL, 
#' by.song = TRUE, ...)
#' @param X object of class 'extended_selection_table' created by the function \code{\link[warbleR]{selection_table}} from the warbleR package. This should be the same data than that used for aligning signals using \code{\link[warbleR]{xcorr}}.
#' @param Y object of class 'data.frame' with the output of \code{\link[warbleR]{find_peaks}} from the warbleR package. 
#' @param output Character vector of length 1 to determine if an extended selection table ('est', default) or a data.frame ("data.frame").
#' @param marker Character string of length 1 with the name of the marker that was used for aligning signals. Default is 'start_marker'.
#' @param path Character string containing the directory path where test (re-recorded) sound files are found. 
#' @param by.song Logical argument to indicate if the extended selection table should be created by song (see 'by.song' \code{\link[warbleR]{selection_table}} argument). Default is \code{TRUE}.
#' @param ...	Additional arguments to be passed to \code{\link[warbleR]{selection_table}} for customizing extended selection table.
#' @return An extended selection table with the aligned signals from test (re-recorded) sound files.
#' @export
#' @name align_test_files
#' @details The function takes the output of \code{\link[warbleR]{find_peaks}} ('Y') and aligns signals found in re-recorded sound files according to a master sound file referenced in 'X'. The function outputs a 'extended selection table'.
#' @seealso \code{\link{spcc_align}}
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
#' # create a matrix that contains the selection/files to be cross-correlated
#' comp_mat <- matrix(c(rep(paste(master.sf$sound.files[1], 
#' master.sf$selec[1], sep = "-"), 2), "example_test1.wav", 
#' "example_test2.wav"), nrow = 2)
#' 
#' # run cross correlation
#' xc <- warbleR::xcorr(master.sf, compare.matrix = comp_mat, wl = 100, 
#' ovlp = 50, path = td, output = "list", pb = FALSE)
#' 
#' # find peaks
#' pks <- warbleR::find_peaks(xc.output = xc, max.peak = TRUE, path = td)
#' 
#' # align signals and output extended selection table
#' alg.tests <- align_test_files(X =  master.sf, Y = pks, path = td)
#' }
#' 
#' @author Marcelo Araya-Salas (\email{marceloa27@@gmail.com}) 
#' @references {
#' Araya-Salas, M. (2020). baRulho: baRulho: quantifying habitat-induced degradation of (animal) acoustic signals in R. R package version 1.0.0
#' }
#last modification on dec-26-2019 (MAS)

align_test_files <- function(X, Y, output = "est", marker = "start_marker", path = NULL, by.song = TRUE, ...){
  
  #check output
  if (!any(output %in% c("est", "data.frame"))) stop("'output' must be either 'est' or 'data.frame'")  
  
  # align each file
 out <- lapply(1:nrow(Y), function(x) {
    
   # name of master sound file
   mstr.nm <- substr(Y$template[x], start = 0, regexpr("\\-[^\\-]*$", Y$template[x]) - 1)
   
    # start on new recording
    start  <- X$start - X$start[X$sound.files == mstr.nm & X$orig.sound.file == marker] + Y$start[x]
    
    # end on new recording
    end  <- X$end - X$end[X$sound.files == mstr.nm & X$orig.sound.file == marker] + Y$end[x]
    
    # make data frame
    W <- data.frame(sound.files = Y$sound.files[x], selec = 1:length(start), start, end, bottom.freq = X$bottom.freq[1:length(start)], top.freq = X$top.freq[1:length(start)], template = X$orig.sound.file[1:length(start)])
    
    return(W)
  })
  
  # put data frames togheter
  sync.sls <- do.call(rbind, out)
  
  
  if (output == "est")
  {
    if (by.song) # if by song add a numeric column to represent sound files
    {
      sync.sls$song <- as.numeric(as.factor(sync.sls$sound.files))
      by.song <- "song"
      } else 
        by.song <- NULL # rewrite by song as null
        
    sync.sls <- selection_table(sync.sls, extended = TRUE, confirm.extended = FALSE, path = path, by.song = by.song, ...)
    }
  
  return(sync.sls)
}
