#' Measure frequency spectral correlation
#'
#' \code{spectrum_correlation} measures frequency spectrum correlation of sounds referenced in an extended selection table.
#' @usage spectrum_correlation(X, parallel = 1, cores = 1, pb = TRUE, method = 1,
#' cor.method = "pearson", output = "est",
#' hop.size = 11.6, wl = NULL, ovlp = 70, path = NULL)
#' @param output Character vector of length 1 to determine if an extended selection table ('est', default) or a data frame ('data.frame').ion(X, parallel = 1, cores = 1, pb = TRUE, method = 1,
#' cor.method = "pearson", output = "est",
#' hop.size = 11.6, wl = NULL, ovlp = 70)
#' @param X Object of class 'data.frame', 'selection_table' or 'extended_selection_table' (the last 2 classes are created by the function \code{\link[warbleR]{selection_table}} from the warbleR package) with the reference to the sounds in the master sound file. Must contain the following columns: 1) "sound.files": name of the .wav files, 2) "selec": unique selection identifier (within a sound file), 3) "start": start time and 4) "end": end time of selections, 5)  "bottom.freq": low frequency for bandpass, 6) "top.freq": high frequency for bandpass and 7) "sound.id": ID of sounds used to identify counterparts across distances. Each sound must have a unique ID within a distance.
#' @param parallel DEPRECATED. Use 'cores' instead.
#' @param cores Numeric vector of length 1. Controls whether parallel computing is applied by specifying the number of cores to be used. Default is 1 (i.e. no parallel computing).
#' If \code{NULL} (default) then the current working directory is used.
#' @param pb Logical argument to control if progress bar is shown. Default is \code{TRUE}.
#' @param method Numeric vector of length 1 to indicate the 'experimental design' to measure frequency spectrum correlation. Two methods are available:
#' \itemize{
#' \item \code{1}: compare all sounds with their counterpart that was recorded at the closest distance to source (e.g. compare a sound recorded at 5m, 10m and 15m with its counterpart recorded at 1m). This is the default method.
#' \item \code{2}: compare all sounds with their counterpart recorded at the distance immediately before (e.g. a sound recorded at 10m compared with the same sound recorded at 5m, then sound recorded at 15m compared with same sound recorded at 10m and so on).
#' }
#' @param cor.method Character string indicating the correlation coefficient to be applied ("pearson", "spearman", or "kendall", see \code{\link[stats]{cor}}).
#' @param output Character vector of length 1 to determine if an extended selection table ('est', default) or a data frame ('data.frame').
#' @param hop.size A numeric vector of length 1 specifying the time window duration (in ms). Default is 11.6 ms, which is equivalent to 512 wl for a 44.1 kHz sampling rate. Ignored if 'wl' is supplied.
#' @param wl A numeric vector of length 1 specifying the window length of the spectrogram, default
#' is NULL. If supplied, 'hop.size' is ignored.
#' @param ovlp Numeric vector of length 1 specifying the percent overlap between two
#'   consecutive windows, as in \code{\link[seewave]{spectro}}. Default is 70.
#' @param path Character string containing the directory path where the sound files are found. Only needed when 'X' is not an extended selection table.
#' @return A data frame, or extended selection table similar to input data (depending on argument 'output'), but also includes a new column ('spectrum.correlation')
#' with the calculated frequency spectrum correlation coefficients.
#' @export
#' @name spectrum_correlation
#' @details spectral correlation measures the similarity of two sounds in the frequency domain. The function measures the spectral correlation coefficients of sounds in which a reference playback has been re-recorded at increasing distances. Values range from 1 (identical frequency spectrum, i.e. no degradation) to 0. The 'sound.id' column must be used to indicate the function to only compare sounds belonging to the same category (e.g. song-types). The function will then compare each sound to the corresponding reference sound. Two methods for calculating spectral correlation are provided (see 'method' argument). Use \code{\link{spectrum_blur_ratio}} to get spectra for plotting.
#' @examples
#' {
#' # load example data
#' data("playback_est")
#'
#' # remove ambient selections
#' pe <- playback_est[playback_est$sound.id != "ambient", ]
#'
#' # method 1
#'spectrum_correlation(X = pe)
#'
#' # method 2
#' spectrum_correlation(X = pe, method = 2)
#' }
#'
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})
#' @seealso \code{\link{envelope_correlation}}, \code{\link{spectrum_blur_ratio}}
#' @references {
#' Araya-Salas, M. (2020). baRulho: baRulho: quantifying habitat-induced degradation of (animal) acoustic signals in R. R package version 1.0.2
#'
#' Apol, C.A., Sturdy, C.B. & Proppe, D.S. (2017). Seasonal variability in habitat structure may have shaped acoustic signals and repertoires in the black-capped and boreal chickadees. Evol Ecol. 32:57-74.
#' }
#last modification on nov-01-2019 (MAS)

spectrum_correlation <-
  function(X,
           parallel = 1,
           cores = 1,
           pb = TRUE,
           method = 1,
           cor.method = "pearson",
           output = "est",
           hop.size = 11.6,
           wl = NULL,
           ovlp = 70,
           path = NULL) {
    # set path if not provided
    if (is.null(path))
      path <- getwd() else
      if (!dir.exists(path))
        stop2("'path' provided does not exist")
    
    # must have the same sampling rate
    if (is_extended_selection_table(X)) {
      if (length(unique(attr(X, "check.results")$sample.rate)) > 1)
        stop2(
          "all wave objects in the extended selection table must have the same sampling rate (they can be homogenized using warbleR::resample_est())"
        )
    } else
      print("assuming all sound files have the same sampling rate")
    
    # If cores is not numeric
    if (!is.numeric(cores))
      stop2("'cores' must be a numeric vector of length 1")
    if (any(!(cores %% 1 == 0), cores < 1))
      stop2("'cores' should be a positive integer")
    
    #check output
    if (!any(output %in% c("est", "data.frame")))
      stop2("'output' must be 'est' or 'data.frame'")
    
    # hopsize
    if (!is.numeric(hop.size) |
        hop.size < 0)
      stop2("'hop.size' must be a positive number")
    
    # adjust wl based on hope.size
    if (is.null(wl))
      wl <-
        round(read_sound_file(
          X,
          index = 1,
          header = TRUE,
          path = path
        )$sample.rate * hop.size  / 1000,
        0)
    
    # make wl even if odd
    if (!(wl %% 2) == 0)
      wl <- wl + 1
    
    # If method is not numeric
    if (!is.character(cor.method))
      stop2("'cor.method' must be a character vector of length 1")
    if (!any(cor.method %in%  c("pearson", "kendall", "spearman")))
      stop2("'method' must be either  'pearson', 'kendall' or 'spearman'")
    
    # check sound.id column
    if (is.null(X$sound.id))
      stop2("'X' must contain a 'sound.id' column")
    
    # add sound file selec column and names to X (weird column name so it does not overwrite user columns)
    if (pb)
      write(file = "",
            x = paste0("Preparing data for analysis (step 1 out of 3):"))
    
    X <-
      prep_X_bRlo_int(X,
                      method = method,
                      cores = cores,
                      pb = pb)
    
    # set clusters for windows OS
    if (Sys.info()[1] == "Windows" & cores > 1)
      cl <-
      parallel::makePSOCKcluster(getOption("cl.cores", cores)) else
      cl <- cores
    
    if (pb)
      write(file = "", x = "Calculating frequency spectrums (step 2 out of 3):")
    
    # calculate all spectra apply function
    spcs <-
      warbleR:::pblapply_wrblr_int(
        pbar = pb,
        X = 1:nrow(X),
        cl = cl,
        FUN = function(y, wle = wl, ovl = ovlp) {
          # load clip
          clp <- warbleR::read_sound_file(X = X,
                                          index = y,
                                          path = path)
          
          # mean spec
          mspc <-
            meanspec(
              wave = clp,
              f = clp@samp.rate,
              plot = FALSE,
              wl = wle,
              ovlp = ovl
            )
          
          return(mspc)
        }
      )
    
    # add sound file selec names to envelopes (weird column name so it does not overwrite user columns)
    names(spcs) <- X$TEMP....sgnl
    
    # function to measure envelope correlation
    # y and z are the sound.files+selec names of the sounds and reference sound (model)
    spctr_cor_FUN <- function(y, z) {
      # if names are the same return NA
      if (y == z)
        cor.spctr <- NA else {
        # extract envelope for sound and model
        sgnl.spctr <- spcs[[which(names(spcs) == y)]]
        mdl.spctr <- spcs[[which(names(spcs) == z)]]
        
        
        ### filter to freq range of sounds and remove freq column
        # get range as lowest bottom and highest top
        frng <-
          c(min(X$bottom.freq[X$TEMP....sgnl %in% c(y, z)]), max(X$top.freq[X$TEMP....sgnl %in% c(y, z)]))
        sgnl.spctr <-
          sgnl.spctr[sgnl.spctr[, 1] > frng[1] &
                       sgnl.spctr[, 1] < frng[2], 2]
        mdl.spctr <-
          mdl.spctr[mdl.spctr[, 1] > frng[1] & mdl.spctr[, 1] < frng[2], 2]
        
        # get correlation assuming they have same length
        cor.spctr <- cor(sgnl.spctr, mdl.spctr, method = cor.method)
      }
      
      return(cor.spctr)
    }
    
    # make a data frame with 2 columns with names of the sounds to be compare
    X$reference <- sapply(1:nrow(X), function(x, meth = method) {
      # extract for single sound and order by distance
      Y <-
        as.data.frame(X[X$sound.id == X$sound.id[X$TEMP....sgnl == X$TEMP....sgnl[x]], , drop = FALSE])
      Y <- Y[order(Y$distance),]
      
      # method 1 compare to closest distance to source
      if (meth == 1)
        z <- Y$TEMP....sgnl[which.min(Y$distance)] else
        # if method 2
        # if not the first row then the previous row
        if (Y$TEMP....sgnl[1] != X$TEMP....sgnl[x])
          z <- X$TEMP....sgnl[x - 1] else
        # else the first row
        z <- Y$TEMP....sgnl[1]
      
      return(z)
    })
    
    if (pb)
      write(file = "", x = "Calculating spectrum correlations (step 3 out of 3):")
    
    # calculate all envelops apply function
    X$spectrum.correlation <-
      pbapply::pbsapply(
        X = 1:nrow(X),
        cl = cl,
        FUN = function(x) {
          spctr_cor_FUN(y = X$TEMP....sgnl[x], z = X$reference[x])
        }
      )
    
    # remove temporal columns
    X$TEMP....sgnl <- NULL
    
    # return data frame
    if (output == "data.frame")
      X <- as.data.frame(X) else
      attributes(X)$call <- base::match.call() # fix call attribute
    
    return(X)
  }
