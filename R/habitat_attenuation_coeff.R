#' Estimate habitat_attenuation of sound pressure level
#' 
#' \code{habitat_attenuation} estimates atmospheric habitat_attenuation and atmospheric absorption.
#' @usage habitat_attenuation(frequency, temp = 20, rh = 60, pa = 101325, dist0, dist, hab.att.coef = 0.02)
#' @param frequency Numeric vector of length 1 with frequency (in Hertz).
#' @param temp Numeric vector of length 1 with frequency (in Celsius). Default is 20.
#' @param rh Numeric vector of length 1 with relative humidity (in percentage). Default is 60.
#' @param pa Numeric vector of length 1 with atmospheric (barometric) pressure in Pa (standard: 101325, default). Used for atmospheric habitat_attenuation.
#' @param dist Numeric vector of length 1 with distance (m) over which a sound propagates. 
#' @param dist0 Numeric vector of length 1 with distance (m) for the reference SPL. 
#' @param hab.att.coef habitat_attenuation coefficient of the habitat (in dB/kHz/m).
#' @return Returns the geometric, atmospheric and habitat habitat_attenuation (in dB) as well as the combined habitat_attenuation.   
#' @export
#' @name habitat_attenuation_coeff
#' @details Calculate the geometric, atmospheric and habitat attenuation and the overall expected attenuation (the sum of the other three) based on temperature, relative humidity, atmospheric pressure and sound frequency. Attenuation values are given in dB.
#' @examples
#' {
#' # measure habitat attenuation coefficient
#' habitat_attenuation_coeff(f = 2000, dist = 50, dist0 = 1)
#' }
#' 
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr}) 
#' @references {
#' Araya-Salas, M. (2020). baRulho: baRulho: quantifying habitat-induced degradation of (animal) acoustic signals in R. R package version 1.0.2
#' }

habitat_attenuation_coeff <-
  function(X,
           cores = 1,
           pb = TRUE,
           output = "est",
           hop.size = 11.6,
           wl = NULL,
           path = NULL,
           temp = 20,
           rh = 60,
           pa = 101325,
           subtract.bgn = TRUE,
           envelope = "abs",
           mar = NULL) {
    
    # error if no mar supplied when subtract.bgn
    if (is.null(mar) & subtract.bgn)
      stop2("'mar' must be supplied when 'subtract.bgn = TRUE'")
    
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
    
    # need mar if subtract.bgn TRUE
    if (subtract.bgn & is.null(mar))
      stop2("'mar' must be supplied if 'subtract.bgn = TRUE'")
    
    # get sampling rate
    sampling_rate <-
      warbleR::read_sound_file(
        X = X,
        index = 1,
        path = path,
        header = TRUE
      )$sample.rate
    
    # If cores is not numeric
    if (!is.numeric(cores))
      stop2("'cores' must be a numeric vector of length 1")
    if (any(!(cores %% 1 == 0), cores < 1))
      stop2("'cores' should be a positive integer")
    
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
    
    # check sound.id column
    if (is.null(X$sound.id))
      stop2("'X' must contain a 'sound.id' column")
    
    #check output
    if (!any(output %in% c("est", "data.frame", "list")))
      stop2("'output' must be 'est', 'data.frame' or 'list'")
    
    # add sound file selec column and names to X (weird column name so it does not overwrite user columns)
    if (pb)
      write(file = "",
            x = paste0("Preparing data for analysis (step 1 out of 3):"))
    
    X <-
      prep_X_bRlo_int(X,
                      method = 1,
                      cores = cores,
                      pb = pb)
    
    # set clusters for windows OS
    if (Sys.info()[1] == "Windows" & cores > 1)
      cl <-
      parallel::makePSOCKcluster(getOption("cl.cores", cores)) else
        cl <- cores
    
    # print message
    if (pb)
      write(file = "", x = "Computing habitat attenuation coefficient (step 2 out of 3):")
    
    # calculate all spectra apply function
    peak_freq_list <-
      warbleR:::pblapply_wrblr_int(
        pbar = pb,
        X = 1:nrow(X),
        cl = cl,
        FUN = function(y, wl)   {
          # load clip
          clp <- warbleR::read_sound_file(X = X,
                                          index = y,
                                          path = path)
          
          # calculate spectrum
          clp.spc <-
            seewave::spec(
              wave = clp,
              f = clp@samp.rate,
              plot = FALSE,
              wl = wl
            )
          
          # get peak frequency
          peak_freq <- seewave::fpeaks(clp.spc, nmax = 1, plot = FALSE)
          
            # get amplitude
            sigamp <- seewave::rms(seewave::env(clp, envt = envelope, plot = FALSE)) 
            
            # convert to dB
            signaldb <- 20 * log10(sigamp)
            
            # remove background SPL
            if (subtract.bgn){
              bg.noise <- read_sound_file(X, index = y, path = path, from = if (X$start[y] - mar < 0) 0 else X$start[y] - mar, to = X$start[y])
              
              noiseamp <- seewave::rms(seewave::env(bg.noise, f = bg.noise@samp.rate, envt = envelope, plot = FALSE))
              noisedb <- 20 * log10(noiseamp)
              
              # remove noise SPL from signal SPL
              signaldb <- warbleR:::lessdB(signal.noise = signaldb, noise = noisedb)
            } 
            
          # put output together
          output <- list(spl = signaldb, peakf = peak_freq[1, 1], dist = X$distance[y])
          return(output)
        }
      )
    
    # add sound file selec names to spectra
    names(peak_freq_list) <- X$TEMP....sgnl
    
    ## function to measure detection distance
    hab_att_coef_FUN <- function(x, temp, rh, pa, ...) {
      
      # get names of sound and reference
      sgnl <- X$TEMP....sgnl[x]
      rfnc <- X$reference[x]
      
      # if sounds are the same or the selection is noise return NA
      if (sgnl == rfnc | any(c(X$sound.id[X$TEMP....sgnl == sgnl], X$sound.id[X$reference == rfnc]) == "ambient"))
        hab_att_coef <- NA else {
          
          # extract data for sound and model
          sgnl.spl <- peak_freq_list[[which(names(peak_freq_list) == sgnl)]]$spl
          rfnc.spl <- peak_freq_list[[which(names(peak_freq_list) == rfnc)]]$spl
          rfnc.pkf <-
          peak_freq_list[[which(names(peak_freq_list) == rfnc)]]$peakf
          dist0 <- peak_freq_list[[which(names(peak_freq_list) == rfnc)]]$dist
          dist <- peak_freq_list[[which(names(peak_freq_list) == sgnl)]]$dist
          
          # get expected attenuations
          exp_att <- attenuation(frequency = rfnc.pkf, temp = temp, rh = rh, pa = pa, dist0 = dist0, dist = dist, hab.att.coef = 0)

          # get observed attenuation
          obs_att <- rfnc.spl - sgnl.spl
                 
          # substract geometric and atmospheric attenuations to observed attenuation to get habitat attenuation
          hab_att <- obs_att - (exp_att$geometric.attenuation + exp_att$atmopheric.attenuation) 
             
          hab_att_coef <- (hab_att * 1000) / (rfnc.pkf * (dist-dist0))
        }
      return(hab_att_coef)
    }
    
    if (pb)
      write(file = "", x = "Computing detection distance (step 3 out of 3):")
    
    # get habitat attenuatio coefficient
    # calculate all spectra apply function
    X$habitat.attenuation.coeff <-
      pbapply::pbsapply(
        X = 1:nrow(X),
        cl = cl,
        FUN = function(x,
                       tempe = temp,
                       rhe = rh,
                       pae = pa,
                       ...)   {
          
          hab_att_coef_FUN(x, tempe, rhe, pae, ...)
        }
      )
    
    # remove temporal columns
    X$TEMP....sgnl <- NULL
    
    # return data frame
    if (output == "data.frame")
      X <- as.data.frame(X) else
        attributes(X)$call <-
      base::match.call() # fix call attribute
    
    return(X)
}


##### modified from http://forum.studiotips.com/viewtopic.php?t=158
## and https://scikit-maad.github.io/generated/maad.spl.attenuation_dB.html#maad.spl.habitat_attenuation_dB