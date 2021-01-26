# internal baRulho function, not to be called by users. It prepares X for comparing signals
# @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})
# last modification on jan-3-2020 (MAS)

prep_X_bRlo_int <- function(X, method = 1, parallel = 1, pb = TRUE) {
  
  # add sound file selec colums to X (weird column name so it does not overwrite user columns)
  X$TEMP....sgnl <- paste(X$sound.files, X$selec, sep = "-")
  
  # set pb options 
  pbapply::pboptions(type = ifelse(as.logical(pb), "timer", "none"))
  
  # set clusters for windows OS
  if (Sys.info()[1] == "Windows" & parallel > 1)
    cl <- parallel::makePSOCKcluster(getOption("cl.cores", parallel)) else cl <- parallel
  
  # add second column with names of the reference signals to be compare against
  X$reference <- pbapply::pbsapply(1:nrow(X), cl = cl, function(x, meth = method){
    
    # extract for single signal and order by distance
    Y <- as.data.frame(X[X$signal.type == X$signal.type[X$TEMP....sgnl == X$TEMP....sgnl[x]], , drop = FALSE])
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
