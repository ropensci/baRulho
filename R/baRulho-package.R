#' baRulho: quantifying habitat-induced acoustic signal degradation
#'
#' `baRulho` is a package intended to quantify habitat-induced degradation of (animal) acoustic signals.
#'
#' The main features of the package are:
#'   \itemize{
#'   \item Loops to apply tasks through acoustic signals referenced in an extended selection table
#'   \item The comparison of playback signals re-recorded at different distances 
#'   }
#' Most functions allow the parallelization of tasks, which distributes the tasks among several processors to improve computational efficiency.
#'   
#' @import pbapply
#' @import tuneR
#' @import seewave
#' @import graphics
#' @import grDevices
#' @import utils
#' @import parallel
#' @import warbleR
#' @importFrom stats cor aggregate approx
#' 
#' @author Marcelo Araya-Salas
#'   
#'   Maintainer: Marcelo Araya-Salas (\email{marceloa27@@gmail.com})
#'   
#' @docType package
#' @name baRulho
#' @details License: GPL (>= 2)  
NULL
#> NULL 
#'
