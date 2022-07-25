#' Estimate atmospheric attenuation and absorption of sound 
#' 
#' \code{attenuation} estimates atmospheric attenuation and atmospheric absorption.
#' @usage soud_attenuation(f, temp, RH, p = 101325, 
#' formula = 1, spi = NULL, dist = NULL)
#' @param f numeric vector of length 1 with frequency (in Hertz).
#' @param temp numeric vector of length 1 with frequency (in Celsius).
#' @param RH numeric vector of length 1 with relative humidity
#' @param p numeric vector of length 1 with ambient pressure in Pa (standard: 101325, default).
#' @param formula 1 = Bazley 1976, 2 = Rossing 2007 (p. 116, see details). 
#' @param spi numeric vector of length 1 with the initial sound pressure in Pa. Required for calculating atmospheric absorption. Default is \code{NULL}. 
#' @param dist numeric vector of length 1 with distance (m) over which a sound propagates. Required for calculating atmospheric absorption. Default is \code{NULL}. 
#' @return Returns atmospheric attenuation (in dB/m) of sound based on supplied parameters. If 'spi' and 'dist' are supplied the function also returns atmospheric absorption (in dB).  
#' @export
#' @name soud_attenuation
#' @details Calculate the atmospheric attenuation based on temperature, relative humidity, pressure and sound frequency. The function can applied to formulae based on:
#' \itemize{
#' \item \code{1}: default. As used by Bazley (1976), Sound absorption in air at frequencies up to 100 kHz. NPL acoustics report Ac 74. 
#' \item \code{2}: as used by Rossing (2007), Handbook of Acoustics, Springer.
#' }
#' If 'spi' and 'dist' are supplied the function also returns the atmospheric absorption (in dB).
#' @examples
#' {
#' # measure atmospheric attenuation formula 1
#' soud_attenuation(f = 20000, temp = 20, RH = 90, p = 88000, formula = 1)
#' }
#' 
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr}) 
#' @references {
#' Araya-Salas, M. (2020). baRulho: baRulho: quantifying habitat-induced degradation of (animal) acoustic signals in R. R package version 1.0.2
#' }
#last modification on nov-01-2019 (MAS)

### NOT SURE IF ABSORPTION WORKS 
soud_attenuation <- function(f, temp, RH, p = 101325, formula = 1, spi = NULL, dist = NULL){
  
  
  
}

##### from http://forum.studiotips.com/viewtopic.php?t=158
