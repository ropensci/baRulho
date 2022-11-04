#' Estimate atmospheric attenuation
#' 
#' \code{atmospheric_attenuation} estimates atmospheric attenuation and atmospheric absorption.
#' @usage atmospheric_attenuation(f, temp, RH, p = 101325, formula = NULL, spi = NULL, dist)
#' @param f numeric vector of length 1 with frequency (in Hertz).
#' @param temp numeric vector of length 1 with frequency (in Celsius).
#' @param RH numeric vector of length 1 with relative humidity
#' @param p numeric vector of length 1 with ambient pressure in Pa (standard: 101325, default).
#' @param formula DEPRECATED.
#' @param spi numeric vector of length 1 with the initial sound pressure in Pa. Default is \code{NULL}. DEPRECATED.
#' @param dist numeric vector of length 1 with distance (m) over which a sound propagates. 
#' @return Returns atmospheric attenuation (in dB) of sound based on supplied parameters.   
#' @export
#' @name atmospheric_attenuation
#' @details Calculate the atmospheric attenuation based on temperature, relative humidity, pressure and sound frequency.
#' @examples
#' {
#' # measure atmospheric attenuation formula 1
#' atmospheric_attenuation(f = 20000, temp = 20, RH = 90, p = 88000, dist = 30)
#' }
#' 
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr}) 
#' @references {
#' Araya-Salas, M. (2020). baRulho: baRulho: quantifying habitat-induced degradation of (animal) acoustic signals in R. R package version 1.0.2
#' }
#last modification on oct-28-2021 (MAS)

atmospheric_attenuation <- function(f, temp, RH, p = 101325, formula = NULL, spi = NULL, dist){
  pr <- 101325
  To1 <- 273.16
  To <- 293.15
  t <- temp + 273.15

  psat <- pr * 10 **(-6.8346 * (To1 / t) ** 1.261 + 4.6151)
  
  h <- RH * (psat / p)
  
  fr0 <- (p / pr) * (24 + 4.04e4 * h * ((0.02 + h) / (0.391 + h)))
  
  frN <- (p / pr) * sqrt(t / To) * (9 + 280 * h * exp(-4.170 * ((t / To) ** (-1 / 3) - 1)))
  
  
  z = 0.1068 * exp(-3352 / t) / (frN + f ** 2 / frN)
  
  y = (t / To) ** (-5/2) * (0.01275 * exp(-2239.1/t) * 1/(fr0 + f ** 2 / fr0) + z)
  
  
  AA_coef <- 8.686 * f ** 2 * ((1.84e-11 * 1/(p / pr) * sqrt(t / To)) + y)
  
  AA <- AA_coef * (dist)
  
  return(AA)  
  
}

##### from http://forum.studiotips.com/viewtopic.php?t=158
