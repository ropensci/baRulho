#' Estimate attenuation of sound pressure level
#' 
#' \code{attenuation} estimates atmospheric attenuation and atmospheric absorption.
#' @usage attenuation(f, temp = 20, RH = 60, pa = 101325, dist0, dist, att.coef = 0.02)
#' @param f numeric vector of length 1 with frequency (in Hertz).
#' @param temp numeric vector of length 1 with frequency (in Celsius). Default is 20.
#' @param RH numeric vector of length 1 with relative humidity (in percentage). Default is 60.
#' @param pa numeric vector of length 1 with ambient pressure in Pa (standard: 101325, default).
#' @param dist numeric vector of length 1 with distance (m) over which a sound propagates. 
#' @param dist0 numeric vector of length 1 with distance (m) for the reference SPL. 
#' @param att.coef attenuation coefficient of the habitat (in dB/kHz/m).
#' @return Returns the geometric, atmospheric and habitat attenuation (in dB) as well as the combined attenuation.   
#' @export
#' @name attenuation
#' @details Calculate the atmospheric attenuation based on temperature, relative humidity, pressure and sound frequency.
#' @examples
#' {
#' # measure attenuation
#' attenuation(f = 2000, dist = 50, dist0 = 1)
#' }
#' 
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr}) 
#' @references {
#' Araya-Salas, M. (2020). baRulho: baRulho: quantifying habitat-induced degradation of (animal) acoustic signals in R. R package version 1.0.2
#' }
#last modification on oct-28-2021 (MAS)

attenuation <- function(f, temp = 20, RH = 60, pa = 101325, dist0, dist, att.coef = 0.02){
    
    # atmospheric attenuation
    pr <- 101325
    To1 <- 273.16
    To <- 293.15
    temp <- temp + 273.15
    
    psat <- pr*10**(-6.8346*(To1/temp)**1.261 + 4.6151)
    
    h <- RH*(psat/pa)
    
    fr0 <- (pa/pr)*(24+4.04e4*h*((0.02+h)/(0.391+h)))
    
    frN <- (pa/pr)*sqrt(temp/To)*(9+280*h*exp(-4.170*((temp/To)**(-1/3)-1)))
    
    
    z <- 0.1068*exp(-3352/temp)/(frN+f**2/frN)
    
    y <- (temp/To)**(-5/2) * (0.01275 * exp(-2239.1/temp) * 1/(fr0+f**2/fr0) + z)
    
    
    atm_att_coef <- 8.686*f**2*((1.84e-11*1/(pa/pr)*sqrt(temp/To))+y)
    atm_att <- atm_att_coef * (dist-dist0)
    
    geom_att <- -20*log10(dist0/dist)
    
    hab_att <- (att.coef*f/1000) * (dist-dist0)
    
    total_att <- geom_att + atm_att + hab_att
    
    outdf <- data.frame(frequency = f, dist = dist, geometric.attenuation = geom_att, atmopheric.attenuation = atm_att, habitat.attenuation = hab_att, combined.attenuation = total_att)
    
    return(outdf)
  }

##### from http://forum.studiotips.com/viewtopic.php?t=158
## and https://scikit-maad.github.io/generated/maad.spl.attenuation_dB.html#maad.spl.attenuation_dB