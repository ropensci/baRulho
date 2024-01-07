#' Estimate attenuation of sound pressure level
#'
#' \code{attenuation} estimates atmospheric attenuation and atmospheric absorption.
#' @usage attenuation(frequency, temp = 20, rh = 60, pa = 101325, dist0, dist, hab.att.coef = 0.02)
#' @param frequency Numeric vector of length 1 with frequency (in Hertz).
#' @param temp Numeric vector of length 1 with frequency (in Celsius). Default is 20.
#' @param rh Numeric vector of length 1 with relative humidity (in percentage). Default is 60.
#' @param pa Numeric vector of length 1 with atmospheric (barometric) pressure in Pa (standard: 101325, default). Used for atmospheric attenuation.
#' @param dist Numeric vector of length 1 with distance (m) over which a sound propagates.
#' @param dist0 Numeric vector of length 1 with distance (m) for the reference SPL.
#' @param hab.att.coef Attenuation coefficient of the habitat (in dB/kHz/m).
#' @return Returns the geometric, atmospheric and habitat attenuation (in dB) as well as the combined attenuation.
#' @export
#' @name attenuation
#' @details Calculate the geometric, atmospheric and habitat attenuation and the overall expected attenuation (the sum of the other three) based on temperature, relative humidity, atmospheric pressure and sound frequency. Attenuation values are given in dB. Modified from http://www.sengpielaudio.com
## and https://scikit-maad.github.io/generated/maad.spl.attenuation_dB.html#maad.spl.attenuation_dB.
#' @examples {
#'   # measure attenuation
#'   attenuation(frequency = 2000, dist = 50, dist0 = 1)
#' }
#' @family miscellaneous
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})
#' @references {
#' Araya-Salas, M. (2020). baRulho: baRulho: quantifying degradation of (animal) acoustic signals in R. R package version 1.0.2
#' }

attenuation <-
  function(frequency,
           temp = 20,
           rh = 60,
           pa = 101325,
           dist0,
           dist,
           hab.att.coef = 0.02) {
    # check arguments
    arguments <- as.list(base::match.call())
    
    # add objects to argument names
    for (i in names(arguments)[-1]) {
      arguments[[i]] <- get(i)
    }
    
    # check each arguments
    check_results <-
      .check_arguments(fun = arguments[[1]], args = arguments)
    
    # report errors
    .report_assertions(check_results)
    # atmospheric attenuation
    pr <- 101325
    To1 <- 273.16
    To <- 293.15
    temp <- temp + 273.15
    
    psat <- pr * 10 ** (-6.8346 * (To1 / temp) ** 1.261 + 4.6151)
    
    h <- rh * (psat / pa)
    
    fr0 <-
      (pa / pr) * (24 + 4.04e4 * h * ((0.02 + h) / (0.391 + h)))
    
    frN <-
      (pa / pr) * sqrt(temp / To) * (9 + 280 * h * exp(-4.170 * ((temp / To) **
                                                                   (-1 / 3) - 1)))
    
    
    z <- 0.1068 * exp(-3352 / temp) / (frN + frequency ** 2 / frN)
    
    y <-
      (temp / To) ** (-5 / 2) * (0.01275 * exp(-2239.1 / temp) * 1 / (fr0 + frequency **
                                                                        2 / fr0) + z)
    
    
    atm_att_coef <-
      8.686 * frequency ** 2 * ((1.84e-11 * 1 / (pa / pr) * sqrt(temp / To)) + y)
    atm_att <- atm_att_coef * (dist - dist0)
    
    geom_att <- -20 * log10(dist0 / dist)
    
    hab_att <-
      hab.att.coef * frequency * (1 / 1000) * (dist - dist0)
    
    total_att <- geom_att + atm_att + hab_att
    
    outdf <-
      data.frame(
        frequency = frequency,
        dist = dist,
        geometric.attenuation = geom_att,
        atmopheric.attenuation = atm_att,
        habitat.attenuation = hab_att,
        combined.attenuation = total_att
      )
    
    return(outdf)
  }

##### modified from  http://www.sengpielaudio.com
## and https://scikit-maad.github.io/generated/maad.spl.attenuation_dB.html#maad.spl.attenuation_dB
