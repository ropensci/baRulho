#' Measure atmospheric attenuation and absorption of sound 
#' 
#' \code{atmospheric_attenuation} measures atmospheric attenuation and atmospheric absorption of signals referenced in an extended selection table.
#' @usage atmospheric_attenuation(f, temp, RH, p = 101325, 
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
#' @name atmospheric_attenuation
#' @details Calculate the atmospheric attenuation based on temperature, relative humidity, pressure and sound frequency. The function can applied to formulae based on:
#' \itemize{
#' \item \code{1}: default. As used by Bazley (1976), Sound absorption in air at frequencies up to 100 kHz. NPL acoustics report Ac 74. 
#' \item \code{2}: as used by Rossing (2007), Handbook of Acoustics, Springer.
#' }
#' If 'spi' and 'dist' are supplied the function also returns the atmospheric absorption (in dB).
#' @examples
#' {
#' # load example data
#' data("playback_est")
#' 
#' #' # remove ambient selections
#' playback_est <- playback_est[playback_est$signal.type != "ambient", ]
#'
#' # measure atmospheric attenuation formula 1
#'atmospheric_attenuation(f = 20000, temp = 20, RH = 90, p = 88000, formula = 1)
#' }
#' 
#' @author Marcelo Araya-Salas (\email{marceloa27@@gmail.com}) 
#' @references {
#' Araya-Salas, M. (2020), baRulho: baRulho: quantifying habitat-induced degradation of (animal) acoustic signals in R. R package version 1.0.0
#' }
#last modification on nov-01-2019 (MAS)

### NOT SURE IF ABSORPTION WORKS 
atmospheric_attenuation <- function(f, temp, RH, p = 101325, formula = 1, spi = NULL, dist = NULL){
  
  if (!is.null(spi) & is.null(dist)) stop("'dist' must also be supplied to calculate atmospheric absorption")

  if (is.null(spi) & !is.null(dist)) stop("'spi' must also be supplied to calculate atmospheric absorption")
  
  Ta <- temp + 273.15      # convert to Kelvin
  Tr <- Ta / 293.15    # convert to relative air temperature (re 20 deg C)
    
  if (formula == 1){ # According to Bazley 1976 (Marc)

    P = p / 101325   #convert to atm
    
    # percentage water molecules:
    h <- ((RH / P)  * Ta ^-4.922 * 10^(20.5318-2939 / Ta))
    
    # molecular attenuation caused by oxygen:
    umax <- (0.0000042425 + 0.000000088168 * temp + 0.00000000054834 * temp^2)
    fo  <- (30560 * P * h^1.3)
    mo <- (2 * umax * (f / (f / fo + fo / f)))
  
    # molecular attenuation caused by nitrogen:
    mn <- 0.0000000171 / sqrt(1 + 0.00366 * T) * (h * P * f^2 / (0.0000277 * f^2 + h^2 * P^2))
    
    # classical and rotational effects:
    mcr <- 0.000000000036 * (1+0.001 * temp) * f^2 / P
    
    atmatt <- 4.343 * (mo + mn + mcr)
  } else
    { # formula 2
    ## According to Rossing 2007 (Peter) - correct, 5.4.2012
    # Saturation Concentration of water vapor.
    # NOTE the *ERROR* in Rossing 2007!! Instead of Tref = 293.15 K (20 C),
    # here the triple-point isotherm temperature (273.15 K + 0.01 K =
    # 273.16 K) has to be used!
    # See http://forum.studiotips.com/viewtopic.php?t=158
    
    P <- p / 101325       # convert to relative pressure
    
    C <- 4.6151 - 6.8346 * ((273.16 / Ta)^1.261)
    
    # percentage molar concentration of water vapor:
    h <- RH * 10^C / P
    
    # relaxation frequencies associated with the vibration of oxygen and nitrogen:  
    frO <- P * (24+4.04e4 * h * (0.02 + h) / (0.391 + h))
    frN <- P *(Tr^(-0.5)) * (9+280 * h * exp(-4.17 * ((Tr^(-1/3))-1)))
    
    # attenuation coefficient (Np/(m*atm)):
    alpha <- f * f * (1.84e-11 * (1 / P) * sqrt(Tr) + (Tr^(-2.5)) * (0.01275 * (exp(-2239.1 / Ta) * 1 / (frO + f * f / frO)) + 0.1068 * (exp(-3352 / Ta) * 1 / (frN + f * f /frN))))
    
    atmatt <- 8.686 * alpha # convert to dB (lg(x/x0)) from Neper (ln(x/x0))
  }
  
  if (!is.null(spi) & !is.null(dist)){
    #        a ........ pure-tone sound attenuation coefficient, in dB/m, for atmospheric absorption
    #        s ........ distance in m through which the sounds propagates
    #        Pi (spi) ....... initial sound pressure amplitude, in Pa
    #        Pt ....... sound pressure amplitude, in Pa
    #        Pa (p / 1000) ...... ambient atmospheric pressure in kPa
    #        Pr ....... reference ambient atmospheric pressure: 101.325 kPa
    #        Psat .. saturation vapor pressure ca equals:
    #          .................. International Meteorological Tables WMO-No.188 TP94
    #        .................. World Meteorological Organization - Geneva Switzerland
    #        T ........ ambient atmospheric temperature in K (Kelvin).
    #        ........... K = 273.15 + Temperature in C (by US known as centigrade, Europe as Celsius)
    #        To ...... reference temperature in K: 293.15 K (20 C)
    #        To1..... triple-point isotherm temp: 273.16 K = 273.15 + 0.01 K (0.01 C)
    #        h ........ molar concentration of water vapor, as a percentage
    #        hr (RH)........ relative humidity as a percentage
    #        f ......... frequency
    #        frO ..... oxygen relaxation frequency
    #        frN ..... nitrogen relaxation frequency
    #        x ........ Just a help factor to shorten formula, improvement on standard by Eric Desart
    #        y ........ Just a help factor to shorten formula
    #        z ........ Just a help factor to shorten formula
    
    To <- 293.15
    To1 <- 273.16
    Pa <- p / 1000
    Pr <- 101.325
    
    Psat <- Pr * 10^(-6.8346 * (To1/Ta)^1.261 + 4.6151)
    
    h <- RH * ((Psat / Pr) / (Pa / Pr))
    
    frO <- (Pa / Pr) * (24 + 4.04 * 10^4 * h * ((0.02 + h) / (0.391 + h)))
    
    frN <- (Pa / Pr) * (Ta / To)^(-1/2) * (9 + 280 * h * exp(-4.170 * ((Ta / To)^(-1/3)-1)))
    
    z <- 0.1068 * exp(-3352 / Ta) * (frN + f ^2 / frN)^-1
    
    y <- (T / To)^(-5/2) * (0.01275 * exp(-2239.1 / T) * (frO + f ^2 / frO)^-1 + z)
    
    a <- 8.686 * f^2 * ((1.84 * 10^-11 * (Pa / Pr)^-1 * (Ta / To)^(1/2)) + y) #[dB/m]
    
    
    x <- 1 / (10 * log((exp(1))^2)) #= ca 0.1151 #(value in norm, formula E. Desart)
    
    as = a * dist #[dB] total absorption at distance s
    
    # Pt <- spi * exp(-x * as) # in Pa
    # 
    # as <- 10 * log( spi^2 / Pt^2 ) #= as [dB] (was Delta.Lt)
    # 
    
    
    reslt <- list(atmospheric.attenuation = atmatt, atmospheric.absorption = as) 
  } else reslt <- list(atmospheric.attenuation = atmatt) 
  
  return(reslt)  
  
}

##### from http://forum.studiotips.com/viewtopic.php?t=158
