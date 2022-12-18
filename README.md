baRulho: quantifying habitat-induced degradation of (animal) acoustic
signals
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![Dependencies](https://tinyverse.netlify.com/badge/baRulho)](https://cran.r-project.org/package=baRulho)
[![Project Status: Active The project has reached a stable, usable state
and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Licence](https://img.shields.io/badge/https://img.shields.io/badge/licence-GPL--2-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
[![minimal R
version](https://img.shields.io/badge/R%3E%3D-%3E=%203.2.1-6666ff.svg)](https://cran.r-project.org/)
[![packageversion](https://img.shields.io/badge/Package%20version-1.0.7-orange.svg?style=flat-square)](commits/develop)
[![Last-changedate](https://img.shields.io/badge/last%20change-2022--12--18-yellowgreen.svg)](/commits/master)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/baRulho)](https://cran.r-project.org/package=baRulho)
[![Total
Downloads](https://cranlogs.r-pkg.org/badges/grand-total/baRulho)](https://cranlogs.r-pkg.org/badges/grand-total/baRulho)

`baRulho` is intended to facilitate acoustic analysis of (animal) sound
transmission experiments, which typically aim to quantify changes in
signal structure when transmitted in a given habitat by broadcasting and
re-recording animal sounds at increasing distances. The package offers a
workflow with functions to prepare the data set for analysis as well as
to calculate and visualize several degradation metrics, including blur
ratio, signal-to-noise ratio, excess attenuation and envelope
correlation among others (Dabelsteen et al 1993).

Install/load the package from CRAN as follows:

``` r
# From CRAN would be
#install.packages("baRulho")

#load package
library(baRulho)
```

To install the latest developmental version from
[github](https://github.com/) you will need the R package
[remotes](https://cran.r-project.org/package=remotes):

``` r
# From github
remotes::install_github("maRce10/baRulho")

#load package
library(baRulho)
```

Please cite [baRulho](https://marce10.github.io/baRulho/) as follows:

Araya-Salas, M. (2020), *baRulho: quantifying habitat-induced
degradation of (animal) acoustic signals in R*. R package version 1.0.0.

# References

1.  Dabelsteen, T., Larsen, O. N., & Pedersen, S. B. (1993).
    *Habitat-induced degradation of sound signals: Quantifying the
    effects of communication sounds and bird location on blur ratio,
    excess attenuation, and signal-to-noise ratio in blackbird song*.
    The Journal of the Acoustical Society of America, 93(4), 2206.

2.  Marten, K., & Marler, P. (1977). *Sound transmission and its
    significance for animal vocalization*. Behavioral Ecology and
    Sociobiology, 2(3), 271-290.

3.  Morton, E. S. (1975). *Ecological sources of selection on avian
    sounds*. The American Naturalist, 109(965), 17-34.
