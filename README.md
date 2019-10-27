# baRulho

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/baRulho)](https://cran.r-project.org/package=baRulho)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/baRulho)](http://www.r-pkg.org/pkg/baRulho)
[![Total Downloads](http://cranlogs.r-pkg.org/badges/grand-total/baRulho)](http://www.r-pkg.org/badges/grand-total/baRulho)
[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)

`baRulho` is a package intended to quantify habitat-induced degradation of (animal) acoustic signals. Most functions are based on the metrics provided in Dabelsteen et al (1993).

Install/load the package from CRAN as follows:

```r

# From CRAN would be
#install.packages("baRulho")

#load package
library(baRulho)

```

To install the latest developmental version from [github](http://github.com/) you will need the R package [devtools](https://cran.r-project.org/package=devtools):

```r
# From CRAN would be
#install.packages("baRulho")

# From github
devtools::install_github("maRce10/baRulho")

#load package
library(baRulho)

```

Please cite [baRulho](https://cran.r-project.org/package=baRulho) as follows:

Araya-Salas, M. (2019), *baRulho: a R package to evaluate habitat-induced degradation of (animal) acoustic signals*. R package version 1.0.0.

# References

Dabelsteen, T., Larsen, O. N., & Pedersen, S. B. (1993). Habitat-induced degradation of sound signals: Quantifying the effects of communication sounds and bird location on blur ratio, excess attenuation, and signal-to-noise ratio in blackbird song. The Journal of the Acoustical Society of America, 93(4), 2206.
