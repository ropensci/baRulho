# baRulho

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/baRulho)](https://cran.r-project.org/package=baRulho)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/baRulho)](http://www.r-pkg.org/pkg/baRulho)
[![Total downloads](https://cranlogs.r-pkg.org/badges/grand-total/baRulho?color=blue)](https://r-pkg.org/pkg/baRulho)
[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)

`baRulho` is intended to facilitate acoustic analysis of (animal) sound transmission experiments, which typically aim to quantify changes in signal structure when transmitted in a given habitat by broadcasting and re-recording animal sounds at increasing distances. The package offers a workflow with functions to prepare the data set for analysis as well as to calculate and visualize several degradation metrics, including blur ratio, signal-to-noise ratio, excess attenuation and envelope correlation among others (Dabelsteen et al 1993).

Install/load the package from CRAN as follows:

```r

# From CRAN would be
#install.packages("baRulho")

#load package
library(baRulho)

```

To install the latest developmental version from [github](http://github.com/) you will need the R package [devtools](https://cran.r-project.org/package=devtools):

```r

# From github
devtools::install_github("maRce10/baRulho")

#load package
library(baRulho)

```

Please cite [baRulho](https://marce10.github.io/baRulho/) as follows:

Araya-Salas, M. (2020), *baRulho: quantifying habitat-induced degradation of (animal) acoustic signals in R*. R package version 1.0.0.

# References

1. Araya-Salas, M. (2017), *Rraven: connecting R and Raven bioacoustic software*. R package version 1.0.0.

1. Araya-Salas M, Smith-Vidaurre G (2017) *warbleR: An R package to streamline analysis of animal acoustic signals*. Methods Ecol Evol 8:184–191.

1. Dabelsteen, T., Larsen, O. N., & Pedersen, S. B. (1993). *Habitat-induced degradation of sound signals: Quantifying the effects of communication sounds and bird location on blur ratio, excess attenuation, and signal-to-noise ratio in blackbird song*. The Journal of the Acoustical Society of America, 93(4), 2206.

1. Marten, K., & Marler, P. (1977). *Sound transmission and its significance for animal vocalization. Behavioral* Ecology and Sociobiology, 2(3), 271-290.

1. Morton, E. S. (1975). *Ecological sources of selection on avian sounds*. The American Naturalist, 109(965), 17-34.

1. Tobias, J. A., Aben, J., Brumfield, R. T., Derryberry, E. P., Halfwerk, W., Slabbekoorn, H., & Seddon, N. (2010). *Song divergence by sensory drive in Amazonian birds*. Evolution, 64(10), 2820-2839.
