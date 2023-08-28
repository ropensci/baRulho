baRulho: quantifying habitat-induced degradation of (animal) sounds
================

<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![Dependencies](https://tinyverse.netlify.com/badge/baRulho)](https://cran.r-project.org/package=baRulho)
[![Project Status: Active The project has reached a stable, usable state
and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Licence](https://img.shields.io/badge/https://img.shields.io/badge/licence-GPL--2-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
[![minimal R
version](https://img.shields.io/badge/R%3E%3D-%3E=%203.2.1-6666ff.svg)](https://cran.r-project.org/)
[![packageversion](https://img.shields.io/badge/Package%20version-2.1.0-orange.svg?style=flat-square)](commits/develop)
[![Last-changedate](https://img.shields.io/badge/last%20change-2023--08--26-yellowgreen.svg)](/commits/master)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/baRulho)](https://cran.r-project.org/package=baRulho)
[![Total
Downloads](https://cranlogs.r-pkg.org/badges/grand-total/baRulho)](https://cranlogs.r-pkg.org/badges/grand-total/baRulho)
[![Codecov test
coverage](https://codecov.io/gh/maRce10/baRulho/branch/master/graph/badge.svg)](https://app.codecov.io/gh/maRce10/baRulho?branch=master)
<!-- badges: end -->

<img src="vignettes/baRulho_sticker.png" alt="baRulho logo" align="right" width = "25%" height="25%"/>

[baRulho](https://cran.r-project.org/package=baRulho) is intended to
facilitate acoustic analysis of (animal) sound transmission experiments,
which typically aim to quantify changes in signal structure when
transmitted in a given habitat by broadcasting and re-recording animal
sounds at increasing distances. A common sequence of steps to
experimentally test hypotheses related to sound transmission is depicted
in the following diagram:

<center>
<img src="vignettes/analysis_workflow.jpg" alt="analysis workflow" width="620">
</center>

[baRulho](https://marce10.github.io/baRulho/) offers functions for
critical steps in this workflow (those in black, including ‘checks’)
that required acoustic data manipulation and analysis.

The main features of the package are:

- The use of loops to apply tasks through sounds referenced in a
  selection table (sensu
  [warbleR](https://cran.r-project.org/package=warbleR))
- The production of image files with graphic representations of sound in
  time and/or frequency that let users verify acoustic analyses
- The use of annotation tables as the object format to input acoustic
  data and annotations and to output results
- The use of parallelization to distribute tasks among several cores to
  improve computational efficiency

[baRulho](https://marce10.github.io/baRulho/) builds upon functions and
data formats from the
[warbleR](https://cran.r-project.org/package=warbleR) and
[seewave](https://cran.r-project.org/package=seewave) packages, so some
experience with these packages is advised.

Install/load the package from CRAN as follows:

``` r
# From CRAN would be
# install.packages("baRulho")

# load package
library(baRulho)
```

To install the latest developmental version from
[github](https://github.com/) you will need the R package
[remotes](https://cran.r-project.org/package=remotes):

``` r
# From github
remotes::install_github("maRce10/baRulho")

# load package
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
