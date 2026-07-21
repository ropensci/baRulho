# baRulho: quantifying degradation of (animal) sounds

![baRulho logo](reference/figures/baRulho_sticker.png)

[baRulho](https://cran.r-project.org/package=baRulho) is intended to
facilitate the implementation of (animal) sound propagation experiments,
which typically aim to quantify changes in signal structure when
transmitted in a given habitat by broadcasting and re-recording animal
sounds at increasing distances.

These experiments aim to answer research questions such as:

- How habitat structure has shaped the propagation properties of animal
  acoustic signals?
- Which acoustic features are shaped by selection for improving
  propagation?
- Which features are more degraded in different habitats?
- How far a acoustic signals can be detected?

A common sequence of steps to experimentally test hypotheses related to
sound propagation is depicted in the following diagram:

![](reference/figures/analysis_workflow.png)

*Diagram depicting a typical workflow for a experiment working on signal
propagation and degradation. Nodes with black font indicate steps that
can be conducted using baRulho functions. Blue nodes denote the
functions that can be used at those steps.*

 

[baRulho](https://docs.ropensci.org/baRulho//) offers functions for
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

[baRulho](https://docs.ropensci.org/baRulho//) builds upon functions and
data formats from the
[warbleR](https://cran.r-project.org/package=warbleR) and
[seewave](https://cran.r-project.org/package=seewave) packages, so some
experience with these packages is advised.

Take a look at the vignettes for an overview of the main features of the
packages:

- [Align test
  sounds](https://docs.ropensci.org/baRulho//articles/align_test_sounds.html)
- [Quantify
  degradation](https://docs.ropensci.org/baRulho//articles/quantify_degradation.html)

## Installing baRulho

Install/load the package from CRAN as follows:

``` r

# From CRAN would be
# install.packages("baRulho")

# load package
library(baRulho)
```

It can also be install from
[R-Universe](https://ropensci.org/blog/2021/06/22/setup-runiverse/) in
this way:

``` r

install.packages("baRulho", repos = "https://ropensci.r-universe.dev")
```

``` R
## Installing package into '/home/marce/R/x86_64-pc-linux-gnu-library/4.3'
## (as 'lib' is unspecified)
```

To install the latest developmental version from
[github](https://github.com/) you will need the R package
[remotes](https://cran.r-project.org/package=remotes):

``` r

# install remotes if not installed
if (!requireNamespace("remotes")) {
  install.packages("remotes")
}

# From github
remotes::install_github("ropensci/baRulho")

# load package
library(baRulho)
```

Further system requirements due to the dependency
[seewave](https://cran.r-project.org/package=seewave) may be needed.
Take a look a [this link](https://rug.mnhn.fr/seewave/inst.html) for
instruction on how to install/troubleshoot these external dependencies.

## Other packages

The packages [seewave](https://cran.r-project.org/package=seewave) and
[tuneR](https://cran.r-project.org/package=tuneR) provide a huge variety
of functions for acoustic analysis and manipulation. They mostly work on
wave objects already imported into the R environment. The package
[warbleR](https://cran.r-project.org/package=warbleR) provides functions
to visualize and measure sounds already referenced in annotation tables,
similar to [baRulho](https://docs.ropensci.org/baRulho//). The package
[Rraven](https://cran.r-project.org/package=Rraven) facilitates the
exchange of data between R and [Raven sound analysis
software](https://www.ravensoundsoftware.com/) ([Cornell Lab of
Ornithology](https://www.birds.cornell.edu/home)) and can be very
helpful for incorporating Raven as the annotating tool into acoustic
analysis workflow in R. The package
[ohun](https://github.com/ropensci/ohun) works on automated detection of
sound events, providing functions to diagnose and optimize detection
routines.

## Citation

Please cite [baRulho](https://docs.ropensci.org/baRulho//) as follows:

Araya-Salas M., E. Grabarczyk, M. Quiroz-Oliva, A. Garcia-Rodriguez, A.
Rico-Guevara. (2023), *baRulho: an R package to quantify degradation in
animal acoustic signals*. bioRxiv 2023.11.22.568305.

## References

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
