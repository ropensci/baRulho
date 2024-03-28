#' baRulho: quantifying acoustic signal degradation
#'
#' `baRulho` is a package intended to quantify habitat-induced degradation of (animal) acoustic signals.
#'
#' The main features of the package are:
#'   \itemize{
#'   \item Loops to apply tasks through sounds referenced in an extended selection table
#'   \item The comparison of playback sounds re-recorded at different distances
#'   }
#' Most functions allow the parallelization of tasks, which distributes the tasks among several processors to improve computational efficiency.
#'
#' @import tuneR
#' @import seewave
#' @import graphics
#' @import grDevices
#' @import utils
#' @import parallel
#' @import warbleR
#' @import fftw
#' @import viridis
#' @import Sim.DiffProc
#' @import png
#' @importFrom rlang arg_match check_installed
#' @importFrom stats cor aggregate approx ave na.omit complete.cases runif rnorm fft
#' @importFrom cli style_bold style_italic make_ansi_style num_ansi_colors
#' @importFrom checkmate assert_logical assert_character assert_integerish assert_numeric assert_multi_class assert_list assert_function assert_class assert_directory assert_data_frame assert_names makeAssertCollection makeAssertionFunction reportAssertions
#' @importFrom methods is
#' @importFrom ohun template_detector template_correlator
#'
#' @author Marcelo Araya-Salas
#'
#'   Maintainer: Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})
#'
#' @docType package
#' @details License: GPL (>= 2)
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
