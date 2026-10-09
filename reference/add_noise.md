# Add synthetic noise to annotations

`add_noise()` adds synthetic noise to sounds referenced in an extended
selection table to decrease the signal-to-noise ratio. This can be
useful, for instance, for evaluating the effect of background noise on
signal structure. Note that the implementation is slow.

## Usage

``` r
add_noise(
  X,
  mar = NULL,
  target.snr = 2,
  precision = 0.1,
  cores = getOption("mc.cores", 1),
  pb = getOption("pb", TRUE),
  max.iterations = 1000,
  kind = c("pink", "white", "brown", "red", "power"),
  alpha = 1,
  seed = 123,
  ...
)
```

## Arguments

- X:

  Object of class `extended_selection_table` (created by
  [`warbleR::selection_table()`](https://marce10.github.io/warbleR/reference/selection_table.html)
  from the **warbleR** package), generated ["by
  element"](https://marce10.github.io/warbleR/articles/b_annotation_data_format.html#by-element-vs-by-song-extended-selection-tables),
  with the test sound files' annotations (which in **baRulho** is
  typically the output of
  [`align_test_files()`](https://marce10.github.io/baRulho/reference/align_test_files.md)).
  Must contain the following columns: 1) `sound.files`: name of the
  `.wav` files, 2) `selec`: unique selection identifier (within a sound
  file), 3) `start`: start time and 4) `end`: end time of selections, 5)
  `bottom.freq`: low frequency for bandpass, and 6) `top.freq`: high
  frequency for bandpass. If the `sound.id` column is supplied, noise is
  only added to those sounds with a `sound.id` different from
  `"ambient"`, `"start_marker"`, or `"end_marker"`.

- mar:

  Numeric vector of length 1. Specifies the margins adjacent to the
  start point of the annotation over which to measure ambient noise.

- target.snr:

  Numeric vector of length 1. Specifies the desired signal-to-noise
  ratio. Must be lower than the current signal-to-noise ratio.
  Annotations showing a signal-to-noise ratio higher than `target.snr`
  will remain unchanged. Must be supplied.

- precision:

  Numeric vector of length 1. Specifies the precision of the adjusted
  signal-to-noise ratio (in dB).

- cores:

  Numeric vector of length 1. Controls whether parallel computing is
  applied by specifying the number of cores to be used. Default `1`
  (i.e. no parallel computing). Can be set globally for the current R
  session via the `"mc.cores"` option (see
  [`options()`](https://rdrr.io/r/base/options.html)).

- pb:

  Logical argument to control if progress bar is shown. Default `TRUE`.
  Can be set globally for the current R session via the `"pb"` option
  (see [`options()`](https://rdrr.io/r/base/options.html)).

- max.iterations:

  Numeric vector of length 1. Specifies the maximum number of iterations
  that the internal signal-to-noise adjusting routine will run before
  stopping. Note that in most cases the default maximum number of
  iterations (1000) is not reached.

- kind:

  Character vector of length 1 indicating the kind of noise: `"white"`,
  `"pink"`, `"power"`, `"brown"`, or `"red"`. Noise is synthesized with
  a modified version of
  [`tuneR::noise()`](https://rdrr.io/pkg/tuneR/man/Waveforms.html).
  Default `"pink"`, which is similar to background noise in natural
  environments.

- alpha:

  Numeric vector of length 1. The power for the power law noise
  (defaults are 1 for pink and 1.5 for red noise). Only used when
  `kind = "power"`.

- seed:

  Numeric vector of length 1. Seed for random number generation. Default
  `123`. If `NULL`, the seed is not set.

- ...:

  Additional arguments to be passed internally to
  [`signal_to_noise_ratio()`](https://marce10.github.io/baRulho/reference/signal_to_noise_ratio.md).
  Note that the `"custom"` noise reference (argument `noise.ref` in
  [`signal_to_noise_ratio()`](https://marce10.github.io/baRulho/reference/signal_to_noise_ratio.md))
  is currently not supported.

## Value

Object `X` in which the wave objects have been modified to match the
target signal-to-noise ratio. It also includes an additional column,
`adjusted.snr`, with the new signal-to-noise ratio values.

## Details

The function adds synthetic noise to sounds referenced in an extended
selection table (class created by
[`warbleR::selection_table()`](https://marce10.github.io/warbleR/reference/selection_table.html)
from the **warbleR** package) by iteratively amplifying the synthesized
noise and mixing it into each sound's waveform until the measured
signal-to-noise ratio reaches `target.snr` (within `precision` dB) or
`max.iterations` is exceeded. Annotations whose signal-to-noise ratio is
already at or below `target.snr` are left unmodified, and a warning
lists how many of these were skipped.

## References

Araya-Salas, M., Grabarczyk, E. E., Quiroz-Oliva, M., Garcia-Rodriguez,
A., & Rico-Guevara, A. (2025). Quantifying degradation in animal
acoustic signals with the R package baRulho. Methods in Ecology and
Evolution, 00, 1-12. https://doi.org/10.1111/2041-210X.14481 Timmer. J
and M. König (1995): On generating power law noise. Astron. Astrophys.
300, 707-710.

## See also

[`signal_to_noise_ratio()`](https://marce10.github.io/baRulho/reference/signal_to_noise_ratio.md),
which this function calls internally to measure the current
signal-to-noise ratio.

Other miscellaneous:
[`attenuation()`](https://marce10.github.io/baRulho/reference/attenuation.md),
[`noise_profile()`](https://marce10.github.io/baRulho/reference/noise_profile.md)

## Author

Marcelo Araya-Salas (<marcelo.araya@ucr.ac.cr>)

## Examples

``` r
if (FALSE) { # \dontrun{
# load example data
data("test_sounds_est")

# make it a 'by element' extended selection table
X <- warbleR::by_element_est(X = test_sounds_est)

# add noise to the first five rows
X_noise <- add_noise(X = X[1:5, ], mar = 0.2, target.snr = 3)
} # }
```
