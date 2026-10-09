# Plot blur ratio

`plot_blur_ratio()` plots time and frequency blur ratio in sounds
referenced in an extended selection table.

## Usage

``` r
plot_blur_ratio(
  X,
  type = c("envelope", "spectrum"),
  cores = getOption("mc.cores", 1),
  pb = getOption("pb", TRUE),
  env.smooth = getOption("env.smooth", 200),
  spec.smooth = getOption("spec.smooth", 5),
  res = 150,
  flim = c("-1", "+1"),
  hop.size = getOption("hop.size", 11.6),
  wl = getOption("wl", NULL),
  ovlp = getOption("ovlp", 70),
  palette = viridis::viridis,
  collevels = seq(-120, 0, 5),
  dest.path = getOption("dest.path", "."),
  path = getOption("sound.files.path", "."),
  colors = viridis::viridis(3),
  n.samples = if (type == "envelope") 100 else 500
)
```

## Arguments

- X:

  The output of
  [`set_reference_sounds()`](https://marce10.github.io/baRulho/reference/set_reference_sounds.md),
  an object of class `data.frame`, `selection_table`, or
  `extended_selection_table` (the last 2 classes are created by
  [`warbleR::selection_table()`](https://marce10.github.io/warbleR/reference/selection_table.html)
  from the **warbleR** package) with the test sound files' annotations.
  Must contain the following columns: 1) `sound.files`: name of the
  `.wav` files, 2) `selec`: unique selection identifier (within a sound
  file), 3) `start`: start time and 4) `end`: end time of selections, 5)
  `bottom.freq`: low frequency for bandpass, 6) `top.freq`: high
  frequency for bandpass, 7) `sound.id`: ID of sounds used to identify
  counterparts across distances, and 8) `reference`: identity of sounds
  to be used as reference for each test sound (row). See
  [`set_reference_sounds()`](https://marce10.github.io/baRulho/reference/set_reference_sounds.md)
  for more details on the structure of `X`.

- type:

  Character vector of length 1 indicating the type of blur ratio to
  plot. The two options are `"envelope"` (for regular blur ratio as in
  [`blur_ratio()`](https://marce10.github.io/baRulho/reference/blur_ratio.md),
  default) and `"spectrum"` (for spectrum blur ratio as in
  [`spectrum_blur_ratio()`](https://marce10.github.io/baRulho/reference/spectrum_blur_ratio.md)).

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

- env.smooth:

  Numeric vector of length 1 determining the length of the sliding
  window (in amplitude samples) used for a sum smooth for amplitude
  envelope calculation (used internally by
  [`seewave::env()`](https://rdrr.io/pkg/seewave/man/env.html)). Default
  `200`.

- spec.smooth:

  Numeric vector of length 1 determining the length of the sliding
  window used for a sum smooth for power spectrum calculation (in kHz).
  Default `5`. Can be set globally for the current R session via the
  `"spec.smooth"` option (see
  [`options()`](https://rdrr.io/r/base/options.html)).

- res:

  Numeric argument of length 1. Controls image resolution. Default `150`
  (faster), although 300-400 is recommended for publication/presentation
  quality.

- flim:

  Numeric vector of length 2 indicating the highest and lowest frequency
  limits (kHz) of the spectrograms, as in
  [`seewave::spectro()`](https://rdrr.io/pkg/seewave/man/spectro.html).
  Default `NULL`. Alternatively, a character vector similar to
  `c("-1", "1")`, in which the first number is the value to be added to
  the minimum bottom frequency in `X` and the second the value to be
  added to the maximum top frequency in `X`. This is computed
  independently for each sound ID, so the frequency limit better fits
  the frequency range of the annotated signals. This is useful when test
  sounds show marked differences in their frequency ranges.

- hop.size:

  Numeric vector of length 1 specifying the time window duration (in
  ms). Default `11.6` ms, which is equivalent to 512 `wl` for a 44.1 kHz
  sampling rate. Ignored if `wl` is supplied. Can be set globally for
  the current R session via the `"hop.size"` option (see
  [`options()`](https://rdrr.io/r/base/options.html)).

- wl:

  A vector with a single even integer number specifying the window
  length of the spectrogram. Default `NULL`. If supplied, `hop.size` is
  ignored. Odd integers will be rounded up to the nearest even number.
  Can be set globally for the current R session via the `"wl"` option
  (see [`options()`](https://rdrr.io/r/base/options.html)).

- ovlp:

  Numeric vector of length 1 specifying the percentage of overlap
  between two consecutive windows, as in
  [`seewave::spectro()`](https://rdrr.io/pkg/seewave/man/spectro.html).
  Only used when plotting. Default `70`. Applied to both spectra and
  spectrograms on image files. Can be set globally for the current R
  session via the `"ovlp"` option (see
  [`options()`](https://rdrr.io/r/base/options.html)).

- palette:

  A color palette function to be used to assign colors in the plot, as
  in
  [`seewave::spectro()`](https://rdrr.io/pkg/seewave/man/spectro.html).
  Default
  [`viridis::viridis()`](https://sjmgarnier.github.io/viridis/reference/reexports.html).

- collevels:

  Numeric vector indicating a set of levels used to partition the
  amplitude range of the spectrogram (in dB), as in
  [`seewave::spectro()`](https://rdrr.io/pkg/seewave/man/spectro.html).
  Default `seq(-120, 0, 5)`.

- dest.path:

  Character string containing the directory path where the image files
  will be saved. If not supplied the current working directory will be
  used instead. Can be set globally for the current R session via the
  `"dest.path"` option (see
  [`options()`](https://rdrr.io/r/base/options.html)).

- path:

  Character string containing the directory path where the sound files
  are found. Only needed when `X` is not an extended selection table. If
  not supplied the current working directory is used. Can be set
  globally for the current R session via the `"sound.files.path"` option
  (see [`options()`](https://rdrr.io/r/base/options.html)).

- colors:

  Character vector of length 4 containing the colors to be used for the
  color to identify the reference sound (element 1), the color to
  identify the test sound (element 2), and the color of the blurred
  region (element 3).

- n.samples:

  Numeric vector of length 1 specifying the number of amplitude samples
  (or frequency bins if `spectrum = TRUE`) to use for representing power
  distributions. Default `100` for `type = "envelope"` and `500` for
  `type = "spectrum"`. If `NULL`, the raw power distribution is used
  (note that this can result in high RAM memory usage for large data
  sets).

## Value

One image file (in `jpeg` format) for each blur ratio estimation,
showing spectrograms of both sounds and the overlaid amplitude envelopes
(or power spectra if `spectrum = TRUE`) as probability mass functions
(PMF). Spectrograms are shown within the frequency range of the
reference sound. It also returns the file path of the images invisibly.

## Details

The function generates image files (in `jpeg` format) for each possible
blur ratio estimation in `X`. The image files show the spectrograms of
both sounds and the overlaid power distribution (either amplitude
envelopes or power spectrum, see the `type` argument) as probability
mass functions (PMF). The output graphs highlight the mismatch between
the compared distributions, which represents the estimated blur ratio
returned by either
[`blur_ratio()`](https://marce10.github.io/baRulho/reference/blur_ratio.md)
or
[`spectrum_blur_ratio()`](https://marce10.github.io/baRulho/reference/spectrum_blur_ratio.md).
Spectrograms are shown within the frequency range of the reference
sound, and also show dotted lines with the time (`type = "envelope"`) or
frequency range (`type = "spectrum"`) in which energy distributions were
computed.

## References

Dabelsteen, T., Larsen, O. N., & Pedersen, S. B. (1993). Habitat-induced
degradation of sound signals: Quantifying the effects of communication
sounds and bird location on blur ratio, excess attenuation, and
signal-to-noise ratio in blackbird song. The Journal of the Acoustical
Society of America, 93(4), 2206.

Araya-Salas, M., Grabarczyk, E. E., Quiroz-Oliva, M., Garcia-Rodriguez,
A., & Rico-Guevara, A. (2025). Quantifying degradation in animal
acoustic signals with the R package baRulho. Methods in Ecology and
Evolution, 00, 1-12. https://doi.org/10.1111/2041-210X.14481

## See also

[`envelope_correlation()`](https://marce10.github.io/baRulho/reference/envelope_correlation.md),
[`spectrum_blur_ratio()`](https://marce10.github.io/baRulho/reference/spectrum_blur_ratio.md),
and
[`blur_ratio()`](https://marce10.github.io/baRulho/reference/blur_ratio.md),
which this function visualizes.

Other quantify degradation:
[`blur_ratio()`](https://marce10.github.io/baRulho/reference/blur_ratio.md),
[`detection_distance()`](https://marce10.github.io/baRulho/reference/detection_distance.md),
[`envelope_correlation()`](https://marce10.github.io/baRulho/reference/envelope_correlation.md),
[`plot_degradation()`](https://marce10.github.io/baRulho/reference/plot_degradation.md),
[`set_reference_sounds()`](https://marce10.github.io/baRulho/reference/set_reference_sounds.md),
[`signal_to_noise_ratio()`](https://marce10.github.io/baRulho/reference/signal_to_noise_ratio.md),
[`spcc()`](https://marce10.github.io/baRulho/reference/spcc.md),
[`spectrum_blur_ratio()`](https://marce10.github.io/baRulho/reference/spectrum_blur_ratio.md),
[`spectrum_correlation()`](https://marce10.github.io/baRulho/reference/spectrum_correlation.md),
[`tail_to_signal_ratio()`](https://marce10.github.io/baRulho/reference/tail_to_signal_ratio.md)

## Author

Marcelo Araya-Salas (<marcelo.araya@ucr.ac.cr>)

## Examples

``` r
{
  # load example data
  data("test_sounds_est")

  # add reference to X
  X <- set_reference_sounds(X = test_sounds_est)

  # create plots
  plot_blur_ratio(X = X, dest.path = tempdir())
}
#> The image files have been saved in the directory path '/tmp/RtmpLKxQlp'
```
