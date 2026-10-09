# Measure frequency spectrum correlation

`spectrum_correlation()` measures frequency spectrum correlation of
sounds referenced in an extended selection table. Spectral correlation
measures the similarity of two sounds in the frequency domain.

## Usage

``` r
spectrum_correlation(
  X,
  cores = getOption("mc.cores", 1),
  pb = getOption("pb", TRUE),
  cor.method = c("pearson", "spearman", "kendall"),
  spec.smooth = getOption("spec.smooth", 5),
  hop.size = getOption("hop.size", 11.6),
  wl = getOption("wl", NULL),
  ovlp = getOption("ovlp", 70),
  path = getOption("sound.files.path", "."),
  n.bins = 100
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

- cor.method:

  Character string indicating the correlation coefficient to be applied
  (`"pearson"`, `"spearman"`, or `"kendall"`, see
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html)).

- spec.smooth:

  Numeric vector of length 1 determining the length of the sliding
  window used for a sum smooth for power spectrum calculation (in kHz).
  Default `5`.

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
  Default `70`. Can be set globally for the current R session via the
  `"ovlp"` option (see
  [`options()`](https://rdrr.io/r/base/options.html)).

- path:

  Character string containing the directory path where the sound files
  are found. Only needed when `X` is not an extended selection table. If
  not supplied the current working directory is used. Can be set
  globally for the current R session via the `"sound.files.path"` option
  (see [`options()`](https://rdrr.io/r/base/options.html)).

- n.bins:

  Numeric vector of length 1 specifying the number of frequency bins to
  use for representing power spectra. Default `100`. If `NULL`, the raw
  power spectrum is used (note that this can result in high RAM memory
  usage for large data sets). Power spectrum values are interpolated
  using [`stats::approx()`](https://rdrr.io/r/stats/approxfun.html).

## Value

Object `X` with an additional column, `spectrum.correlation`, containing
the computed frequency spectrum correlation coefficients.

## Details

The function measures the spectral correlation coefficients of sounds in
which a reference playback has been re-recorded at increasing distances.
Values range from 1 (identical frequency spectrum, i.e. no degradation)
to 0. The `sound.id` column must be used to tell the function to only
compare sounds belonging to the same category (e.g. song-types). The
function will then compare each sound to the corresponding reference
sound. Two methods for computing spectral correlation are provided (see
the `method` argument). The function uses
[`seewave::meanspec()`](https://rdrr.io/pkg/seewave/man/meanspec.html)
internally to compute power spectra. Use
[`spectrum_blur_ratio()`](https://marce10.github.io/baRulho/reference/spectrum_blur_ratio.md)
to extract raw spectra values. `NA` is returned if at least one of the
power spectra cannot be computed.

## References

Araya-Salas, M., Grabarczyk, E. E., Quiroz-Oliva, M., Garcia-Rodriguez,
A., & Rico-Guevara, A. (2025). Quantifying degradation in animal
acoustic signals with the R package baRulho. Methods in Ecology and
Evolution, 00, 1-12. https://doi.org/10.1111/2041-210X.14481 Apol, C.A.,
Sturdy, C.B. & Proppe, D.S. (2017). Seasonal variability in habitat
structure may have shaped acoustic signals and repertoires in the
black-capped and boreal chickadees. Evol Ecol. 32:57-74.

## See also

[`envelope_correlation()`](https://marce10.github.io/baRulho/reference/envelope_correlation.md)
and
[`spectrum_blur_ratio()`](https://marce10.github.io/baRulho/reference/spectrum_blur_ratio.md).

Other quantify degradation:
[`blur_ratio()`](https://marce10.github.io/baRulho/reference/blur_ratio.md),
[`detection_distance()`](https://marce10.github.io/baRulho/reference/detection_distance.md),
[`envelope_correlation()`](https://marce10.github.io/baRulho/reference/envelope_correlation.md),
[`plot_blur_ratio()`](https://marce10.github.io/baRulho/reference/plot_blur_ratio.md),
[`plot_degradation()`](https://marce10.github.io/baRulho/reference/plot_degradation.md),
[`set_reference_sounds()`](https://marce10.github.io/baRulho/reference/set_reference_sounds.md),
[`signal_to_noise_ratio()`](https://marce10.github.io/baRulho/reference/signal_to_noise_ratio.md),
[`spcc()`](https://marce10.github.io/baRulho/reference/spcc.md),
[`spectrum_blur_ratio()`](https://marce10.github.io/baRulho/reference/spectrum_blur_ratio.md),
[`tail_to_signal_ratio()`](https://marce10.github.io/baRulho/reference/tail_to_signal_ratio.md)

## Author

Marcelo Araya-Salas (<marcelo.araya@ucr.ac.cr>)

## Examples

``` r
{
  # load example data
  data("test_sounds_est")

  # method 1
  # add reference column
  Y <- set_reference_sounds(X = test_sounds_est)

  # run spectrum correlation
  spectrum_correlation(X = Y)

  # method 2
  Y <- set_reference_sounds(X = test_sounds_est, method = 2)
  # spectrum_correlation(X = Y)
}
```
