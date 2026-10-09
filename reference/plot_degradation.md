# Save multipanel plots with reference and test sounds

`plot_degradation()` creates multipanel plots (as image files) with
reference and test sounds by distance and transect.

## Usage

``` r
plot_degradation(
  X,
  nrow = 4,
  env.smooth = getOption("env.smooth", 200),
  hop.size = getOption("hop.size", 11.6),
  wl = getOption("wl", NULL),
  ovlp = getOption("ovlp", 70),
  path = getOption("sound.files.path", "."),
  dest.path = getOption("dest.path", "."),
  cores = getOption("mc.cores", 1),
  pb = getOption("pb", TRUE),
  collevels = seq(-120, 0, 5),
  palette = viridis::viridis,
  flim = c("-1", "+1"),
  envelope = TRUE,
  spectrum = TRUE,
  heights = c(4, 1),
  widths = c(5, 1),
  margins = c(2, 1),
  row.height = 2,
  col.width = 2,
  cols = viridis::mako(4, alpha = 0.3),
  res = 120,
  ...
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

- nrow:

  Numeric vector of length 1 with the number of rows per image file.
  Default `4`. This will be dynamically adjusted if more rows than
  needed are set.

- env.smooth:

  Numeric vector of length 1 determining the length of the sliding
  window (in amplitude samples) used for a sum smooth for amplitude
  envelope and power spectrum calculations (used internally by
  [`seewave::env()`](https://rdrr.io/pkg/seewave/man/env.html)). Default
  `200`.

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

- path:

  Character string containing the directory path where the sound files
  are found. Only needed when `X` is not an extended selection table. If
  not supplied the current working directory is used. Can be set
  globally for the current R session via the `"sound.files.path"` option
  (see [`options()`](https://rdrr.io/r/base/options.html)).

- dest.path:

  Character string containing the directory path where the image files
  will be saved. If not supplied the current working directory will be
  used instead. Can be set globally for the current R session via the
  `"dest.path"` option (see
  [`options()`](https://rdrr.io/r/base/options.html)).

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

- collevels:

  Numeric vector indicating a set of levels used to partition the
  amplitude range of the spectrogram (in dB), as in
  [`seewave::spectro()`](https://rdrr.io/pkg/seewave/man/spectro.html).
  Default `seq(-120, 0, 5)`.

- palette:

  A color palette function to be used to assign colors in the plot, as
  in
  [`seewave::spectro()`](https://rdrr.io/pkg/seewave/man/spectro.html).
  Default
  [`viridis::viridis()`](https://sjmgarnier.github.io/viridis/reference/reexports.html).

- flim:

  Numeric vector of length 2 indicating the highest and lowest frequency
  limits (kHz) of the spectrogram, as in
  [`seewave::spectro()`](https://rdrr.io/pkg/seewave/man/spectro.html).
  Default `NULL`. Alternatively, a character vector similar to
  `c("-1", "1")`, in which the first number is the value to be added to
  the minimum bottom frequency in `X` and the second the value to be
  added to the maximum top frequency in `X`. This is computed
  independently for each sound ID, so the frequency limit better fits
  the frequency range of the annotated signals. This is useful when test
  sounds show marked differences in their frequency ranges.

- envelope:

  Logical to control if envelopes are plotted. Default `TRUE`.

- spectrum:

  Logical to control if power spectra are plotted. Default `TRUE`.

- heights:

  Numeric vector of length 2 to control the relative heights of the
  spectrogram (first number) and amplitude envelope (second number) when
  `envelope = TRUE`. Default `c(4, 1)`.

- widths:

  Numeric vector of length 2 to control the relative widths of the
  spectrogram (first number) and power spectrum (second number) when
  `spectrum = TRUE`. Default `c(5, 1)`.

- margins:

  Numeric vector of length 2 to control the relative time of the test
  sound (first number) and adjacent margins (i.e. adjacent background
  noise, second number) to be included in the spectrogram when
  `spectrum = TRUE`. Default `c(2, 1)`, which means that each margin
  next to the sound is half the duration of the sound. Note that all
  spectrograms will have the same time length, so margins will be
  calculated to ensure all spectrograms match the duration of the
  spectrogram of the longest sound. As such, this argument controls the
  margin on the longest sound.

- row.height:

  Numeric vector of length 1 controlling the height (in inches) of sound
  panels in the output image file. Default `2`.

- col.width:

  Numeric vector of length 1 controlling the width (in inches) of sound
  panels in the output image file. Default `2`.

- cols:

  Character vector of length 4 containing the colors to be used for the
  background of column and row title panels (element 1), the color of
  amplitude envelopes (element 2), the color of power spectra (element
  3), and the background color of envelopes and spectra (element 4).

- res:

  Numeric argument of length 1. Controls image resolution. Default `120`
  (faster), although 300-400 is recommended for publication/presentation
  quality.

- ...:

  Additional arguments to be passed to the internal spectrogram-creating
  function for customizing graphical output. The function is a modified
  version of
  [`seewave::spectro()`](https://rdrr.io/pkg/seewave/man/spectro.html),
  so it takes the same arguments.

## Value

One or more image files with a multipanel figure of spectrograms of test
sounds by distance, sound ID, and transect. It also returns the file
path of the images invisibly.

## Details

The function aims to simplify the visual inspection of sound degradation
by producing multipanel figures (saved in `dest.path`) containing
visualizations of each test sound and its reference. Sounds are sorted
by distance (columns) and transect (if more than 1). Visualizations
include spectrograms, amplitude envelopes, and power spectra (the last 2
are optional). Each row includes all the copies of a sound ID for a
given transect (the row label includes the sound ID in the first line
and transect in the second line), also including its reference if it
comes from another transect. Ambient noise annotations (`sound.id`
`"ambient"`) are excluded. Amplitude envelopes and power spectra are
computed using
[`warbleR::envelope()`](https://marce10.github.io/warbleR/reference/envelope.html)
and [`seewave::spec()`](https://rdrr.io/pkg/seewave/man/spec.html)
respectively. These two visualizations show the power distribution in
time and frequency between the minimum and maximum power values for each
sound. Therefore, scales are not necessarily comparable across panels.
If a transect ID is not supplied in `X` (i.e. no `transect` column), the
function assumes that all sounds are from the same transect. Only a
single copy of a test sound (`sound.id` label) is allowed per
transect/distance combination. The function uses
[`seewave::spectro()`](https://rdrr.io/pkg/seewave/man/spectro.html)
internally to create the spectrograms.

## References

Araya-Salas, M., Grabarczyk, E. E., Quiroz-Oliva, M., Garcia-Rodriguez,
A., & Rico-Guevara, A. (2025). Quantifying degradation in animal
acoustic signals with the R package baRulho. Methods in Ecology and
Evolution, 00, 1-12. https://doi.org/10.1111/2041-210X.14481

## See also

[`blur_ratio()`](https://marce10.github.io/baRulho/reference/blur_ratio.md)
and
[`plot_aligned_sounds()`](https://marce10.github.io/baRulho/reference/plot_aligned_sounds.md);
also
[`plot_blur_ratio()`](https://marce10.github.io/baRulho/reference/plot_blur_ratio.md),
for the analogous multipanel plot of blur ratio estimations.

Other quantify degradation:
[`blur_ratio()`](https://marce10.github.io/baRulho/reference/blur_ratio.md),
[`detection_distance()`](https://marce10.github.io/baRulho/reference/detection_distance.md),
[`envelope_correlation()`](https://marce10.github.io/baRulho/reference/envelope_correlation.md),
[`plot_blur_ratio()`](https://marce10.github.io/baRulho/reference/plot_blur_ratio.md),
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
# \donttest{
  # load example data
  data("test_sounds_est")

  # order so spectrograms from same sound id as close in the graph
  test_sounds_est <- test_sounds_est[order(test_sounds_est$sound.id), ]

  # set directory to save image files
  options(dest.path = tempdir())

  # method 1
  Y <- set_reference_sounds(X = test_sounds_est)

  # plot degradation spectrograms
  plot_degradation(
    X = Y, nrow = 3, ovlp = 95
  )
#> The image files have been saved in the directory path '/tmp/Rtmp0XZ3jd'

  # using other color palettes
  plot_degradation(
    X = Y, nrow = 3, ovlp = 95,
    cols = viridis::magma(4, alpha = 0.3),
    palette = viridis::magma
  )
#> The image files have been saved in the directory path '/tmp/Rtmp0XZ3jd'

  # missing some data, 2 rows
  plot_degradation(
    X = Y[-3, ], nrow = 2, ovlp = 95,
    cols = viridis::mako(4, alpha = 0.4), palette = viridis::mako, wl = 200
  )
#> The image files have been saved in the directory path '/tmp/Rtmp0XZ3jd'

  # changing marging and high overlap
  plot_degradation(X = Y, margins = c(5, 1), nrow = 6, ovlp = 95)
#> The image files have been saved in the directory path '/tmp/Rtmp0XZ3jd'

  # more rows than needed (will adjust it automatically)
  plot_degradation(X = Y, nrow = 10, ovlp = 90)
#> The image files have been saved in the directory path '/tmp/Rtmp0XZ3jd'
# }
```
