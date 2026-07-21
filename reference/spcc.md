# Measure spectrographic cross-correlation as a measure of sound distortion

`spcc` measures spectrographic cross-correlation as a measure of sound
distortion in sounds referenced in an extended selection table.

## Usage

``` r
spcc(
  X,
  cores = getOption("mc.cores", 1),
  pb = getOption("pb", TRUE),
  cor.method = c("pearson", "spearman", "kendall"),
  hop.size = getOption("hop.size", 11.6),
  wl = getOption("wl", NULL),
  ovlp = getOption("ovlp", 90),
  wn = "hanning",
  path = getOption("sound.files.path", ".")
)
```

## Arguments

- X:

  The output of
  [`set_reference_sounds`](https://marce10.github.io/baRulho/reference/set_reference_sounds.md)
  which is an object of class 'data.frame', 'selection_table' or
  'extended_selection_table' (the last 2 classes are created by the
  function
  [`selection_table`](https://marce10.github.io/warbleR/reference/selection_table.html)
  from the warbleR package) with the test sound files' annotations .
  Must contain the following columns: 1) "sound.files": name of the .wav
  files, 2) "selec": unique selection identifier (within a sound
  file), 3) "start": start time and 4) "end": end time of selections, 5)
  "bottom.freq": low frequency for bandpass, 6) "top.freq": high
  frequency for bandpass, 7) "sound.id": ID of sounds used to identify
  counterparts across distances and 8) "reference": identity of sounds
  to be used as reference for each test sound (row). See
  [`set_reference_sounds`](https://marce10.github.io/baRulho/reference/set_reference_sounds.md)
  for more details on the structure of 'X'.

- cores:

  Numeric vector of length 1. Controls whether parallel computing is
  applied by specifying the number of cores to be used. Default is 1
  (i.e. no parallel computing). Can be set globally for the current R
  session via the "mc.cores" option (see
  [`options`](https://rdrr.io/r/base/options.html)).

- pb:

  Logical argument to control if progress bar is shown. Default is
  `TRUE`. Can be set globally for the current R session via the "pb"
  option (see [`options`](https://rdrr.io/r/base/options.html)).

- cor.method:

  Character string indicating the correlation coefficient to be applied
  ("pearson", "spearman", or "kendall", see
  [`cor`](https://rdrr.io/r/stats/cor.html)).

- hop.size:

  A numeric vector of length 1 specifying the time window duration (in
  ms). Default is 11.6 ms, which is equivalent to 512 wl for a 44.1 kHz
  sampling rate. Ignored if 'wl' is supplied. Can be set globally for
  the current R session via the "hop.size" option (see
  [`options`](https://rdrr.io/r/base/options.html)).

- wl:

  a vector with a single even integer number specifying the window
  length of the spectrogram, default is `NULL`. If supplied, 'hop.size'
  is ignored. Odd integers will be rounded up to the nearest even
  number. Can be set globally for the current R session via the "wl"
  option (see [`options`](https://rdrr.io/r/base/options.html)).

- ovlp:

  Numeric vector of length 1 specifying % of overlap between two
  consecutive windows, as in
  [`spectro`](https://rdrr.io/pkg/seewave/man/spectro.html). Default
  is 90. High values of ovlp slow down the function but produce more
  accurate results. Can be set globally for the current R session via
  the "ovlp" option (see
  [`options`](https://rdrr.io/r/base/options.html)).

- wn:

  A character vector of length 1 specifying the window name as in
  [`ftwindow`](https://rdrr.io/pkg/seewave/man/ftwindow.html).

- path:

  Character string containing the directory path where the sound files
  are found. Only needed when 'X' is not an extended selection table. If
  not supplied the current working directory is used. Can be set
  globally for the current R session via the "sound.files.path" option
  (see [`options`](https://rdrr.io/r/base/options.html)).

## Value

Object 'X' with an additional column, 'cross.correlation', containing
the computed spectrogram cross-correlation coefficients.

## Details

Spectrographic cross-correlation measures frequency distortion of sounds
as a similarity metric. Values close to 1 means very similar
spectrograms (i.e. little sound distortion has occurred).
Cross-correlation is measured of sounds in which a reference playback
has been re-recorded at increasing distances. The 'sound.id' column must
be used to indicate the function to only compare sounds belonging to the
same category (e.g. song-types). The function compares each sound to the
corresponding reference sound within the supplied frequency range (e.g.
bandpass) of the reference sound ('bottom.freq' and 'top.freq' columns
in 'X'). Two methods for computing cross-correlation are provided (see
'method' argument). The function is a wrapper on warbleR's
[`cross_correlation`](https://marce10.github.io/warbleR/reference/cross_correlation.html)
function.

## References

Araya-Salas, M., Grabarczyk, E. E., Quiroz-Oliva, M., Garcia-Rodriguez,
A., & Rico-Guevara, A. (2025). Quantifying degradation in animal
acoustic signals with the R package baRulho. Methods in Ecology and
Evolution, 00, 1-12. https://doi.org/10.1111/2041-210X.14481 Clark,
C.W., Marler, P. & Beeman K. (1987). Quantitative analysis of animal
vocal phonology: an application to Swamp Sparrow song. Ethology.
76:101-115.

## See also

[`blur_ratio`](https://marce10.github.io/baRulho/reference/blur_ratio.md),
[`manual_realign`](https://marce10.github.io/baRulho/reference/manual_realign.md),
[`cross_correlation`](https://marce10.github.io/warbleR/reference/cross_correlation.html)

Other quantify degradation:
[`blur_ratio()`](https://marce10.github.io/baRulho/reference/blur_ratio.md),
[`detection_distance()`](https://marce10.github.io/baRulho/reference/detection_distance.md),
[`envelope_correlation()`](https://marce10.github.io/baRulho/reference/envelope_correlation.md),
[`plot_blur_ratio()`](https://marce10.github.io/baRulho/reference/plot_blur_ratio.md),
[`plot_degradation()`](https://marce10.github.io/baRulho/reference/plot_degradation.md),
[`set_reference_sounds()`](https://marce10.github.io/baRulho/reference/set_reference_sounds.md),
[`signal_to_noise_ratio()`](https://marce10.github.io/baRulho/reference/signal_to_noise_ratio.md),
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

  # get spcc
  spcc(X = X)
}
#> Object of class 'extended_selection_table'
#> * The output of the following call:
#> `[.extended_selection_table`(X = test_sounds_est, i = order(test_sounds_est$sound.id,  test_sounds_est$transect, decreasing = FALSE))
#> 
#> Contains: 
#> *  A selection table data frame with 25 row(s) and 11 columns:
#>      sound.files       selec   start      end   bottom.freq   top.freq
#> ---  ---------------  ------  ------  -------  ------------  ---------
#> 8    10m_closed.wav        1    0.05   0.2000        1.3333     2.6667
#> 23   30m_closed.wav        1    0.05   0.2000        1.3333     2.6667
#> 13   10m_open.wav          1    0.05   0.2000        1.3333     2.6667
#> 18   1m_open.wav           1    0.05   0.2000        1.3333     2.6667
#> 28   30m_open.wav          1    0.05   0.2000        1.3333     2.6667
#> 11   10m_closed.wav        4    1.80   2.0001        0.4220     1.2230
#> ... 5 more column(s) (sound.id, transect, distance, reference, cross.correlation)
#>  and 19 more row(s)
#> 
#> * 5 wave object(s) (as attributes): 
#> 10m_closed.wav, 10m_open.wav, 1m_open.wav, 30m_closed.wav, 30m_open.wav
#> 
#> * A data frame (check.results) with 25 rows generated by check_sels() (as an attribute)
#> 
#> Additional information:
#> * The selection table was created by song (see description in '?selection_table')
#> * 1 sampling rate(s) (in kHz): 22.05
#> * 1 bit depth(s): 16
#> * Created by warbleR < 1.1.21
```
