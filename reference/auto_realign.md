# Fix small misalignments in the time position test sounds

`auto_realign` fixes small misalignments in the time position of test
sounds in an extended selection table using spectrographic
cross-correlation

## Usage

``` r
auto_realign(
  X,
  Y,
  cores = getOption("mc.cores", 1),
  pb = getOption("pb", TRUE),
  hop.size = getOption("hop.size", 11.6),
  wl = getOption("wl", NULL),
  ovlp = getOption("ovlp", 90),
  wn = c("hanning", "hamming", "bartlett", "blackman", "flattop", "rectangle"),
  bp = NULL
)
```

## Arguments

- X:

  object of class 'extended_selection_table' (created by the function
  [`selection_table`](https://marce10.github.io/warbleR/reference/selection_table.html)
  from the warbleR package) with the test sound files' annotations to be
  aligned. Must contain the following columns: 1) "sound.files": name of
  the .wav files, 2) "selec": unique selection identifier (within a
  sound file), 3) "start": start time and 4) "end": end time of
  selections, 5) "bottom.freq": low frequency for bandpass, 6)
  "top.freq": high frequency for bandpass and 7) "sound.id": ID of
  sounds used to identify counterparts across distances. Each sound must
  have a unique ID within a given distance.. The object must include the
  following additional columns: 'sound.id', 'bottom.freq' and
  'top.freq'.

- Y:

  object of class 'extended_selection_table' (a class created by the
  function
  [`selection_table`](https://marce10.github.io/warbleR/reference/selection_table.html)
  from the warbleR package) with the master sound file annotations. This
  should be the same data than that was used for finding the position of
  markers in
  [`find_markers`](https://marce10.github.io/baRulho/reference/find_markers.md).
  It should also contain a 'sound.id' column.

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

  Numeric vector of length 1 specifying the percentage of overlap
  between two consecutive windows, as in
  [`spectro`](https://rdrr.io/pkg/seewave/man/spectro.html). Default
  is 90. High values slow down the function but produce more accurate
  results. Can be set globally for the current R session via the "ovlp"
  option (see [`options`](https://rdrr.io/r/base/options.html)).

- wn:

  A character vector of length 1 specifying the window name as in
  [`ftwindow`](https://rdrr.io/pkg/seewave/man/ftwindow.html).

- bp:

  Numeric vector of length 2 giving the lower and upper limits of a
  frequency bandpass filter (in kHz). Default is `NULL`.

## Value

Object 'X' in which time parameters (columns 'start' and 'end') have
been tailored to more closely match the start and end of the reference
sound.

## Details

Precise alignment is crucial for downstream measures of sound
degradation. This function uses spectrogram cross-correlation to improve
the time position alignment of test sounds. The master sound file is
used as reference. The function calls warbleR's
[`cross_correlation`](https://marce10.github.io/warbleR/reference/cross_correlation.html)
internally to align sounds using cross-correlation. The output extended
selection table contains the new start and end values after alignment.
**Note that 1) this function only works to further improve alignments if
the estimated position of the test sound is already close to the actual
position and 2) both 'X' and 'Y' must be extended selection tables sensu
[`selection_table`](https://marce10.github.io/warbleR/reference/selection_table.html)**.
The function might not work properly with annotations with a small
frequency range (e.g. pure tones).

## References

Araya-Salas, M., Grabarczyk, E. E., Quiroz-Oliva, M., Garcia-Rodriguez,
A., & Rico-Guevara, A. (2025). Quantifying degradation in animal
acoustic signals with the R package baRulho. Methods in Ecology and
Evolution, 00, 1-12. https://doi.org/10.1111/2041-210X.14481

Clark, C.W., Marler, P. & Beeman K. (1987). Quantitative analysis of
animal vocal phonology: an application to Swamp Sparrow song. Ethology.
76:101-115.

## See also

[`blur_ratio`](https://marce10.github.io/baRulho/reference/blur_ratio.md),
[`cross_correlation`](https://marce10.github.io/warbleR/reference/cross_correlation.html)

Other test sound alignment:
[`align_test_files()`](https://marce10.github.io/baRulho/reference/align_test_files.md),
[`find_markers()`](https://marce10.github.io/baRulho/reference/find_markers.md),
[`manual_realign()`](https://marce10.github.io/baRulho/reference/manual_realign.md),
[`plot_aligned_sounds()`](https://marce10.github.io/baRulho/reference/plot_aligned_sounds.md)

## Author

Marcelo Araya-Salas (<marcelo.araya@ucr.ac.cr>)

## Examples

``` r
{
  # load example data
  data("test_sounds_est")
  data("master_est")
  
  # create "unaligned_test_sounds_est" by
  # adding error to "test_sounds_est" start and end
  unaligned_test_sounds_est <- test_sounds_est
  set.seed(123)
  noise_time <- sample(c(0.009, -0.01, 0.03, -0.03, 0, 0.07, -0.007),
  nrow(unaligned_test_sounds_est),
  replace = TRUE)
  
  attr(unaligned_test_sounds_est, "check.res")$start <- 
  unaligned_test_sounds_est$start <- 
  unaligned_test_sounds_est$start + noise_time
  attr(unaligned_test_sounds_est, "check.res")$end <- 
  unaligned_test_sounds_est$end  <- 
  unaligned_test_sounds_est$end + noise_time

# re align
realigned_est <- auto_realign(X = unaligned_test_sounds_est, Y = master_est)
}
```
