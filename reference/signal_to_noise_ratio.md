# Measure attenuation as signal-to-noise ratio

`signal_to_noise_ratio` measures attenuation as signal-to-noise ratio of
sounds referenced in an extended selection table.

## Usage

``` r
signal_to_noise_ratio(
  X,
  mar = NULL,
  cores = getOption("mc.cores", 1),
  pb = getOption("pb", TRUE),
  eq.dur = FALSE,
  noise.ref = c("adjacent", "custom"),
  snr.formula = 1,
  bp = "freq.range",
  hop.size = getOption("hop.size", 1),
  wl = getOption("wl", NULL),
  ovlp = getOption("ovlp", 0),
  path = getOption("sound.files.path", ".")
)
```

## Arguments

- X:

  Object of class 'data.frame', 'selection_table' or
  'extended_selection_table' (the last 2 classes are created by the
  function
  [`selection_table`](https://marce10.github.io/warbleR/reference/selection_table.html)
  from the warbleR package) with the test sound files' annotations
  (typically the output of
  [`align_test_files`](https://marce10.github.io/baRulho/reference/align_test_files.md)).
  Must contain the following columns: 1) "sound.files": name of the .wav
  files, 2) "selec": unique selection identifier (within a sound
  file), 3) "start": start time and 4) "end": end time of selections, 5)
  "bottom.freq": low frequency for bandpass, 6) "top.freq": high
  frequency for bandpass and 7) "sound.id": ID of sounds used to
  identify counterparts across distances (only needed for "custom" noise
  reference, see "noise.ref" argument). If the "sound.id" column is
  supplied then SNR is only computed for those rows with a sound.id
  different from "ambient", "start_marker" or "end_marker".

- mar:

  numeric vector of length 1. Specifies the margins adjacent to the
  start point of the annotation over which to measure ambient noise.

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

- eq.dur:

  Logical. Controls whether the ambient noise segment that is measured
  has the same duration to that of the sound (if `TRUE`. Default is
  `FALSE`). If `TRUE` then 'mar' and 'noise.ref' arguments are ignored.

- noise.ref:

  Character vector of length 1 to determined which noise segment must be
  used for measuring ambient noise. Two options are available:

  - `adjacent`: measure ambient noise right before test sounds (using
    argument 'mar' to define duration of ambient noise segments).

  - `custom`: measure ambient noise segments referenced in the selection
    table (labeled as 'ambient' in the 'sound.id' column). Those
    segments will be used to apply the same ambient noise reference to
    all sounds in a sound file. Therefore, at least one 'ambient'
    selection for each sound file must be provided. If several 'ambient'
    selections by sound file are supplied, then the root mean square of
    the amplitude envelope will be averaged across those selections.

- snr.formula:

  Integer vector of length 1. Selects the formula to be used to
  calculate the signal-to-noise ratio (S = signal , N = background
  noise):

  - `1`: ratio of S amplitude envelope root mean square to N amplitude
    envelope root mean square (`20 * log10(rms(env(S))/rms(env(N)))`) as
    described by Darden (2008).

  - `2`: ratio of the difference between S amplitude envelope root mean
    square and N amplitude envelope root mean square to N amplitude
    envelope root mean square
    (`20 * log10((rms(env(S)) - rms(env(N)))/rms(env(N)))`, as described
    by Dabelsteen et al (1993).

- bp:

  Numeric vector of length 2 giving the lower and upper limits of a
  frequency bandpass filter (in kHz). Alternatively, when set to
  'freq.range' (default), the function will use the 'bottom.freq' and
  'top.freq' for each sound as the bandpass range.

- hop.size:

  A numeric vector of length 1 specifying the time window duration (in
  ms). Default is 1 ms, which is equivalent to ~45 wl for a 44.1 kHz
  sampling rate. Ignored if 'wl' is supplied. Can be set globally for
  the current R session via the "hop.size" option (see
  [`options`](https://rdrr.io/r/base/options.html)).

- wl:

  A numeric vector of length 1 specifying the window length of the
  spectrogram, default is NULL. Ignored if `bp = NULL`. If supplied,
  'hop.size' is ignored. Note that lower values will increase time
  resolution, which is more important for amplitude ratios calculations.

- ovlp:

  Numeric vector of length 1 specifying the percentage of overlap
  between two consecutive windows, as in
  [`spectro`](https://rdrr.io/pkg/seewave/man/spectro.html). Default
  is 0. Only used for bandpass filtering. Can be set globally for the
  current R session via the "ovlp" option (see
  [`options`](https://rdrr.io/r/base/options.html)).

- path:

  Character string containing the directory path where the sound files
  are found. Only needed when 'X' is not an extended selection table. If
  not supplied the current working directory is used. Can be set
  globally for the current R session via the "sound.files.path" option
  (see [`options`](https://rdrr.io/r/base/options.html)).

## Value

Object 'X' with an additional column, 'signal.to.noise.ratio', with the
signal-to-noise ratio values (in dB).

## Details

Signal-to-noise ratio (SNR) measures sound amplitude level in relation
to ambient noise. Noise is measured on the background noise immediately
before the test sound. A general margin in which ambient noise will be
measured must be specified. Alternatively, a selection of ambient noise
can be used as reference (see 'noise.ref' argument). When margins
overlap with another sound nearby, SNR will be inaccurate, so margin
length should be carefully considered. Any SNR less than or equal to one
suggests background noise is equal to or overpowering the sound. The
function will measure signal-to-noise ratio within the supplied
frequency range (e.g. bandpass) of the reference signal ('bottom.freq'
and 'top.freq' columns in 'X') by default (that is, when
`bp = 'freq.range'`. SNR can be ~0 when both tail and signal have very
low amplitude.

## References

Araya-Salas, M., Grabarczyk, E. E., Quiroz-Oliva, M., Garcia-Rodriguez,
A., & Rico-Guevara, A. (2025). Quantifying degradation in animal
acoustic signals with the R package baRulho. Methods in Ecology and
Evolution, 00, 1-12. https://doi.org/10.1111/2041-210X.14481 Holland J,
Dabelsteen T, Pedersen SB, Paris AL (2001) Potential ranging cues
contained within the energetic pauses of transmitted wren song.
Bioacoustics 12(1):3-20. Darden, SK, Pedersen SB, Larsen ON, &
Dabelsteen T. (2008). Sound transmission at ground level in a
short-grass prairie habitat and its implications for long-range
communication in the swift fox \*Vulpes velox\*. The Journal of the
Acoustical Society of America, 124(2), 758-766.

## See also

[`excess_attenuation`](https://marce10.github.io/baRulho/reference/excess_attenuation.md)

Other quantify degradation:
[`blur_ratio()`](https://marce10.github.io/baRulho/reference/blur_ratio.md),
[`detection_distance()`](https://marce10.github.io/baRulho/reference/detection_distance.md),
[`envelope_correlation()`](https://marce10.github.io/baRulho/reference/envelope_correlation.md),
[`plot_blur_ratio()`](https://marce10.github.io/baRulho/reference/plot_blur_ratio.md),
[`plot_degradation()`](https://marce10.github.io/baRulho/reference/plot_degradation.md),
[`set_reference_sounds()`](https://marce10.github.io/baRulho/reference/set_reference_sounds.md),
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

  # using measure ambient noise reference selections
  signal_to_noise_ratio(X = test_sounds_est, mar = 0.05, noise.ref = "custom")

  # using margin for ambient noise of 0.05 and adjacent measure ambient noise reference
  signal_to_noise_ratio(X = test_sounds_est, mar = 0.05, noise.ref = "adjacent")
}
#> Object of class 'extended_selection_table'
#> * The output of the following call:
#> `[.extended_selection_table`(X = test_sounds_est, i = order(test_sounds_est$sound.id,  test_sounds_est$transect, decreasing = FALSE))
#> 
#> Contains: 
#> *  A selection table data frame with 25 row(s) and 10 columns:
#>      sound.files       selec   start      end   bottom.freq   top.freq
#> ---  ---------------  ------  ------  -------  ------------  ---------
#> 8    10m_closed.wav        1    0.05   0.2000        1.3333     2.6667
#> 23   30m_closed.wav        1    0.05   0.2000        1.3333     2.6667
#> 13   10m_open.wav          1    0.05   0.2000        1.3333     2.6667
#> 18   1m_open.wav           1    0.05   0.2000        1.3333     2.6667
#> 28   30m_open.wav          1    0.05   0.2000        1.3333     2.6667
#> 11   10m_closed.wav        4    1.80   2.0001        0.4220     1.2230
#> ... 4 more column(s) (sound.id, transect, distance, signal.to.noise.ratio)
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
