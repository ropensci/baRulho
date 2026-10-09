# Measure reverberations as tail-to-signal ratio

`tail_to_signal_ratio()` measures reverberations as the tail-to-signal
ratio of sounds referenced in an extended selection table.

## Usage

``` r
tail_to_signal_ratio(
  X,
  mar,
  cores = getOption("mc.cores", 1),
  pb = getOption("pb", TRUE),
  tsr.formula = 1,
  bp = "freq.range",
  hop.size = getOption("hop.size", 1),
  wl = getOption("wl", NULL),
  ovlp = getOption("ovlp", 0),
  path = getOption("sound.files.path", ".")
)
```

## Arguments

- X:

  Object of class `data.frame`, `selection_table`, or
  `extended_selection_table` (the last 2 classes are created by
  [`warbleR::selection_table()`](https://marce10.github.io/warbleR/reference/selection_table.html)
  from the **warbleR** package) with the test sound files' annotations.
  Must contain the following columns: 1) `sound.files`: name of the
  `.wav` files, 2) `selec`: unique selection identifier (within a sound
  file), 3) `start`: start time and 4) `end`: end time of selections, 5)
  `bottom.freq`: low frequency for bandpass, 6) `top.freq`: high
  frequency for bandpass, and 7) `sound.id`: ID of sounds used to
  identify counterparts across distances.

- mar:

  Numeric vector of length 1. Specifies the margins adjacent to the end
  of the sound over which to measure tail power.

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

- tsr.formula:

  Integer vector of length 1. Determines the formula to be used to
  calculate the tail-to-signal ratio (S = signal, T = tail, N =
  background noise):

  - **`1`**: ratio of T amplitude envelope root mean square to S
    amplitude envelope root mean square
    (`20 * log10(rms(env(T))/rms(env(S)))`) as described by Dabelsteen
    et al. (1993).

  - **`2`**: ratio of T amplitude envelope root mean square to N
    amplitude envelope root mean square
    (`20 * log10(rms(env(T))/rms(env(N)))`). N is measured in the margin
    right before the sound. So `tsr.formula = 2` actually measures
    tail-to-noise ratio.

- bp:

  Numeric vector of length 2 giving the lower and upper limits of a
  frequency bandpass filter (in kHz). Alternatively, when set to
  `"freq.range"` (default), the function will use the `bottom.freq` and
  `top.freq` for each sound as the bandpass range.

- hop.size:

  Numeric vector of length 1 specifying the time window duration (in
  ms). Default `1` ms, which is equivalent to ~45 `wl` for a 44.1 kHz
  sampling rate. Ignored if `wl` is supplied. Can be set globally for
  the current R session via the `"hop.size"` option (see
  [`options()`](https://rdrr.io/r/base/options.html)). Note that this
  might be internally adjusted if the number of samples in the tail is
  lower than `hop.size`.

- wl:

  Numeric vector of length 1 specifying the window length of the
  spectrogram. Default `NULL`. Ignored if `bp = NULL`. If supplied,
  `hop.size` is ignored. Note that lower values will increase time
  resolution, which is more important for amplitude calculations.

- ovlp:

  Numeric vector of length 1 specifying the percentage of overlap
  between two consecutive windows, as in
  [`seewave::spectro()`](https://rdrr.io/pkg/seewave/man/spectro.html).
  Default `0`. Only used for bandpass filtering. Can be set globally for
  the current R session via the `"ovlp"` option (see
  [`options()`](https://rdrr.io/r/base/options.html)).

- path:

  Character string containing the directory path where the sound files
  are found. Only needed when `X` is not an extended selection table. If
  not supplied the current working directory is used. Can be set
  globally for the current R session via the `"sound.files.path"` option
  (see [`options()`](https://rdrr.io/r/base/options.html)).

## Value

Object `X` with an additional column, `tail.to.signal.ratio`, with the
tail-to-signal ratio values (in dB).

## Details

Tail-to-signal ratio (TSR) measures the ratio of power in the tail of
reverberations to that in the test sound. A general margin in which the
reverberation tail will be measured must be specified. The function will
measure TSR within the supplied frequency range (e.g. bandpass) of the
reference sound (`bottom.freq` and `top.freq` columns in `X`). Two
methods for computing reverberations are provided (see the `tsr.formula`
argument). Note that `tsr.formula = 2` is not equivalent to the original
description of TSR in Dabelsteen et al. (1993), and is better referred
to as tail-to-noise ratio. Tail-to-signal ratio values are typically
negative, as signals tend to have higher power than that in the
reverberating tail. TSR can be ~0 when both tail and signal have very
low amplitude.

## References

Araya-Salas, M., Grabarczyk, E. E., Quiroz-Oliva, M., Garcia-Rodriguez,
A., & Rico-Guevara, A. (2025). Quantifying degradation in animal
acoustic signals with the R package baRulho. Methods in Ecology and
Evolution, 00, 1-12. https://doi.org/10.1111/2041-210X.14481 Darden, SK,
Pedersen SB, Larsen ON, & Dabelsteen T. (2008). Sound transmission at
ground level in a short-grass prairie habitat and its implications for
long-range communication in the swift fox *Vulpes velox*. The Journal of
the Acoustical Society of America, 124(2), 758-766. Mathevon, N.,
Dabelsteen, T., & Blumenrath, S. H. (2005). Are high perches in the
blackcap Sylvia atricapilla song or listening posts? A sound
transmission study. The Journal of the Acoustical Society of America,
117(1), 442-449.

## See also

[`excess_attenuation()`](https://marce10.github.io/baRulho/reference/excess_attenuation.md),
for a related degradation metric.

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
[`spectrum_correlation()`](https://marce10.github.io/baRulho/reference/spectrum_correlation.md)

## Author

Marcelo Araya-Salas (<marcelo.araya@ucr.ac.cr>)

## Examples

``` r
{
  # load example data

  data("test_sounds_est")

  # set global options
  options(pb = FALSE)

  # using margin for noise of 0.01
  tsr <- tail_to_signal_ratio(X = test_sounds_est, mar = 0.01)

  # use tsr.formula 2 which is equivalent to tail-to-noise ratio
  tsr <- tail_to_signal_ratio(X = test_sounds_est, mar = 0.01, tsr.formula = 2)
}
```
