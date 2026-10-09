# Measure full spectrum sound noise profiles

`noise_profile()` measures full spectrum sound pressure levels (i.e.
noise profiles) in sound files or extended selection tables.

## Usage

``` r
noise_profile(
  X = NULL,
  files = NULL,
  mar = NULL,
  noise.ref = c("adjacent", "custom"),
  cores = getOption("mc.cores", 1),
  pb = getOption("pb", TRUE),
  path = getOption("sound.files.path", "."),
  bp = NULL,
  hop.size = getOption("hop.size", 1),
  wl = getOption("wl", NULL),
  PSD = FALSE,
  norm = TRUE,
  dB = c("A", "B", "C", "D", "ITU", "max0"),
  averaged = TRUE
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
  identify counterparts across distances (needed for `"custom"` noise
  reference, see the `noise.ref` argument). Default `NULL`.

- files:

  Character vector with names of wave files to be analyzed. Files must
  be found in the supplied `path` (or in the working directory if `path`
  is not supplied). Default `NULL`.

- mar:

  Numeric vector of length 1. Specifies the margins adjacent to the
  start and end points of a selection over which to measure ambient
  noise. Required if `X` is supplied, and ignored if `X` is not
  supplied. Default `NULL`.

- noise.ref:

  Character vector of length 1 determining which noise segment must be
  used for measuring ambient noise. Ignored if `X` is not supplied. Two
  options are available:

  - **`adjacent`**: measure ambient noise right before the sound (using
    the `mar` argument to define duration of ambient noise segments).

  - **`custom`**: measure ambient noise segments referenced in the
    selection table (labeled as `"ambient"` in the `sound.id` column).

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

- path:

  Character string containing the directory path where the sound files
  are found. Only needed when `X` is not an extended selection table. If
  not supplied the current working directory is used. Can be set
  globally for the current R session via the `"sound.files.path"` option
  (see [`options()`](https://rdrr.io/r/base/options.html)).

- bp:

  Numeric vector of length 2 giving the lower and upper limits of a
  frequency bandpass filter (in kHz). Default `NULL`.

- hop.size:

  Numeric vector of length 1 specifying the time window duration (in
  ms). Default `1` ms, which is equivalent to ~45 `wl` for a 44.1 kHz
  sampling rate. Ignored if `wl` is supplied. Can be set globally for
  the current R session via the `"hop.size"` option (see
  [`options()`](https://rdrr.io/r/base/options.html)).

- wl:

  Numeric vector of length 1 specifying the window length of the
  spectrogram. Default `NULL`. Ignored if `bp = NULL`. If supplied,
  `hop.size` is ignored. Note that lower values will increase time
  resolution, which is more important for amplitude ratio calculations.

- PSD:

  Logical to control whether the Probability Mass Function (the
  probability distribution of frequencies) is returned. See
  [`seewave::meanspec()`](https://rdrr.io/pkg/seewave/man/meanspec.html).
  Default `FALSE`.

- norm:

  Logical to control whether amplitude values are normalized (divided by
  the maximum) so the highest value is 1. See
  [`seewave::meanspec()`](https://rdrr.io/pkg/seewave/man/meanspec.html).
  Default `TRUE`.

- dB:

  Character string of length 1 specifying the type of dB to return:
  `"max0"` for a maximum dB value at 0, `"A"`, `"B"`, `"C"`, `"D"`, and
  `"ITU"` for common dB weights. See
  [`seewave::meanspec()`](https://rdrr.io/pkg/seewave/man/meanspec.html).
  Default `"A"`.

- averaged:

  Logical to control if frequency spectra are averaged within a sound
  file. Default `TRUE`.

## Value

A `data.frame` containing the frequency spectra for each sound file or
wave object (if `X` is supplied and is of class
`extended_selection_table`).

## Details

The function estimates full spectrum sound pressure levels (i.e. noise
profiles) of ambient noise. This can be done on data frames/(extended)
selection tables (using the segments containing no target sound or the
`"ambient"` sound ID) or over complete sound files in the working
directory (or supplied `path`). The function uses
[`seewave::meanspec()`](https://rdrr.io/pkg/seewave/man/meanspec.html)
internally to calculate frequency spectra.

## References

Araya-Salas, M., Grabarczyk, E. E., Quiroz-Oliva, M., Garcia-Rodriguez,
A., & Rico-Guevara, A. (2025). Quantifying degradation in animal
acoustic signals with the R package baRulho. Methods in Ecology and
Evolution, 00, 1-12. https://doi.org/10.1111/2041-210X.14481.

## See also

[`excess_attenuation`](https://marce10.github.io/baRulho/reference/excess_attenuation.md)

Other miscellaneous:
[`add_noise()`](https://marce10.github.io/baRulho/reference/add_noise.md),
[`attenuation()`](https://marce10.github.io/baRulho/reference/attenuation.md)

## Author

Marcelo Araya-Salas (<marcelo.araya@ucr.ac.cr>)

## Examples

``` r
{
  # load example data
  data("test_sounds_est")

  # measure on custom noise reference
  noise_profile(X = test_sounds_est, mar = 0.01, pb = FALSE, noise.ref = "custom")

  # remove noise selections so noise is measured right before the signals
  pe <- test_sounds_est[test_sounds_est$sound.id != "ambient", ]

  noise_profile(X = pe, mar = 0.01, pb = FALSE, noise.ref = "adjacent")
}
#>       sound.files      freq           amp
#> 1  10m_closed.wav  1.002273  -1.478756370
#> 2  10m_closed.wav  2.004545  -8.554912560
#> 3  10m_closed.wav  3.006818 -15.653879395
#> 4  10m_closed.wav  4.009091 -20.594564397
#> 5  10m_closed.wav  5.011364 -24.596826317
#> 6  10m_closed.wav  6.013636 -27.741536121
#> 7  10m_closed.wav  7.015909 -28.631141207
#> 8  10m_closed.wav  8.018182 -29.343274598
#> 9  10m_closed.wav  9.020455 -31.581649783
#> 10 10m_closed.wav 10.022727 -37.142756979
#> 11   10m_open.wav  1.002273  -1.470944432
#> 12   10m_open.wav  2.004545  -5.716210091
#> 13   10m_open.wav  3.006818 -14.060376878
#> 14   10m_open.wav  4.009091 -21.637019767
#> 15   10m_open.wav  5.011364 -28.007595476
#> 16   10m_open.wav  6.013636 -31.998254047
#> 17   10m_open.wav  7.015909 -33.040067144
#> 18   10m_open.wav  8.018182 -33.446761775
#> 19   10m_open.wav  9.020455 -36.054384166
#> 20   10m_open.wav 10.022727 -40.262096085
#> 21    1m_open.wav  1.002273   0.001384514
#> 22    1m_open.wav  2.004545  -4.328099010
#> 23    1m_open.wav  3.006818 -12.731830547
#> 24    1m_open.wav  4.009091 -20.519891620
#> 25    1m_open.wav  5.011364 -24.786032150
#> 26    1m_open.wav  6.013636 -27.331145913
#> 27    1m_open.wav  7.015909 -27.522403062
#> 28    1m_open.wav  8.018182 -27.601169568
#> 29    1m_open.wav  9.020455 -30.958903505
#> 30    1m_open.wav 10.022727 -36.213758740
#> 31 30m_closed.wav  1.002273  -5.018620398
#> 32 30m_closed.wav  2.004545 -19.660800804
#> 33 30m_closed.wav  3.006818 -26.467263845
#> 34 30m_closed.wav  4.009091 -29.190019017
#> 35 30m_closed.wav  5.011364 -33.368022266
#> 36 30m_closed.wav  6.013636 -37.485873353
#> 37 30m_closed.wav  7.015909 -39.725404754
#> 38 30m_closed.wav  8.018182 -42.088478865
#> 39 30m_closed.wav  9.020455 -44.455196777
#> 40 30m_closed.wav 10.022727 -47.253627732
#> 41   30m_open.wav  1.002273  -0.794205765
#> 42   30m_open.wav  2.004545  -5.394132145
#> 43   30m_open.wav  3.006818 -16.362837354
#> 44   30m_open.wav  4.009091 -21.711879855
#> 45   30m_open.wav  5.011364 -26.204882614
#> 46   30m_open.wav  6.013636 -26.312227870
#> 47   30m_open.wav  7.015909 -28.104324225
#> 48   30m_open.wav  8.018182 -28.673914281
#> 49   30m_open.wav  9.020455 -32.218532962
#> 50   30m_open.wav 10.022727 -34.076725570
```
