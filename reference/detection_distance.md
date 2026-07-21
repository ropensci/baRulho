# Measure detection distance of sound

`detection_distance` detection distance of sounds.

## Usage

``` r
detection_distance(
  X,
  cores = getOption("mc.cores", 1),
  pb = getOption("pb", TRUE),
  hop.size = getOption("hop.size", 11.6),
  wl = getOption("wl", NULL),
  path = getOption("sound.files.path", "."),
  spl = NULL,
  spl.cutoff = NULL,
  temp = 20,
  rh = 60,
  pa = 101325,
  hab.att.coef = 0.02,
  max.distance = 1000,
  resolution = 0.1,
  subtract.bgn = TRUE,
  envelope = c("abs", "hil"),
  mar = NULL
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

- path:

  Character string containing the directory path where the sound files
  are found. Only needed when 'X' is not an extended selection table. If
  not supplied the current working directory is used. Can be set
  globally for the current R session via the "sound.files.path" option
  (see [`options`](https://rdrr.io/r/base/options.html)).

- spl:

  A numeric vector of length 1 specifying the sound pressure level of
  sounds. If not supplied then it will be measured from the sounds
  themselves.

- spl.cutoff:

  A numeric vector of length 1 specifying the sound pressure level
  cutoff to define if the sound is no longer detected. Ideally it should
  be estimated based on the sound detection threshold of the species.

- temp:

  Numeric vector of length 1 with frequency (in Celsius). Default is 20.

- rh:

  Numeric vector of length 1 with relative humidity (in percentage).
  Default is 60.

- pa:

  Numeric vector of length 1 with ambient pressure in Pa (standard:
  101325, default). Used for Atmospheric attenuation.

- hab.att.coef:

  Attenuation coefficient of the habitat (in dB/kHz/m).

- max.distance:

  Numeric vector of length 1 with the maximum distance (in m) at which
  detection would be evaluated. Note that the function calculates the
  expected sound pressure level values along a vector of distances to
  find the distance at which the expected sound pressure level equates
  'spl.cutoff'. Default is 1000 m.

- resolution:

  Numeric vector of length 1 with the distance resolution (in m) for
  estimated detection distance. Higher resolutions take longer to
  estimate. Default is 0.1 m.

- subtract.bgn:

  Logical argument to control if SPL from background noise is excluded
  from the measured signal SPL. Default is `FALSE`.

- envelope:

  Character string vector with the method to calculate amplitude
  envelopes (in which SPL is measured, only used if 'spl' is not
  supplied), as in [`env`](https://rdrr.io/pkg/seewave/man/env.html).
  Must be either 'abs' (absolute envelope, default) or 'hil' (Hilbert
  transformation).

- mar:

  numeric vector of length 1. Specifies the margins adjacent to the
  start and end points of selection over which to measure background
  noise. This is required to subtract background noise sound pressure
  level (so only needed when 'subtract.bgn = TRUE').

## Value

Object 'X' with an additional column, 'detection.distance', containing
the computed detection distances (in m).

## Details

The function computes the maximum distance at which a sound would be
detected, which is calculated as the distance in which the sound
pressure level (SPL) goes below the specified SPL cutoff
('spl.cutoff')). This is returned as an additional column
'detection.distance' (in m). The function uses internally
[`attenuation`](https://marce10.github.io/baRulho/reference/attenuation.md)
to estimate SPL at increasing values until it reaches the defined
cutoff. The peak frequency (calculated on the power spectrum of the
reference sound) of the reference sound for each sound ID is used as the
carrier frequency for distance estimation. The sound recorded at the
lowest distance is used as reference. **This function assumes that all
recordings have been made at the same recording volume**.

## References

Araya-Salas, M., Grabarczyk, E. E., Quiroz-Oliva, M., Garcia-Rodriguez,
A., & Rico-Guevara, A. (2025). Quantifying degradation in animal
acoustic signals with the R package baRulho. Methods in Ecology and
Evolution, 00, 1-12. https://doi.org/10.1111/2041-210X.14481

Clark, C.W., Marler, P. & Beeman K. (1987). Quantitative analysis of
animal vocal phonology: an application to Swamp Sparrow song. Ethology.
76:101-115.

## See also

[`attenuation`](https://marce10.github.io/baRulho/reference/attenuation.md)

Other quantify degradation:
[`blur_ratio()`](https://marce10.github.io/baRulho/reference/blur_ratio.md),
[`envelope_correlation()`](https://marce10.github.io/baRulho/reference/envelope_correlation.md),
[`plot_blur_ratio()`](https://marce10.github.io/baRulho/reference/plot_blur_ratio.md),
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
if (FALSE) { # \dontrun{
# load example data
data("test_sounds_est")

# add reference to X
X <- set_reference_sounds(X = test_sounds_est)

detection_distance(X = X[X$distance %in% c(1, 10), ], spl.cutoff = 5, mar = 0.05)
} # }
```
